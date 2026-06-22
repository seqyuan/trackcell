"""
Jupyter-native interactive ROI selection for TrackCell.

This module intentionally does **not** use napari/Qt.  It relies on matplotlib's
widget selectors, which are safer in Jupyter notebooks when used with an
interactive backend such as ipympl (``%matplotlib widget``).

Key features:

* **No blocking ``input()``** — ROIs are auto-named (``ROI_1``, ``ROI_2``, …)
  with a configurable prefix, so the interactive workflow is never interrupted.
* **Inline toolbar** — clickable **Rect** / **Ellipse** / **Lasso** buttons
  plus **Clear** and **Undo** — no need to remember keyboard shortcuts.
* **Keyboard shortcuts** still work (``r``/``e``/``l``) for power users.
* **Selected cells are highlighted** in real-time with coloured scatter points.
* **Post-hoc rename** via ``selector.rename_roi(old_name, new_name)``.
"""

from __future__ import annotations

from typing import Optional, Dict, List, Literal, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.widgets import (
    RectangleSelector,
    EllipseSelector,
    LassoSelector,
    Button,
)
from shapely.geometry import Polygon, Point
from shapely import wkt

from .plot import spatial_cell, spatial_squarebin, _resolve_library_id


ModeType = Literal["auto", "cellbin", "squarebin"]
SelectorMode = Literal["rectangle", "ellipse", "lasso"]


# ---------------------------------------------------------------------------
#  helpers
# ---------------------------------------------------------------------------


def _detect_mode(adata, library_id: str, mode: Optional[str], basis: str) -> str:
    """Detect cellbin vs squarebin from TrackCell/AnnData spatial fields."""
    if mode is not None and mode not in ("auto", "cellbin", "squarebin"):
        raise ValueError("`mode` must be 'auto', 'cellbin', or 'squarebin'.")
    if mode in ("cellbin", "squarebin"):
        return mode

    spatial_info = adata.uns.get("spatial", {}).get(library_id, {})
    if spatial_info.get("geometries") is not None or "geometry" in adata.obs.columns:
        return "cellbin"
    if basis in adata.obsm:
        return "squarebin"
    raise ValueError(
        "Cannot auto-detect mode. Expected cell geometries for cellbin or "
        f"adata.obsm['{basis}'] for squarebin."
    )


def _valid_ids(adata, spatial_info: dict, mode: str, basis: str) -> pd.Index:
    if mode == "cellbin":
        geometries = spatial_info.get("geometries", None)
        if geometries is not None:
            return geometries.index.intersection(adata.obs_names)
        if "geometry" in adata.obs.columns:
            return adata.obs_names[adata.obs["geometry"].notna()]
        raise ValueError(
            "No cell geometries found. Expected adata.uns['spatial'][library_id]"
            "['geometries'] or adata.obs['geometry'] WKT strings."
        )
    if basis not in adata.obsm:
        raise ValueError(f"`adata.obsm['{basis}']` not found for squarebin mode.")
    return adata.obs_names


def _get_geometry(adata, cid: str, geometries, use_wkt: bool):
    try:
        if use_wkt:
            raw = adata.obs.loc[cid, "geometry"]
            if pd.isna(raw):
                return None
            return wkt.loads(raw) if isinstance(raw, str) else raw
        if cid not in geometries.index:
            return None
        return geometries.loc[cid, "geometry"]
    except Exception:
        return None


def _ellipse_to_polygon(center, width, height, angle=0.0, n=64):
    """Sample an ellipse into a polygon of ``n`` vertices (for intersection tests)."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = center[0] + (width / 2) * np.cos(t)
    y = center[1] + (height / 2) * np.sin(t)
    if angle != 0:
        theta = np.deg2rad(angle)
        cos_a, sin_a = np.cos(theta), np.sin(theta)
        x -= center[0]
        y -= center[1]
        xr = x * cos_a - y * sin_a
        yr = x * sin_a + y * cos_a
        x = xr + center[0]
        y = yr + center[1]
    return np.column_stack([x, y])


# ---------------------------------------------------------------------------
#  colour palette for ROI highlights
# ---------------------------------------------------------------------------

_HIGHLIGHT_COLORS = [
    "#e6194b",  # red
    "#3b75af",  # blue
    "#44aa44",  # green
    "#ff8c00",  # orange
    "#911eb4",  # purple
    "#46f0f0",  # cyan
    "#f032e6",  # magenta
    "#d2f53c",  # lime
    "#fabebe",  # pink
    "#008080",  # teal
]


# =============================================================================
#  RegionSelector
# =============================================================================


class RegionSelector:
    """
    Controller returned by :func:`select_regions`.

    All shape modes (rectangle, ellipse, lasso) are available via both the
    **inline toolbar** (clickable buttons at the bottom of the figure) and
    **keyboard shortcuts** (``r``/``e``/``l`` while the figure has focus).

    ROIs are auto-named (``ROI_1``, ``ROI_2``, …).  Use
    :meth:`rename_roi` to change names after the fact.

    Attributes
    ----------
    rois : dict
        Mapping of ROI names to lists of observation IDs.
    polygons : dict
        Mapping of ROI names to arrays of ROI boundary vertices in data coords.
    ax, fig
        Matplotlib axes and figure.
    """

    # --- visual defaults ---------------------------------------------------
    _EDGE_COLOR = "red"
    _FACE_COLOR = "none"
    _LINEWIDTH = 2.0

    def __init__(
        self,
        adata,
        ax: plt.Axes,
        mode: str,
        basis: str,
        valid_ids: pd.Index,
        geometries=None,
        use_wkt: bool = False,
        key_added: str = "ROI",
        inplace: bool = True,
        roi_prefix: str = "ROI",
    ):
        self.adata = adata
        self.ax = ax
        self.fig = ax.figure
        self.mode = mode
        self.basis = basis
        self.valid_ids = valid_ids
        self.geometries = geometries
        self.use_wkt = use_wkt
        self.key_added = key_added
        self.inplace = inplace
        self.roi_prefix = roi_prefix

        # state
        self._current_mode: SelectorMode = "rectangle"
        self.rois: Dict[str, List] = {}
        self.polygons: Dict[str, np.ndarray] = {}
        self.patches: List[Any] = []
        self._roi_counter = 0
        self._highlight_collections: List[Any] = []  # scatter highlight artists

        # selectors (all created, only one active at a time)
        self._selectors: Dict[str, Any] = {}
        self._connect_all_selectors()

        # inline toolbar
        self._build_toolbar()

    # ------------------------------------------------------------------
    #  selector wiring
    # ------------------------------------------------------------------

    def _connect_all_selectors(self) -> None:
        """Create all three selectors and then activate the default."""
        # rectangle
        self._selectors["rectangle"] = RectangleSelector(
            self.ax,
            self._on_rectangle,
            useblit=True,
            button=[1],
            minspanx=1,
            minspany=1,
            spancoords="data",
            interactive=True,
        )
        # ellipse
        self._selectors["ellipse"] = EllipseSelector(
            self.ax,
            self._on_ellipse,
            useblit=True,
            button=[1],
            minspanx=1,
            minspany=1,
            spancoords="data",
            interactive=True,
        )
        # lasso
        self._selectors["lasso"] = LassoSelector(
            self.ax,
            self._on_lasso,
            useblit=True,
            props={"color": self._EDGE_COLOR, "linewidth": self._LINEWIDTH},
        )

        # all start inactive except the default
        self._set_active("rectangle")

        # keyboard toggle
        self._key_cid = self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    def _set_active(self, mode: SelectorMode) -> None:
        """Activate one selector and deactivate all others."""
        for name, sel in self._selectors.items():
            sel.set_active(name == mode)
        self._current_mode = mode
        self._update_title()
        self._update_toolbar_highlight()

    def _update_title(self) -> None:
        suffix = f" | inplace={self.inplace}" if not self.inplace else ""
        self.ax.set_title(
            f"Mode: [{self._current_mode}]  |  {len(self.rois)} ROI(s){suffix}"
        )

    def _on_key(self, event) -> None:
        if event.key in ("r", "e", "l"):
            mode_map = {"r": "rectangle", "e": "ellipse", "l": "lasso"}
            self._set_active(mode_map[event.key])

    # ------------------------------------------------------------------
    #  inline toolbar
    # ------------------------------------------------------------------

    def _build_toolbar(self) -> None:
        """Build clickable buttons at the bottom of the figure."""
        self.fig.subplots_adjust(bottom=0.13)

        y0, bw, bh = 0.02, 0.09, 0.055
        spacing = 0.125

        # --- Mode toggles ---
        mode_defs = [
            ("■ Rect", "rectangle"),
            ("● Ellipse", "ellipse"),
            ("✎ Lasso", "lasso"),
        ]
        self._mode_buttons: Dict[str, Any] = {}
        for i, (label, mode_name) in enumerate(mode_defs):
            ax_btn = self.fig.add_axes([0.06 + i * spacing, y0, bw, bh])
            btn = Button(ax_btn, label)
            # capture both mode_name and self by value via default args
            btn.on_clicked(lambda e, m=mode_name, s=self: s._toolbar_set_mode(m))
            self._mode_buttons[mode_name] = btn

        # --- Clear ---
        clear_ax = self.fig.add_axes([0.06 + 3 * spacing, y0, bw, bh])
        self._clear_btn = Button(clear_ax, "✕ Clear", color="lightcoral", hovercolor="coral")
        self._clear_btn.on_clicked(lambda e: self.clear())

        # --- Undo ---
        undo_ax = self.fig.add_axes([0.06 + 4 * spacing, y0, bw, bh])
        self._undo_btn = Button(undo_ax, "↩ Undo")
        self._undo_btn.on_clicked(lambda e: self.undo())

        self._update_toolbar_highlight()

    def _toolbar_set_mode(self, mode_name: str) -> None:
        """Called when a mode button is clicked."""
        self._set_active(mode_name)

    def _update_toolbar_highlight(self) -> None:
        """Highlight the currently active mode button."""
        active_color = "#2563eb"  # blue
        inactive_color = "#e0e0e0"
        for mode_name, btn in self._mode_buttons.items():
            btn.color = active_color if mode_name == self._current_mode else inactive_color

    # ------------------------------------------------------------------
    #  shape callbacks
    # ------------------------------------------------------------------

    def _on_rectangle(self, eclick, erelease) -> None:
        if eclick.xdata is None or erelease.xdata is None:
            return
        x0, x1 = sorted([eclick.xdata, erelease.xdata])
        y0, y1 = sorted([eclick.ydata, erelease.ydata])
        vertices = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
        self._finish_roi(vertices, shape="rectangle")

    def _on_ellipse(self, eclick, erelease) -> None:
        if eclick.xdata is None or erelease.xdata is None:
            return
        center = (
            (eclick.xdata + erelease.xdata) / 2,
            (eclick.ydata + erelease.ydata) / 2,
        )
        width = abs(erelease.xdata - eclick.xdata)
        height = abs(erelease.ydata - eclick.ydata)
        vertices = _ellipse_to_polygon(center, width, height)
        self._finish_roi(vertices, shape="ellipse")

    def _on_lasso(self, vertices) -> None:
        vertices = np.asarray(vertices, dtype=float)
        if vertices.shape[0] < 3:
            return
        self._finish_roi(vertices, shape="lasso")

    # ------------------------------------------------------------------
    #  common finish path  (extract → name → save → draw)
    # ------------------------------------------------------------------

    def _finish_roi(self, vertices: np.ndarray, shape: str) -> None:
        """Shared pipeline: extract IDs, auto-name, optionally save, draw patch + highlight."""
        # --- extract --------------------------------------------------------
        ids = self._extract_ids(vertices)

        # --- auto-name (no input() blocking) --------------------------------
        self._roi_counter += 1
        roi_name = f"{self.roi_prefix}_{self._roi_counter}"

        # --- store ----------------------------------------------------------
        self.rois[roi_name] = ids
        self.polygons[roi_name] = vertices

        if self.inplace:
            self._save_to_adata()

        # --- draw boundary patch --------------------------------------------
        patch = MplPolygon(
            vertices,
            closed=True,
            fill=(self._FACE_COLOR != "none"),
            facecolor=self._FACE_COLOR,
            edgecolor=self._EDGE_COLOR,
            linewidth=self._LINEWIDTH,
            alpha=0.9,
            zorder=1000,
        )
        self.ax.add_patch(patch)
        self.patches.append(patch)

        # --- highlight selected cells ---------------------------------------
        color_idx = (self._roi_counter - 1) % len(_HIGHLIGHT_COLORS)
        self._highlight_ids(ids, color=_HIGHLIGHT_COLORS[color_idx])

        self._update_title()
        self.fig.canvas.draw_idle()

        print(f"[{shape}] {roi_name}: selected {len(ids)} observations")

    # ------------------------------------------------------------------
    #  visual feedback: highlight selected cells
    # ------------------------------------------------------------------

    def _highlight_ids(self, ids: List, color: str = "#e6194b") -> None:
        """Add a semi-transparent scatter highlight over selected cells/bins."""
        # Use obsm coordinates: self.basis for squarebin, "spatial" for cellbin
        coord_key = self.basis if self.mode == "squarebin" else "spatial"
        if coord_key not in self.adata.obsm:
            return  # no coordinates available for highlighting
        coords = self.adata.obsm[coord_key]
        idxs = [self.adata.obs_names.get_loc(cid) for cid in ids]
        pts = coords[idxs]

        if len(pts) == 0:
            return

        sc = self.ax.scatter(
            pts[:, 0],
            pts[:, 1],
            s=12,
            color=color,
            alpha=0.5,
            edgecolors="none",
            zorder=999,
            rasterized=True,
        )
        self._highlight_collections.append(sc)

    def _clear_highlights(self) -> None:
        """Remove all highlight scatter collections."""
        for sc in self._highlight_collections:
            try:
                sc.remove()
            except Exception:
                pass
        self._highlight_collections.clear()
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    #  naming (no input() — always auto-name)
    # ------------------------------------------------------------------

    def _prompt_roi_name(self) -> str:
        """Auto-generate ROI name — never calls ``input()`` in this version."""
        self._roi_counter += 1
        return f"{self.roi_prefix}_{self._roi_counter}"

    def _save_to_adata(self) -> None:
        """Write current ROI labels to ``adata.obs`` (inplace mode)."""
        key = self.key_added
        if key not in self.adata.obs.columns:
            self.adata.obs[key] = None
            self.adata.obs[key] = self.adata.obs[key].astype(object)
        else:
            self.adata.obs[key] = self.adata.obs[key].astype(object)
            self.adata.obs[key] = None
        for roi_name, ids in self.rois.items():
            self.adata.obs.loc[ids, key] = roi_name

    # ------------------------------------------------------------------
    #  geometry intersection
    # ------------------------------------------------------------------

    def _extract_ids(self, vertices: np.ndarray) -> list:
        """Return observation IDs whose geometries intersect the ROI polygon."""
        roi_poly = Polygon(vertices)
        minx, miny, maxx, maxy = roi_poly.bounds

        if self.mode == "squarebin":
            coords = self.adata.obsm[self.basis]
            mask = (
                (coords[:, 0] >= minx)
                & (coords[:, 0] <= maxx)
                & (coords[:, 1] >= miny)
                & (coords[:, 1] <= maxy)
            )
            candidates = self.valid_ids[mask]
            selected = []
            for cid in candidates:
                idx = self.adata.obs_names.get_loc(cid)
                pt = coords[idx]
                if roi_poly.contains(Point(pt[0], pt[1])):
                    selected.append(cid)
            return selected

        # cellbin
        selected = []
        for cid in self.valid_ids:
            geom = _get_geometry(self.adata, cid, self.geometries, self.use_wkt)
            if geom is None:
                continue
            try:
                if geom.intersects(roi_poly):
                    selected.append(cid)
            except Exception:
                continue
        return selected

    # ------------------------------------------------------------------
    #  public API
    # ------------------------------------------------------------------

    def add_roi(self, vertices: np.ndarray, name: Optional[str] = None) -> str:
        """Manually add a ROI (programmatic use).

        Parameters
        ----------
        vertices : np.ndarray
            N×2 array of polygon vertices.
        name : str or None
            ROI name.  Auto-generated (``<prefix>_N``) when ``None``.

        Returns
        -------
        str
            The assigned ROI name.
        """
        ids = self._extract_ids(vertices)
        if name is None:
            self._roi_counter += 1
            name = f"{self.roi_prefix}_{self._roi_counter}"
        self.rois[name] = ids
        self.polygons[name] = np.asarray(vertices, dtype=float)

        if self.inplace:
            self._save_to_adata()

        patch = MplPolygon(
            vertices,
            closed=True,
            fill=(self._FACE_COLOR != "none"),
            facecolor=self._FACE_COLOR,
            edgecolor=self._EDGE_COLOR,
            linewidth=self._LINEWIDTH,
            alpha=0.9,
            zorder=1000,
        )
        self.ax.add_patch(patch)
        self.patches.append(patch)

        # highlight
        color_idx = len(self.rois) % len(_HIGHLIGHT_COLORS)
        self._highlight_ids(ids, color=_HIGHLIGHT_COLORS[color_idx])

        self._update_title()
        self.fig.canvas.draw_idle()
        print(f"[manual] {name}: selected {len(ids)} observations")
        return name

    def rename_roi(self, old_name: str, new_name: str) -> None:
        """Rename a previously created ROI.

        Parameters
        ----------
        old_name : str
            Current ROI name.
        new_name : str
            New ROI name (must not already exist).
        """
        if old_name not in self.rois:
            raise KeyError(
                f"ROI '{old_name}' not found. "
                f"Existing ROIs: {list(self.rois.keys())}"
            )
        if new_name in self.rois:
            raise KeyError(f"ROI '{new_name}' already exists.")
        self.rois[new_name] = self.rois.pop(old_name)
        self.polygons[new_name] = self.polygons.pop(old_name)
        if self.inplace:
            self._save_to_adata()
        self._update_title()
        print(f"Renamed: '{old_name}' → '{new_name}'")

    def undo(self) -> None:
        """Remove the last added ROI (boundary, highlight, and data)."""
        if not self.rois:
            print("Nothing to undo.")
            return
        last_name = list(self.rois.keys())[-1]
        del self.rois[last_name]
        del self.polygons[last_name]
        # remove boundary patch
        if self.patches:
            patch = self.patches.pop()
            try:
                patch.remove()
            except Exception:
                pass
        # remove last highlight
        if self._highlight_collections:
            sc = self._highlight_collections.pop()
            try:
                sc.remove()
            except Exception:
                pass
        if self.inplace:
            self._save_to_adata()
        self._update_title()
        self.fig.canvas.draw_idle()
        print(f"Undone: '{last_name}'")

    def save(self, key_added: Optional[str] = None) -> None:
        """Force-write current ROI labels to ``adata.obs`` (useful with ``inplace=False``)."""
        key = key_added or self.key_added
        self._save_to_adata()
        print(f"Stored {len(self.rois)} ROI(s) in adata.obs['{key}']")

    def to_adata(self):
        """
        Return a **copy** of the AnnData object with the ROI column added.

        This is the primary way to retrieve results when ``inplace=False``.
        """
        import copy as _copy

        ad = _copy.deepcopy(self.adata)
        key = self.key_added
        ad.obs[key] = None
        ad.obs[key] = ad.obs[key].astype(object)
        for roi_name, ids in self.rois.items():
            ad.obs.loc[ids, key] = roi_name
        return ad

    def clear(self) -> None:
        """Remove all ROIs, visual patches, and highlights."""
        self.rois.clear()
        self.polygons.clear()
        for patch in self.patches:
            try:
                patch.remove()
            except Exception:
                pass
        self.patches.clear()
        self._clear_highlights()
        self._roi_counter = 0
        if self.inplace and self.key_added in self.adata.obs.columns:
            self.adata.obs[self.key_added] = None
        self._update_title()
        self.fig.canvas.draw_idle()
        print("Cleared all ROIs.")

    def disconnect(self) -> None:
        """Disconnect all selectors, toolbar buttons, and keyboard callback."""
        for sel in self._selectors.values():
            sel.disconnect_events()
        if hasattr(self, "_key_cid"):
            self.fig.canvas.mpl_disconnect(self._key_cid)
        print("Disconnected ROI selector.")


# =============================================================================
#  Public entry point
# =============================================================================


def select_regions(
    adata,
    color: Optional[str] = None,
    library_id: Optional[str] = None,
    img_key: Optional[str] = None,
    mode: ModeType = "auto",
    basis: str = "spatial",
    key_added: str = "ROI",
    inplace: bool = True,
    roi_prefix: str = "ROI",
    ax: Optional[plt.Axes] = None,
    show: bool = True,
    figsize: Optional[tuple] = None,
    invert_y: bool = True,
    **kwargs,
) -> RegionSelector:
    """
    Select spatial regions interactively in Jupyter using matplotlib widgets.

    .. important::
       For the best notebook experience run ``%matplotlib widget`` before
       calling this function.

    **No more blocking ``input()``** — ROIs are auto-named with a configurable
    prefix (default ``ROI_1``, ``ROI_2``, …).  Use
    ``selector.rename_roi("ROI_1", "my_region")`` to rename afterwards.

    An **inline toolbar** at the bottom of the figure provides clickable buttons:

    * ``■ Rect`` — rectangle selector
    * ``● Ellipse`` — ellipse / circle selector
    * ``✎ Lasso`` — lasso / freehand selector
    * ``✕ Clear`` — remove all ROIs
    * ``↩ Undo`` — remove the last ROI

    Keyboard shortcuts also work while the figure has focus:

    * ``r`` — rectangle
    * ``e`` — ellipse / circle
    * ``l`` — lasso / freehand

    Selected cells are highlighted with coloured scatter points in real-time.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    color : str or None
        Column or gene used to colour the background plot.
    library_id : str or None
        Spatial library id (auto-detected when omitted).
    img_key : str or None
        Key for the background image in ``uns['spatial']``.
    mode : {'auto', 'cellbin', 'squarebin'}
        Spatial data mode.  Auto-detected when ``'auto'``.
    basis : str
        Key in ``adata.obsm`` for squarebin coordinates (default ``'spatial'``).
    key_added : str
        Column name in ``adata.obs`` where ROI labels are stored.
    inplace : bool
        If ``True`` (default) ROI labels are written to ``adata.obs[key_added]``
        after every selection.  If ``False``, results are only stored on the
        returned ``RegionSelector`` — use ``selector.to_adata()`` to get a copy.
    roi_prefix : str
        Prefix for auto-generated ROI names (default ``'ROI'`` → ``ROI_1``, ``ROI_2``, …).
    ax : matplotlib.axes.Axes or None
        Pre-existing axes to draw on.
    show : bool
        Whether to call ``plt.show()`` at the end.
    figsize : tuple or None
        Figure size passed to ``plt.subplots()`` when ``ax`` is not provided.
    invert_y : bool
        Match image coordinates (y increases top-to-bottom).  Default ``True``.
    **kwargs
        Forwarded to the underlying plotting function (``spatial_cell`` or
        ``spatial_squarebin``).

    Returns
    -------
    RegionSelector
        Controller with ``.rois``, ``.polygons``, ``.rename_roi()``,
        ``.undo()``, ``.save()``, ``.clear()``, ``.to_adata()``, and
        ``.disconnect()``.
    """

    library_id = _resolve_library_id(adata, library_id)
    spatial_info = adata.uns["spatial"][library_id]
    mode = _detect_mode(adata, library_id, mode, basis)
    valid_ids = _valid_ids(adata, spatial_info, mode, basis)
    if len(valid_ids) == 0:
        raise ValueError("No valid observations found for ROI selection.")

    geometries = spatial_info.get("geometries", None)
    use_wkt = (mode == "cellbin") and geometries is None and "geometry" in adata.obs.columns

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (8, 8))

    # Draw background
    if mode == "squarebin":
        spatial_squarebin(
            adata,
            color=color,
            library_id=library_id,
            img_key=img_key,
            basis=basis,
            ax=ax,
            show=False,
            invert_y=invert_y,
            **kwargs,
        )
    else:
        spatial_cell(
            adata,
            color=color,
            library_id=library_id,
            img_key=img_key,
            ax=ax,
            show=False,
            invert_y=invert_y,
            **kwargs,
        )

    selector = RegionSelector(
        adata=adata,
        ax=ax,
        mode=mode,
        basis=basis,
        valid_ids=valid_ids,
        geometries=geometries,
        use_wkt=use_wkt,
        key_added=key_added,
        inplace=inplace,
        roi_prefix=roi_prefix,
    )

    print(
        "Interactive ROI selector ready.\n"
        "  Click toolbar buttons at the bottom:  ■ Rect  ● Ellipse  ✎ Lasso  ✕ Clear  ↩ Undo\n"
        "  Or use keyboard:  r=rect  e=ellipse  l=lasso  (figure must have focus)\n"
        "  ROIs are auto-named — use selector.rename_roi(old, new) to rename.\n"
        "  Run '%matplotlib widget' in Jupyter for best results."
    )
    if show:
        plt.show()
    return selector
