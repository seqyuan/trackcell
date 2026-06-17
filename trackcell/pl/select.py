"""
Jupyter-native interactive ROI selection for TrackCell.

This module intentionally does **not** use napari/Qt.  It relies on matplotlib's
widget selectors, which are safer in Jupyter notebooks when used with an
interactive backend such as ipympl (``%matplotlib widget``).

Key features:

* Keyboard-toggle between rectangle (``r``), ellipse (``e``), and lasso/freehand
  (``l``) modes — all available without a ``shape_type`` parameter.
* Interactive ROI naming via ``input()`` after each selection.
* ``inplace=True`` (default) writes ROI labels directly to ``adata.obs[key_added]``;
  ``inplace=False`` stores results only on the ``RegionSelector`` object.
"""

from __future__ import annotations

from typing import Optional, Dict, List, Literal, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.widgets import RectangleSelector, EllipseSelector, LassoSelector
from shapely.geometry import Polygon, Point
from shapely import wkt

from .plot import spatial_cell, spatial_squarebin, _resolve_library_id


ModeType = Literal["auto", "cellbin", "squarebin"]
SelectorMode = Literal["rectangle", "ellipse", "lasso"]


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


class RegionSelector:
    """
    Controller returned by :func:`select_regions`.

    All shape modes (rectangle, ellipse, lasso) are active simultaneously and
    toggled via keyboard shortcuts while the figure has focus.

    Keyboard shortcuts
    ------------------
    * ``r`` — rectangle selector
    * ``e`` — ellipse selector
    * ``l`` — lasso / freehand selector

    Attributes
    ----------
    rois : dict
        Mapping of ROI names to lists of observation IDs.
    polygons : dict
        Mapping of ROI names to arrays of ROI boundary vertices in data coords.
    ax, fig
        Matplotlib axes and figure.
    """

    # --- visual defaults (no longer exposed to the user) -----
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

        # state
        self._current_mode: SelectorMode = "rectangle"
        self.rois: Dict[str, List] = {}
        self.polygons: Dict[str, np.ndarray] = {}
        self.patches: List[Any] = []
        self._roi_counter = 0

        # selectors (all created, only one active at a time)
        self._selectors: Dict[str, Any] = {}
        self._connect_all_selectors()

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

    def _update_title(self) -> None:
        suffix = f" | inplace={self.inplace}" if not self.inplace else ""
        self.ax.set_title(
            f"Mode: [{self._current_mode}]  (r=rect  e=ellipse  l=lasso)"
            f"  —  {len(self.rois)} ROI(s){suffix}"
        )

    def _on_key(self, event) -> None:
        if event.key in ("r", "e", "l"):
            mode_map = {"r": "rectangle", "e": "ellipse", "l": "lasso"}
            self._set_active(mode_map[event.key])

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
        """Shared pipeline: extract IDs, name the ROI, optionally save, draw patch."""

        # --- extract --------------------------------------------------------
        ids = self._extract_ids(vertices)

        # --- name -----------------------------------------------------------
        roi_name = self._prompt_roi_name()

        # --- store ----------------------------------------------------------
        self.rois[roi_name] = ids
        self.polygons[roi_name] = vertices

        if self.inplace:
            self._save_to_adata()

        # --- draw -----------------------------------------------------------
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
        self._update_title()
        self.fig.canvas.draw_idle()

        print(f"[{shape}] {roi_name}: selected {len(ids)} observations")

    def _prompt_roi_name(self) -> str:
        """Ask the user for a ROI name via ``input()``; auto-generate on empty."""
        self._roi_counter += 1
        default = f"ROI_{self._roi_counter}"
        try:
            name = input(f"ROI name (press Enter for '{default}'): ").strip()
        except (EOFError, KeyboardInterrupt):
            # fallback for non-interactive contexts
            name = ""
        if not name:
            name = default
        return name

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
        """Manually add a ROI (programmatic use)."""
        ids = self._extract_ids(vertices)
        if name is None:
            self._roi_counter += 1
            name = f"ROI_{self._roi_counter}"
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
        self._update_title()
        self.fig.canvas.draw_idle()
        print(f"[manual] {name}: selected {len(ids)} observations")
        return name

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
        """Remove all ROIs and visual patches."""
        self.rois.clear()
        self.polygons.clear()
        for patch in self.patches:
            patch.remove()
        self.patches.clear()
        self._roi_counter = 0
        if self.inplace and self.key_added in self.adata.obs.columns:
            self.adata.obs[self.key_added] = None
        self._update_title()
        self.fig.canvas.draw_idle()

    def disconnect(self) -> None:
        """Disconnect all selectors and keyboard callback."""
        for sel in self._selectors.values():
            sel.disconnect_events()
        if hasattr(self, "_key_cid"):
            self.fig.canvas.mpl_disconnect(self._key_cid)


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

    All shape modes are available simultaneously.  While the figure has focus
    press:

    * ``r`` — rectangle
    * ``e`` — ellipse / circle
    * ``l`` — lasso / freehand

    After each drawn ROI you will be prompted for a name (press Enter to accept
    an auto-generated name such as ``ROI_1``).

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
        Controller with ``.rois``, ``.polygons``, ``.to_adata()``, ``.save()``,
        ``.clear()``, and ``.disconnect()``.
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
    )

    print(
        "Interactive ROI selector ready.  "
        "Press r=rect  e=ellipse  l=lasso while the figure has focus.\n"
        "Use `%matplotlib widget` in Jupyter for best results."
    )
    if show:
        plt.show()
    return selector
