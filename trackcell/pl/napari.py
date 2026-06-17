"""
Napari-based interactive spatial visualization and ROI selection for TrackCell.

Supports both **cellbin** (polygon geometries) and **squarebin** (regular grid
coordinates) data.  Auto-detects the mode, or accepts an explicit ``mode``
parameter.

Provides functions to open napari viewers with H&E images and annotations,
and to interactively draw regions of interest (ROIs) for cell/bin extraction.

Requires ``napari`` to be installed: ``pip install 'trackcell[napari]'``
"""

from __future__ import annotations

import warnings
from typing import Optional, Union, Dict, List, Tuple, Any, Literal

import numpy as np
import pandas as pd
from shapely.geometry import Polygon, Point

from .plot import _resolve_library_id


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _detect_mode(
    adata,
    library_id: str,
    mode: Optional[str],
) -> str:
    """
    Determine whether the adata is ``'cellbin'`` or ``'squarebin'``.

    If *mode* is ``'auto'`` or ``None``, auto-detect:
      1. Geometries (GeoDataFrame or WKT) -> ``'cellbin'``
      2. ``adata.obsm['spatial']`` present  -> ``'squarebin'``
      3. Otherwise raise.
    """
    if mode is not None and mode not in ("auto", "cellbin", "squarebin"):
        raise ValueError(
            f"`mode` must be 'auto', 'cellbin', or 'squarebin', got '{mode}'."
        )

    if mode in ("cellbin", "squarebin"):
        return mode

    # auto-detect
    spatial_info = adata.uns["spatial"].get(library_id, {})
    has_geometries = spatial_info.get("geometries") is not None
    has_wkt = "geometry" in adata.obs.columns

    if has_geometries or has_wkt:
        return "cellbin"

    if "spatial" in adata.obsm:
        return "squarebin"

    raise ValueError(
        "Cannot auto-detect data mode.  Expected either:\n"
        "  - cellbin:  adata.uns['spatial'][lib]['geometries'] (GeoDataFrame) "
        "or adata.obs['geometry'] (WKT), or\n"
        "  - squarebin: adata.obsm['spatial'] with coordinate array.\n"
        "Set `mode='cellbin'` or `mode='squarebin'` explicitly."
    )


def _get_geometry(adata, cell_id: str, geometries, use_wkt: bool = False):
    """Safely retrieve a shapely geometry for a given cell_id."""
    try:
        if use_wkt:
            from shapely import wkt
            raw = adata.obs.loc[cell_id, "geometry"]
            if pd.isna(raw):
                return None
            return wkt.loads(raw) if isinstance(raw, str) else raw
        else:
            if cell_id not in geometries.index:
                return None
            return geometries.loc[cell_id, "geometry"]
    except Exception:
        return None


def _he_image_to_napari(
    viewer,
    spatial_info: dict,
    img_key: Optional[str] = None,
) -> Optional[Any]:
    """
    Add the H&E background image as a napari Image layer.

    Returns the layer or None if no image is available.
    """
    if img_key is None:
        img_key = "hires" if "hires" in spatial_info.get("images", {}) else None

    if not img_key or "images" not in spatial_info or img_key not in spatial_info["images"]:
        return None

    img = spatial_info["images"][img_key]
    scalefactors = spatial_info.get("scalefactors", {})

    scale_key = f"tissue_{img_key}_scalef"
    scale_factor = scalefactors.get(scale_key, 1.0)

    # In napari the leading axis is y (rows), second axis is x (cols).
    # scale=(1/sf, 1/sf) maps image pixels -> data coordinates.
    layer = viewer.add_image(
        img,
        name=f"H&E ({img_key})",
        scale=(1.0 / scale_factor, 1.0 / scale_factor),
        blending="translucent",
        opacity=0.6,
        visible=True,
    )
    return layer


def _get_centroids_for_mode(
    adata,
    mode: str,
    geometries,
    valid_ids: pd.Index,
    use_wkt: bool,
    basis: str,
) -> Optional[Tuple[np.ndarray, pd.Index]]:
    """
    Return (coords, valid_index) where *coords* is an (N,2) array of
    [y, x] centroids (napari convention) and *valid_index* is the
    subset of *valid_ids* that produced a valid centroid.

    *mode* is ``'cellbin'`` or ``'squarebin'``.
    """
    centroids_x = []
    centroids_y = []
    keep = []

    if mode == "cellbin":
        for cid in valid_ids:
            geom = _get_geometry(adata, cid, geometries, use_wkt)
            if geom is None or not hasattr(geom, "centroid"):
                continue
            c = geom.centroid
            if c is None:
                continue
            centroids_x.append(c.x)
            centroids_y.append(c.y)
            keep.append(cid)
    else:
        # squarebin: coordinates come from adata.obsm[basis]
        if basis not in adata.obsm:
            raise ValueError(f"`adata.obsm['{basis}']` not found for squarebin mode.")
        coords_all = np.asarray(adata.obsm[basis])
        if coords_all.ndim != 2 or coords_all.shape[1] < 2:
            raise ValueError(
                f"`adata.obsm['{basis}']` must be n_obs x 2 for squarebin mode."
            )
        for cid in valid_ids:
            loc = adata.obs_names.get_loc(cid)
            centroids_x.append(coords_all[loc, 0])
            centroids_y.append(coords_all[loc, 1])
            keep.append(cid)

    if not keep:
        return None

    # napari Points data: (N, 2) -> [[y, x], ...]
    coords = np.column_stack([np.array(centroids_y), np.array(centroids_x)])
    return coords, pd.Index(keep)


def _cell_centroids_as_points(
    viewer,
    adata,
    mode: str,
    geometries,
    valid_ids: pd.Index,
    color: Optional[str],
    palette: Optional[Union[dict, list, np.ndarray]],
    cmap: str,
    point_size: float,
    use_wkt: bool,
    basis: str,
) -> Optional[Any]:
    """
    Add centroids as a napari Points layer, optionally coloured by *color*.
    Works for both ``'cellbin'`` and ``'squarebin'`` modes.
    """
    result = _get_centroids_for_mode(adata, mode, geometries, valid_ids, use_wkt, basis)
    if result is None:
        return None

    coords, valid_idx = result

    # --- colour handling ---
    face_color = "white"
    properties = None
    metadata = {"valid_cells": valid_idx}

    if color is not None:
        if color in adata.var_names:
            expr = adata[valid_idx, color].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = np.asarray(expr).flatten()

            from matplotlib.cm import get_cmap as _get_cmap
            from matplotlib.colors import Normalize

            cmap_obj = _get_cmap(cmap)
            norm = Normalize(vmin=np.nanmin(expr), vmax=np.nanmax(expr))
            rgba = cmap_obj(norm(expr))
            face_color = rgba

        elif color in adata.obs.columns:
            vals = adata.obs.loc[valid_idx, color]
            if pd.api.types.is_numeric_dtype(vals):
                arr = vals.values.astype(float)
                from matplotlib.cm import get_cmap as _get_cmap
                from matplotlib.colors import Normalize

                cmap_obj = _get_cmap(cmap)
                norm = Normalize(vmin=np.nanmin(arr), vmax=np.nanmax(arr))
                rgba = cmap_obj(norm(arr))
                face_color = rgba
            else:
                cats = vals.astype(str)
                unique_cats = sorted(cats.dropna().unique())
                cat_colors = _build_categorical_colors(unique_cats, palette)
                rgba = np.array([cat_colors.get(c, "gray") for c in cats])
                from matplotlib.colors import to_rgba
                rgba = np.array([to_rgba(c) for c in rgba])
                face_color = rgba
        else:
            warnings.warn(
                f"`color` key '{color}' not found in adata.var_names or adata.obs.columns."
            )

    label = "Bins" if mode == "squarebin" else "Cells"
    layer = viewer.add_points(
        coords,
        name=f"{label} ({color})" if color else label,
        size=point_size,
        face_color=face_color,
        edge_color="black",
        edge_width=0.0,
        opacity=0.8,
        properties=properties,
        metadata=metadata,
    )
    return layer


def _build_categorical_colors(
    categories: List[str],
    palette: Optional[Union[dict, list, np.ndarray]] = None,
) -> Dict[str, str]:
    """Build a mapping of category -> color string."""
    if palette is None:
        from matplotlib import colormaps
        cmap = colormaps.get("tab10") if len(categories) <= 10 else colormaps.get("tab20")
        from matplotlib.colors import to_hex
        return {cat: to_hex(cmap(i / max(len(categories) - 1, 1)))
                for i, cat in enumerate(categories)}

    if isinstance(palette, dict):
        result = {}
        for i, cat in enumerate(categories):
            result[cat] = palette.get(cat, "gray")
        return result

    arr = np.asarray(palette)
    return {cat: arr[i % len(arr)] for i, cat in enumerate(categories)}


def _extract_cells_in_shapes(
    adata,
    mode: str,
    shapes_data: List[np.ndarray],
    geometries,
    valid_ids: pd.Index,
    use_wkt: bool,
    basis: str,
    shape_types: Optional[List[str]] = None,
) -> Dict[str, List]:
    """
    For each user-drawn shape, find all cells/bins whose centroid falls
    inside the shape polygon.

    - *mode='cellbin'*: uses ``geom.intersects(roi_poly)``.
    - *mode='squarebin'*: uses ``Point(x,y).within(roi_poly)``.

    Returns dict: ``{"ROI_1": [cell_id, ...], "ROI_2": [...]}``
    """
    results = {}

    if mode == "squarebin":
        coords_all = np.asarray(adata.obsm[basis])

    for i, vertices in enumerate(shapes_data):
        if len(vertices) < 3:
            continue

        try:
            roi_poly = Polygon(vertices)
            if not roi_poly.is_valid:
                roi_poly = roi_poly.buffer(0)
            if roi_poly.is_empty:
                continue
        except Exception:
            continue

        cells_in_roi = []

        if mode == "cellbin":
            for cid in valid_ids:
                geom = _get_geometry(adata, cid, geometries, use_wkt)
                if geom is None:
                    continue
                try:
                    if geom.intersects(roi_poly):
                        cells_in_roi.append(cid)
                except Exception:
                    continue
        else:
            for cid in valid_ids:
                loc = adata.obs_names.get_loc(cid)
                x = coords_all[loc, 0]
                y = coords_all[loc, 1]
                if Point(x, y).within(roi_poly):
                    cells_in_roi.append(cid)

        key = f"ROI_{i + 1}"
        results[key] = cells_in_roi

    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def napari_view(
    adata,
    color: Optional[str] = None,
    library_id: Optional[str] = None,
    img_key: Optional[str] = None,
    mode: Literal["auto", "cellbin", "squarebin"] = "auto",
    basis: str = "spatial",
    palette: Optional[Union[dict, list, np.ndarray]] = None,
    cmap: str = "viridis",
    point_size: float = 5.0,
) -> "napari.Viewer":
    """
    Open a napari viewer with spatial transcriptomics data.

    Supports both **cellbin** (polygon geometries) and **squarebin** (regular
    grid coordinates).  The mode is auto-detected by default.

    Layers added:
    - H&E background image (if available)
    - Cell / bin centroids as a Points layer (optionally coloured by *color*)

    Use this function for interactive exploration.  After closing napari
    you can call :func:`napari_extract` to extract cells/bins within user-drawn
    ROIs.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial information.
    color : str, optional
        Column in ``adata.obs`` or gene in ``adata.var_names`` to colour
        centroids by.  If ``None``, shown in white.
    library_id : str, optional
        Key in ``adata.uns["spatial"]``.  Auto-detected if ``None``.
    img_key : str, optional
        Image key in ``adata.uns["spatial"][library_id]["images"]``.
        Defaults to ``"hires"``.
    mode : str, default "auto"
        Data mode: ``"auto"``, ``"cellbin"``, or ``"squarebin"``.
    basis : str, default "spatial"
        Key in ``adata.obsm`` for squarebin coordinates.  Ignored for cellbin.
    palette : dict, list, or array, optional
        Colour palette for categorical *color*.
    cmap : str, default "viridis"
        Matplotlib colormap for continuous *color*.
    point_size : float, default 5.0
        Size of centroid dots in napari.

    Returns
    -------
    napari.Viewer
        The napari Viewer instance with layers already added.
    """
    try:
        import napari
    except ImportError:
        raise ImportError(
            "napari is required for interactive ROI selection. "
            "Install with:  pip install 'trackcell[napari]'"
        )

    library_id = _resolve_library_id(adata, library_id)
    spatial_info = adata.uns["spatial"][library_id]
    mode = _detect_mode(adata, library_id, mode)

    geometries = spatial_info.get("geometries", None)
    use_wkt = (mode == "cellbin") and geometries is None and "geometry" in adata.obs.columns

    if mode == "cellbin" and geometries is None and not use_wkt:
        raise ValueError(
            "No cell geometries found. Expected:\n"
            "  - adata.uns['spatial'][library_id]['geometries'] (GeoDataFrame), or\n"
            "  - adata.obs['geometry'] (WKT strings)."
        )

    if mode == "cellbin":
        if use_wkt:
            valid_ids = adata.obs_names[adata.obs["geometry"].notna()]
        else:
            valid_ids = geometries.index.intersection(adata.obs_names)
    else:
        if basis not in adata.obsm:
            raise ValueError(
                f"`adata.obsm['{basis}']` not found for squarebin mode."
            )
        valid_ids = adata.obs_names

    if len(valid_ids) == 0:
        raise ValueError("No valid observations found.")

    label = "Bin" if mode == "squarebin" else "Cell"
    viewer = napari.Viewer(title=f"TrackCell - {library_id} ({label}s)")

    _he_image_to_napari(viewer, spatial_info, img_key)
    _cell_centroids_as_points(
        viewer, adata, mode, geometries, valid_ids,
        color, palette, cmap, point_size, use_wkt, basis,
    )

    return viewer


def select_regions(
    adata,
    color: Optional[str] = None,
    library_id: Optional[str] = None,
    img_key: Optional[str] = None,
    mode: Literal["auto", "cellbin", "squarebin"] = "auto",
    basis: str = "spatial",
    shape_type: str = "rectangle",
    palette: Optional[Union[dict, list, np.ndarray]] = None,
    cmap: str = "viridis",
    point_size: float = 5.0,
    key_added: str = "ROI",
    copy: bool = False,
) -> Optional[Dict[str, List]]:
    """
    Interactive ROI selection using napari.

    Opens a napari viewer with the H&E image and cell/bin centroids.  The user
    draws rectangles / polygons on the **"ROI regions"** Shapes layer, then
    presses **Enter** to finish.  Cells/bins whose centroids fall inside each
    ROI are collected.

    Supports both **cellbin** and **squarebin** data.  Mode is auto-detected
    by default.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial information.
    color : str, optional
        Column in ``adata.obs`` or gene name in ``adata.var_names`` to colour
        centroids by.
    library_id : str, optional
        Key in ``adata.uns["spatial"]``.  Auto-detected if ``None``.
    img_key : str, optional
        Image key.  Defaults to ``"hires"``.
    mode : str, default "auto"
        Data mode: ``"auto"``, ``"cellbin"``, or ``"squarebin"``.
    basis : str, default "spatial"
        Key in ``adata.obsm`` for squarebin coordinates.  Ignored for cellbin.
    shape_type : str, default "rectangle"
        Allowed shape types for ROI drawing.  One of:
        ``"rectangle"``, ``"polygon"``, ``"free"`` (freehand path),
        or ``"any"`` (all four types).
    palette : dict, list, or array, optional
        Colour palette for categorical *color*.
    cmap : str, default "viridis"
        Colormap for continuous *color*.
    point_size : float, default 5.0
        Size of centroid dots in napari.
    key_added : str, default "ROI"
        Key in ``adata.obs`` where ROI annotations are stored (when
        ``copy=False``).  Each cell/bin is assigned the ROI label
        (``"ROI_1"``, ``"ROI_2"``, ...), or ``None`` if not in any ROI.
    copy : bool, default False
        If ``False``, ROI annotations are written into ``adata.obs[key_added]``
        and the function returns ``None``.  If ``True``, ``adata`` is not
        modified; instead a ``dict`` mapping ROI labels -> ID lists is
        returned.

    Returns
    -------
    dict or None
        If ``copy=True``: ``{"ROI_1": [cell_id, ...], "ROI_2": [...]}``.
        If ``copy=False``: ``None`` (annotations stored in ``adata.obs``).

    Examples
    --------
    Cellbin data:

    >>> import trackcell as tcl
    >>> rois = tcl.pl.select_regions(adata, color="CellType", copy=True)
    >>> print(rois.keys())

    Squarebin data:

    >>> rois = tcl.pl.select_regions(adata, mode="squarebin", copy=True)
    >>> # or let it auto-detect:
    >>> tcl.pl.select_regions(adata, key_added="ROI_bins")
    """
    try:
        import napari
    except ImportError:
        raise ImportError(
            "napari is required for interactive ROI selection. "
            "Install with:  pip install 'trackcell[napari]'"
        )

    library_id = _resolve_library_id(adata, library_id)
    spatial_info = adata.uns["spatial"][library_id]
    mode = _detect_mode(adata, library_id, mode)

    _shape_map = {
        "rectangle": "rectangle",
        "polygon": "polygon",
        "free": "path",
        "any": ["rectangle", "polygon", "path", "ellipse"],
    }
    napari_shape = _shape_map.get(shape_type)
    if napari_shape is None:
        raise ValueError(
            f"Unknown shape_type '{shape_type}'. "
            f"Choose from: {list(_shape_map.keys())}"
        )

    geometries = spatial_info.get("geometries", None)
    use_wkt = (mode == "cellbin") and geometries is None and "geometry" in adata.obs.columns

    if mode == "cellbin" and geometries is None and not use_wkt:
        raise ValueError(
            "No cell geometries found. Expected:\n"
            "  - adata.uns['spatial'][library_id]['geometries'] (GeoDataFrame), or\n"
            "  - adata.obs['geometry'] (WKT strings)."
        )

    if mode == "cellbin":
        if use_wkt:
            valid_ids = adata.obs_names[adata.obs["geometry"].notna()]
        else:
            valid_ids = geometries.index.intersection(adata.obs_names)
    else:
        if basis not in adata.obsm:
            raise ValueError(f"`adata.obsm['{basis}']` not found for squarebin mode.")
        valid_ids = adata.obs_names

    if len(valid_ids) == 0:
        raise ValueError("No valid observations found.")

    # --- Build viewer ---
    label = "Bin" if mode == "squarebin" else "Cell"
    viewer = napari.Viewer(
        title=f"TrackCell ROI - {library_id} ({label}s)  [Press Enter when done]"
    )

    _he_image_to_napari(viewer, spatial_info, img_key)
    _cell_centroids_as_points(
        viewer, adata, mode, geometries, valid_ids,
        color, palette, cmap, point_size, use_wkt, basis,
    )

    roi_layer = viewer.add_shapes(
        [],
        shape_type=napari_shape,
        edge_width=3,
        edge_color="red",
        face_color="royalblue",
        opacity=0.3,
        name="ROI regions",
    )

    _result: Dict[str, Any] = {"rois": None}

    @viewer.bind_key("Enter")
    def _finish(v: napari.Viewer):
        shapes_data = list(roi_layer.data)
        shape_types = list(roi_layer.shape_type)
        _result["rois"] = _extract_cells_in_shapes(
            adata, mode, shapes_data, geometries, valid_ids, use_wkt, basis,
            shape_types=shape_types,
        )
        v.close()

    @viewer.window._qt_window.destroyed.connect
    def _on_close():
        if _result["rois"] is None:
            shapes_data = list(roi_layer.data)
            shape_types = list(roi_layer.shape_type)
            _result["rois"] = _extract_cells_in_shapes(
                adata, mode, shapes_data, geometries, valid_ids, use_wkt, basis,
                shape_types=shape_types,
            )

    print(
        f"\n{'='*60}\n"
        f"  Mode: {mode}  |  Shape type: {shape_type}\n"
        f"  Draw ROI shapes on the 'ROI regions' layer.\n"
        f"  Press ENTER to finish and extract {label.lower()}s.\n"
        f"{'='*60}\n"
    )

    napari.run()

    rois = _result.get("rois") or {}

    if copy:
        return rois

    if key_added not in adata.obs.columns:
        adata.obs[key_added] = None
        adata.obs[key_added] = adata.obs[key_added].astype(object)

    for roi_name, cell_ids in rois.items():
        adata.obs.loc[cell_ids, key_added] = roi_name

    print(f"Stored {len(rois)} ROI(s) in adata.obs['{key_added}']")
    return None


def napari_extract(
    adata,
    viewer: "napari.Viewer",
    mode: Literal["auto", "cellbin", "squarebin"] = "auto",
    basis: str = "spatial",
    shapes_layer: str = "ROI regions",
    key_added: str = "ROI",
    copy: bool = False,
) -> Optional[Dict[str, List]]:
    """
    Extract cells/bins within user-drawn shapes from a napari viewer.

    Use this after :func:`napari_view` when you have added a Shapes layer
    for ROIs manually.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    viewer : napari.Viewer
        The napari Viewer containing a Shapes layer with ROIs.
    mode : str, default "auto"
        Data mode: ``"auto"``, ``"cellbin"``, or ``"squarebin"``.
    basis : str, default "spatial"
        Key in ``adata.obsm`` for squarebin coordinates.  Ignored for cellbin.
    shapes_layer : str, default "ROI regions"
        Name of the Shapes layer containing the ROIs.
    key_added : str, default "ROI"
        Key in ``adata.obs`` where annotations are stored (``copy=False``).
    copy : bool, default False
        If ``True``, return a dict without modifying ``adata``.

    Returns
    -------
    dict or None
    """
    from .plot import _resolve_library_id

    library_id = _resolve_library_id(adata, None)
    spatial_info = adata.uns["spatial"][library_id]
    mode = _detect_mode(adata, library_id, mode)

    geometries = spatial_info.get("geometries", None)
    use_wkt = (mode == "cellbin") and geometries is None and "geometry" in adata.obs.columns

    if mode == "cellbin":
        if use_wkt:
            valid_ids = adata.obs_names[adata.obs["geometry"].notna()]
        else:
            valid_ids = geometries.index.intersection(adata.obs_names)
    else:
        if basis not in adata.obsm:
            raise ValueError(f"`adata.obsm['{basis}']` not found for squarebin mode.")
        valid_ids = adata.obs_names

    if shapes_layer not in viewer.layers:
        available = [l.name for l in viewer.layers]
        raise ValueError(
            f"Shapes layer '{shapes_layer}' not found. "
            f"Available layers: {available}"
        )

    sl = viewer.layers[shapes_layer]
    shapes_data = list(sl.data)
    shape_types = list(sl.shape_type) if hasattr(sl, "shape_type") else None

    rois = _extract_cells_in_shapes(
        adata, mode, shapes_data, geometries, valid_ids, use_wkt, basis,
        shape_types=shape_types,
    )

    if copy:
        return rois

    if key_added not in adata.obs.columns:
        adata.obs[key_added] = None
        adata.obs[key_added] = adata.obs[key_added].astype(object)

    for roi_name, cell_ids in rois.items():
        adata.obs.loc[cell_ids, key_added] = roi_name

    print(f"Stored {len(rois)} ROI(s) in adata.obs['{key_added}']")
    return None
