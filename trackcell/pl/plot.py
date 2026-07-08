"""
Plotting functions for TrackCell package.

This module provides functions for visualizing spatial transcriptomics data,
including cell polygon visualization.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap, to_rgba
from matplotlib.patches import Rectangle, Circle
from matplotlib.collections import PatchCollection
from matplotlib.cm import ScalarMappable
from typing import Optional, Union, List, Dict
import warnings

try:
    import geopandas as gpd
    from shapely import wkt
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    warnings.warn("geopandas and shapely are required for spatial_cell function")


def _process_background_image(spatial_info, img_key, data_coords_range=None):
    """
    Process background image similar to scanpy's approach.
    
    Parameters
    ----------
    spatial_info : dict
        Spatial information dictionary from adata.uns['spatial'][library_id]
    img_key : str or None
        Key for the image to use
    data_coords_range : tuple, optional
        Tuple of (x_min, y_min, x_max, y_max) for data coordinate range
        
    Returns
    -------
    img : numpy.ndarray or None
        Background image array
    img_extent : list or None
        Image extent [left, right, bottom, top] in data coordinates
    """
    if img_key is None:
        img_key = "hires" if "hires" in spatial_info.get("images", {}) else None
    
    if not img_key or "images" not in spatial_info or img_key not in spatial_info["images"]:
        return None, None
    
    img = spatial_info["images"][img_key]
    scalefactors = spatial_info.get("scalefactors", {})
    
    # Get scale factor for the image (scanpy way)
    # For hires: tissue_hires_scalef (~0.05-0.2)
    # For lowres: tissue_lowres_scalef (~0.03)
    # For fullres: scale_factor = 1.0
    scale_key = f"tissue_{img_key}_scalef"
    scale_factor = scalefactors.get(scale_key, 1.0)
    
    # Calculate image extent in data coordinates
    # Image shape is (height, width) in pixels
    img_height, img_width = img.shape[:2]
    
    if scale_factor < 1.0:
        # Image is downscaled, need to scale up the extent
        # The image pixels represent a smaller region in full-res coordinates
        img_extent = [0, img_width / scale_factor, img_height / scale_factor, 0]
    else:
        # Full resolution image
        img_extent = [0, img_width, img_height, 0]
    
    return img, img_extent


def _draw_categorical_edges(
    ax: plt.Axes,
    temp_gdf: gpd.GeoDataFrame,
    adata,
    valid_cells: list,
    edge_color_col: str,
    edge_palette: Optional[Union[dict, list, np.ndarray]] = None,
    linewidth: float = 0.5,
    alpha: float = 0.8,
    legend: bool = True,
):
    """Draw categorical cell boundaries on top of filled polygons."""
    from matplotlib.patches import Patch
    
    if len(valid_cells) == 0:
        return
    
    # Get edge category data
    edge_vals = adata.obs.loc[valid_cells, edge_color_col]
    temp_gdf['_edge_cat'] = edge_vals.values
    
    # Determine categories and colors
    if pd.api.types.is_categorical_dtype(edge_vals):
        edge_categories = edge_vals.cat.categories.tolist()
    else:
        edge_categories = sorted(edge_vals.dropna().unique())
    
    n_cats = len(edge_categories)
    if n_cats == 0:
        return
    
    # Build edge color map
    if edge_palette is not None:
        if isinstance(edge_palette, dict):
            edge_colors = [edge_palette.get(cat, 'gray') for cat in edge_categories]
        elif isinstance(edge_palette, (list, np.ndarray)):
            pa = np.asarray(edge_palette)
            edge_colors = [pa[i % len(pa)] for i in range(n_cats)]
        else:
            raise ValueError(f"Unsupported edge_palette type: {type(edge_palette)}")
    else:
        default_cmap = plt.get_cmap('tab10' if n_cats <= 10 else 'tab20')
        edge_colors = [default_cmap(i / n_cats) for i in range(n_cats)]
    
    # Draw edges per category
    for cat, cat_color in zip(edge_categories, edge_colors):
        mask = adata.obs.loc[valid_cells, edge_color_col] == cat
        cat_cells = [cid for cid, m in zip(valid_cells, mask) if m]
        if not cat_cells:
            continue
        cat_gdf = temp_gdf.loc[temp_gdf.index.isin(cat_cells)]
        if len(cat_gdf) == 0:
            continue
        try:
            cat_gdf.boundary.plot(
                ax=ax, color=cat_color, linewidth=linewidth, alpha=alpha
            )
        except Exception:
            pass
    
    # Edge legend
    if legend and n_cats <= 30:
        legend_elements = [
            Patch(facecolor='none', edgecolor=ec, linewidth=linewidth * 2,
                  label=str(cat))
            for cat, ec in zip(edge_categories, edge_colors)
        ]
        # Place edge legend at bottom of plot
        ax.legend(
            handles=legend_elements,
            title=edge_color_col,
            bbox_to_anchor=(0.5, -0.12),
            loc='upper center',
            ncol=min(n_cats, 8),
            frameon=True,
            fontsize='small',
        )
    
    # Clean up temp column
    if '_edge_cat' in temp_gdf.columns:
        del temp_gdf['_edge_cat']


def _is_hex_color_series(series: pd.Series) -> bool:
    """Detect if a Series contains raw hex color strings (e.g., from multigene_blend)."""
    non_null = series.dropna()
    if len(non_null) == 0:
        return False
    sample = non_null.iloc[:100] if len(non_null) > 100 else non_null
    try:
        return all(
            isinstance(v, str) and v.startswith('#') and len(v) == 7
            for v in sample
        )
    except Exception:
        return False


def _plot_raw_hex_colors(
    ax: plt.Axes,
    temp_gdf: gpd.GeoDataFrame,
    plot_column: str,
    edgecolor: str = 'none',
    linewidth: float = 0,
    alpha: float = 1.0,
    use_equal_aspect: bool = False,
):
    """Plot polygons with per-cell raw hex colors (from multigene_blend)."""
    grouped = temp_gdf.groupby(plot_column, sort=False)
    n_groups = len(grouped)
    
    for hex_color, group in grouped:
        if not isinstance(hex_color, str) or not hex_color.startswith('#'):
            continue
        plot_kw = {
            'ax': ax,
            'color': hex_color,
            'edgecolor': edgecolor,
            'linewidth': linewidth,
            'alpha': alpha,
        }
        if use_equal_aspect:
            plot_kw['aspect'] = 'equal'
        try:
            group.plot(**plot_kw)
        except ValueError as e:
            if 'aspect must be finite and positive' in str(e) and not use_equal_aspect:
                plot_kw['aspect'] = 'equal'
                try:
                    group.plot(**plot_kw)
                except Exception:
                    pass
            else:
                pass




def _resolve_library_id(adata, library_id: Optional[str] = None) -> str:
    """Resolve spatial library id similar to scanpy/scalpy style helpers."""
    if "spatial" not in adata.uns:
        raise ValueError("`adata.uns['spatial']` is required but missing.")

    if library_id is None:
        available_library_ids = list(adata.uns["spatial"].keys())
        if len(available_library_ids) == 0:
            raise ValueError("No library_id found in `adata.uns['spatial']`.")
        library_id = available_library_ids[0]
        if len(available_library_ids) > 1:
            warnings.warn(
                f"Multiple library_ids found: {available_library_ids}. "
                f"Using '{library_id}'. Specify `library_id` explicitly to use a different one."
            )

    if library_id not in adata.uns["spatial"]:
        raise ValueError(
            f"`library_id` '{library_id}' not found in `adata.uns['spatial']`. "
            f"Available library_ids: {list(adata.uns['spatial'].keys())}"
        )
    return library_id


def _has_valid_geometry_bounds(geom) -> bool:
    if geom is None:
        return False
    try:
        if pd.isna(geom):
            return False
    except (TypeError, ValueError):
        pass
    if not hasattr(geom, 'bounds'):
        return False
    if hasattr(geom, 'is_empty') and geom.is_empty:
        return False
    if hasattr(geom, 'is_valid') and not geom.is_valid:
        return False
    try:
        bounds = geom.bounds
        return all(np.isfinite(bounds)) and bounds[2] > bounds[0] and bounds[3] > bounds[1]
    except Exception:
        return False


def _bounds_are_finite(bounds) -> bool:
    try:
        return all(np.isfinite(bounds)) and bounds[2] > bounds[0] and bounds[3] > bounds[1]
    except Exception:
        return False


def _filter_gdf_to_valid_bounds(temp_gdf: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, bool]:
    """Return a GeoDataFrame with invalid geometry bounds removed."""
    if len(temp_gdf) == 0:
        return temp_gdf, False
    try:
        if _bounds_are_finite(temp_gdf.total_bounds):
            return temp_gdf, False
    except Exception:
        pass

    valid_mask = temp_gdf.geometry.apply(_has_valid_geometry_bounds)
    return temp_gdf.loc[valid_mask].copy(), bool((~valid_mask).any())


def _sync_geometries_to_obs(adata, library_id: str, warn: bool = True):
    """Filter/reorder stored geometries to match current adata.obs_names."""
    spatial_info = adata.uns.get("spatial", {}).get(library_id, {})
    geometries = spatial_info.get("geometries")
    if geometries is None or not HAS_GEOPANDAS:
        return geometries
    if not isinstance(geometries, (gpd.GeoDataFrame, gpd.GeoSeries)):
        return geometries

    obs_index = pd.Index(adata.obs_names)
    geom_index = pd.Index(geometries.index)
    if geom_index.equals(obs_index):
        return geometries

    geom_ids = set(geom_index)
    common_ids = [cid for cid in obs_index if cid in geom_ids]
    if len(common_ids) == 0:
        if warn:
            warnings.warn(
                "No observation IDs match the geometry index. Plotting may fall back "
                "to adata.obs['geometry'] if available."
            )
        return geometries

    missing_ids = obs_index.difference(geom_index)
    if len(missing_ids) > 0 and warn:
        warnings.warn(
            f"{len(missing_ids)} observations have no stored geometry and will be skipped. "
            f"Examples: {missing_ids[:5].tolist()}"
        )

    synced = geometries.loc[common_ids].copy()
    spatial_info["geometries"] = synced
    return synced


def _normalize_color_argument(color):
    if color is None:
        return [None]
    if isinstance(color, str):
        return [color]
    return list(color)


def _validate_color_keys(adata, colors_to_plot):
    for color_key in colors_to_plot:
        if color_key is None:
            continue
        if color_key not in adata.obs.columns and color_key not in adata.var_names:
            if hasattr(adata, 'layers') and color_key in adata.layers:
                raise ValueError(
                    f"`color` key '{color_key}' found in `adata.layers`, but layer-based plotting "
                    f"is not yet supported. Please use gene names from `adata.var_names` "
                    f"or metadata from `adata.obs.columns`."
                )
            raise ValueError(
                f"`color` key '{color_key}' not found in `adata.obs.columns` or `adata.var_names`. "
                f"Available obs keys: {list(adata.obs.columns[:10])}... "
                f"Available var names (genes): {list(adata.var_names[:10])}..."
            )


def _create_axes(colors_to_plot, figsize=None, ax=None):
    if ax is None:
        if figsize is None:
            figsize = (5 * len(colors_to_plot), 5) if len(colors_to_plot) > 1 else (10, 10)
        if len(colors_to_plot) > 1:
            fig, axes = plt.subplots(1, len(colors_to_plot), figsize=figsize, sharex=True, sharey=True)
            if len(colors_to_plot) == 1:
                axes = [axes]
            else:
                axes = list(np.ravel(axes))
        else:
            fig, single_ax = plt.subplots(1, 1, figsize=figsize)
            axes = [single_ax]
    else:
        fig = ax.figure
        axes = [ax]
        if len(colors_to_plot) > 1:
            warnings.warn("Multiple colors specified but single ax provided. Only first color will be plotted.")
            colors_to_plot = [colors_to_plot[0]]
    return fig, axes, colors_to_plot


def _filter_obs_mask(adata, color_key=None, groups=None, groupby=None):
    if groups is None:
        return np.ones(adata.n_obs, dtype=bool)

    filter_column = None
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise ValueError(f"`groupby` column '{groupby}' not found in `adata.obs.columns`.")
        filter_column = groupby
    elif color_key is not None and color_key in adata.obs.columns:
        filter_column = color_key
    else:
        raise ValueError(
            "`groups` requires either `groupby` to be specified or `color` to be a column in `adata.obs`."
        )

    return adata.obs[filter_column].isin(groups).to_numpy()


def _extract_color_values(adata, color_key, mask):
    if color_key is None:
        return None
    if color_key in adata.obs.columns:
        return adata.obs.loc[mask, color_key]
    gene_idx = adata.var_names.get_loc(color_key)
    values = adata.X[mask, gene_idx]
    if hasattr(values, 'toarray'):
        values = values.toarray().ravel()
    else:
        values = np.asarray(values).ravel()
    return pd.Series(values, index=adata.obs_names[mask], name=color_key)


def _compute_data_extent_from_coords(coords, pad_fraction=0.05, extra_pad=0.0):
    if coords is None or len(coords) == 0:
        x_min = y_min = 0.0
        x_max = y_max = 1.0
    else:
        x_min = float(np.nanmin(coords[:, 0]))
        y_min = float(np.nanmin(coords[:, 1]))
        x_max = float(np.nanmax(coords[:, 0]))
        y_max = float(np.nanmax(coords[:, 1]))

    x_range = x_max - x_min
    y_range = y_max - y_min
    x_padding = (x_range * pad_fraction if x_range > 0 else 1.0) + extra_pad
    y_padding = (y_range * pad_fraction if y_range > 0 else 1.0) + extra_pad
    return (x_min, y_min, x_max, y_max), (x_padding, y_padding)


def _apply_spatial_axis_formatting(ax, x_min, y_min, x_max, y_max, x_padding, y_padding, xlabel, ylabel, show_ticks, force_show_ticks=False, invert_y=True):
    ax.set_aspect('equal')
    if invert_y:
        ax.invert_yaxis()
        ax.set_ylim(y_max + y_padding, y_min - y_padding)
    else:
        ax.set_ylim(y_min - y_padding, y_max + y_padding)
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if force_show_ticks:
        ax.tick_params(axis='both', which='major', labelsize=10)
    elif not show_ticks:
        ax.set_xticks([])
        ax.set_yticks([])


def _infer_square_size(adata, library_id, coords, binsize=None):
    spatial_info = adata.uns.get('spatial', {}).get(library_id, {})
    if binsize is None:
        binsize = spatial_info.get('binsize')

    scalefactors = spatial_info.get('scalefactors', {})
    microns_per_pixel = scalefactors.get('microns_per_pixel')
    if binsize is not None and microns_per_pixel not in (None, 0):
        try:
            size = float(binsize) / float(microns_per_pixel)
            if np.isfinite(size) and size > 0:
                return size
        except Exception:
            pass

    if coords is not None and len(coords) > 1:
        diffs = []
        for axis in (0, 1):
            uniq = np.unique(np.asarray(coords[:, axis], dtype=float))
            if len(uniq) > 1:
                d = np.diff(np.sort(uniq))
                d = d[np.isfinite(d) & (d > 0)]
                if len(d) > 0:
                    diffs.append(np.median(d))
        if diffs:
            size = float(np.min(diffs))
            if np.isfinite(size) and size > 0:
                return size

    return 1.0


def _draw_background_only(ax, spatial_info, img_key, x_min, y_min, x_max, y_max, xlabel, ylabel, invert_y=True):
    img, img_extent = _process_background_image(spatial_info, img_key)
    if img is not None and img_extent is not None:
        ax.imshow(img, extent=img_extent, origin='upper', alpha=1.0)
        ax.set_xlim(img_extent[0], img_extent[1])
        if invert_y:
            ax.set_ylim(img_extent[2], img_extent[3])
        else:
            ax.set_ylim(img_extent[3], img_extent[2])
    else:
        ax.set_xlim(x_min, x_max)
        if invert_y:
            ax.set_ylim(y_max, y_min)
        else:
            ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')
    if invert_y:
        ax.invert_yaxis()
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.tick_params(axis='both', which='major', labelsize=10)


def _plot_squarebin_values(ax, coords, values, square_size, cmap, palette, vmin, vmax, alpha, legend, edges_width, edges_color, na_color, rasterized=False, shape='circle'):
    if len(coords) == 0:
        return None

    half = square_size / 2.0
    if shape == 'square':
        patches = [Rectangle((x - half, y - half), square_size, square_size) for x, y in coords]
    else:
        patches = [Circle((x, y), radius=half) for x, y in coords]

    edgecolor = 'none' if edges_width == 0 else edges_color

    if values is None:
        pc = PatchCollection(patches, facecolor='none', edgecolor=edgecolor, linewidth=edges_width, alpha=alpha, rasterized=rasterized)
        ax.add_collection(pc)
        return None

    if isinstance(values, pd.Series):
        values_arr = values.to_numpy()
    else:
        values_arr = np.asarray(values)

    is_categorical = (
        pd.api.types.is_categorical_dtype(values)
        or pd.api.types.is_object_dtype(values)
        or pd.api.types.is_bool_dtype(values)
    )

    if is_categorical:
        values_ser = pd.Series(values_arr, index=np.arange(len(values_arr)))
        categories = pd.Index(pd.unique(values_ser.dropna()))
        if len(categories) == 0:
            facecolors = [to_rgba(na_color, alpha=alpha)] * len(values_ser)
            pc = PatchCollection(patches, facecolor=facecolors, edgecolor=edgecolor, linewidth=edges_width, rasterized=rasterized)
            ax.add_collection(pc)
            return None

        if palette is not None:
            if isinstance(palette, dict):
                cat_to_color = {cat: palette.get(cat, na_color) for cat in categories}
            else:
                pal = list(palette)
                cat_to_color = {cat: pal[i % len(pal)] for i, cat in enumerate(categories)}
        else:
            default_cmap = plt.get_cmap('tab20' if len(categories) <= 20 else 'tab20b')
            cat_to_color = {cat: default_cmap(i / max(len(categories), 1)) for i, cat in enumerate(categories)}

        facecolors = []
        for v in values_ser:
            if pd.isna(v):
                facecolors.append(to_rgba(na_color, alpha=alpha))
            else:
                facecolors.append(to_rgba(cat_to_color.get(v, na_color), alpha=alpha))

        pc = PatchCollection(patches, facecolor=facecolors, edgecolor=edgecolor, linewidth=edges_width, rasterized=rasterized)
        ax.add_collection(pc)

        if legend:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=cat_to_color[c], edgecolor='none', label=str(c)) for c in categories]
            if handles:
                ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
        return None

    values_num = np.asarray(values_arr, dtype=float)
    finite_mask = np.isfinite(values_num)
    if not finite_mask.any():
        facecolors = [to_rgba(na_color, alpha=alpha)] * len(values_num)
        pc = PatchCollection(patches, facecolor=facecolors, edgecolor=edgecolor, linewidth=edges_width, rasterized=rasterized)
        ax.add_collection(pc)
        return None
    data_vmin = float(np.nanmin(values_num[finite_mask])) if vmin is None else vmin
    data_vmax = float(np.nanmax(values_num[finite_mask])) if vmax is None else vmax
    if data_vmax == data_vmin:
        data_vmax = data_vmin + 1e-12
    norm = Normalize(vmin=data_vmin, vmax=data_vmax)
    cmap_obj = plt.get_cmap(cmap)
    facecolors = []
    for v in values_num:
        if np.isnan(v):
            facecolors.append(to_rgba(na_color, alpha=alpha))
        else:
            facecolors.append(to_rgba(cmap_obj(norm(v)), alpha=alpha))
    pc = PatchCollection(patches, facecolor=facecolors, edgecolor=edgecolor, linewidth=edges_width, rasterized=rasterized)
    ax.add_collection(pc)
    if legend:
        sm = ScalarMappable(norm=norm, cmap=cmap_obj)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    return None

def spatial_cell(
    adata,
    color: Optional[Union[str, List[str]]] = None,
    groups: Optional[List[str]] = None,
    groupby: Optional[str] = None,
    library_id: Optional[str] = None,
    size: float = 1.0,
    figsize: Optional[tuple] = None,
    cmap: str = "viridis",
    palette: Optional[Union[dict, list, np.ndarray]] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    img_key: Optional[str] = None,
    basis: str = "spatial",
    edges_width: float = 0.5,
    edges_color: str = "black",
    edge_color: Optional[str] = None,
    edge_palette: Optional[Union[dict, list, np.ndarray]] = None,
    alpha: float = 0.8,
    alpha_img: float = 0.5,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    legend: bool = True,
    xlabel: Optional[str] = "spatial 1",
    ylabel: Optional[str] = "spatial 2",
    show_ticks: bool = False,
    invert_y: bool = True,
    **kwargs
):
    """
    Plot spatial transcriptomics data with cell polygons instead of points.
    
    This function visualizes cells as polygons (from cell segmentation) rather than
    simple points, providing a more accurate representation of cell boundaries.
    Uses GeoDataFrame.plot() for efficient rendering and automatic legend generation.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial information and cell geometries.
    color : str or list of str, optional
        Keys for observation/categorical or continuous variables to color cells.
        Can be a single key or a list of keys for multiple plots.
        Can be:
        - A column name in `adata.obs` (metadata)
        - A gene name in `adata.var_names` (gene expression)
        - None: Only display the H&E background image without cell polygons.
          When None, axis ticks and labels are automatically shown.
    groups : list of str, optional
        Subset of groups to plot. If None, plots all groups.
        Requires either `color` to be a categorical column in `adata.obs` or `groupby` to be specified.
    groupby : str, optional
        Column name in `adata.obs` to use for filtering with `groups` parameter.
        If None and `groups` is specified, will use `color` if it's a categorical column in `adata.obs`.
        This is useful when `color` is a continuous variable (e.g., gene expression) but you want
        to filter by a categorical column (e.g., 'classification').
    library_id : str, optional
        Key in adata.uns["spatial"] containing spatial information.
        If None, uses the first available library_id (similar to sc.pl.spatial).
        Default is None, which will auto-select the first library_id.
    size : float, default 1.0
        Size scaling factor for cells (not used for polygons, kept for API compatibility).
        This parameter is currently not implemented for polygon-based visualization.
    figsize : tuple, optional
        Figure size (width, height) in inches.
    cmap : str, default "viridis"
        Colormap for continuous values.
    vmin : float, optional
        Minimum value for colormap normalization. If None, uses the minimum value in the data.
    vmax : float, optional
        Maximum value for colormap normalization. If None, uses the maximum value in the data.
    palette : dict, list, or array, optional
        Color palette for categorical variables. Can be:
        - A dictionary mapping category names to colors (e.g., {'A': 'red', 'B': 'blue'})
        - A list/array of colors that will be assigned to categories in sorted order
          (e.g., ['red', 'blue', 'green'] will assign colors to categories alphabetically)
    img_key : str, optional
        Key in adata.uns["spatial"][library_id]["images"] for background image.
        If None, uses "hires" if available.
    basis : str, default "spatial"
        Key in adata.obsm containing spatial coordinates (for fallback).
    edges_width : float, default 0.5
        Width of cell polygon edges.
    edges_color : str, default "black"
        Color of cell polygon edges (when `edge_color` is not specified).
    edge_color : str, optional
        Column name in ``adata.obs`` for categorical coloring of cell edges.
        When specified, cell edges are colored by category (e.g., cell type),
        while ``color`` controls the fill (typically gene expression).
        Overrides ``edges_color`` when set.
    edge_palette : dict, list, or array, optional
        Color palette for ``edge_color`` categories. Same format as ``palette``.
    alpha : float, default 0.8
        Transparency of cell polygons.
    alpha_img : float, default 0.5
        Transparency of background image.
    show : bool, default True
        Whether to display the plot.
    ax : matplotlib.Axes, optional
        Axes object to plot on. If None, creates a new figure.
    legend : bool, default True
        Whether to show legend for categorical values or colorbar for continuous values.
    xlabel : str, optional, default "spatial 1"
        Label for the x-axis. Set to None to hide the label.
    ylabel : str, optional, default "spatial 2"
        Label for the y-axis. Set to None to hide the label.
    show_ticks : bool, default False
        Whether to show axis ticks and tick labels. If False, ticks are hidden.
        Note: When `color=None`, ticks are automatically shown regardless of this setting.
    invert_y : bool, default True
        Whether to invert the y-axis so spatial coordinates increase from top to bottom
        (image convention, matching the H&E background). Set to ``False`` to use
        Cartesian convention (y increases from bottom to top).
    **kwargs
        Additional arguments passed to GeoDataFrame.plot().
    
    Returns
    -------
    matplotlib.Axes or list of matplotlib.Axes
        Axes object(s) containing the plot.
    
    Examples
    --------
    >>> import trackcell as tcl
    >>> adata = tcl.io.read_hd_cellseg("path/to/data", sample="sample1")
    >>> # Plot by metadata (categorical)
    >>> tcl.pl.spatial_cell(adata, color="classification")
    >>> # Plot by metadata (continuous)
    >>> tcl.pl.spatial_cell(adata, color="Cluster-2_dist", cmap="Reds")
    >>> # Plot by gene expression
    >>> tcl.pl.spatial_cell(adata, color="PDPN", cmap="viridis")
    >>> # Plot with groups filter
    >>> tcl.pl.spatial_cell(adata, color="classification", groups=["Cluster-1", "Cluster-2"])
    >>> # Plot only H&E image (no cell polygons)
    >>> tcl.pl.spatial_cell(adata, color=None)
    """
    if not HAS_GEOPANDAS:
        raise ImportError("geopandas and shapely are required for spatial_cell function")
    
    # Check if geometries are available
    if "spatial" not in adata.uns:
        raise ValueError("`adata.uns['spatial']` is required but missing.")
    
    # Auto-select library_id if not provided (similar to sc.pl.spatial)
    if library_id is None:
        available_library_ids = list(adata.uns["spatial"].keys())
        if len(available_library_ids) == 0:
            raise ValueError("No library_id found in `adata.uns['spatial']`.")
        library_id = available_library_ids[0]
        if len(available_library_ids) > 1:
            warnings.warn(
                f"Multiple library_ids found: {available_library_ids}. "
                f"Using '{library_id}'. Specify `library_id` explicitly to use a different one."
            )
    
    if library_id not in adata.uns["spatial"]:
        raise ValueError(
            f"`library_id` '{library_id}' not found in `adata.uns['spatial']`. "
            f"Available library_ids: {list(adata.uns['spatial'].keys())}"
        )
    
    spatial_info = adata.uns["spatial"][library_id]
    
    # Get geometries - try GeoDataFrame first, then fallback to WKT strings in obs
    geometries = None
    use_wkt = False
    
    if "geometries" in spatial_info:
        geometries = _sync_geometries_to_obs(adata, library_id)
    elif "geometry" in adata.obs.columns:
        # Fallback: use WKT strings from obs
        use_wkt = True
        warnings.warn(
            "Using geometry from adata.obs['geometry'] (WKT strings). "
            "For better performance, use adata.uns['spatial']['geometries'] (GeoDataFrame)."
        )
    else:
        raise ValueError(
            f"Cell geometries not found. Expected either:\n"
            f"  - adata.uns['spatial']['{library_id}']['geometries'] (GeoDataFrame), or\n"
            f"  - adata.obs['geometry'] (WKT strings).\n"
            "Please use `tcl.io.read_hd_cellseg()` to load data with geometries."
        )
    
    # Handle color parameter (can be single string or list)
    if color is None:
        colors_to_plot = [None]
    elif isinstance(color, str):
        colors_to_plot = [color]
    else:
        colors_to_plot = color
    
    # Validate all color keys upfront (before creating figures)
    # This avoids repeated checks later and provides early error feedback
    for color_key in colors_to_plot:
        if color_key is not None:
            # Check if color_key exists in obs.columns or var_names
            if color_key not in adata.obs.columns and color_key not in adata.var_names:
                # Check if it's in layers (not supported yet)
                if hasattr(adata, 'layers') and color_key in adata.layers:
                    raise ValueError(
                        f"`color` key '{color_key}' found in `adata.layers`, but gene expression "
                        f"from layers is not yet supported. Please use gene names from `adata.var_names` "
                        f"or metadata from `adata.obs.columns`."
                    )
                else:
                    raise ValueError(
                        f"`color` key '{color_key}' not found in `adata.obs.columns` or `adata.var_names`. "
                        f"Available obs keys: {list(adata.obs.columns[:10])}... "
                        f"Available var names (genes): {list(adata.var_names[:10])}..."
                    )
    
    # Create figure if needed
    if ax is None:
        if figsize is None:
            if len(colors_to_plot) > 1:
                figsize = (5 * len(colors_to_plot), 5)
            else:
                figsize = (10, 10)
        
        if len(colors_to_plot) > 1:
            fig, axes = plt.subplots(1, len(colors_to_plot), figsize=figsize, sharex=True, sharey=True)
            if len(colors_to_plot) == 1:
                axes = [axes]
        else:
            fig, axes = plt.subplots(1, 1, figsize=figsize)
            axes = [axes]
    else:
        fig = ax.figure
        axes = [ax]
        if len(colors_to_plot) > 1:
            warnings.warn("Multiple colors specified but single ax provided. Only first color will be plotted.")
            colors_to_plot = [colors_to_plot[0]]
    
    axes_list = []
    
    for idx, color_key in enumerate(colors_to_plot):
        if idx < len(axes):
            current_ax = axes[idx]
        else:
            # Should not happen, but handle gracefully
            continue
        
        # Filter cells if groups is specified
        if groups is not None:
            # Determine which column to use for filtering
            filter_column = None
            if groupby is not None:
                # Use explicitly specified groupby column
                if groupby not in adata.obs.columns:
                    raise ValueError(
                        f"`groupby` column '{groupby}' not found in `adata.obs.columns`. "
                        f"Available columns: {list(adata.obs.columns[:10])}..."
                    )
                filter_column = groupby
            elif color_key is not None and color_key in adata.obs.columns:
                # Use color_key if it's a categorical column in obs
                filter_column = color_key
            else:
                # Cannot determine filter column
                raise ValueError(
                    "`groups` parameter requires either:\n"
                    "  - `color` to be a categorical column in `adata.obs`, or\n"
                    "  - `groupby` to specify the column name for filtering.\n"
                    f"Current `color` value: {color_key}"
                )
            
            # Apply groups filter
            mask = adata.obs[filter_column].isin(groups)
        else:
            mask = pd.Series(True, index=adata.obs_names)
        
        # Get cell indices to plot
        cells_to_plot = adata.obs_names[mask]
        
        if len(cells_to_plot) == 0:
            warnings.warn("No cells to plot after filtering.")
            axes_list.append(current_ax)
            continue
        
        # Calculate coordinate range from actual data to be plotted
        # This ensures image extent matches the data range, especially for subset data
        coords_list = []
        for cell_id in cells_to_plot:
            if use_wkt:
                if cell_id in adata.obs.index and pd.notna(adata.obs.loc[cell_id, "geometry"]):
                    try:
                        geom = wkt.loads(adata.obs.loc[cell_id, "geometry"])
                        if hasattr(geom, 'bounds'):
                            coords_list.append(geom.bounds)  # (minx, miny, maxx, maxy)
                    except Exception:
                        continue
            else:
                if cell_id in geometries.index:
                    geom = geometries.loc[cell_id, "geometry"]
                    if geom is not None and hasattr(geom, 'bounds'):
                        coords_list.append(geom.bounds)
        
        if coords_list:
            # Calculate overall bounds from all geometries
            all_bounds = np.array(coords_list)
            x_min = all_bounds[:, 0].min()
            y_min = all_bounds[:, 1].min()
            x_max = all_bounds[:, 2].max()
            y_max = all_bounds[:, 3].max()
        else:
            # Fallback to spatial coordinates if available
            if basis in adata.obsm and len(adata.obsm[basis]) > 0:
                spatial_coords = adata.obsm[basis][mask]
                if len(spatial_coords) > 0:
                    x_min, y_min = spatial_coords.min(axis=0)
                    x_max, y_max = spatial_coords.max(axis=0)
                else:
                    x_min = y_min = 0
                    x_max = y_max = 1
            else:
                x_min = y_min = 0
                x_max = y_max = 1
        
        data_coords_range = (x_min, y_min, x_max, y_max)
        
        # Process and draw background image (scanpy way)
        img, img_extent = _process_background_image(spatial_info, img_key, data_coords_range)
        if img is not None and img_extent is not None:
            if color_key is None:
                if invert_y:
                    current_ax.imshow(img, extent=img_extent, origin='upper', alpha=1.0)
                else:
                    current_ax.imshow(img[::-1], extent=img_extent, origin='upper', alpha=1.0)
            else:
                if invert_y:
                    current_ax.imshow(img, extent=img_extent, origin='upper', alpha=alpha_img)
                else:
                    current_ax.imshow(img[::-1], extent=img_extent, origin='upper', alpha=alpha_img)
        
        # If color is None, skip geometry plotting and only show HE image
        if color_key is None:
            # Set axis limits based on image extent or data coordinates
            if img_extent is not None:
                current_ax.set_xlim(img_extent[0], img_extent[1])
                if invert_y:
                    current_ax.set_ylim(img_extent[2], img_extent[3])
                else:
                    current_ax.set_ylim(img_extent[3], img_extent[2])
            else:
                # Fallback to data coordinates
                current_ax.set_xlim(x_min, x_max)
                if invert_y:
                    current_ax.set_ylim(y_max, y_min)
                else:
                    current_ax.set_ylim(y_min, y_max)
            
            # Set axis properties
            current_ax.set_aspect('equal')
            # invert_yaxis() is intentionally NOT called here:
            # set_ylim(h, 0) with origin='upper' already produces image convention;
            # for Cartesian (invert_y=False) the image was pre-flipped via img[::-1].
            
            # Set axis labels
            if xlabel is not None:
                current_ax.set_xlabel(xlabel)
            if ylabel is not None:
                current_ax.set_ylabel(ylabel)
            
            # Always show ticks when color is None
            current_ax.tick_params(axis='both', which='major', labelsize=10)
            
            axes_list.append(current_ax)
            continue  # Skip to next color or finish
        
        # Create temporary GeoDataFrame for plotting
        # This combines geometry and color data in one structure
        if use_wkt:
            # Convert WKT strings to geometries
            geom_list = []
            valid_cells = []
            for cell_id in cells_to_plot:
                if cell_id in adata.obs.index and pd.notna(adata.obs.loc[cell_id, "geometry"]):
                    try:
                        geom = wkt.loads(adata.obs.loc[cell_id, "geometry"])
                        geom_list.append(geom)
                        valid_cells.append(cell_id)
                    except Exception:
                        continue
            temp_geometries = gpd.GeoSeries(geom_list, index=valid_cells)
        else:
            # Use geometries from GeoDataFrame
            # Filter to only include cells that exist in geometries index
            cells_in_geometries = cells_to_plot[cells_to_plot.isin(geometries.index)]
            if len(cells_in_geometries) == 0:
                warnings.warn(
                    "No cells found in geometries index. "
                    "Geometries may not be synchronized after subset. "
                    "Consider using tcl.io.sync_geometries_after_subset() after subsetting."
                )
                axes_list.append(current_ax)
                continue
            
            # Get geometries for cells that exist in the geometries index
            temp_geometries_raw = geometries.loc[cells_in_geometries, "geometry"]
            
            # Filter out invalid geometries (None, NaN, or invalid geometry objects)
            valid_cells = []
            for cell_id in cells_in_geometries:
                if cell_id not in temp_geometries_raw.index:
                    continue
                    
                geom = temp_geometries_raw.loc[cell_id]
                
                # Check if geometry is valid
                if geom is None or pd.isna(geom):
                    continue
                
                if not hasattr(geom, 'bounds'):
                    continue
                
                # Check if geometry is valid shapely object
                if hasattr(geom, 'is_valid') and not geom.is_valid:
                    continue
                
                # Verify bounds are finite
                try:
                    bounds = geom.bounds
                    if not all(np.isfinite(bounds)):
                        continue
                    # Check bounds are valid (max > min)
                    if bounds[2] <= bounds[0] or bounds[3] <= bounds[1]:
                        continue
                except Exception:
                    continue
                
                valid_cells.append(cell_id)
            
            if len(valid_cells) == 0:
                warnings.warn(
                    "No valid geometries found after filtering invalid geometries. "
                    "This may be due to geometries not being synchronized after subset. "
                    "Consider using tcl.io.sync_geometries_after_subset() after subsetting."
                )
                axes_list.append(current_ax)
                continue
            
            temp_geometries = temp_geometries_raw.loc[valid_cells]
        
        if len(temp_geometries) == 0:
            warnings.warn("No valid geometries found for plotting.")
            axes_list.append(current_ax)
            continue
        
        # Additional validation: check if total_bounds would be valid
        try:
            test_gdf = gpd.GeoDataFrame(geometry=temp_geometries)
            bounds = test_gdf.total_bounds
            if not all(np.isfinite(bounds)) or bounds[2] <= bounds[0] or bounds[3] <= bounds[1]:
                warnings.warn(
                    f"Invalid geometry bounds detected. This may cause plotting errors. "
                    f"Bounds: {bounds}. Consider using sync_geometries_after_subset() after subsetting."
                )
        except Exception as e:
            warnings.warn(
                f"Could not validate geometry bounds: {e}. "
                f"This may cause plotting errors. Consider using sync_geometries_after_subset() after subsetting."
            )
        
        # Create GeoDataFrame with color data
        # color_key has already been validated upfront, so we can safely proceed
        if color_key is None:
            # No coloring, use default
            temp_gdf = gpd.GeoDataFrame(
                geometry=temp_geometries,
                index=valid_cells
            )
            plot_column = None
        elif color_key in adata.obs.columns:
            # Get color data from obs (metadata)
            color_data = adata.obs.loc[valid_cells, color_key]
            temp_gdf = gpd.GeoDataFrame(
                {color_key: color_data},
                geometry=temp_geometries,
                index=valid_cells
            )
            plot_column = color_key
        else:
            # color_key must be in var_names (already validated upfront)
            # Get color data from gene expression (var)
            gene_idx = adata.var_names.get_loc(color_key)
            
            # Get expression values for valid cells
            # Handle both sparse and dense matrices
            if hasattr(adata.X, 'toarray'):
                # Sparse matrix
                expression_values = adata.X[adata.obs_names.get_indexer(valid_cells), gene_idx].toarray().flatten()
            else:
                # Dense matrix
                expression_values = adata.X[adata.obs_names.get_indexer(valid_cells), gene_idx]
            
            # Create Series with valid cells as index
            color_data = pd.Series(expression_values, index=valid_cells, name=color_key)
            
            temp_gdf = gpd.GeoDataFrame(
                {color_key: color_data},
                geometry=temp_geometries,
                index=valid_cells
            )
            plot_column = color_key
        
        temp_gdf, filtered_invalid_bounds = _filter_gdf_to_valid_bounds(temp_gdf)
        if filtered_invalid_bounds:
            warnings.warn(
                "Filtered out geometries with invalid or non-finite bounds before plotting."
            )
        if len(temp_gdf) == 0:
            warnings.warn("No valid geometries remaining after filtering. Skipping plot.")
            axes_list.append(current_ax)
            continue

        valid_cells = temp_gdf.index.tolist()
        try:
            use_equal_aspect = not _bounds_are_finite(temp_gdf.total_bounds)
        except Exception as e:
            warnings.warn(
                f"Could not validate geometry bounds after filtering: {e}. "
                "Using aspect='equal' to avoid GeoPandas aspect errors."
            )
            use_equal_aspect = True
        
        # Determine if continuous or categorical immediately after creating temp_gdf
        # This allows all subsequent logic to use is_categorical directly
        is_categorical = False
        use_custom_palette = False
        custom_cmap = None
        categories = None
        color_list_for_legend = None  # Store color list for legend creation
        
        if color_key is not None:
            # plot_column is guaranteed to be color_key at this point (already validated)
            is_categorical = not pd.api.types.is_numeric_dtype(temp_gdf[plot_column])
            # Detect raw hex color mode (from multigene_blend etc.)
            raw_color = _is_hex_color_series(temp_gdf[plot_column])
        else:
            raw_color = False
        
        # Process categorical data (palette, categories, etc.)
        if color_key is not None and is_categorical:
            # Get categories (prioritize Categorical's original order)

            if pd.api.types.is_categorical_dtype(temp_gdf[plot_column]):
                # Categorical type: use original order
                categories = temp_gdf[plot_column].cat.categories.tolist()
            else:
                # Regular column: use sorted unique values
                categories = sorted(temp_gdf[plot_column].dropna().unique())
            
            color_list = []
            if palette is not None:
                use_custom_palette = True
                # Convert palette to color list (in categories order)
                if isinstance(palette, dict):
                    # Dictionary: map categories to colors
                    # Check for missing categories and warn
                    missing_cats = [cat for cat in categories if cat not in palette]
                    if missing_cats:
                        warnings.warn(
                            f"Palette dictionary is missing colors for {len(missing_cats)} categories: "
                            f"{missing_cats[:5]}{'...' if len(missing_cats) > 5 else ''}. "
                            f"Using 'gray' as default color for missing categories."
                        )
                    color_list = [palette.get(cat, 'gray') for cat in categories]
                elif isinstance(palette, (list, np.ndarray)):
                    # List/array: assign colors sequentially
                    palette_array = np.asarray(palette)
                    if len(palette_array) < len(categories):
                        warnings.warn(
                            f"Palette has {len(palette_array)} colors but there are "
                            f"{len(categories)} categories. Colors will be cycled."
                        )
                    color_list = [
                        palette_array[i % len(palette_array)]
                        for i in range(len(categories))
                    ]
                else:
                    raise ValueError(f"Unsupported palette type: {type(palette)}")
                
                # Create custom colormap from color list
                custom_cmap = ListedColormap(color_list)
                color_list_for_legend = color_list  # Store for legend
                
                # Convert column to Categorical type with specified categories order
                # This ensures GeoPandas uses the correct order for color mapping
                # GeoPandas will automatically use cat.categories for color assignment
                
        # Prepare plot arguments for GeoDataFrame.plot()
        # ── Dual-color mode: edge_color as a column name for categorical edge coloring ──
        dual_color = (edge_color is not None and edge_color in adata.obs.columns)
        if dual_color and color_key is not None:
            edges_draw = edges_color  # kept as fallback name
            edges_color = 'none'      # fill pass draws no edges
        
        plot_kwargs = {
            'ax': current_ax,
            'edgecolor': edges_color,
            'linewidth': edges_width,
            'alpha': alpha,
            **kwargs
        }
        
        if color_key is not None:
            # plot_column is guaranteed to be color_key at this point (already validated)
            plot_kwargs['column'] = plot_column
            
            if raw_color:
                # ── Raw hex color mode (multi-gene blend) ──
                # Group by unique hex colors and plot each group
                _plot_raw_hex_colors(
                    ax=current_ax,
                    temp_gdf=temp_gdf,
                    plot_column=plot_column,
                    edgecolor=plot_kwargs.get('edgecolor', 'none'),
                    linewidth=plot_kwargs.get('linewidth', 0),
                    alpha=plot_kwargs.get('alpha', 1.0),
                    use_equal_aspect=use_equal_aspect,
                )
                
            elif is_categorical:
                # Categorical values
                plot_kwargs['legend'] = False  # Disable automatic legend, we'll add it manually
                
                # Only set categorical=True if column is not already Categorical type
                # GeoPandas can automatically detect Categorical columns, so we don't need
                # to set categorical=True for Categorical columns (as shown in the example)
                #if not pd.api.types.is_categorical_dtype(temp_gdf[plot_column]):
                #    plot_kwargs['categorical'] = True
                
                # If using custom palette, use custom colormap
                # GeoPandas will automatically use the Categorical column's cat.categories
                # for color assignment, so we don't need to set categories parameter
                # Reference: https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.plot.html
                if use_custom_palette:
                    # Use custom colormap
                    # The column is already converted to Categorical with correct order
                    # GeoPandas will use cat.categories automatically
                    plot_kwargs['cmap'] = custom_cmap
                # else: use default GeoPandas categorical plotting with column
                
                # Set aspect='equal' if bounds are invalid to avoid calculation errors
                if use_equal_aspect:
                    plot_kwargs['aspect'] = 'equal'
                
                # Plot using GeoDataFrame.plot() - uses column + categorical + cmap
                # Use try-except to catch aspect calculation errors and retry with aspect='equal'
                try:
                    temp_gdf.plot(**plot_kwargs)
                except ValueError as e:
                    if "aspect must be finite and positive" in str(e):
                        # Retry with explicit aspect='equal'
                        plot_kwargs['aspect'] = 'equal'
                        warnings.warn(
                            f"Aspect calculation failed. Using aspect='equal' instead. "
                            f"Consider using tcl.io.sync_geometries_after_subset() after subsetting."
                        )
                        temp_gdf.plot(**plot_kwargs)
                    else:
                        raise
                
                # Create legend manually from categories
                if legend:
                    from matplotlib.patches import Patch
                    
                    if use_custom_palette:
                        # Use colors from stored color list (same order as categories)
                        legend_elements = [
                            Patch(facecolor=color_list_for_legend[i], label=str(cat))
                            for i, cat in enumerate(categories)
                        ]
                    else:
                        # Generate colors matching GeoPandas default
                        n_cats = len(categories)
                        default_cmap = plt.get_cmap('tab20' if n_cats <= 20 else 'tab20b')
                        legend_elements = [
                            Patch(facecolor=default_cmap(i / n_cats), label=str(cat))
                            for i, cat in enumerate(categories)
                        ]
                    
                    if legend_elements:
                        current_ax.legend(handles=legend_elements, 
                                        bbox_to_anchor=(1.05, 1), 
                                        loc='upper left',
                                        frameon=True)
            else:
                # Continuous values
                plot_kwargs['cmap'] = cmap
                # Use GeoPandas automatic colorbar for continuous values
                plot_kwargs['legend'] = legend
                
                # Handle vmin and vmax for continuous values
                # If not provided, GeoPandas will use data min/max automatically
                if vmin is not None:
                    plot_kwargs['vmin'] = vmin
                if vmax is not None:
                    plot_kwargs['vmax'] = vmax
                
                # Set aspect='equal' if bounds are invalid to avoid calculation errors
                if use_equal_aspect:
                    plot_kwargs['aspect'] = 'equal'
                
                # Plot using GeoDataFrame.plot()
                # GeoPandas will automatically create colorbar if legend=True
                # Use try-except to catch aspect calculation errors and retry with aspect='equal'
                try:
                    temp_gdf.plot(**plot_kwargs)
                except ValueError as e:
                    if "aspect must be finite and positive" in str(e):
                        # Retry with explicit aspect='equal'
                        plot_kwargs['aspect'] = 'equal'
                        warnings.warn(
                            f"Aspect calculation failed. Using aspect='equal' instead. "
                            f"Consider using tcl.io.sync_geometries_after_subset() after subsetting."
                        )
                        temp_gdf.plot(**plot_kwargs)
                    else:
                        raise
        else:
            # No coloring - only show HE image (background image), no cell polygons
            # This is useful for viewing just the tissue image with coordinates
            pass  # Skip geometry plotting, only background image will be shown
        
        # ── Dual-color edge pass: draw categorical cell boundaries on top of fill ──
        if dual_color and color_key is not None:
            _draw_categorical_edges(
                ax=current_ax,
                temp_gdf=temp_gdf,
                adata=adata,
                valid_cells=valid_cells,
                edge_color_col=edge_color,
                edge_palette=edge_palette,
                linewidth=edges_width,
                alpha=alpha,
                legend=legend,
            )
        
        # Set axis properties
        current_ax.set_aspect('equal')
        if invert_y:
            current_ax.invert_yaxis()  # Match image coordinates
        
        # Set axis limits based on actual data range
        # Add small padding (5% of range) for better visualization
        x_range = x_max - x_min
        y_range = y_max - y_min
        x_padding = x_range * 0.05 if x_range > 0 else 1
        y_padding = y_range * 0.05 if y_range > 0 else 1
        
        current_ax.set_xlim(x_min - x_padding, x_max + x_padding)
        if invert_y:
            current_ax.set_ylim(y_max + y_padding, y_min - y_padding)  # Inverted for y-axis
        else:
            current_ax.set_ylim(y_min - y_padding, y_max + y_padding)
        
        # Set axis labels
        if xlabel is not None:
            current_ax.set_xlabel(xlabel)
        if ylabel is not None:
            current_ax.set_ylabel(ylabel)
        
        # Control ticks visibility
        # If color is None, always show ticks to display coordinates
        if color_key is None:
            # When color=None, show ticks and labels by default
            current_ax.tick_params(axis='both', which='major', labelsize=10)
        elif not show_ticks:
            current_ax.set_xticks([])
            current_ax.set_yticks([])
        
        if color_key:
            current_ax.set_title(color_key)
        
        axes_list.append(current_ax)
    
    if show:
        # Only adjust layout for figures created inside this function. User-provided
        # axes may belong to a larger subplot layout that should not be changed here.
        if ax is None:
            fig.tight_layout(rect=[0, 0, 0.95, 1])
        plt.show()
    
    if len(axes_list) == 1:
        return axes_list[0]
    else:
        return axes_list




def spatial_squarebin(
    adata,
    color: Optional[Union[str, List[str]]] = None,
    groups: Optional[List[str]] = None,
    groupby: Optional[str] = None,
    library_id: Optional[str] = None,
    binsize: Optional[float] = None,
    figsize: Optional[tuple] = None,
    cmap: str = "viridis",
    palette: Optional[Union[dict, list, np.ndarray]] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    img_key: Optional[str] = None,
    basis: str = "spatial",
    edges_width: float = 0.0,
    edges_color: str = "none",
    alpha: float = 0.8,
    alpha_img: float = 0.5,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    legend: bool = True,
    xlabel: Optional[str] = "spatial 1",
    ylabel: Optional[str] = "spatial 2",
    show_ticks: bool = False,
    crop_coord: Optional[tuple] = None,
    na_color: str = "#d3d3d3",
    rasterized: bool = False,
    invert_y: bool = True,
    shape: str = 'circle',
    **kwargs
):
    """
    Plot Visium HD square-bin data as circles (default) or squares over an optional H&E image.

    This function is designed for outputs loaded by :func:`trackcell.io.read_hd_bin`.
    It mirrors the user-facing behavior of :func:`spatial_cell` where possible, but
    renders regular square bins instead of cell polygons.

    Parameters
    ----------
    adata
        AnnData object with bin-level spatial coordinates in ``adata.obsm[basis]`` and
        image metadata in ``adata.uns['spatial'][library_id]``.
    color
        Observation column, gene name, or list of either. If ``None``, only the H&E
        image and coordinate range are displayed.
    groups, groupby
        Optional filtering of bins, following the same semantics as ``spatial_cell``.
    library_id
        Spatial library identifier. If ``None``, the first entry in
        ``adata.uns['spatial']`` is used.
    binsize
        Bin size in micrometers. If ``None``, uses
        ``adata.uns['spatial'][library_id]['binsize']`` when available.
    figsize, cmap, palette, vmin, vmax, img_key, basis, alpha, alpha_img, show, ax, legend,
    xlabel, ylabel, show_ticks
        Behave similarly to ``spatial_cell``.
    edges_width, edges_color
        Styling for square boundaries. Defaults favor performance by disabling edges.
    crop_coord
        Optional ``(x_min, x_max, y_min, y_max)`` crop in spatial coordinates.
    na_color
        Color used for missing values.
    rasterized
        Whether to rasterize square patches for smaller vector outputs on large datasets.
    invert_y
        Whether to invert the y-axis so spatial coordinates increase from top to bottom
        (image convention, matching the H&E background). Default ``True``.
        Set to ``False`` to use Cartesian convention (y increases from bottom to top).
    shape
        Shape of each bin marker. ``'circle'`` (default) or ``'square'``.
    **kwargs
        Reserved for future extensions. Currently unused.
    """
    if basis not in adata.obsm:
        raise ValueError(f"`adata.obsm['{basis}']` is required but missing.")

    if shape not in ('circle', 'square'):
        raise ValueError("`shape` must be 'circle' or 'square'.")

    library_id = _resolve_library_id(adata, library_id)
    spatial_info = adata.uns['spatial'][library_id]

    colors_to_plot = _normalize_color_argument(color)
    _validate_color_keys(adata, colors_to_plot)
    fig, axes, colors_to_plot = _create_axes(colors_to_plot, figsize=figsize, ax=ax)

    coords_all = np.asarray(adata.obsm[basis])
    if coords_all.ndim != 2 or coords_all.shape[1] < 2:
        raise ValueError(f"`adata.obsm['{basis}']` must be an n_obs × 2 array of spatial coordinates.")
    coords_all = coords_all[:, :2]

    axes_list = []

    for idx, color_key in enumerate(colors_to_plot):
        current_ax = axes[idx]
        mask = _filter_obs_mask(adata, color_key=color_key, groups=groups, groupby=groupby)

        if crop_coord is not None:
            if len(crop_coord) != 4:
                raise ValueError("`crop_coord` must be a 4-tuple: (x_min, x_max, y_min, y_max).")
            x0, x1, y0, y1 = crop_coord
            crop_mask = (
                (coords_all[:, 0] >= x0) & (coords_all[:, 0] <= x1) &
                (coords_all[:, 1] >= y0) & (coords_all[:, 1] <= y1)
            )
            mask = mask & crop_mask

        coords = coords_all[mask]
        square_size = _infer_square_size(adata, library_id, coords, binsize=binsize)
        (x_min, y_min, x_max, y_max), (x_padding, y_padding) = _compute_data_extent_from_coords(coords, extra_pad=square_size / 2.0)

        img, img_extent = _process_background_image(spatial_info, img_key)
        if img is not None and img_extent is not None:
            if color_key is None:
                if invert_y:
                    current_ax.imshow(img, extent=img_extent, origin='upper', alpha=1.0)
                else:
                    current_ax.imshow(img[::-1], extent=img_extent, origin='upper', alpha=1.0)
            else:
                if invert_y:
                    current_ax.imshow(img, extent=img_extent, origin='upper', alpha=alpha_img)
                else:
                    current_ax.imshow(img[::-1], extent=img_extent, origin='upper', alpha=alpha_img)

        if color_key is None:
            if img is not None and img_extent is not None:
                current_ax.set_xlim(img_extent[0], img_extent[1])
                current_ax.set_aspect('equal')
                if invert_y:
                    current_ax.set_ylim(img_extent[2], img_extent[3])
                else:
                    current_ax.set_ylim(img_extent[3], img_extent[2])
                if xlabel is not None:
                    current_ax.set_xlabel(xlabel)
                if ylabel is not None:
                    current_ax.set_ylabel(ylabel)
                current_ax.tick_params(axis='both', which='major', labelsize=10)
            else:
                _apply_spatial_axis_formatting(
                    current_ax, x_min, y_min, x_max, y_max, x_padding, y_padding,
                    xlabel, ylabel, show_ticks=True, force_show_ticks=True,
                    invert_y=invert_y
                )
            axes_list.append(current_ax)
            continue

        values = _extract_color_values(adata, color_key, mask)

        _plot_squarebin_values(
            current_ax,
            coords=coords,
            values=values,
            square_size=square_size,
            cmap=cmap,
            palette=palette,
            vmin=vmin,
            vmax=vmax,
            alpha=alpha,
            legend=legend,
            edges_width=edges_width,
            edges_color=edges_color,
            na_color=na_color,
            rasterized=rasterized,
            shape=shape,
        )

        _apply_spatial_axis_formatting(
            current_ax, x_min, y_min, x_max, y_max, x_padding, y_padding,
            xlabel, ylabel, show_ticks=show_ticks, force_show_ticks=False,
            invert_y=invert_y
        )
        current_ax.set_title(color_key)
        axes_list.append(current_ax)

    if show:
        if ax is None:
            fig.tight_layout(rect=[0, 0, 0.95, 1])
        plt.show()

    return axes_list[0] if len(axes_list) == 1 else axes_list


# Short alias for consistency with bin-level workflows
spatial_bin = spatial_squarebin

def mark_region(
    ax: plt.Axes,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    edges_color: str = 'red',
    edges_width: float = 2.0,
    fill_color: Optional[str] = None,
    fill_alpha: float = 0.15,
    zorder: int = 100,
    refresh: bool = True,
    show: bool = True,
):
    """
    Mark a rectangular region on a spatial plot by drawing a rectangle.

    This function draws a rectangle on the given axes to highlight a specific
    spatial region. It works with any spatial plot (``spatial_cell``,
    ``spatial_squarebin``, ``sc.pl.spatial``) regardless of the ``invert_y``
    setting.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to draw the rectangle on.
    xlim : tuple, optional
        Tuple of (x_min, x_max) to define the x-range of the region.
        If None, uses the current x-axis limits.
    ylim : tuple, optional
        Tuple of (y_min, y_max) to define the y-range of the region.
        If None, uses the current y-axis limits.
    edges_color : str, default 'red'
        Color of the rectangle edges.
    edges_width : float, default 2.0
        Width of the rectangle edges.
    fill_color : str, optional
        If provided, fills the rectangle with this color at ``fill_alpha``
        opacity, making the region more visible against dense H&E backgrounds.
    fill_alpha : float, default 0.15
        Opacity of the fill when ``fill_color`` is set.
    zorder : int, default 100
        Z-order for rendering. High value ensures the rectangle is drawn on
        top of all other plot elements.
    refresh : bool, default True
        If True, calls ``ax.figure.canvas.draw_idle()`` to update the display
        after adding the rectangle. Set to False when adding multiple regions
        before a single refresh.
    show : bool, default True
        Whether to display the figure via ``plt.show()``. Set to ``False``
        when you want to defer display (e.g., when adding multiple regions
        or calling ``plt.show()`` manually).

    Returns
    -------
    matplotlib.patches.Rectangle
        The rectangle patch object that was added to the axes.

    Examples
    --------
    >>> import trackcell as tcl
    >>>
    >>> # Simplest usage: mark_region auto-shows by default
    >>> # (spatial_cell returns ax even when show=True)
    >>> ax = tcl.pl.spatial_cell(adata, color="CellType")
    >>> tcl.pl.mark_region(ax, xlim=(54500, 56000), ylim=(15000, 16000))
    >>>
    >>> # For full control, use show=False in both functions
    >>> fig, ax = plt.subplots(figsize=(10, 10))
    >>> tcl.pl.spatial_cell(adata, color="CellType", ax=ax, show=False)
    >>> tcl.pl.mark_region(
    ...     ax, xlim=(54500, 56000), ylim=(15000, 16000),
    ...     fill_color='red',
    ...     edges_width=3.0
    ... )
    >>> plt.show()
    >>>
    >>> # Mark multiple regions efficiently
    >>> fig, ax = plt.subplots(figsize=(10, 10))
    >>> tcl.pl.spatial_cell(adata, color="CellType", ax=ax, show=False)
    >>> tcl.pl.mark_region(ax, xlim=(40000, 42000), ylim=(5000, 7000),
    ...                    edges_color='cyan', fill_color='cyan',
    ...                    refresh=False, show=False)
    >>> tcl.pl.mark_region(ax, xlim=(55000, 57000), ylim=(15000, 17000),
    ...                    edges_color='yellow', fill_color='yellow',
    ...                    refresh=False, show=False)
    >>> tcl.pl.mark_region(ax, xlim=(60000, 62000), ylim=(10000, 12000),
    ...                    edges_color='magenta', fill_color='magenta')
    >>> # plt.show() is called by the last mark_region (show=True)
    """
    from matplotlib.patches import Rectangle

    # Get current axis limits if xlim/ylim are not provided
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()

    x_min, x_max = xlim
    y_min, y_max = ylim

    # Normalize: ensure positive width / height for robust rendering.
    # When the y-axis is inverted (invert_y=True), ax.get_ylim() returns
    # (bottom, top) with bottom > top, giving negative height.  Normalizing
    # guarantees the Rectangle patch is always well-formed.
    if x_min > x_max:
        x_min, x_max = x_max, x_min
    if y_min > y_max:
        y_min, y_max = y_max, y_min

    width = x_max - x_min
    height = y_max - y_min

    # Create rectangle
    rect = Rectangle(
        (x_min, y_min),
        width,
        height,
        linewidth=edges_width,
        edgecolor=edges_color,
        facecolor=fill_color if fill_color else 'none',
        alpha=fill_alpha if fill_color else None,
        zorder=zorder,
    )

    # Add rectangle to axes
    ax.add_patch(rect)

    # Refresh the display (for interactive backends)
    if refresh and ax.figure is not None:
        ax.figure.canvas.draw_idle()

    # Auto-show the figure so the rectangle is visible
    if show and ax.figure is not None:
        plt.show()

    return rect
