"""
Plotting functions for TrackCell package.

This module provides functions for visualizing spatial transcriptomics data,
including cell polygon visualization.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize, ListedColormap
from typing import Optional, Union, List
import warnings

try:
    import geopandas as gpd
    from shapely import wkt
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    warnings.warn("geopandas and shapely are required for spatial_cell function")


def spatial_cell(
    adata,
    color: Optional[Union[str, List[str]]] = None,
    groups: Optional[List[str]] = None,
    library_id: str = "spatial",
    size: float = 1.0,
    figsize: Optional[tuple] = None,
    cmap: str = "viridis",
    palette: Optional[dict] = None,
    img_key: Optional[str] = None,
    basis: str = "spatial",
    edges_width: float = 0.5,
    edges_color: str = "black",
    alpha: float = 0.8,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    **kwargs
):
    """
    Plot spatial transcriptomics data with cell polygons instead of points.
    
    This function visualizes cells as polygons (from cell segmentation) rather than
    simple points, providing a more accurate representation of cell boundaries.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial information and cell geometries.
    color : str or list of str, optional
        Keys for observation/categorical or continuous variables to color cells.
        Can be a single key or a list of keys for multiple plots.
    groups : list of str, optional
        Subset of groups to plot. If None, plots all groups.
    library_id : str, default "spatial"
        Key in adata.uns["spatial"] containing spatial information.
    size : float, default 1.0
        Size scaling factor for cells (not used for polygons, kept for API compatibility).
    figsize : tuple, optional
        Figure size (width, height) in inches.
    cmap : str, default "viridis"
        Colormap for continuous values.
    palette : dict, optional
        Dictionary mapping category names to colors for categorical variables.
    img_key : str, optional
        Key in adata.uns["spatial"][library_id]["images"] for background image.
        If None, uses "hires" if available.
    basis : str, default "spatial"
        Key in adata.obsm containing spatial coordinates (for fallback).
    edges_width : float, default 0.5
        Width of cell polygon edges.
    edges_color : str, default "black"
        Color of cell polygon edges.
    alpha : float, default 0.8
        Transparency of cell polygons.
    show : bool, default True
        Whether to display the plot.
    ax : matplotlib.Axes, optional
        Axes object to plot on. If None, creates a new figure.
    **kwargs
        Additional arguments passed to matplotlib plotting functions.
    
    Returns
    -------
    matplotlib.Axes or list of matplotlib.Axes
        Axes object(s) containing the plot.
    
    Examples
    --------
    >>> import trackcell as tcl
    >>> adata = tcl.io.read_hd_cellseg("path/to/data", sample="sample1")
    >>> tcl.pl.spatial_cell(adata, color="classification")
    >>> tcl.pl.spatial_cell(adata, color="Cluster-2_dist", cmap="Reds")
    >>> tcl.pl.spatial_cell(adata, color="classification", groups=["Cluster-1", "Cluster-2"])
    """
    if not HAS_GEOPANDAS:
        raise ImportError("geopandas and shapely are required for spatial_cell function")
    
    # Check if geometries are available
    if "spatial" not in adata.uns:
        raise ValueError("`adata.uns['spatial']` is required but missing.")
    
    if library_id not in adata.uns["spatial"]:
        raise ValueError(f"`library_id` '{library_id}' not found in `adata.uns['spatial']`.")
    
    spatial_info = adata.uns["spatial"][library_id]
    
    # Get geometries - try GeoDataFrame first, then fallback to WKT strings in obs
    geometries = None
    use_wkt = False
    
    if "geometries" in spatial_info:
        geometries = spatial_info["geometries"]
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
        
        # Get background image if available
        if img_key is None:
            img_key = "hires" if "hires" in spatial_info.get("images", {}) else None
        
        if img_key and "images" in spatial_info and img_key in spatial_info["images"]:
            img = spatial_info["images"][img_key]
            current_ax.imshow(img, extent=[0, img.shape[1], img.shape[0], 0], origin='upper', alpha=0.5)
        
        # Filter cells if groups is specified
        if groups is not None:
            if color_key is not None and color_key in adata.obs.columns:
                mask = adata.obs[color_key].isin(groups)
            else:
                # If no color key, need another way to filter
                raise ValueError("`groups` parameter requires a valid `color` parameter.")
        else:
            mask = pd.Series(True, index=adata.obs_names)
        
        # Get cell indices to plot
        cells_to_plot = adata.obs_names[mask]
        
        # Prepare color values
        if color_key is None:
            # No coloring, use default color
            color_values = None
            face_colors = ['lightblue'] * len(cells_to_plot)
        elif color_key in adata.obs.columns:
            color_data = adata.obs.loc[cells_to_plot, color_key]
            
            # Check if continuous or categorical
            if pd.api.types.is_numeric_dtype(color_data):
                # Continuous values
                norm = Normalize(vmin=color_data.min(), vmax=color_data.max())
                cmap_obj = plt.get_cmap(cmap)
                face_colors = [cmap_obj(norm(v)) for v in color_data]
                
                # Add colorbar
                sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
                sm.set_array([])
                plt.colorbar(sm, ax=current_ax, label=color_key)
            else:
                # Categorical values
                unique_cats = color_data.unique()
                
                if palette is not None:
                    # Use provided palette
                    face_colors = [palette.get(cat, 'gray') for cat in color_data]
                else:
                    # Generate colors
                    if hasattr(adata.uns, 'classification_colors') and color_key == 'classification':
                        # Use stored classification colors if available
                        cat_palette = adata.uns['classification_colors']
                    else:
                        # Generate default colors
                        n_cats = len(unique_cats)
                        default_cmap = plt.get_cmap('tab20' if n_cats <= 20 else 'tab20b')
                        cat_palette = {cat: default_cmap(i / n_cats) for i, cat in enumerate(unique_cats)}
                    
                    face_colors = [cat_palette.get(cat, 'gray') for cat in color_data]
        else:
            raise ValueError(f"`color` key '{color_key}' not found in `adata.obs`.")
        
        # Create patches for polygons
        patches = []
        valid_colors = []
        
        for cell_id in cells_to_plot:
            # Get geometry
            if use_wkt:
                # Read from WKT string in obs
                if cell_id in adata.obs.index and pd.notna(adata.obs.loc[cell_id, "geometry"]):
                    geom_str = adata.obs.loc[cell_id, "geometry"]
                    try:
                        geom = wkt.loads(geom_str)
                    except Exception:
                        continue
                else:
                    continue
            else:
                # Read from GeoDataFrame
                if cell_id in geometries.index:
                    geom = geometries.loc[cell_id, "geometry"]
                else:
                    continue
            
            if geom is not None:
                # Convert shapely geometry to matplotlib polygon
                if hasattr(geom, 'exterior'):
                    # Polygon
                    coords = np.array(geom.exterior.coords)
                elif hasattr(geom, 'coords'):
                    # Point or LineString (fallback to point)
                    coords = np.array(geom.coords)
                else:
                    continue
                
                # Create polygon patch
                polygon = mpatches.Polygon(coords, closed=True)
                patches.append(polygon)
                
                # Get corresponding color
                if color_key is None:
                    valid_colors.append('lightblue')
                else:
                    color_idx = list(cells_to_plot).index(cell_id)
                    valid_colors.append(face_colors[color_idx])
        
        # Create patch collection for efficient rendering
        if patches:
            p = PatchCollection(
                patches,
                facecolors=valid_colors,
                edgecolors=edges_color,
                linewidths=edges_width,
                alpha=alpha,
                **kwargs
            )
            current_ax.add_collection(p)
        
        # Set axis properties
        current_ax.set_aspect('equal')
        current_ax.invert_yaxis()  # Match image coordinates
        current_ax.set_xlabel('X coordinate')
        current_ax.set_ylabel('Y coordinate')
        
        if color_key:
            current_ax.set_title(color_key)
        
        axes_list.append(current_ax)
    
    if show:
        plt.tight_layout()
        plt.show()
    
    if len(axes_list) == 1:
        return axes_list[0]
    else:
        return axes_list

