"""
Plotting functions for TrackCell package.

This module provides functions for visualizing spatial transcriptomics data,
including cell polygon visualization.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from typing import Optional, Union, List
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


def spatial_cell(
    adata,
    color: Optional[Union[str, List[str]]] = None,
    groups: Optional[List[str]] = None,
    library_id: Optional[str] = None,
    size: float = 1.0,
    figsize: Optional[tuple] = None,
    cmap: str = "viridis",
    palette: Optional[dict] = None,
    img_key: Optional[str] = None,
    basis: str = "spatial",
    edges_width: float = 0.5,
    edges_color: str = "black",
    alpha: float = 0.8,
    alpha_img: float = 0.5,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    legend: bool = True,
    xlabel: Optional[str] = "spatial 1",
    ylabel: Optional[str] = "spatial 2",
    show_ticks: bool = False,
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
    groups : list of str, optional
        Subset of groups to plot. If None, plots all groups.
    library_id : str, optional
        Key in adata.uns["spatial"] containing spatial information.
        If None, uses the first available library_id (similar to sc.pl.spatial).
        Default is None, which will auto-select the first library_id.
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
            current_ax.imshow(img, extent=img_extent, origin='upper', alpha=alpha_img)
        
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
            temp_geometries = geometries.loc[cells_to_plot, "geometry"]
            valid_cells = list(cells_to_plot)
        
        if len(temp_geometries) == 0:
            warnings.warn("No valid geometries found for plotting.")
            axes_list.append(current_ax)
            continue
        
        # Create GeoDataFrame with color data
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
        elif color_key in adata.var_names:
            # Get color data from gene expression (var)
            # Find the gene index
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
        else:
            # Check if it's in layers
            if hasattr(adata, 'layers') and color_key in adata.layers:
                # Try to get from layers (assuming it's a gene name)
                # This is less common, but handle it if needed
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
        
        # Determine if continuous or categorical
        is_categorical = False
        use_custom_palette = False
        custom_palette = None
        
        if plot_column is not None:
            is_categorical = not pd.api.types.is_numeric_dtype(temp_gdf[plot_column])
            
            # Check if we need to use custom palette
            if is_categorical:
                if palette is not None:
                    use_custom_palette = True
                    custom_palette = palette
                elif hasattr(adata.uns, 'classification_colors') and color_key == 'classification':
                    use_custom_palette = True
                    custom_palette = adata.uns['classification_colors']
        
        # Prepare plot arguments for GeoDataFrame.plot()
        # Filter out parameters that GeoPandas doesn't support or that we handle separately
        # GeoPandas plot() doesn't support 'palette' parameter directly
        # Note: 'cmap' is supported for continuous values, but we handle it explicitly
        excluded_params = {'palette', 'palete'}  # Handle palette separately
        filtered_kwargs = {k: v for k, v in kwargs.items() if k not in excluded_params}
        
        plot_kwargs = {
            'ax': current_ax,
            'edgecolor': edges_color,
            'linewidth': edges_width,
            'alpha': alpha,
            **filtered_kwargs
        }
        
        if plot_column is not None:
            plot_kwargs['column'] = plot_column
            
            if is_categorical:
                # Categorical values
                plot_kwargs['categorical'] = True
                plot_kwargs['legend'] = False  # Disable automatic legend, we'll add it manually
                
                # If using custom palette, map categories to colors and use color parameter
                # This is more efficient than manually plotting each category
                # GeoPandas.plot() supports color parameter as array/Series (see documentation)
                if use_custom_palette:
                    # Map each category to its color from palette
                    # Create a color Series with same index as temp_gdf
                    color_series = temp_gdf[plot_column].map(
                        lambda x: custom_palette.get(x, 'gray') if pd.notna(x) else 'gray'
                    )
                    # Use color parameter instead of column for custom palette
                    plot_kwargs['color'] = color_series
                    # Remove 'column' from plot_kwargs when using color
                    plot_kwargs.pop('column', None)
                # else: use default GeoPandas categorical plotting with column
                # (plot_kwargs already has 'column' set)
                
                # Plot using GeoDataFrame.plot() - it handles both column and color
                temp_gdf.plot(**plot_kwargs)
                
                # Create legend manually from unique categories
                if legend:
                    unique_cats = sorted(temp_gdf[plot_column].dropna().unique())
                    from matplotlib.patches import Patch
                    
                    if use_custom_palette:
                        # Use colors from custom palette
                        legend_elements = [
                            Patch(facecolor=custom_palette.get(cat, 'gray'), 
                                  label=str(cat))
                            for cat in unique_cats
                        ]
                    else:
                        # Generate colors matching GeoPandas default
                        n_cats = len(unique_cats)
                        default_cmap = plt.get_cmap('tab20' if n_cats <= 20 else 'tab20b')
                        legend_elements = [
                            Patch(facecolor=default_cmap(i / n_cats), label=str(cat))
                            for i, cat in enumerate(unique_cats)
                        ]
                    
                    if legend_elements:
                        current_ax.legend(handles=legend_elements, 
                                        bbox_to_anchor=(1.05, 1), 
                                        loc='upper left',
                                        frameon=True)
            else:
                # Continuous values
                plot_kwargs['cmap'] = cmap
                plot_kwargs['legend'] = False  # Disable automatic colorbar, we'll add it manually
                # Plot using GeoDataFrame.plot()
                temp_gdf.plot(**plot_kwargs)
                
                # Add colorbar outside the plot if requested
                if legend:
                    # Get the mappable from the plot
                    # GeoPandas plot returns a collection, we need to get the colormap
                    import matplotlib.cm as cm
                    from matplotlib.colors import Normalize
                    
                    # Get value range for normalization
                    values = temp_gdf[plot_column]
                    vmin, vmax = values.min(), values.max()
                    norm = Normalize(vmin=vmin, vmax=vmax)
                    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)
                    sm.set_array([])
                    
                    # Add colorbar outside the plot
                    cbar = plt.colorbar(sm, ax=current_ax, 
                                       label=color_key,
                                       shrink=0.8,
                                       pad=0.02)
        else:
            # No coloring, just plot geometries
            temp_gdf.plot(ax=current_ax, color='lightblue',
                         edgecolor=edges_color, linewidth=edges_width,
                         alpha=alpha, **kwargs)
        
        # Set axis properties
        current_ax.set_aspect('equal')
        current_ax.invert_yaxis()  # Match image coordinates
        
        # Set axis limits based on actual data range
        # Add small padding (5% of range) for better visualization
        x_range = x_max - x_min
        y_range = y_max - y_min
        x_padding = x_range * 0.05 if x_range > 0 else 1
        y_padding = y_range * 0.05 if y_range > 0 else 1
        
        current_ax.set_xlim(x_min - x_padding, x_max + x_padding)
        current_ax.set_ylim(y_max + y_padding, y_min - y_padding)  # Inverted for y-axis
        
        # Set axis labels
        if xlabel is not None:
            current_ax.set_xlabel(xlabel)
        if ylabel is not None:
            current_ax.set_ylabel(ylabel)
        
        # Control ticks visibility
        if not show_ticks:
            current_ax.set_xticks([])
            current_ax.set_yticks([])
        
        if color_key:
            current_ax.set_title(color_key)
        
        axes_list.append(current_ax)
    
    if show:
        # Adjust layout to make room for legend/colorbar on the right
        # This is similar to scanpy's approach
        if ax is None:  # Only adjust if we created the figure
            fig.tight_layout(rect=[0, 0, 0.95, 1])  # Leave 5% space on the right
        else:
            fig.tight_layout()
        plt.show()
    
    if len(axes_list) == 1:
        return axes_list[0]
    else:
        return axes_list
