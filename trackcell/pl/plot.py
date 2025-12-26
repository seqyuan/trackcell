"""
Plotting functions for TrackCell package.

This module provides functions for visualizing spatial transcriptomics data,
including cell polygon visualization.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
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
        Note: When `color=None`, ticks are automatically shown regardless of this setting.
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
            # If color is None, use full opacity for the image
            img_alpha = 1.0 if color_key is None else alpha_img
            current_ax.imshow(img, extent=img_extent, origin='upper', alpha=img_alpha)
        
        # If color is None, skip geometry plotting and only show HE image
        if color_key is None:
            # Set axis limits based on image extent or data coordinates
            if img_extent is not None:
                current_ax.set_xlim(img_extent[0], img_extent[1])
                current_ax.set_ylim(img_extent[2], img_extent[3])
            else:
                # Fallback to data coordinates
                current_ax.set_xlim(x_min, x_max)
                current_ax.set_ylim(y_max, y_min)  # Inverted for y-axis
            
            # Set axis properties
            current_ax.set_aspect('equal')
            current_ax.invert_yaxis()  # Match image coordinates
            
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
        
        # Validate and fix total_bounds after creating temp_gdf
        # This is critical for subset data where geometries may not be properly synchronized
        max_retries = 2
        retry_count = 0
        while retry_count <= max_retries:
            try:
                # Validate total_bounds
                bounds = temp_gdf.total_bounds
                if not all(np.isfinite(bounds)) or bounds[2] <= bounds[0] or bounds[3] <= bounds[1]:
                    # Invalid bounds detected - filter out problematic geometries
                    if retry_count == 0:
                        warnings.warn(
                            f"Invalid geometry bounds detected (bounds: {bounds}). "
                            f"Filtering out problematic geometries. "
                            f"Consider using tcl.io.sync_geometries_after_subset() after subsetting."
                        )
                    
                    # Filter geometries with valid bounds
                    valid_mask = pd.Series(True, index=temp_gdf.index)
                    for idx in temp_gdf.index:
                        try:
                            geom = temp_gdf.loc[idx, 'geometry']
                            if geom is None or pd.isna(geom):
                                valid_mask.loc[idx] = False
                                continue
                            geom_bounds = geom.bounds
                            if not all(np.isfinite(geom_bounds)) or geom_bounds[2] <= geom_bounds[0] or geom_bounds[3] <= geom_bounds[1]:
                                valid_mask.loc[idx] = False
                        except Exception:
                            valid_mask.loc[idx] = False
                    
                    temp_gdf = temp_gdf[valid_mask]
                    
                    if len(temp_gdf) == 0:
                        warnings.warn("No valid geometries remaining after filtering. Skipping plot.")
                        axes_list.append(current_ax)
                        break  # Exit the retry loop and continue to next color
                    
                    # Update valid_cells to match filtered temp_gdf
                    valid_cells = temp_gdf.index.tolist()
                    
                    # Re-validate bounds
                    bounds = temp_gdf.total_bounds
                    if not all(np.isfinite(bounds)) or bounds[2] <= bounds[0] or bounds[3] <= bounds[1]:
                        # Still invalid after filtering - will use aspect='equal' as fallback
                        if retry_count == 0:
                            warnings.warn(
                                f"Still invalid bounds after filtering (bounds: {bounds}). "
                                f"Will use aspect='equal' to avoid errors."
                            )
                        retry_count += 1
                        continue
                    else:
                        # Bounds are now valid
                        break
                else:
                    # Bounds are valid
                    break
            except Exception as e:
                if retry_count == 0:
                    warnings.warn(
                        f"Error validating geometry bounds: {e}. "
                        f"Will use aspect='equal' to avoid errors."
                    )
                retry_count += 1
                continue
        
        # If we still have invalid bounds after retries, we'll set aspect='equal' later
        use_equal_aspect = False
        if retry_count > max_retries:
            try:
                bounds = temp_gdf.total_bounds
                if not all(np.isfinite(bounds)) or bounds[2] <= bounds[0] or bounds[3] <= bounds[1]:
                    use_equal_aspect = True
            except Exception:
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
            
            if is_categorical:
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


def mark_region(
    ax: plt.Axes,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    edges_color: str = 'red',
    edges_width: float = 1.0
):
    """
    Mark a rectangular region on a spatial plot by drawing a rectangle.
    
    This function draws a rectangle on the given axes to highlight a specific
    spatial region. It can be used with any spatial plot.
    
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
    edges_width : float, default 1.0
        Width of the rectangle edges.
    
    Returns
    -------
    matplotlib.patches.Rectangle
        The rectangle patch object that was added to the axes.
    
    Examples
    --------
    >>> import trackcell as tcl
    >>> import matplotlib.pyplot as plt
    >>> 
    >>> # Plot with spatial_cell and mark a region
    >>> fig, ax = plt.subplots(figsize=(10, 10))
    >>> tcl.pl.spatial_cell(adata, color="CellType", ax=ax)
    >>> tcl.pl.mark_region(ax, xlim=(54500, 56000), ylim=(15000, 16000))
    >>> 
    >>> # Mark a region on any plot
    >>> fig, ax = plt.subplots(figsize=(10, 10))
    >>> # ... create your plot on ax ...
    >>> tcl.pl.mark_region(ax, xlim=(54500, 56000), ylim=(15000, 16000), 
    ...                    edges_color='blue', edges_width=2.0)
    """
    from matplotlib.patches import Rectangle
    
    # Get current axis limits if xlim/ylim are not provided
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()
    
    x_min, x_max = xlim
    y_min, y_max = ylim
    
    # Create rectangle
    rect = Rectangle(
        (x_min, y_min),
        x_max - x_min,
        y_max - y_min,
        linewidth=edges_width,
        edgecolor=edges_color,
        facecolor='none'
    )
    
    # Add rectangle to axes
    ax.add_patch(rect)
    
    return rect
