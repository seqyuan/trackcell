Visualization
=============

This section covers visualization of spatial transcriptomics data using TrackCell.

Basic Plotting with Cell Polygons
----------------------------------

TrackCell provides a specialized plotting function that visualizes cells as polygons
instead of points, providing a more accurate representation of cell boundaries:

.. code-block:: python

   # Plot cells as polygons (requires data loaded with read_hd_cellseg)
   tcl.pl.spatial_cell(
       adata, 
       color="classification",  # Color by cell type
       groups=['Cluster-2', 'Cluster-3'],  # Optional: filter specific groups
       figsize=(6, 6),
       edges_width=0.5,
       edges_color="black",
       alpha=0.8
   )

   # Plot continuous values (e.g., distance to a label)
   tcl.pl.spatial_cell(
       adata,
       color="Cluster-2_dist",  # Distance to Cluster-2
       cmap="Reds",
       figsize=(6, 6)
   )

   # Plot gene expression
   tcl.pl.spatial_cell(
       adata,
       color="PDPN",  # Gene name
       cmap="viridis",
       figsize=(6, 6)
   )

   # Plot multiple variables in subplots
   tcl.pl.spatial_cell(
       adata,
       color=["classification", "Cluster-2_dist"],  # Two subplots
       figsize=(12, 6)
   )

Performance Optimization for Large Datasets
--------------------------------------------

When working with large datasets (e.g., 200,000+ cells), visualization can be slow. Here are recommended strategies to improve performance:

Recommended Combination Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For datasets with 200,000+ cells, we recommend using a combination of optimization techniques:

**Option 1: Optimized Polygon Plotting**

Use ``groups`` parameter to filter cells, disable edges, and use low-resolution images:

.. code-block:: python

   import trackcell as tcl
   
   # Optimized polygon plotting for large datasets
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       groups=['Cluster-2', 'Cluster-3'],  # Only plot cells of interest
       edges_width=0,                      # Disable edges for better performance
       img_key="lowres",                   # Use low-resolution image
       figsize=(10, 10)
   )

**Option 2: Fast Point-based Preview**

For quick exploration, use point-based visualization which is much faster:

.. code-block:: python

   import scanpy as sc
   
   # Fast point-based preview
   sc.pl.spatial(
       adata,
       color='classification',
       spot_size=0.5,
       size=0.3,
       groups=['Cluster-2', 'Cluster-3']
   )

Performance Comparison
~~~~~~~~~~~~~~~~~~~~~~

Expected performance for different approaches with ~230,000 cells:

* **Full polygon plotting**: Several minutes to tens of minutes
* **Using groups** (filtering to ~10% of cells): ~10-30 seconds
* **Point-based plotting**: ~1-5 seconds
* **Downsampling to 10,000 cells**: ~5-15 seconds

Optimization Strategies
~~~~~~~~~~~~~~~~~~~~~~~

**Strategy 1: Filter by Cell Groups**

The most effective optimization is to plot only cells of interest using the ``groups`` parameter:

.. code-block:: python

   # Plot only specific cell types
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       groups=['Cluster-2', 'Cluster-3', 'Cluster-5']  # Only these cell types
   )

**Strategy 2: Spatial Region Cropping**

Crop to a specific spatial region of interest:

.. code-block:: python

   import numpy as np
   
   # Define region of interest
   x_min, x_max = 1000, 5000
   y_min, y_max = 1000, 5000
   
   # Create mask for spatial coordinates
   spatial_coords = adata.obsm['spatial']
   mask = ((spatial_coords[:, 0] >= x_min) & (spatial_coords[:, 0] <= x_max) &
           (spatial_coords[:, 1] >= y_min) & (spatial_coords[:, 1] <= y_max))
   
   # Create subset
   adata_subset = adata[mask].copy()
   
   # Plot subset
   tcl.pl.spatial_cell(adata_subset, color="classification")


**Strategy 3: Use Point-based Visualization**

For large datasets, point-based visualization is much faster:

.. code-block:: python

   import scanpy as sc
   
   # Point-based visualization (much faster)
   sc.pl.spatial(
       adata,
       color='classification',
       spot_size=1,      # Small spots
       size=0.5,         # Further reduce size
       groups=['Cluster-2', 'Cluster-3']
   )

**Strategy 4: Disable Edge Rendering**

Disable polygon edges to improve rendering performance:

.. code-block:: python

   tcl.pl.spatial_cell(
       adata,
       color="classification",
       edges_width=0,  # Disable edges
       groups=['Cluster-2', 'Cluster-3']
   )

**Strategy 5: Use Low-Resolution Images**

Use low-resolution background images when available:

.. code-block:: python

   tcl.pl.spatial_cell(
       adata,
       color="classification",
       img_key="lowres",  # Use low-resolution image
       groups=['Cluster-2', 'Cluster-3']
   )

Best Practices
~~~~~~~~~~~~~~

1. **Always use GeoDataFrame format**: Ensure your data uses ``adata.uns['spatial'][sample]['geometries']`` (GeoDataFrame) rather than WKT strings for better performance.

2. **Start with point-based visualization**: Use ``sc.pl.spatial()`` for initial exploration, then switch to polygon-based visualization for detailed analysis.

3. **Filter before plotting**: Always use ``groups`` or spatial cropping to reduce the number of cells before plotting.

4. **Combine strategies**: Use multiple optimization strategies together for best results.

5. **Save intermediate results**: For repeated visualization, consider saving filtered subsets to avoid repeated filtering operations.

Example Workflow
~~~~~~~~~~~~~~~~

Here's a recommended workflow for visualizing large datasets:

.. code-block:: python

   import trackcell as tcl
   import scanpy as sc
   
   # Step 1: Quick overview with point-based plot
   sc.pl.spatial(
       adata,
       color='classification',
       spot_size=0.5,
       size=0.3
   )
   
   # Step 2: Detailed view of specific cell types with polygons
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       groups=['Cluster-2', 'Cluster-3'],  # Focus on specific types
       edges_width=0,                      # Optimize performance
       img_key="lowres",                   # Use low-res image
       figsize=(6, 6)
   )
   
   # Step 3: High-resolution view of specific region
   # (Use spatial cropping as shown in Strategy 2)

