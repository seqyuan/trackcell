Visualization Guide
===================

This guide provides best practices and optimization strategies for visualizing spatial transcriptomics data, especially for large datasets with hundreds of thousands of cells.

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
-----------------------

Strategy 1: Filter by Cell Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most effective optimization is to plot only cells of interest using the ``groups`` parameter:

.. code-block:: python

   # Plot only specific cell types
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       groups=['Cluster-2', 'Cluster-3', 'Cluster-5']  # Only these cell types
   )

**Benefits**: Dramatically reduces the number of cells to render
**Use case**: When you're interested in specific cell types

Strategy 2: Spatial Region Cropping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

**Benefits**: Only renders cells in the region of interest
**Use case**: When focusing on specific spatial regions

Strategy 3: Downsampling
~~~~~~~~~~~~~~~~~~~~~~~~

For quick exploration, randomly sample a subset of cells:

.. code-block:: python

   import numpy as np
   
   # Random sampling
   np.random.seed(42)
   n_sample = 10000  # Number of cells to sample
   sample_idx = np.random.choice(len(adata), size=n_sample, replace=False)
   adata_sample = adata[sample_idx].copy()
   
   # Plot sampled data
   tcl.pl.spatial_cell(adata_sample, color="classification")

**Benefits**: Fast preview of overall distribution
**Use case**: Quick exploration of data patterns

Strategy 4: Use Point-based Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

**Benefits**: Very fast, suitable for large datasets
**Use case**: When precise cell boundaries are not required

Strategy 5: Disable Edge Rendering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Disable polygon edges to improve rendering performance:

.. code-block:: python

   tcl.pl.spatial_cell(
       adata,
       color="classification",
       edges_width=0,  # Disable edges
       groups=['Cluster-2', 'Cluster-3']
   )

**Benefits**: Reduces rendering overhead
**Use case**: When cell boundaries are not critical

Strategy 6: Use Low-Resolution Images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use low-resolution background images when available:

.. code-block:: python

   tcl.pl.spatial_cell(
       adata,
       color="classification",
       img_key="lowres",  # Use low-resolution image
       groups=['Cluster-2', 'Cluster-3']
   )

**Benefits**: Faster image loading and rendering
**Use case**: When high-resolution images are not necessary

Strategy 7: Tile-based Plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For saving high-resolution images, divide the plot into tiles:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   # Define regions
   spatial_coords = adata.obsm['spatial']
   
   regions = [
       {'x': (0, 5000), 'y': (0, 5000), 'name': 'region1'},
       {'x': (5000, 10000), 'y': (0, 5000), 'name': 'region2'},
       # Add more regions as needed
   ]
   
   for region in regions:
       mask = ((spatial_coords[:, 0] >= region['x'][0]) & 
               (spatial_coords[:, 0] < region['x'][1]) &
               (spatial_coords[:, 1] >= region['y'][0]) & 
               (spatial_coords[:, 1] < region['y'][1]))
       
       adata_region = adata[mask].copy()
       ax = tcl.pl.spatial_cell(adata_region, color="classification", show=False)
       plt.savefig(f"{region['name']}.png", dpi=300, bbox_inches='tight')
       plt.close()

**Benefits**: Allows high-resolution output for large datasets
**Use case**: When you need to save publication-quality images

Best Practices
---------------

1. **Always use GeoDataFrame format**: Ensure your data uses ``adata.uns['spatial'][sample]['geometries']`` (GeoDataFrame) rather than WKT strings for better performance.

2. **Start with point-based visualization**: Use ``sc.pl.spatial()`` for initial exploration, then switch to polygon-based visualization for detailed analysis.

3. **Filter before plotting**: Always use ``groups`` or spatial cropping to reduce the number of cells before plotting.

4. **Combine strategies**: Use multiple optimization strategies together for best results.

5. **Save intermediate results**: For repeated visualization, consider saving filtered subsets to avoid repeated filtering operations.

Example Workflow
----------------

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
       figsize=(10, 10)
   )
   
   # Step 3: High-resolution view of specific region
   # (Use spatial cropping as shown in Strategy 2)

GPU Acceleration
----------------

**Note**: GPU acceleration is currently **not available** for polygon-based visualization because:

* The underlying libraries (geopandas, shapely, matplotlib) are CPU-based
* The main bottlenecks are in geometry conversion and rendering, not computation

For large datasets, the optimization strategies above are more effective than GPU acceleration.

