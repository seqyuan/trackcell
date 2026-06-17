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

Custom Color Palettes
----------------------

You can customize colors for categorical variables using the ``palette`` parameter.
The ``palette`` parameter accepts either a dictionary or a list/array of colors.

**Using a Dictionary (Category-to-Color Mapping)**

When using a dictionary, you explicitly map each category to a color:

.. code-block:: python

   # Define custom color palette as dictionary
   custom_palette = {
       'Cluster-1': 'red',
       'Cluster-2': 'blue',
       'Cluster-3': 'green',
       'Cluster-4': 'orange'
   }
   
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       palette=custom_palette,
       figsize=(6, 6)
   )

**Using a List/Array (Sequential Color Assignment)**

When using a list or array, colors are assigned to categories in alphabetical order:

.. code-block:: python

   # Define custom color palette as list
   # Colors will be assigned to categories in sorted order
   color_list = ['#FF0000', '#0000FF', '#00FF00', '#FFA500']
   
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       palette=color_list,
       figsize=(6, 6)
   )

   # Or using numpy array
   import numpy as np
   color_array = np.array(['red', 'blue', 'green', 'yellow'])
   
   tcl.pl.spatial_cell(
       adata,
       color="classification",
       palette=color_array,
       figsize=(6, 6)
   )

**Note**: If the palette has fewer colors than categories, colors will be cycled.
A warning will be issued if this occurs.


Dual-Color Visualization (Fill + Edge)
------------------------------------------

Use ``edge_color`` to color cell boundaries by a categorical column (e.g., cell type),
while ``color`` controls the fill (gene expression or continuous value).
This creates a two-layer visualization where:

* **Fill color** maps a continuous variable (gene expression, UMI counts, etc.)
* **Edge color** maps a categorical variable (cell type, cluster, etc.)

Fill = Gene Expression, Edge = Cell Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   # Fill: gene expression (continuous), Edge: cell type (categorical)
   tcl.pl.spatial_cell(
       adata,
       color='EPCAM',              # Fill: gene expression
       cmap='Reds',
       edge_color='cell_type',     # Edge: cell type
       edge_palette={
           'T cell': '#e41a1c',
           'B cell': '#377eb8',
           'Myeloid': '#4daf4a',
       },
       edges_width=1.2,
       alpha=0.7,
   )

Fill = Continuous obs, Edge = Categorical obs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Fill: continuous observation, Edge: categorical observation
   tcl.pl.spatial_cell(
       adata,
       color='total_counts',       # Fill: UMIs per cell
       cmap='YlOrRd',
       edge_color='cell_type',     # Edge: cell type
       edges_width=1.5,
       alpha=0.7,
   )

Using edge_palette
~~~~~~~~~~~~~~~~~~~

The ``edge_palette`` parameter accepts the same formats as ``palette``:

* **Dictionary**: Map categories to specific colors
* **List/Array**: Assign colors to categories in sorted order

If ``edge_palette`` is not specified, a default colormap (``tab10`` or ``tab20``) is used.

Both legends are automatically generated:

* **Fill**: colorbar (continuous) or legend (categorical), positioned on the right
* **Edge**: categorical legend, positioned below the plot

.. important::

   **Boundary overlap behavior**: Adjacent cells share edges.  When ``edge_color``
   uses a categorical column, overlapping edges are drawn by category order —
   the **last-drawn** category determines the visible edge colour at shared
   boundaries.  This is a cosmetic rendering artefact and does not affect the
   underlying data.


Multi-Gene Visualization (cell2location-style)
-----------------------------------------------

TrackCell provides ``tl.multigene_blend()`` for co-expression visualization
in two modes, controlled by the ``mode`` parameter.

Mode 1: Blended Composite (``mode='blend'``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps each gene to a base color, blends by expression level using weighted
RGB averaging. The result is a hex color per cell, stored in ``adata.obs``:

.. code-block:: python

   import trackcell as tcl

   # Compute blended colors
   tcl.tl.multigene_blend(
       adata,
       genes=['EPCAM', 'PECAM1', 'VWF'],
       mode='blend',
   )

   # Visualize as composite
   tcl.pl.spatial_cell(adata, color='multigene_blend')

Custom colors and gamma correction:

.. code-block:: python

   tcl.tl.multigene_blend(
       adata,
       genes=['EPCAM', 'PECAM1', 'VWF', 'ACTA2', 'PTPRC'],
       colors=['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3'],
       vmax_percentile=98,
       gamma=0.8,
   )

**How blending works**: Each gene gets a base RGB color. Per-cell expression
is normalized (percentile-clipped to [0, 1]), then blended as:

|color| = Σ (normalized_expr_i × RGB_i) / Σ normalized_expr_i

Zero-expression cells get the ``background`` color (default ``#e0e0e0``).

Mode 2: Faceted Subplots (``mode='facet'``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each gene gets its own panel with a single-hue colormap (white → gene color).
This matches the cell2location paper style:

.. code-block:: python

   tcl.tl.multigene_blend(
       adata,
       genes=['EPCAM', 'PECAM1', 'VWF', 'ACTA2'],
       mode='facet',
       ncols=2,           # 2 columns of subplots
       edges_width=0.3,
   )

Facet-specific parameters:

* ``ncols``: Number of columns (default 3)
* ``figsize``: Figure size (auto-computed if None)
* ``edges_width``: Cell boundary width (default 0.3)
* ``edges_color``: Cell boundary color (default ``'#cccccc30'``)

Square-bin visualization
------------------------

For Visium HD square-bin outputs loaded with ``read_hd_bin()``, use
``spatial_squarebin()``:

.. code-block:: python

   adata_bin = tcl.io.read_hd_bin(
       datapath="SpaceRanger4.0/Cse1/binned_outputs/square_016um",
       sample="Cse1",
       binsize=16,
   )

   # Show H&E image only, with coordinate range
   tcl.pl.spatial_squarebin(adata_bin, color=None)

   # Equivalent alias
   tcl.pl.spatial_bin(adata_bin, color=None)

   # Plot gene expression on square bins
   tcl.pl.spatial_squarebin(
       adata_bin,
       color="EPCAM",
       cmap="Reds",
       alpha=0.8,
       alpha_img=0.4,
   )

Key ``spatial_squarebin()`` parameters:

* ``color``: Key in ``adata.obs``, a gene name, ``None`` (HE-only), or a list of keys
* ``library_id``: Spatial library key (auto-detected if ``None``)
* ``binsize``: Bin size in micrometers (optional, auto-detected)
* ``shape``: Marker shape — ``'circle'`` (default) or ``'square'``
* ``invert_y``: If ``True`` (default), y increases top-to-bottom (image coordinates);
  set ``False`` for Cartesian convention
* ``crop_coord``: Optional ``(x_min, x_max, y_min, y_max)`` crop in spatial coordinates
* ``na_color``: Color for missing values (default ``'#d3d3d3'``)
* ``rasterized``: Whether to rasterize patches (default ``False``)

Compatibility with sc.pl.spatial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``sc.pl.spatial`` (point-based): use ``mode='facet'`` or manually create
individual gene subplots. The blend mode hex colors work best with trackcell's
``spatial_cell`` polygon rendering.


Marking Regions
---------------

The ``mark_region()`` function draws a rectangular highlight on any spatial plot,
useful for annotating regions of interest (ROI) in figures.  It works with
``spatial_cell``, ``spatial_squarebin``, ``sc.pl.spatial``, and any
matplotlib ``Axes``.

By default, ``mark_region`` calls ``plt.show()`` automatically (``show=True``),
so the simplest usage works out of the box:

.. code-block:: python

   import trackcell as tcl

   # spatial_cell returns ax even when show=True
   ax = tcl.pl.spatial_cell(adata, color="CellType")
   tcl.pl.mark_region(ax, xlim=(54500, 56000), ylim=(15000, 16000))

For full control over display timing, use ``show=False`` and call
``plt.show()`` manually:

.. code-block:: python

   import matplotlib.pyplot as plt
   import trackcell as tcl

   fig, ax = plt.subplots(figsize=(10, 10))
   tcl.pl.spatial_cell(adata, color="CellType", ax=ax, show=False)
   tcl.pl.mark_region(
       ax,
       xlim=(54500, 56000),
       ylim=(15000, 16000),
       edges_color='red',
       edges_width=2.0,
       show=False           # defer display
   )
   plt.show()

With ``spatial_squarebin()`` (Visium HD):

.. code-block:: python

   fig, ax = plt.subplots(figsize=(10, 10))
   tcl.pl.spatial_squarebin(adata, color="GeneA", ax=ax, show=False)

   tcl.pl.mark_region(
       ax,
       xlim=(3000, 5000),
       ylim=(2000, 4000),
       edges_color='cyan',
       edges_width=2.0,
       show=False           # defer display
   )
   plt.show()

Filled regions for better visibility against H&E backgrounds:

.. code-block:: python

   tcl.pl.mark_region(
       ax,
       xlim=(54500, 56000),
       ylim=(15000, 16000),
       fill_color='red',      # semi-transparent fill
       fill_alpha=0.15,        # fill opacity
       edges_width=3.0,
       show=False              # defer display
   )

Marking multiple regions efficiently (suppress intermediate refreshes):

.. code-block:: python

   tcl.pl.mark_region(ax, xlim=(40000, 42000), ylim=(5000, 7000),
                      edges_color='cyan', fill_color='cyan',
                      refresh=False, show=False)
   tcl.pl.mark_region(ax, xlim=(55000, 57000), ylim=(15000, 17000),
                      edges_color='yellow', fill_color='yellow',
                      refresh=False, show=False)
   tcl.pl.mark_region(ax, xlim=(60000, 62000), ylim=(10000, 12000),
                      edges_color='magenta', fill_color='magenta')
   # plt.show() is called by the last mark_region (show=True by default)

Key parameters:

* ``ax``: The matplotlib ``Axes`` to annotate (required).
* ``xlim``: ``(x_min, x_max)`` tuple for the region. If ``None``, uses the current x-axis limits.
* ``ylim``: ``(y_min, y_max)`` tuple for the region. If ``None``, uses the current y-axis limits.
* ``edges_color``: Rectangle edge color (default ``'red'``).
* ``edges_width``: Rectangle edge line width (default ``2.0``).
* ``fill_color``: Optional fill color with ``fill_alpha`` opacity (default ``None``).
* ``fill_alpha``: Fill opacity when ``fill_color`` is set (default ``0.15``).
* ``zorder``: Z-order for rendering (default ``100``). Ensures the rectangle is
  drawn on top of all other plot elements.
* ``refresh``: Whether to call ``ax.figure.canvas.draw_idle()`` after adding
  the rectangle (default ``True``). Set to ``False`` when adding multiple
  regions before a single refresh.
* ``show``: Whether to call ``plt.show()`` after adding the rectangle
  (default ``True``). Set to ``False`` when chaining multiple ``mark_region``
  calls or when you want to call ``plt.show()`` manually.

.. note::

   When the y-axis is inverted (``invert_y=True``, the default), ``ax.get_ylim()``
   returns ``(bottom, top)`` with ``bottom > top``.  ``mark_region`` normalizes
   these limits so the ``Rectangle`` always has positive dimensions, ensuring
   robust rendering across all matplotlib backends.

The function returns the ``matplotlib.patches.Rectangle`` object, which can be
further customized or removed later:

.. code-block:: python

   rect = tcl.pl.mark_region(ax, xlim=(1000, 2000), ylim=(3000, 4000))
   rect.set_linestyle('--')  # Make it dashed
   rect.set_linewidth(3.0)


Interactive ROI Selection (Jupyter / matplotlib)
--------------------------------------------------

TrackCell provides notebook-native interactive region-of-interest selection via
:func:`tcl.pl.select_regions`.  The backend uses **matplotlib widgets / ipympl**
instead of napari, avoiding Qt event-loop crashes in Jupyter.

.. important::

   In Jupyter Notebook / JupyterLab, run this **before** calling the selector:

   .. code-block:: python

      %matplotlib widget

   ``ipympl`` is installed with TrackCell.  If your environment lacks it, run
   ``pip install ipympl``.

Quick start
^^^^^^^^^^^

.. code-block:: python

   import trackcell as tcl

   selector = tcl.pl.select_regions(
       adata,
       color="CellType",
       key_added="ROI",
       inplace=True,           # write adata.obs["ROI"] after each ROI
   )

   # Draw ROIs interactively.  While the figure has focus:
   #   r = rectangle    e = ellipse    l = lasso
   #
   # After each ROI you'll be prompted for a name (press Enter for auto-name).
   #
   # Results are available immediately:
   selector.rois
   adata.obs["ROI"].value_counts(dropna=False)

Keyboard shortcuts
^^^^^^^^^^^^^^^^^^

While the figure has focus, press:

* ``r`` — rectangle selector
* ``e`` — ellipse / circle selector
* ``l`` — lasso / freehand selector

ROI naming
^^^^^^^^^^

After drawing each ROI you will see an ``input()`` prompt:

.. code-block:: text

   ROI name (press Enter for 'ROI_1'):

* Type a custom name (e.g. ``tumor_region``) and press Enter.
* Press Enter without typing to accept the auto-generated name (``ROI_1``,
  ``ROI_2``, …).

Returned object
^^^^^^^^^^^^^^^

:func:`tcl.pl.select_regions` returns a ``RegionSelector`` controller:

.. code-block:: python

   selector.rois        # dict: ROI name -> observation IDs
   selector.polygons    # dict: ROI name -> ROI vertices
   selector.save()      # write / rewrite adata.obs[key_added]
   selector.clear()     # remove all ROIs and visual patches
   selector.disconnect()

With ``inplace=False``, use ``selector.to_adata()`` to get an ``AnnData`` copy:

.. code-block:: python

   selector = tcl.pl.select_regions(adata, inplace=False)
   adata_with_rois = selector.to_adata()
   adata_with_rois.obs["ROI"].value_counts()

Squarebin data
^^^^^^^^^^^^^^

For squarebin data (e.g. Visium HD bin-level), set ``mode="squarebin"`` or let
TrackCell auto-detect it.  Bin centroids from ``adata.obsm[basis]`` are tested
against each ROI polygon.

.. code-block:: python

   selector = tcl.pl.select_regions(
       adata_bins,
       mode="squarebin",
       basis="spatial",
       color="gene_of_interest",
       key_added="ROI_bins",
   )

Cellbin data
^^^^^^^^^^^^

For cellbin data, cell geometries are used and cells are selected when their
geometry intersects the ROI polygon.

.. code-block:: python

   selector = tcl.pl.select_regions(
       adata_cell,
       mode="cellbin",
       color="CellType",
       key_added="ROI_cell",
   )

Key parameters for :func:`tcl.pl.select_regions`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``color`` — column in ``adata.obs`` or gene in ``adata.var_names`` to colour
  the spatial plot.
* ``mode`` — ``"auto"`` (default), ``"cellbin"``, or ``"squarebin"``.
* ``basis`` — key in ``adata.obsm`` for squarebin coordinates (default
  ``"spatial"``).
* ``key_added`` — column name in ``adata.obs`` for ROI labels.
* ``inplace`` — ``True`` (default) writes labels to ``adata.obs[key_added]``
  after each selection; ``False`` stores results only on the ``RegionSelector``
  — use ``selector.to_adata()`` to get a copy.
* ``invert_y`` — default ``True`` to match image/spatial coordinates where y
  increases from top to bottom.
* Extra keyword arguments are forwarded to ``spatial_cell`` or
  ``spatial_squarebin``.

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
   import trackcell as tcl
   
   # Define region of interest
   x_min, x_max = 1000, 5000
   y_min, y_max = 1000, 5000
   
   # Create mask for spatial coordinates
   spatial_coords = adata.obsm['spatial']
   mask = ((spatial_coords[:, 0] >= x_min) & (spatial_coords[:, 0] <= x_max) &
           (spatial_coords[:, 1] >= y_min) & (spatial_coords[:, 1] <= y_max))
   
   # Create subset
   adata_subset = adata[mask].copy()
   
   # IMPORTANT: Synchronize geometries after subsetting
   # This is required when data was loaded with read_hd_cellseg()
   tcl.io.sync_geometries_after_subset(adata_subset, sample="Cse1")
   
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
       figsize=(6, 6)
   )
   
   # Step 3: High-resolution view of specific region
   # (Use spatial cropping as shown in Strategy 2)

