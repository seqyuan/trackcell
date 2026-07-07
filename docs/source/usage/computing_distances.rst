Computing Distances to a Label
================================

TrackCell can compute the distance from every cell to the nearest cell annotated with
a specific label. This is useful for analyzing spatial relationships between different
cell types or regions.

Basic Usage
-----------

.. code-block:: python

   # Compute distance to a specific annotation label
   tcl.tl.hd_labeldist(
       adata,
       groupby="classification",    # obs column containing cell type annotations
       label="Cluster-2",           # target label to measure distances from
       inplace=True                 # add "{label}_px" and "{label}_dist" to adata.obs
   )

When ``inplace=True``, the function adds two columns to ``adata.obs``:

* ``{label}_px``: Distance in pixels (on the hires/registered image)
* ``{label}_dist``: Physical distance in microns

The ``method`` parameter controls the distance computation backend:

* ``"kdtree"`` (default): Uses ``scipy.spatial.cKDTree`` for efficient O(n log m)
  queries on large datasets.
* ``"cdist"``: Uses ``scipy.spatial.distance.cdist`` for fast all-pairs
  computation on **small** datasets (warning: O(n×m) memory).

Returning Results as DataFrame
-------------------------------

When ``inplace=False``, the function returns a DataFrame with the distance columns:

.. code-block:: python

   # When inplace=False, returns a DataFrame
   dist_df = tcl.tl.hd_labeldist(
       adata, 
       groupby="classification", 
       label="Cluster-2", 
       inplace=False
   )
   
   # The DataFrame contains two columns:
   # - Cluster-2_px: pixel distance
   # - Cluster-2_dist: physical distance in microns

Visualizing Distances
---------------------

After computing distances, you can visualize them using ``tcl.pl.spatial_cell()``:

.. code-block:: python

   # Visualize the distance
   tcl.pl.spatial_cell(adata, color='Cluster-2_dist', cmap='Reds', figsize=(10, 10))

Coordinate Resolution Detection
-------------------------------

The function automatically detects whether coordinates are in hires or full-res resolution
by comparing the computed tissue size with the expected chip size (6.5mm for 10X HD),
ensuring compatibility with both SpaceRanger output and bin2cell-processed data.

* **Standard SpaceRanger output**: When ``hires_scale`` is 0.05-0.2, correctly identifies as hires coordinates
* **bin2cell processed data**: When ``hires_scale`` is close to 1.0, correctly identifies as full-res coordinates
* **Edge cases**: When between 0.9-0.99, uses size comparison to determine


Distance to Colony Centroids (GC / Tumor Nests)
================================================

After spatial colony clustering (see :doc:`xenium_slice_separation`), you can
identify colony centroids and compute per-cell distances to the nearest centroid.
This workflow is used in CCHD for GC B cell / tumor nest distance analysis.

Mark Colony Centroids
---------------------

:func:`~trackcell.tl.mark_colony_centroids` finds the cell nearest to the geometric
centroid of each colony and labels it (e.g. ``"GCC"`` for germinal center centroids):

.. code-block:: python

   import trackcell as tcl

   # After spatial_colony_cluster on a filtered subset
   tcl.tl.spatial_colony_cluster(gc, eps=50, min_samples=5, min_cluster_size=40)
   tcl.tl.mark_colony_centroids(
       gc,
       colony_key="spatial_colony",
       centroid_key="cell_centroid_type",
       centroid_label="GCC",         # label for centroid cells
   )

This adds:

* ``obs['cell_centroid_type']`` — ``"GCC"`` for centroid cells, ``NaN`` for others
* ``uns['colony_centroids']`` — list of dicts with colony ID, size, and centroid coordinates

Distance to Nearest Centroid
----------------------------

:func:`~trackcell.tl.distance_to_nearest_centroids` uses a **cKDTree** to compute
the distance from every cell in the full section to the nearest marked centroid:

.. code-block:: python

   # Compute distance on the full slice, not just the filtered subset
   tcl.tl.distance_to_nearest_centroids(
       adata,                           # full-section AnnData
       centroid_key="cell_centroid_type",
       distance_key="distance_to_nearest_gcc",
       centroid_label="GCC",
   )

The resulting ``obs['distance_to_nearest_gcc']`` is in the same units as
``obsm['spatial']`` (µm for Xenium centroids, pixel coords for Visium HD).

Visualize Centroid Distances
-----------------------------

.. code-block:: python

   # Distance heatmap on the full tissue section
   tcl.pl.spatial_cell(adata, color="distance_to_nearest_gcc", cmap="Reds")

   # Colony identification
   tcl.pl.spatial_cell(gc, color="spatial_colony")


Multi-Gene Color Blending
==========================

TrackCell can compute blended colors for multi-gene co-expression visualization.
See the full documentation in :doc:`visualization`.

.. code-block:: python

   import trackcell as tcl

   # Blend mode: single composite image
   tcl.tl.multigene_blend(adata, genes=['EPCAM', 'PECAM1', 'VWF'], mode='blend')
   tcl.pl.spatial_cell(adata, color='multigene_blend')

   # Facet mode: subplots per gene
   tcl.tl.multigene_blend(
       adata,
       genes=['EPCAM', 'PECAM1', 'VWF', 'ACTA2'],
       mode='facet', ncols=2,
   )

