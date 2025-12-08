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

