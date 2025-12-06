Usage Guide
===========

This guide provides examples of how to use TrackCell for reading and analyzing spatial transcriptomics data.

Reading SpaceRanger Output
---------------------------

Reading Cell Segmentation Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TrackCell can read 10X HD SpaceRanger cell segmentation output:

.. code-block:: python

   import trackcell as tcl

   # Read SpaceRanger cell segmentation output
   adata = tcl.io.read_hd_cellseg(
       datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
       sample="Cse1"
   )

The resulting AnnData object contains:

* Expression matrix in ``.X``
* Cell metadata in ``.obs``
* Gene metadata in ``.var``
* Spatial coordinates in ``.obsm["spatial"]``
* Tissue images in ``.uns["spatial"][sample]["images"]``
* Scalefactors in ``.uns["spatial"][sample]["scalefactors"]``
* Cell geometries in ``.uns["spatial"][sample]["geometries"]`` (GeoDataFrame)
* Cell geometries in ``.obs["geometry"]`` (WKT strings for serialization)

Reading Bin-Level Data (2um/8um/16um)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For bin-level data:

.. code-block:: python

   import trackcell as tcl

   # Read SpaceRanger bin-level output (2um/8um/16um bins)
   adata = tcl.io.read_hd_bin(
       datapath="SpaceRanger4.0/Cse1/binned_outputs",
       sample="Cse1",
       binsize=16  # Bin size in micrometers (default: 16, common values: 2, 8, or 16)
   )

   # Access the bin size information
   print(f"Bin size: {adata.uns['spatial']['Cse1']['binsize']} um")

Visualization
-------------

Plotting with Cell Polygons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TrackCell provides a specialized plotting function that visualizes cells as polygons
instead of points, providing a more accurate representation of cell boundaries:

.. code-block:: python

   # Plot cells as polygons (requires data loaded with read_hd_cellseg)
   tcl.pl.spatial_cell(
       adata, 
       color="classification",  # Color by cell type
       groups=['Cluster-2', 'Cluster-3'],  # Optional: filter specific groups
       figsize=(10, 10),
       edges_width=0.5,
       edges_color="black",
       alpha=0.8
   )

   # Plot continuous values (e.g., distance to a label)
   tcl.pl.spatial_cell(
       adata,
       color="Cluster-2_dist",  # Distance to Cluster-2
       cmap="Reds",
       figsize=(10, 10)
   )

   # Plot multiple variables in subplots
   tcl.pl.spatial_cell(
       adata,
       color=["classification", "Cluster-2_dist"],  # Two subplots
       figsize=(20, 10)
   )

Computing Distances to a Label
-------------------------------

TrackCell can compute the distance from every cell to the nearest cell annotated with
a specific label:

.. code-block:: python

   # Compute distance to a specific annotation label
   tcl.tl.hd_labeldist(
       adata,
       groupby="classification",    # obs column containing cell type annotations
       label="Cluster-2",           # target label to measure distances from
       inplace=True                 # add "{label}_px" and "{label}_dist" to adata.obs
   )

   # When inplace=False, returns a DataFrame
   dist_df = tcl.tl.hd_labeldist(
       adata, 
       groupby="classification", 
       label="Cluster-2", 
       inplace=False
   )

   # Visualize the distance
   tcl.pl.spatial_cell(adata, color='Cluster-2_dist', cmap='Reds', figsize=(10, 10))

The function automatically detects whether coordinates are in hires or full-res resolution
by comparing the computed tissue size with the expected chip size (6.5mm for 10X HD),
ensuring compatibility with both SpaceRanger output and bin2cell-processed data.

