Reading SpaceRanger Output
============================

This section covers how to read spatial transcriptomics data from SpaceRanger output.

Reading Cell Segmentation Data
-------------------------------

TrackCell can read 10X HD SpaceRanger cell segmentation output.

Required Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function expects the SpaceRanger output directory to have the following structure:

::

   segmented_outputs/
   ├── graphclust_annotated_cell_segmentations.geojson  # Cell segmentation polygons (default)
   │   OR
   ├── cell_segmentations.geojson                       # Alternative filename (auto-detected)
   ├── filtered_feature_cell_matrix.h5                  # Expression matrix
   └── spatial/
       ├── tissue_hires_image.png                        # High-resolution tissue image
       ├── tissue_lowres_image.png                       # Low-resolution tissue image
       └── scalefactors_json.json                        # Image scaling factors

**Note**: The function will automatically try alternative filenames if the default 
``graphclust_annotated_cell_segmentations.geojson`` is not found. Supported alternatives include:
``cell_segmentations.geojson``, ``cell_segmentations_annotated.geojson``, and 
``annotated_cell_segmentations.geojson``. You can also explicitly specify the filename using 
the ``cell_segmentations_file`` parameter.

Usage
~~~~~

.. code-block:: python

   import trackcell as tcl

   # Read SpaceRanger cell segmentation output
   # The function will auto-detect the segmentation file if default name is not found
   adata = tcl.io.read_hd_cellseg(
       datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
       sample="Cse1",
       # cell_segmentations_file is optional, default is "graphclust_annotated_cell_segmentations.geojson"
   )
   
   # Or explicitly specify the segmentation file name
   adata = tcl.io.read_hd_cellseg(
       datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
       sample="Cse1",
       cell_segmentations_file="cell_segmentations.geojson"  # Custom filename
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
--------------------------------------

TrackCell can read 10X HD SpaceRanger bin-level output for 2um, 8um, or 16um bins.

Required Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function expects the SpaceRanger output directory to have the following structure:

::

   square_016um/
   ├── filtered_feature_bc_matrix.h5          # Expression matrix (H5 format, preferred)
   │   OR
   ├── filtered_feature_bc_matrix/            # Expression matrix (directory format)
   │   ├── barcodes.tsv.gz
   │   ├── features.tsv.gz
   │   └── matrix.mtx.gz
   └── spatial/
       ├── tissue_positions.parquet           # Spatial coordinates (parquet format, preferred)
       │   OR
       ├── tissue_positions.csv               # Spatial coordinates (CSV format)
       ├── tissue_hires_image.png             # High-resolution tissue image
       ├── tissue_lowres_image.png            # Low-resolution tissue image
       └── scalefactors_json.json             # Image scaling factors

Usage
~~~~~

.. code-block:: python

   import trackcell as tcl

   # Read SpaceRanger bin-level output (2um/8um/16um bins)
   adata = tcl.io.read_hd_bin(
       datapath="SpaceRanger4.0/Cse1/binned_outputs/square_016um",
       sample="Cse1",
       binsize=16  # Bin size in micrometers (default: 16, common values: 2, 8, or 16)
   )

   # Access the bin size information
   print(f"Bin size: {adata.uns['spatial']['Cse1']['binsize']} um")

Subsetting Data and Synchronizing Geometries
---------------------------------------------

When you subset an AnnData object loaded with ``read_hd_cellseg()``, the cell geometries 
stored in ``adata.uns["spatial"][sample]["geometries"]`` are **not automatically updated**. 
This can cause errors when plotting subsetted data.

**Important**: Always call ``sync_geometries_after_subset()`` after subsetting data loaded 
with ``read_hd_cellseg()`` to ensure geometries are synchronized.

Usage
~~~~~

.. code-block:: python

   import trackcell as tcl
   import numpy as np
   
   # Read data
   adata = tcl.io.read_hd_cellseg(
       datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
       sample="Cse1"
   )
   
   # Method 1: Subset by spatial region
   x_min, x_max = 16000, 18000
   y_min, y_max = 14000, 18000
   
   spatial_coords = adata.obsm['spatial']
   mask = ((spatial_coords[:, 0] >= x_min) & (spatial_coords[:, 0] <= x_max) &
           (spatial_coords[:, 1] >= y_min) & (spatial_coords[:, 1] <= y_max))
   
   adata_subset = adata[mask].copy()
   
   # IMPORTANT: Synchronize geometries after subsetting
   tcl.io.sync_geometries_after_subset(adata_subset, sample="Cse1")
   
   # Now you can safely plot the subset
   tcl.pl.spatial_cell(adata_subset, color="classification")
   
   # Method 2: Subset by cell metadata
   adata_subset2 = adata[adata.obs['classification'] == 'Cluster-1'].copy()
   
   # IMPORTANT: Synchronize geometries after subsetting
   tcl.io.sync_geometries_after_subset(adata_subset2, sample="Cse1")
   
   # Now you can safely plot the subset
   tcl.pl.spatial_cell(adata_subset2, color="classification")

What Gets Synchronized
~~~~~~~~~~~~~~~~~~~~~~

The ``sync_geometries_after_subset()`` function:

* Filters ``adata.uns["spatial"][sample]["geometries"]`` (GeoDataFrame) to only include 
  cells present in the subsetted ``adata.obs_names``
* Ensures the geometries match the subsetted data

**Note**: ``adata.obs["geometry"]`` (WKT strings) is automatically subset when you subset 
the AnnData object, so it doesn't need manual synchronization. However, the plotting function 
prefers using the GeoDataFrame format for better performance.

Why This Is Necessary
~~~~~~~~~~~~~~~~~~~~~

When you subset an AnnData object:

* ``adata.obs`` and ``adata.obsm`` are automatically subset (they are indexed by cell IDs)
* ``adata.uns["spatial"][sample]["geometries"]`` is **NOT** automatically subset (it's a 
  separate GeoDataFrame object)

If you try to plot without synchronizing, the plotting function may:

* Fail with errors like ``ValueError: aspect must be finite and positive``
* Attempt to access geometries for cells that no longer exist in the subset
* Produce incorrect visualizations

Always call ``sync_geometries_after_subset()`` after subsetting to avoid these issues.

