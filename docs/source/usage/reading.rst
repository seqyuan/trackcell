Reading SpaceRanger Output
============================

This section covers how to read spatial transcriptomics data from SpaceRanger output.

Reading Xenium Output
---------------------

TrackCell can read 10x Xenium Analyzer cell segmentation output.

Required Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function expects the Xenium output directory to have the following structure:

::

   output-XETG00070__0085811__yuanfa__20260529__052421/
   ├── cell_feature_matrix.h5      # Expression matrix (10x HDF5, CSC format)
   ├── cells.parquet               # Cell metadata (centroids, area, counts)
   ├── cell_boundaries.parquet     # Cell boundary vertices (long-table)
   ├── nucleus_boundaries.parquet  # Nucleus boundary vertices (optional)
   ├── experiment.xenium           # Experiment metadata (JSON)
   └── gene_panel.json             # Gene panel information (JSON)

**Key format differences from Visium HD**:

* **Cell boundaries**: Stored as long-table parquet (one row per vertex) instead of GeoJSON.
  Automatically converted to Shapely polygons.
* **Expression matrix**: CSC (genes × cells) format, automatically transposed to CSR (cells × genes).
* **Cell IDs**: Include the ``-1`` suffix natively (no stripping needed).

Usage
~~~~~

.. code-block:: python

   import trackcell as tcl

   # Read Xenium Analyzer output
   adata = tcl.io.read_xenium_cellseg(
       datapath="/path/to/xenium/output",
       sample="sample1"
   )

The resulting AnnData object contains:

* Expression matrix in ``.X`` (CSR, cells × genes)
* Cell metadata in ``.obs`` (centroids, area, counts, segmentation_method, etc.)
* Gene metadata in ``.var`` (gene_ids, feature_types)
* Spatial coordinates in ``.obsm["spatial"]``
* Cell polygons in ``.uns["spatial"][sample]["geometries"]`` (GeoDataFrame)
* Cell boundary arrays in ``.uns["cell_boundaries"]`` (compact vertex arrays)
* Nucleus boundary arrays in ``.uns["nucleus_boundaries"]`` (if available)
* WKT geometry strings in ``.obs["geometry"]`` (for serialization)
* Experiment metadata in ``.uns["experiment"]``

After loading, the AnnData is fully compatible with ``tcl.pl.spatial_cell()`` for
cell polygon visualization.

TMA multi-core slice separation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Xenium TMA data (multiple tissue cores in one region), enable DBSCAN slice
separation during load or as a separate step. See :doc:`xenium_slice_separation`
for the full workflow.

.. code-block:: python

   adata = tcl.io.read_xenium_cellseg(
       datapath="/path/to/xenium/region",
       sample="85811_S",
       slice_separate=True,
       slice_eps=80,
   )
   print(adata.obs["slice_id"].value_counts())


Reading Cell Segmentation Data (Visium HD)
-------------------------------------------

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

   # Show H&E image only with coordinate range
   tcl.pl.spatial_squarebin(adata, color=None)

   # Equivalent alias
   tcl.pl.spatial_bin(adata, color=None)

   # Plot a gene on square bins
   tcl.pl.spatial_squarebin(
       adata,
       color="EPCAM",
       cmap="Reds",
       alpha=0.8,
       alpha_img=0.4,
   )


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

Next Steps
~~~~~~~~~~

After loading cell segmentation data, you can run spatial clustering with
:doc:`spatial_clustering` (``tcl.tl.spatial_cluster``) or visualize
Space Ranger's pre-computed ``classification`` labels directly with
:doc:`visualization`.

Always call ``sync_geometries_after_subset()`` after subsetting to avoid these issues.


Reading STOmics (Stereo-seq) Data
---------------------------------

TrackCell supports reading STOmics GEF format data for both cellbin (cell
segmentation) and squarebin (tissue bin) analysis.

.. _sto_data_structure:

Data Structure
~~~~~~~~~~~~~~

A typical STOmics pipeline output looks like:

::

   C57_7/
   ├── 03.register/
   │   ├── ssDNA_B02825B2_regist.tif      # 注册后的 ssDNA 组织底图
   │   ├── B02825B2.rpi                   # 金字塔图像 (HDF5)
   │   └── B02825B2.reregist.ipr          # 配准参数 (HDF5)
   ├── 04.tissuecut/
   │   ├── B02825B2.tissue.gef            # squarebin 表达矩阵 (HDF5)
   │   └── B02825B2.tissue.gem.gz         # squarebin 表达矩阵 (文本)
   └── 041.cellcut/
       ├── B02825B2.cellbin.gef           # cellbin 表达矩阵 + 细胞轮廓 (HDF5)
       └── B02825B2.cellbin.gem           # cellbin 表达矩阵 (文本)

**物理分辨率**: 1 GEF 坐标单位 = 500nm = 0.5μm。图像像素与 GEF 坐标 1:1 对齐。

Read Cellbin Data
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   # 支持文件夹路径（自动发现）或文件路径
   adata = tcl.io.read_sto_cellbin(
       "path/to/C57_7/041.cellcut/",   # 或直接传 ".cellbin.gef" 文件
       sample="C57_7"
   )

返回的 AnnData 包含:

* ``.X``: CSR 表达矩阵 (cells × genes)
* ``.obs``: cell_id, area, dnb_count, gene_count, exp_count, geometry (WKT)
* ``.obsm['spatial']``: 细胞质心坐标
* ``.obs['geometry']``: WKT 多边形字符串
* ``.uns['spatial'][sample]['geometries']``: GeoDataFrame 细胞多边形
* ``.uns['spatial'][sample]['images']``: 自动加载的 ssDNA 底图（如存在）
* ``.uns['spatial'][sample]['metadata']``: resolution_nm, pixel_size_um, coordinate_range

Read Squarebin Data
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # squarebin 数据用于 bin 级别的组织分析
   adata = tcl.io.read_sto_bin(
       "path/to/C57_7/04.tissuecut/",   # 或直接传 ".tissue.gef" 文件
       bin_size=50,                      # 合并 50 个 500nm 单位 = 25μm per bin
       sample="C57_7"
   )

See the full workflow in the :doc:`STOmics notebook </notebooks/STOmics_mouse_brain_demo>`.

.. note::

   GEF 文件本身 **不包含** ssDNA 组织底图。底图存储在 ``03.register/`` 目录中。
   ``read_sto_cellbin()`` 和 ``read_sto_bin()`` 会自动发现并加载底图。

Using the ssDNA Image with scanpy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

STOmics ssDNA images are single-channel grayscale.  When loaded by
``read_sto_cellbin()`` / ``read_sto_bin()`` they are automatically converted
to 3-channel (H, W, 3) RGB so that ``sc.pl.spatial`` and
``tcl.pl.spatial_cell`` display them correctly.

.. code-block:: python

   # Load with image — automatic grayscale→RGB conversion
   adata = tcl.io.read_sto_cellbin(
       "path/to/cellbin.gef",
       image_path="path/to/regist.tif",
       sample="sample1"
   )

   # Direct plot — no purple background
   sc.pl.spatial(adata, color="leiden", spot_size=10, alpha_img=0.8)

   # trackcell polygon plot with ssDNA background
   tcl.pl.spatial_cell(adata, color="leiden", alpha_img=0.5)

Background visibility tips
   * ``alpha_img=0`` — hide the image entirely (white background)
   * ``alpha_img=0.3`` — dimmed background, cells pop out
   * ``alpha_img=0.8`` — darker image, better for gene expression overlay


Converting annohdcell Output
-----------------------------

TrackCell provides two functions for converting annohdcell (bin2cell) output
to trackcell-compatible format with polygon geometries:

* :func:`~trackcell.io.convert_annohdcell_to_trackcell` — create a new
  cell-level h5ad from a 2μm bin h5ad (simple count summation).
* :func:`~trackcell.io.add_geometries_to_annohdcell_output` — add polygon
  geometries to annohdcell's final cell h5ad (preserves exact aggregation).

After reading an h5ad that contains WKT geometry strings, restore the
GeoDataFrame for visualization with :func:`~trackcell.io.restore_geometries`.

See :doc:`annohdcell_conversion` for the full tutorial with examples,
parameter details, and troubleshooting.

