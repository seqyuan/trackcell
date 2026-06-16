Changelog
=========

Version 0.3.19
--------------

* **Added dedicated HD squarebin plotting**:
  * Added ``tcl.pl.spatial_squarebin()`` for squarebin/bin-level HD visualization
  * Added ``tcl.pl.spatial_bin()`` alias for a shorter, intuitive API
  * Supports plotting obs columns, genes, categorical values, and continuous values
  * Supports ``crop_coord`` for coordinate-based regional zooming

* **Improved no-color spatial rendering for HD plots**:
  * When ``color=None``, squarebin plots now display the H&E image and spatial coordinate extent
  * Added axis/tick control for coordinate-aware inspection, matching the expected HD browsing workflow
  * Designed to complement ``tcl.pl.spatial_cell()`` for users working with non-cellseg squarebin data

* **Documentation and examples updated**:
  * Added dedicated squarebin usage examples to README and Sphinx docs
  * Added example script for squarebin visualization workflows
  * Reduced Sphinx autodoc build noise via mocked heavy imports for more stable online docs builds

Version 0.3.18
--------------

* **Added ``read_xenium_cellseg`` function**:
  * Reads 10x Xenium Analyzer output (cell_feature_matrix.h5, cells.parquet, cell_boundaries.parquet)
  * Automatically transposes CSC (genes × cells) to CSR (cells × genes)
  * Converts long-table boundary parquet to Shapely polygons + compact vertex arrays in a single pass
  * Reads optional ``nucleus_boundaries.parquet``, ``experiment.xenium``, and ``gene_panel.json``
  * Stores GeoDataFrame geometries (compatible with ``spatial_cell()``) and WKT strings for serialization
  * Squidpy-compatible ``uns['spatial']`` structure

* **Dual-color visualization in ``spatial_cell``**:
  * Added ``edge_color`` parameter to color cell boundaries by a categorical column
  * Added ``edge_palette`` parameter for custom edge color mapping
  * Fill (``color``) + Edge (``edge_color``) two-pass rendering for rich spatial views
  * Automatic dual legend: colorbar for fill, categorical legend for edges

* **Multi-gene visualization**:
  * ``tl.multigene_blend()`` with two modes:
    * ``mode='blend'``: Weighted RGB blending, cell2location-style composite
    * ``mode='facet'``: Faceted subplots with single-hue colormaps per gene
  * Raw hex color auto-detection in ``spatial_cell`` for blend mode output
  * Supports up to 7 genes with customizable colors, percentile clipping, and gamma

Version 0.2.8
-------------

* **Major refactoring of ``spatial_cell`` function**:
  * Refactored to use GeoDataFrame.plot() for efficient rendering and automatic legend generation
  * Integrated scanpy-style background image processing with scale_factor support
  * Fixed coordinate alignment issue when plotting subset data (image extent now matches data range)
  * Added support for gene expression visualization (color parameter can be a gene name from adata.var_names)
  * Legend and colorbar now positioned outside the plot area (right side) similar to sc.pl.spatial
  * Added ``xlabel``, ``ylabel``, and ``show_ticks`` parameters for axis customization
  * ``library_id`` parameter now auto-selects first available library_id when not specified (similar to sc.pl.spatial)
  * Added ``alpha_img`` parameter to control background image transparency
  * Improved categorical legend generation with custom palette support

* **Enhanced ``read_hd_cellseg`` function**:
  * Added automatic detection of alternative segmentation file names (e.g., ``cell_segmentations.geojson``)
  * Automatically adds ``spot_diameter_fullres`` to scalefactors if missing (required for sc.pl.spatial compatibility)
  * Improved error messages with suggestions for alternative filenames

* **Documentation improvements**:
  * Reorganized usage documentation into separate files (reading.rst, visualization.rst, computing_distances.rst)
  * Added directory structure requirements to reading documentation
  * Added comprehensive visualization optimization guide for large datasets
  * Updated examples with gene expression visualization

Version 0.2.5
-------------

* Fixed ``read_hd_cellseg`` bug where ``cellid`` column was lost after ``reset_index()``
* Added robust handling for cases where ``reset_index()`` creates ``'index'`` column instead of ``'cellid'``
* Improved error messages for better debugging

Version 0.2.4
-------------

* Fixed ``read_hd_cellseg`` bug where ``cellid`` column was missing when creating GeoDataFrame
* Added field name detection for ``cell_id`` column (supports multiple naming variants)

Version 0.2.3
-------------

* Fixed ``read_hd_cellseg`` to properly store cell geometries in GeoDataFrame format

Version 0.2.2
-------------

* Added ``spatial_cell`` plotting function for visualizing cells as polygons
* Modified ``read_hd_cellseg`` to store cell geometries in both GeoDataFrame and WKT format
* Renamed ``trackcell/io/spatial.py`` to ``trackcell/io/read_data.py``
* Renamed ``trackcell/pl/spatial.py`` to ``trackcell/pl/plot.py``
* Improved documentation and examples

Version 0.2.1
-------------

* Added ``read_hd_bin`` function for reading bin-level data (2um/8um/16um)
* Enhanced ``hd_labeldist`` with automatic coordinate resolution detection
* Improved memory efficiency with cKDTree method

Version 0.2.0
-------------

* Initial release with core functionality
* ``read_hd_cellseg`` for reading cell segmentation data
* ``hd_labeldist`` for computing distances to labels

