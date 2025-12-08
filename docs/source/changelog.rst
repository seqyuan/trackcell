Changelog
=========

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

