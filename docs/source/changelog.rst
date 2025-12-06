Changelog
=========

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

