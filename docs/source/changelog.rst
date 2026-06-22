Changelog

Version 0.3.29
--------------

* **Fix notebook ROI selector UX — no more blocking ``input()``**:

  * ``tcl.pl.select_regions`` no longer calls ``input()`` after each ROI.
    ROIs are now **auto-named** (``ROI_1``, ``ROI_2``, …) with a configurable
    ``roi_prefix`` parameter — the interactive workflow is never interrupted.
  * **Inline toolbar** added at the bottom of the figure with clickable
    buttons: ``■ Rect``, ``● Ellipse``, ``✎ Lasso``, ``✕ Clear``, ``↩ Undo``.
    The active mode is highlighted in blue.
  * **Selected cells are now highlighted** in real-time with coloured scatter
    points (cycling through 10 distinct colours).
  * New method ``selector.rename_roi(old, new)`` for post-hoc renaming.
  * New method ``selector.undo()`` to remove the last ROI (boundary, highlight,
    and data).
  * Keyboard shortcuts ``r``/``e``/``l`` kept for power users.



Version 0.3.28
--------------

* **Replace napari ROI selection with Jupyter-native matplotlib selector**:

  * ``tcl.pl.select_regions`` keeps the original function name but now uses
    matplotlib widgets / ipympl instead of napari/Qt.
  * Keyboard-toggle between all shape modes: ``r``=rectangle, ``e``=ellipse,
    ``l``=lasso — no ``shape_type`` parameter needed.
  * Interactive ROI naming via ``input()`` after each selection
    (auto-generates ``ROI_1``, ``ROI_2``, … on empty input).
  * ``inplace`` parameter replaces ``copy`` (scanpy convention):
    ``inplace=True`` (default) writes to ``adata.obs``;
    ``inplace=False`` → use ``selector.to_adata()``.
  * Supports both cellbin geometry intersection and squarebin centroid selection.
  * Returns a ``RegionSelector`` controller with ``rois``, ``polygons``,
    ``save()``, ``clear()``, ``to_adata()``, and ``disconnect()``.
  * Removes napari from public plotting API and package dependencies.



Version 0.3.27
--------------

* **Fix Jupyter kernel crash in select_regions**: auto-detect Jupyter
  environment and enable Qt GUI integration (%gui qt) before calling
  napari.run() to prevent Qt/IPython event loop conflict.



Version 0.3.26
--------------

* **napari ROI selection now supports squarebin data**:

  * New ``mode`` parameter (``"auto"`` / ``"cellbin"`` / ``"squarebin"``).
  * New ``basis`` parameter for squarebin coordinate key in ``adata.obsm``.
  * Auto-detection: cellbin if geometries found, squarebin otherwise.
  * Squarebin extraction uses point-in-polygon (``Point.within``) on bin centroids.



Version 0.3.25
--------------

* **Fix PyPI auto-publish**: update pyproject.toml version to match __init__.py
  so poetry build produces correct package version.
* **Fix workflow**: username → user for pypa/gh-action-pypi-publish >= v1.12.


=========

Version 0.3.24
--------------

* **Added napari-based interactive ROI selection**:

  * New module ``trackcell.pl.napari`` with three public functions.
  * Supports rectangle, polygon, freehand, and any shape.
  * ``copy`` parameter for in-place vs. dict return.
  * napari is an optional extra (``pip install 'trackcell[napari]'``).

Version 0.3.23
--------------

* **Added ``show`` parameter to ``mark_region``**:

  * ``show=True`` (default) now calls ``plt.show()`` automatically after adding
    the rectangle, so ``mark_region`` works out of the box without requiring
    ``show=False`` in the plotting function and a manual ``plt.show()``.

  * Set ``show=False`` to defer display when chaining multiple ``mark_region``
    calls or when calling ``plt.show()`` manually.

  * **Motivation**: The previous workflow required ``show=False``,
    ``ax=ax``, and ``plt.show()`` — a pattern many users overlooked, causing
    the rectangle to appear invisible.  Now the simplest usage works:

    .. code-block:: python

       ax = tcl.pl.spatial_cell(adata, color="CellType")
       tcl.pl.mark_region(ax, xlim=(54500, 56000), ylim=(15000, 16000))

* **Updated documentation**:

  * Simplified ``mark_region`` examples in the user guide to showcase the new
    default-behavior workflow
  * Documented the new ``show`` parameter in both the docstring and the
    visualization guide

Version 0.3.22
--------------

* **Rewrote ``mark_region`` for robustness**:
  * Fixed invisible rectangle bug on inverted y-axis: normalized negative
    width/height from ``get_ylim()`` returning ``(bottom, top)`` with ``bottom > top``
  * Added ``zorder=100`` default to ensure the rectangle always renders on top of cells/bins
  * Added ``fill_color`` and ``fill_alpha`` parameters for optional semi-transparent fill
  * Added ``refresh=True`` parameter — calls ``canvas.draw_idle()`` after adding the patch
    (set ``refresh=False`` for batch marking then manually refresh at end)
  * Increased ``edges_width`` default from 1.0 to 2.0 for better visibility

* **Updated documentation**:
  * Added ``.. important::`` admonition about ``ax=ax, show=False`` pattern for ``mark_region``
  * Added examples for ``spatial_cell``, ``spatial_squarebin``, filled regions, and multi-region batching
  * Documented all new parameters: ``fill_color``, ``fill_alpha``, ``zorder``, ``refresh``
  * Noted y-axis normalization behavior with ``invert_y``

Version 0.3.21
--------------

* **Fixed ``invert_y`` not working when ``color=None``**:
  * Root cause: ``set_ylim(h, 0)`` already produces image convention with ``origin='upper'``;
    the subsequent ``invert_yaxis()`` call double-inverted the y-axis, making both
    ``invert_y=True`` and ``invert_y=False`` visually identical.
  * Fix: removed the redundant ``invert_yaxis()`` from both ``spatial_cell`` and
    ``spatial_squarebin`` color-None code paths.
  * For ``invert_y=False`` (Cartesian), the H&E image is now pre-flipped via ``img[::-1]``
    before ``imshow()`` to keep it right-side-up.

* **Cleaned up ``mark_region``**:
  * Removed unused ``Circle`` import inside the function body.
  * Added comprehensive documentation to the user guide with usage examples
    for both ``spatial_cell`` and ``spatial_squarebin``.
  * Documented coordinate convention behavior with ``invert_y``.

Version 0.3.20
--------------

* **Added ``shape`` parameter to ``spatial_squarebin``**:
  * ``shape='circle'`` (default): renders bins as centered circles (radius = bin_size / 2)
  * ``shape='square'``: renders bins as filled rectangles (original behavior)
  * Invalid shape values raise ``ValueError``
  * Also aliased in ``tcl.pl.spatial_bin()``

* **Added ``invert_y`` parameter to ``spatial_cell`` and ``spatial_squarebin``**:
  * ``invert_y=True`` (default): y-axis increases top-to-bottom (image coordinates, matches H&E)
  * ``invert_y=False``: y-axis increases bottom-to-top (Cartesian convention)

* **Documentation updates**:
  * Documented ``edge_color`` boundary overlap behavior in ``spatial_cell``:
    adjacent shared boundaries are drawn by category; the last-drawn category determines the visible color

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

