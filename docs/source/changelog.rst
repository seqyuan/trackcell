Changelog


Version 0.3.39
--------------

* **Fix background image direction when ``invert_y=False``**: ``spatial_cell``
  and ``spatial_squarebin`` now respect ``invert_y=False`` for colored
  background image overlays by vertically flipping the image, matching the
  existing ``color=None`` behavior. Previously, the background image was
  rendered unflipped while axis limits were inverted, causing ROI views to
  appear inconsistent when compared to manually-cropped images.


Version 0.3.38
--------------

* **Fix ssDNA grayscale image → RGB** in ``io.read_sto``: STOmics ssDNA images
  are single-channel grayscale. ``sc.pl.spatial`` applies viridis colormap to
  2D arrays, producing misleading purple backgrounds. The reader now converts
  to (H, W, 3) RGB at load time.

* **``trackcell[clustering]`` extra**: ``leidenalg`` and ``python-igraph`` are
  now an optional extra (``pip install trackcell[clustering]``) rather than
  hard dependencies.


Version 0.3.37
--------------

* **Relax Python version constraint** from ``>=3.10,<3.12`` to ``>=3.10,<3.13``
  to allow installation on Python 3.12 environments.

* **Docs fixes** (:pr:`#—`):

  * STOmics notebook now has embedded figures (PCA, UMAP, spatial cell polygons
    with ssDNA background, ROI zoom with cell boundaries).
  * Gallery thumbnails for STOmics and YardCluster notebooks.
  * ``usage.rst`` updated to mention STOmics support.


Version 0.3.36
--------------

* **STOmics (Stereo-seq) GEF/GEM data support** (:pr:`#—`):

  **New IO module** (``io/read_sto.py``):

  * ``read_sto_cellbin()`` — reads STOmics cellbin GEF (HDF5) and creates an
    AnnData with cell polygon geometries (``cellBorder`` → Shapely polygons).
    Supports both file and folder path input with auto-discovery.
  * ``read_sto_bin()`` — reads STOmics squarebin/tissue GEF files for bin-level
    analysis. Supports automatic bin merging and GEM.gzip fallback with
    auto-redirect to matching GEF.
  * Auto-discovers ssDNA tissue images from ``03.register/`` directory
    (sibling or parent of GEF path).
  * Reads ``resolution`` attribute (500nm/unit) from GEF metadata and stores
    proper physical scalefactors (``pixel_size_um``, ``resolution_nm``,
    ``tissue_hires_scalef``).
  * Handles ``cellBorder`` sentinel stripping (32767 padding for cells with
    fewer than 32 vertices).

  **New example notebooks**:

  * ``STOmics_mouse_brain_demo`` — full STOmics cellbin read→QC→Leiden
    clustering→spatial visualization→subset→cell boundary zoom-in workflow.
  * ``YardCluster_STOmics_demo`` — YardCluster spatial clustering on STOmics
    mouse brain (auto mode: celltype + domain output).

  **Documentation updates**:

  * ``usage/reading.rst`` — added "Reading STOmics (Stereo-seq) Data" section.
  * ``examples.rst`` — linked both new notebooks.

* **Image auto-discovery for STOmics readers** — ``image_path`` parameter
  accepts file, register directory, or None (auto-discover from
  ``03.register/`` sibling of GEF).


Version 0.3.35
--------------

* **Defensive fixes from code review** (:pr:`#—`):

  **Spatial weights NaN / zero-division protection**
  (``tl._spatial_graph.build_spatial_weights``):

  * Handle duplicate spatial coordinates correctly (use neighbour index, not
    zero-distance, to exclude self).
  * Filter out non-finite distances before computing the Gaussian bandwidth.
  * Guard against zero or non-finite weight-sum with uniform fallback.
  * Validate input shapes (``coords.ndim == 2``, ``batch_labels`` length).

  **Annohdcell polygon generation robustness**
  (``io.convert_annohdcell.bins_to_cell_polygon``):

  * Colinear bins now produce a valid buffered Polygon instead of a
    LineString.
  * Exception fallback uses ``Point(…).buffer(…)`` instead of a degenerate
    3-identical-point polygon.

  **Input alignment checks in SpaceRanger readers**
  (``io.read_data.read_hd_bin`` / ``read_hd_cellseg``):

  * Reject duplicate barcodes or missing tissue positions early with clear
    error messages.
  * Warn when expression barcodes lack segmentation geometries, and vice
    versa.
  * Use ``reindex`` instead of ``.loc[]`` for safer alignment.
  * Protect against ``IndexError`` on empty classification / geometry
    columns.

  **YardCluster sparse-matrix scaling** (``tl.spatial_cluster``):

  * ``_scale_for_yardcluster`` now uses ``zero_center=False`` when
    ``adata.X`` is sparse, avoiding accidental densification / OOM.

  **DE-merge guard against expensive pair enumeration**
  (``tl._cluster_merge.merge_clusters_de``):

  * New ``max_de_tests`` parameter (default 2000). When the number of
    potential cluster-pair comparisons exceeds the cap, the merge is skipped
    with a warning instead of running a large number of Wilcoxon tests.
  * Parameter exposed through ``yard_cluster`` and ``spatial_cluster`` as
    ``merge_max_de_tests``.

  **Bug fix: ``integrate='separate'`` batch label corruption**
  (``tl.spatial_cluster``):

  * Fixed a bug where batch-prefixed labels were written as a single
    string containing the entire array instead of per-cell labels.

  **Classification color-parsing robustness**
  (``io.read_data.convert_classification_to_color_dict``):

  * Rewritten to tolerate mixed dict / JSON-string / list formats without
    relying on ``DataFrame.explode().unique()``.

  **``hd_labeldist`` explicit coordinate-system control**
  (``tl.spatial.hd_labeldist``):

  * New ``coordinate_system`` parameter: ``'fullres'``, ``'hires'``, or
    ``'auto'`` (default, preserves old heuristic but warns when ambiguous).
  * New ``library_id``, ``microns_per_pixel``, ``hires_scale`` overrides.
  * Resolved parameters stored in ``adata.uns`` / ``DataFrame.attrs`` for
    reproducibility.

  **Plotting geometry auto-sync after subset**
  (``pl.plot.spatial_cell`` / ``pl.select.select_regions``):

  * ``_sync_geometries_to_obs`` automatically filters and re-orders
    ``uns['spatial'][sample]['geometries']`` to match current
    ``adata.obs_names`` on first use after subsetting.
  * Simplified ``spatial_cell`` geometry-bounds validation: replaced nested
    retry loop with a single pass + ``aspect='equal'`` fallback.
  * ``spatial_cell`` and ``spatial_squarebin`` no longer call
    ``fig.tight_layout()`` on user-provided axes.
  * Added ``from __future__ import annotations`` so geopandas type hints
    don't break module import when geopandas is absent.

  * The following files were changed as part of this work:

    * ``trackcell/tl/_spatial_graph.py``
    * ``trackcell/tl/_cluster_merge.py``
    * ``trackcell/tl/spatial_cluster.py``
    * ``trackcell/tl/spatial.py``
    * ``trackcell/io/convert_annohdcell.py``
    * ``trackcell/io/read_data.py``
    * ``trackcell/pl/plot.py``
    * ``trackcell/pl/select.py``
    * ``tests/test_yardcluster.py``
    * ``tests/test_convert_annohdcell.py`` (new)
    * ``tests/test_read_data_helpers.py`` (new)
    * ``tests/test_spatial.py`` (new)
    * ``tests/test_plot_helpers.py`` (new)


Version 0.3.34
--------------

* **Documentation overhaul** — synchronized all usage tutorials with actual
  function signatures and parameters:

  * Fixed ``spatial_cell(..., sample=...)`` → ``library_id=...`` in all
    code examples (12 places across docstrings, RST, Markdown).
  * **New usage tutorial**: ``docs/source/usage/annohdcell_conversion.rst``
    — annohdcell (bin2cell) conversion, now integrated into Sphinx.
  * ``README.md``: added YardCluster, DBSCAN slice separation, interactive ROI
    selection, and mark_region sections.
  * ``spatial_clustering.rst``: documented missing advanced parameters
    (``preprocess``, ``layer``, ``use_raw``, ``hvg_by_batch``,
    ``scale_by_batch``, ``harmony_integrate``, ``copy``).
  * ``computing_distances.rst``: added ``hd_labeldist`` ``method`` parameter
    (kdtree vs cdist) and colony centroid distance section
    (``mark_colony_centroids``, ``distance_to_nearest_centroids``).
  * ``visualization.rst``: fixed garbled key-parameters section for
    ``select_regions``, added missing ``alpha_facet``/``show`` for facet mode,
    ``library_id``/``figsize`` for ROI select.
  * ``examples.rst``: replaced duplicate toctree with ``:doc:`` links.
  * ``annohdcell_conversion.rst``: added ``create_polygons`` and ``n_jobs``
    method-specific parameter docs.


Version 0.3.33
--------------

* **DBSCAN spatial slice separation** (``tcl.tl.spatial_slice_cluster``):

  * Split Xenium TMA regions into physical tissue sections (``S001``–``SNNN``)
    via density clustering on cell centroids (CCHD defaults: ``eps=80`` µm,
    ``min_samples=10``, ``min_cells=1000``).
  * ``read_xenium_cellseg(..., slice_separate=True)`` for one-step load + annotate.
  * Helpers: ``split_by_slice``, ``slice_cluster_summary``, ``write_slice_annotation``.

* **Micro-region colony clustering** (GC / tumor nests):

  * ``spatial_colony_cluster`` — DBSCAN on BANKSY-filtered subsets
    (defaults: ``eps=50``, ``min_samples=5``, ``min_cluster_size=40``).
  * ``mark_colony_centroids`` and ``distance_to_nearest_centroids`` for GCC/TNC
    distance analysis (CCHD ``54_1_GC_dist.ipynb`` workflow).
  * New tutorial: ``docs/source/usage/xenium_slice_separation.rst``.



Version 0.3.32
--------------

* **YardCluster spatial clustering** (``tcl.tl.spatial_cluster``):

  * Lightweight CPU pipeline combining BANKSY-style neighborhood features,
    dual-channel identity/context embeddings, and scanpy Leiden clustering.
  * Function group: ``spatial_neighbors``, ``neighborhood_features``,
    ``yard_embed``, ``yard_cluster``.
  * Multi-sample support: ``integrate='separate'`` / ``'joint'`` (Harmony).
  * Optional DE-guided cluster merging, sketch mode for large datasets, and
    ``gradient_mode='cosphi'|'regression'``.
  * New usage tutorial: ``docs/source/usage/spatial_clustering.rst``; Colon Cancer
    demo notebook updated.



Version 0.3.31
--------------

* **Fix selectors not responding in ipympl**:

  * Changed all selectors from ``useblit=True`` to ``useblit=False``
    (blitting breaks RectangleSelector/EllipseSelector/LassoSelector in the
    ipympl backend).
  * Added ``fig.canvas.draw()`` after selector creation to ensure widgets
    initialize properly.



Version 0.3.30
--------------

* **Hotfix: fix ``_mode_buttons`` AttributeError** — ``_build_toolbar()`` is
  now called before ``_connect_all_selectors()`` so that
  ``_update_toolbar_highlight()`` has access to the toolbar button dict.



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

