Atera
=====

TrackCell reads 10x Atera bundles in both flavours:

* **official** (Atera 1.x): ``cell_feature_matrix.zarr.zip`` +
  ``cells.zarr.zip`` + ``transcripts.zarr.zip`` + ``experiment.spatial``
* **preview**: Xenium-v4 compatible ``cell_feature_matrix.h5`` +
  ``cells.parquet`` + ``cell_boundaries.parquet`` + ``experiment.xenium``

The reader auto-detects the format.  The resulting ``AnnData`` has the same
schema for both formats, so downstream code works unchanged.

.. seealso::

   :doc:`../notebooks/Atera_demo` — end-to-end notebook for an Atera kidney
   ROI: official Zarr, polygon visualization, cell QC, Q-score, transcript
   density, nuclear fraction, single-molecule point cloud and a preview bundle.

First look at a bundle
----------------------

Atera bundles can be hundreds of GB.  Inspect one before loading anything:

.. code-block:: python

   import trackcell as tcl

   tcl.io.atera_info("/path/to/ROI_bundle")

This prints cells × features, Gene Expression feature count, polygon
coverage, transcript count / grid tiles and the software version, without
loading the expression matrix.

Reading
-------

.. code-block:: python

   # auto-detect preview vs official
   adata = tcl.io.read_atera("/path/to/ROI_bundle")

   # explicit entry points
   adata = tcl.io.read_atera_cellseg("/path/to/ROI_bundle")
   adata = tcl.io.read_atera_preview("/path/to/preview_bundle")

Common options:

.. code-block:: python

   adata = tcl.io.read_atera(
       "/path/to/ROI_bundle",
       sample="C001_ROI1",       # default: experiment region_name
       gene_name_index=True,     # var_names = gene symbols
       build_polygons=True,      # ~1 min for 400k cells
       include_nucleus=None,     # default: same as build_polygons
       store_wkt=False,          # WKT costs memory; add later on demand
       keep_de=False,            # drop hundreds of official de_* columns
       load_panel=False,         # panel_config.json is ~11 MB (customer data)
       region=(xmin, ymin, xmax, ymax),   # µm, or a shapely geometry
       n_cells=None,             # random subset
   )

.. note::
   ``region`` is ``(xmin, ymin, xmax, ymax)``, matching shapely / GeoPandas.
   Passing ``(xmin, xmax, ymin, ymax)`` raises an error instead of silently
   selecting the wrong area.

What the AnnData contains
-------------------------

``.X``
   Raw counts (CSR, cells × features).  **Not** normalised or log-transformed.

``.obs``
   ``barcode`` (string, same as ``obs_names``),
   ``centroid_column`` / ``centroid_row`` (µm),
   ``cell_area`` / ``nucleus_area`` (µm²),
   ``transcript_counts`` (gene-only), ``total_counts`` (**includes control
   probes/codewords**), official cluster / Azimuth columns.
   Official data also has ``cell_id`` (uint32), ``cell_prefix`` and
   ``dataset_suffix``.

``.obsm['spatial']``
   ``(x, y)`` in **micrometres**.  Atera centroids are already µm; do not
   multiply them by ``pixel_size`` (that is for morphology images only).

``.var``
   ``gene_ids`` (Ensembl), ``feature_types`` (``Gene Expression`` /
   ``Negative Control Codeword`` / ...) and ``is_gene_expression``.

``.uns['atera']``
   ``coord_unit='um'``, ``roi`` / ``slide`` / ``run``, ``pixel_size``,
   ``z_step_size``, ``chemistry_version``, ``analysis_version``,
   ``panel_path``, ``transcripts_path``.

``.uns['spatial'][sample]['geometries']``
   Cell polygons (GeoDataFrame).  ``nucleus_geometries`` is added when
   ``include_nucleus`` is enabled (the default follows ``build_polygons``).

``.uns['cell_boundaries']`` / ``['nucleus_boundaries']``
   Compact vertex arrays: ``cell_idx`` maps each ``obs`` row to an entry in
   ``n_vertices`` (``-1`` = no polygon); ``vertex_x`` / ``vertex_y`` are
   concatenated in that order.  Note that this is a per-cell layout and
   differs from the offset-based layout written by
   :func:`trackcell.io.read_xenium_cellseg`.

Filter control features before analysis
---------------------------------------

Atera panels contain thousands of negative-control codewords.  They must be
removed before HVG / PCA / clustering:

.. code-block:: python

   adata = tcl.io.filter_atera_genes(adata)
   print(adata.var["feature_types"].value_counts())

Cell polygons
-------------

.. code-block:: python

   tcl.pl.spatial_cell(adata, color="total_counts")

Cells without a polygon (a small fraction) are skipped with a warning.
Use :func:`trackcell.io.sync_geometries_after_subset` after subsetting.

If you need WKT strings for serialization:

.. code-block:: python

   tcl.io.add_wkt_geometry(adata)          # obs['geometry']
   tcl.io.add_wkt_geometry(adata, nucleus=True)

Transcript point cloud
----------------------

``transcripts.zarr.zip`` contains single molecules and is **not** loaded into
``adata.X``.  A kidney ROI has ~1.19 billion molecules, so always stream:

.. code-block:: python

   for chunk in tcl.io.iter_atera_transcripts(
           roi, genes=["NPHS1", "UMOD"], min_qv=20,
           region=(0, 0, 2000, 2000)):
       print(chunk.head())

Columns: ``x``, ``y``, ``z`` (µm), ``qv``, ``cell_index`` (row in
``adata.obs``, ``-1`` = unassigned), ``cell_prefix``, ``gene_index``,
``gene`` (categorical), ``codeword``, ``overlaps_nucleus``.

Pass ``adata=adata`` when the AnnData was subset, so ``cell_index`` refers to
the subset's ``obs`` order:

.. code-block:: python

   sub = tcl.io.read_atera(roi, region=(0, 0, 2000, 2000))
   df = tcl.io.read_atera_transcripts(
       roi, genes=["NPHS1"], adata=sub
   )

* ``min_qv=20`` matches the official cell matrix.  The point cloud itself
  also contains lower-quality molecules — pass ``min_qv=None`` to keep them.
* ``genes=`` uses the per-tile ``gene_offset`` index, so only the requested
  gene ranges are read.  Targeting one gene on a full kidney ROI takes about
  a second instead of ~5 minutes.
* ``read_atera_transcripts(...)`` returns a single DataFrame for small
  queries.  It refuses to run without ``genes`` / ``region`` /
  ``max_rows`` (use ``force=True`` to override).

QC and plots
------------

.. code-block:: python

   # cell-level QC (from obs / X only)
   qc = tcl.io.atera_cell_qc(adata)

   # transcript QC (streams the point cloud; ~5 min for 1.19B molecules)
   summary = tcl.io.atera_transcript_qc(roi, adata=adata, min_qv=20)

   tcl.pl.atera_qc_summary(adata)            # 2×2 cell QC panel
   tcl.pl.atera_transcript_density(adata)    # stored density map
   tcl.pl.atera_nuclear_fraction(
       adata, genes=["UMOD", "DDX17"]
   )

``atera_transcript_qc`` writes back to ``adata``:

* ``obs['tx_q20_counts']``, ``obs['tx_nuclear_counts']``,
  ``obs['nuclear_fraction']`` — per-cell Q≥20 transcript count, nuclear
  count and fraction.
* ``uns['atera']['transcript_qc']``, ``['qscore']``, ``['density']``,
  ``['gene_qc']`` — summary, histogram, density map and per-gene nuclear
  fractions.

Because all results live in ``adata``, the plotting functions can be called
repeatedly (they never consume a one-shot iterator).

Multi-ROI and multi-customer delivery
-------------------------------------

Atera has **no customer barcodes**; the delivery unit is an ROI.  A run
directory contains one sub-directory per ROI:

.. code-block:: python

   rois = tcl.io.list_atera_rois("/path/to/run")
   print(rois[["roi", "n_cells", "has_transcripts"]])

   run = "/path/to/run"
   for name, adata in tcl.io.iter_atera_rois(run, rois=["C001_ROI1"]):
       ...

Merging the same customer's ROIs requires explicit offsets in µm:

.. code-block:: python

   merged = tcl.io.concat_atera_rois(
       {"C001_ROI1": a, "C001_ROI2": b},
       offsets={"C001_ROI1": (0, 0), "C001_ROI2": (12000, 0)},
       index_unique="-",
   )

There is no reliable way to derive offsets from
``overview_scan_fov_locations.json`` (5–10 µm error), so the function refuses
to guess.  ``index_unique`` suffixes barcodes per ROI (Atera barcodes can
repeat across ROIs) and keeps geometry indices consistent.

Supported analyses and limits
-----------------------------

Supported:

* gene-level expression, cell polygons, cell QC
* Q-score, transcript density, nuclear localisation
  (``overlaps_nucleus``)
* spatial pseudotime / trajectory (Monocle3, Slingshot, CellRank
  expression/connectivity kernels) — **not** RNA velocity
* spatial masking to split shared ROIs downstream

Not supported (Atera does not provide the data):

* isoform / transcript-level quantification (the WTA panel is gene-level)
* splice junctions / alternative splicing
* spliced / unspliced counts and RNA velocity
  (``overlaps_nucleus`` is **not** unspliced)
* customer demultiplexing (there are no customer barcodes)

Performance notes
-----------------

Measured on the public kidney WTA dataset (399,120 cells,
1.19B transcripts):

* ``read_atera(..., build_polygons=False)``: ~3 s
* ``build_polygons=True``: ~60 s (398k cell + 398k nucleus polygons)
* ``atera_transcript_qc`` full scan: ~320 s, ~400 MB RSS
* ``iter_atera_transcripts(genes=["NPHS1"])``: ~1 s (gene-range fast path)
* reading a 717k-cell preview bundle: ~70 s

The transcript reader skips background grid tiles (negative ``grid`` keys)
that contain only placeholder rows outside the ROI.

Panel configuration
-------------------

``panel_config.json`` is ~11 MB and may contain customer / designer
information, so it is not loaded by default.  Only ``uns['atera']['panel_path']``
is recorded.  Load it explicitly when needed:

.. code-block:: python

   panel = tcl.io.read_atera_panel("/path/to/ROI_bundle")
   adata = tcl.io.read_atera(roi, load_panel=True)   # uns['gene_panel']
