Changelog


Version 0.4.2
-------------

**Upgrade note**: ``vst(do_regularize=False)`` now returns real residuals
where 0.4.1 returned NaN, and ``integrate_pair_xgb`` applies a larger (and
correct) correction. Re-run affected analyses rather than mixing versions.

**Fixed**:

* ``tl._sctransform.vst(do_regularize=False)`` ignored R's rule that ``n_genes``
  is dropped when regularisation is off (``sctransform::vst``: "will use all
  genes"). Step 1 stayed subsampled, the extrapolation step never ran, and
  ``model_pars_fit.reindex(genes)`` left the un-subsampled genes NaN — 80% of
  the residual matrix in a 500-gene/``n_genes=100`` run. Step 1 now uses every
  gene, and any gene still without a model raises instead of yielding NaN.
* ``tl.integration._xgboost_integration.integrate_pair_xgb`` computed
  ``delta = y_pred @ pca.components_``, which omits ``pca.mean_``. The mean *is*
  the bulk of a batch shift, so the correction was systematically too small by
  exactly that term. Now uses ``pca.inverse_transform``.
* ``tl.slice_separation.dbscan_slice_labels`` failed silently in three ways:
  non-finite coordinates surfaced as sklearn's opaque "Input contains NaN";
  every cluster below ``min_cells`` marked 100% of cells ``debris`` without a
  word; and a pixel-coordinate slide given a µm ``eps`` collapsed to one slice.
  Each now raises or warns with the actual cause and the relevant numbers,
  including a units hint when ``eps`` is more than 25x the median
  nearest-neighbour spacing of a single surviving cluster.
* ``tl._gene_utils.build_gene_alias_lookup`` let the first gene to claim a
  ``.N``/``-N`` alias keep it, so with both ``A.1`` and ``A-1`` present the real
  ``A-1`` resolved to ``A.1`` — a different gene, silently. Exact names now win
  over aliases, and ``align_gene_list`` deduplicates with a set instead of an
  O(n^2) membership scan.
* ``tl._nhood.diff_nhood_enrichment(pvalues=True)`` produced anti-conservative
  p-values: a category absent from one condition became NaN after reindexing,
  and ``NaN >= x`` is always False, so the count never incremented for those
  cells. Values now use the ``(count+1)/(n_perms+1)`` estimator (p can no longer
  be 0) and a Benjamini-Hochberg matrix is returned as ``pvalue_adj``.
* ``io.read_sto_*`` left duplicate GEF gene symbols in ``var_names``, which makes
  ``get_loc`` return a *slice*: ``X[mask, slice]`` silently becomes an (n, k)
  matrix and plotting fails somewhere unrelated. Duplicates are renamed with a
  numeric suffix, warned about, and the original symbols kept in
  ``var['gene_symbol']``.
* ``io.read_sto_bin`` packed bin ids as ``(x << 32) | y`` after casting to
  ``uint64``. GEF coordinates are relative to an origin that is not always the
  array corner, so a negative bin index wrapped to ~4.29e9 and produced
  nonsense bin centres. The packing now shifts to a non-negative origin and
  unshifts when decoding, so negative coordinates stay exact.
* ``pl.spatial_cell`` classified every integer ``obs`` column as continuous, so
  integer-coded cluster labels (``leiden`` as int8, cluster ids as int64) were
  drawn as a viridis gradient with no legend. Columns are now classified by
  cardinality (``Categorical``/object/bool always, integers with at most 30
  levels), and a new ``categorical=`` argument overrides the inference either
  way.
* ``pl.select`` produced arbitrary ROI membership for a self-intersecting lasso
  (tracing a region and retracing a diagonal is a normal mouse artefact): an
  invalid polygon's ``contains`` results are undefined, and one test case
  selected only one of the two cells the user had drawn around. The ring is now
  repaired with ``make_valid`` before testing membership.
* ``pl.select.RegionSelector.save`` re-initialised the whole ``obs`` column to
  ``None`` before writing the current ROIs, so adding one ROI to an object
  restored from a previous h5ad discarded every earlier label. Only cells inside
  a ROI recorded in this session are overwritten now.
* ``io.read_xenium_cellseg`` only printed a warning when the expression matrix
  and ``cells.parquet`` disagreed, and ``verbose=False`` silenced it — a
  partial-matrix run therefore "succeeded" with far fewer cells than expected.
  Zero overlap raises with both id lists; below 90% raises; the remainder logs
  at ``warning`` through the logger.
* ``io.read_sto_cellbin`` divided by ``img.max()`` when normalising the tissue
  image, emitting ``invalid value encountered in divide`` on a blank capture.
  It also stored ``hires`` and ``lowres`` as the *same array object*, so the
  h5ad carried the full image twice and mutating one changed the other.
  Blank captures render black, and ``lowres`` is now a genuinely downsampled
  copy.
* The wheel no longer bundles ``trackcell/benchmark`` — 2592 lines of
  development scaffolding (including two ``.R`` scripts) with zero imports from
  library code. It was 20 of 64 wheel entries and 12.5% of the package, and it
  skewed every coverage number. ``scripts/check_package.py`` now fails if it
  reappears. Tests keep importing it from the checkout.
* Optional-dependency upper bounds were guesses that the maintainer's own
  environment already violated: ``python-igraph<1.0`` against an installed
  1.0.0, ``harmonypy<1.0`` against 2.1.0. CI therefore installed versions the
  maintainer never ran. Bounds widened to the majors actually in use.
* ``trackcell/py.typed`` contained ``partial``, which has meaning only for
  stub-only distributions; mypy only checks that the file exists. Emptied, which
  is the correct marker for an inline-typed package.
* Coverage ratchets added for ``io/serialize.py`` (>= 70%, currently 87%) and
  ``io/convert_annohdcell.py`` (>= 12%, currently 12% — a floor that can only
  move up). The previously ratcheted files were ``read_data`` and ``read_sto``
  only, so the third IO module with a README tutorial was unguarded.
* ``io.read_atera`` with ``region=`` / ``n_cells=`` (the preview path) subset the
  object *after* building polygons and left ``uns['cell_boundaries']['cell_idx']``
  at its pre-subset length, so every vertex offset addressed the wrong cell —
  silently, since the geometry table itself was trimmed correctly.
* ``io.read_sto_*`` resolved an ambiguous GEF directory by taking
  ``candidates[0]`` in ``glob`` order, which is filesystem-dependent: the same
  directory could yield a different chip on a different machine. Candidates are
  sorted by name and the warning names all of them.
* ``tl._autores`` / ``tl.spatial_cluster`` fell back from the igraph Leiden
  backend on ``(ImportError, TypeError, ValueError)``. A parameter error was
  therefore retried on a second backend instead of surfacing, and the fallback
  ran with leidenalg's default ``n_iterations`` (-1 vs 2) — so a resolution
  sweep could compare clusters from two different algorithms and report the
  difference as stability. Only ``ImportError`` triggers a fallback now, both
  branches use ``n_iterations=2``, and the switch warns.

**Changed**:

* Development tooling is now pinned to the exact versions used locally, in all
  three places that run it. ``black`` was pinned ``>=26.0,<27.0`` in
  ``pyproject.toml`` but ``24.10.0`` in ``.pre-commit-config.yaml``; ``flake8``
  was ``^6.0.0`` in one place and ``7.1.1`` in the other, while the maintainer
  ran ``7.4.1``; ``pytest`` was ``^7.0.0`` against a local ``9.0.2``. Two
  formatter versions in one repository make "pre-commit reformats what CI just
  accepted" an expected outcome rather than a risk. ``pyproject.toml``,
  ``.pre-commit-config.yaml`` and CONTRIBUTING now name the same versions.


Version 0.4.1
-------------

**Upgrade note**: every ``integrate_rpca`` / ``integrate_sct_rpca`` /
``integrate_pair`` result changes. Anchor weighting was inverted (see *Fixed*),
so integrated embeddings, UMAPs and cluster labels produced by 0.4.0 are not
comparable to 0.4.1. Re-run integration rather than mixing outputs from the two
versions.

**Fixed**:

* ``tl.integration`` anchor weighting was **inverted**. Seurat's
  ``IntegrateData`` converts kNN distances to similarities before calling
  ``FindWeightsC`` (Seurat 4.4.0, ``R/integration.R:3820``:
  ``distances <- 1 - (distances / distances[, ncol(x = distances)])``);
  ``integrate_data.integrate_pair`` passed raw Annoy distances into the same
  formula (``src/integration.cpp:49``), so the *nearest* anchor received the
  *smallest* weight. Every ``integrate_rpca`` / ``integrate_sct_rpca`` /
  ``integrate_pair`` run was affected, and no test could catch it because the
  only integration test compared two copies of the same formula. The
  conversion now lives in ``_distance_to_similarity`` and is covered by
  ``tests/test_integrate_data.py::test_weights_decrease_with_distance``.
* ``integrate_pair`` now emits a ``warning`` when fewer query cells carry
  anchors than ``k_weight``, instead of silently narrowing the weight kernel.
* ``tl._sctransform.get_residuals`` labelled residual rows with the *requested*
  gene order but filled them in ``var_names`` order, so whenever
  ``genes != umi_genes`` — i.e. any ``sctransform(residual_features=[...])``
  call with a user-ordered marker panel — every residual was attached to the
  wrong gene, silently. Rows are now selected by position from ``umi_genes``;
  duplicate gene names raise instead of mismatching, and
  ``_xgboost_integration`` indexes residuals by label rather than position.
* ``io.concat_atera_rois`` stored per-ROI offsets as a dict of ``tuple`` in
  ``uns['atera']['offsets']``. anndata has no serializer for ``tuple``, so the
  merged object could never be written with ``save_h5ad``/``write_h5ad``.
  Offsets are now ``float64`` arrays of length 2.
* ``io.serialize`` geometry round trips: a missing polygon no longer aborts the
  write (``None`` in an object column raised ``TypeError``; it is now stored as
  an empty string and read back as ``None``), every non-geometry column is
  preserved (attribute columns and the active column order used to be dropped),
  and ``GeoDataFrame.crs`` survives via a ``<key>__crs`` sidecar. The WKT codec
  now lives in one place — ``io.restore_geometries`` delegates to
  ``serialize.restore_geometries`` instead of keeping a second, lossier copy.

**Changed**:

* ``tests/conftest.py`` skips every ``@pytest.mark.slow`` test when the
  benchmark data root is absent. A clean checkout previously produced ten
  ``FileNotFoundError`` failures for ``pytest -m slow`` despite README and
  CONTRIBUTING both promising skips.
* CI now runs the *same* ``flake8 trackcell tests scripts`` command documented
  in CONTRIBUTING (previously it ran two narrow ``--select`` subsets, so the
  documented command was permanently red with 122 findings). ``.flake8``
  ignores ``E501`` (black owns wrapping) and grants ``E402`` only to the three
  files where the matplotlib backend must be configured before importing
  pyplot.
* ``scripts/check_package.py`` gains ``--expect-version``. The release workflow
  now passes the triggering git tag, which is the only value that can disagree
  with ``pyproject.toml`` — the previous checks compared artifacts built from
  the same source and therefore could never fail.


Version 0.4.0
-------------

**Behavior change**:

* Library code no longer prints progress to stdout. Readers and
  ``multigene_blend`` emit ``INFO`` records under the ``trackcell`` logger
  (silent by default; opt in via :mod:`logging`). Non-fatal data problems
  (missing tissue image or ``scalefactors_json.json``) now raise
  ``RuntimeWarning`` via :func:`warnings.warn`. Intentionally displayed report
  output (``atera_info``, the interactive ROI selector) still prints.  See the
  *Logging* usage page.

**Added**:

* :func:`trackcell._logging.get_logger` — structured logging bootstrap with a
  ``NullHandler`` so ``import trackcell`` stays silent.
* IO test coverage for ``read_hd_bin`` / ``read_hd_cellseg`` / ``read_sto``
  driven by synthetic but structurally faithful vendor fixtures
  (``tests/conftest.py``). CI enforces coverage minimums
  (``read_data`` ≥ 75%, ``read_sto`` ≥ 60%).
* ``scripts/check_package.py`` now validates the built wheel's ``Version``
  against the source version, closing the stale-wheel false positive.

**Changed**:

* ``pl.plot`` and ``pl.select`` no longer silently swallow exceptions;
  optional rendering steps log at ``debug`` and a visible draw failure logs a
  ``warning``.
* ``read_data``'s H5-fallback notice is logged at ``INFO``.

**Fixed**:

* ``read_sto_cellbin`` now reads scalar ``offsetX``/``offsetY`` GEF attributes
  correctly (a parsing bug previously forced them to ``0``) and records them in
  ``uns["spatial"][sample]["metadata"]`` as ``offset_x``/``offset_y``/``offset_um``.
  Coordinates remain in the GEF/image frame (offsets are not added).
* Removed dead code in ``pl.plot.spatial_cell`` (``edges_draw``; dual-color
  edges are already drawn by ``_draw_categorical_edges``) and in
  ``_sctransform_v2.build_poisson_gene_tables`` (unused Poisson-set locals);
  no behavior change.

**Hygiene**:

* Stop tracking ``.annovibe/im.json`` (contained bot credentials); ignore
  ``.annovibe/``.
* Moved benchmark/debug helper scripts out of ``tests/`` into ``dev/``; renamed
  ``tests/sct_benchmark_data.py`` → ``tests/_sct_benchmark_data.py``.
* Removed root-level ``test_import.py``; moved ``promotion_article.md`` to
  ``docs/marketing/``.


Version 0.3.47
--------------

* **Bug fixes**:

  * ``sctransform(copy=True)`` now returns the modified copy (previously it
    silently returned ``None`` and discarded the result).
  * ``store_residuals_in='layer'`` residuals are now aligned to
    ``adata.var_names`` instead of being written in residual-variance order.
    The batch path now honours ``store_residuals_in`` as well.
  * ``read_xenium_cellseg`` no longer raises ``UnboundLocalError`` when
    ``build_polygons=False`` or when ``cell_boundaries.parquet`` is missing but
    ``nucleus_boundaries.parquet`` exists.
  * New :func:`trackcell.io.save_h5ad` / :func:`trackcell.io.flatten_geometries`
    serialize polygon GeoDataFrames to WKT (``adata.write_h5ad`` fails
    otherwise); ``restore_geometries`` now also handles
    ``nucleus_geometries``.
  * ``split_by_slice`` now reindexes ``uns['cell_boundaries']`` using the
    original observation positions (previously every slice received the first
    N cells' vertices).
  * ``spatial_cluster`` validates ``mode`` / ``method`` / ``integrate`` /
    ``graph_mode`` (invalid values used to silently do nothing), and its
    preprocessing no longer requires ``scikit-misc``: ``seurat_v3`` is used on
    raw counts when available, otherwise it falls back to ``seurat`` on
    log-normalized data.  ``layer``/``use_raw`` are part of the neighborhood
    cache key.
  * ``prep_sct_integration`` now defaults to ``key_added='sct_prep'``, matching
    ``integrate_rpca`` / ``run_sct_integration``.
  * ``score_anchors_seurat`` no longer mixes batch1 and batch2 index spaces
    (shared-neighbour counts were inflated by spurious cross-dataset matches).
  * ``multigene_blend`` supports more than seven genes without an explicit
    palette, honours ``layer``/percentile clipping in facet mode, and no
    longer writes ``uns`` when ``inplace=False``.
  * ``RegionSelector.save(key_added=...)`` writes the requested column.
  * ``read_hd_bin``/``read_hd_cellseg`` infer sample names correctly from
    SpaceRanger layouts; ``read_hd_cellseg`` handles list-valued
    classifications.
  * ``hd_labeldist`` falls back to ``metadata['pixel_size']`` / Xenium
    ``experiment['pixel_size']`` when no scalefactor is available.
  * STOmics: circular cell diameter from area, ``gene_name_index`` honoured for
    GEM input, and ``image_path`` supported for GEM input.
  * ``diff_nhood_enrichment`` restricts the permutation null to the two
    compared conditions.

* **Tests**: regression coverage for the fixes above (SCT copy/layer order,
  h5ad geometry round-trip, Xenium boundaries, slice boundary alignment,
  anchor scoring, multigene blends, ROI saving, preprocessing fallback).


Version 0.3.46
--------------

* Pin ``pyarrow`` to ``<22``: pyarrow 22+ no longer works with the
  NumPy 1.x ABI used by trackcell, which broke ``pd.read_parquet`` for
  Xenium / SpaceRanger tables on Python 3.11–3.12.
* Skip the RPCA cache-release test when the optional ``annoy`` package is not
  installed (it is covered by the CI integration job).


Version 0.3.45
--------------

* **Packaging fixes**:

  * ``trackcell[clustering]`` now actually installs ``leidenalg`` and
    ``python-igraph`` (the extra was previously empty and ``poetry check``
    failed).
  * New ``integration`` extra (``annoy``, ``harmonypy``, ``bbknn``) and
    ``interactive`` extra (``ipympl``); both were previously either mandatory
    or documented as manual installs.
  * Direct dependencies (``anndata``, ``scipy``, ``scikit-learn``,
    ``statsmodels``, ``patsy``, ``joblib``, ``matplotlib``, ``pyarrow``) are now
    declared explicitly instead of relying on scanpy transitives.
    ``pyarrow`` is required to read Xenium / SpaceRanger parquet tables.
  * Added the MIT license (``LICENSE`` + package metadata).
  * Removed the machine-specific ``/Volumes/process/tmp/tcl_test`` benchmark
    path; use ``TRACKCELL_BENCHMARK_ROOT`` (default
    ``~/trackcell_benchmark``).

* **Testing**: added a CI workflow (Python 3.10–3.12) with lint, wheel
  metadata validation, and the fast test suite. Heavy R/Seurat parity tests are
  marked ``slow`` and require external benchmark data.

* **Fixes**:

  * The numpy fallback RNG used when the optional ``rrng`` is not installed now
    normalizes sampling probabilities like R's ``sample(prob=)``, fixing
    ``ValueError: probabilities do not sum to 1`` in SCT v2.
  * Integration tests that require the optional ``annoy`` package are skipped
    when it is not installed.

* **Performance**: ``_row_gmean_sparse`` is now vectorized (~30x faster on
  realistic sparse matrices, bit-identical results). Removed ineffective
  object-mode numba caching/``nogil`` from SCT residuals and GLM helpers,
  eliminating import-time ``NumbaWarning`` noise.


Version 0.3.44
--------------

* **Atera support** (``trackcell.io.read_atera``):

  * ``read_atera`` auto-detects official (Atera 1.x Zarr) and preview
    (Xenium-v4 compatible) bundles; ``read_atera_cellseg`` /
    ``read_atera_preview`` are the explicit entry points.
  * ``atera_info`` reports cells/genes/polygons/transcripts without loading
    the matrix.
  * Unified schema across formats: ``var_names`` are gene symbols,
    ``var['gene_ids']`` / ``var['feature_types']`` / ``var['is_gene_expression']``,
    ``obs['barcode']``, ``obsm['spatial']`` in µm, polygons in
    ``uns['spatial'][sample]``.
  * ``filter_atera_genes`` removes negative-control codewords before
    HVG/PCA.
  * ``iter_atera_transcripts`` / ``read_atera_transcripts`` stream the
    single-molecule point cloud with ``genes`` / ``region`` / ``min_qv``
    filters; the gene-range fast path reads only the requested genes.
  * ``atera_cell_qc`` / ``atera_transcript_qc`` (Q-score, assignment rate,
    nuclear fraction, per-gene nuclear fraction, density).
  * ``list_atera_rois`` / ``iter_atera_rois`` / ``concat_atera_rois`` for
    run directories; offsets are required unless ROIs are already aligned.
  * ``pl.atera_qc_summary``, ``pl.atera_transcript_density``,
    ``pl.atera_nuclear_fraction`` read everything from ``adata``.
  * ``read_xenium_cellseg`` gained ``build_polygons``, ``include_nucleus``,
    ``store_wkt`` and ``verbose`` (backward compatible).
  * Documented limits: no isoform, splice-junction, spliced/unspliced or
    RNA-velocity support; no customer barcodes.
  * New docs: :doc:`usage/atera`; new tests: ``tests/test_atera.py``.

* **YardCluster CellCharter-inspired extensions**:

  * Multi-hop neighborhood aggregation (``neighborhood_features(n_layers=...)``,
    ``spatial_cluster(n_layers=...)``).
  * Automatic Leiden resolution selection by clustering stability
    (``resolution="auto"``, ``tcl.tl.sweep_resolutions``,
    ``tcl.tl.select_resolution``).
  * Gaussian-mixture clustering backend (``yard_cluster(method="gmm")``,
    BIC-based ``n_clusters`` selection).
  * Post-analysis tools: ``tcl.tl.connected_components``,
    ``tcl.tl.nhood_enrichment``, ``tcl.tl.diff_nhood_enrichment``, and
    ``tcl.tl.remove_long_links``.
  * Enrichment heatmaps: ``tcl.pl.nhood_enrichment``,
    ``tcl.pl.diff_nhood_enrichment`` (diverging colormap, p-value stars).

* **SCT runtime backend cleanup**:

  * Remove the R subprocess bridge and the public ``method='r'`` /
    ``r_vst_export_*`` paths. ``resolve_sct_integration_method`` is gone with
    ``sct_method``; the backend is selected explicitly via ``sct_v2_backend``.
  * Add explicit ``v2_backend='python'|'pyglmGamPoi'`` selection. Pure Python is
    the deterministic default; ``pyglmGamPoi`` is optional and never selected
    automatically.
  * Add a vectorized ``MASS::theta.ml``-compatible NB offset fallback and separate
    the ``sqrt(N)`` HVG-ranking clip from the tighter ``sqrt(N/30)`` stored assay
    clip.
  * Move readers and R export scripts into ``trackcell.benchmark`` so offline
    Seurat parity checks remain available without a runtime R dependency.

* **SCT + RPCA benchmark on project AA/A04** (3 × 10x samples, 27,129 cells):

  * Pure Python ``v2_backend='python'`` vs Seurat R v4 ``SCTransform(vst.flavor='v2')``.
  * SCT HVG Jaccard median **0.939**, IntegrationFeatures Jaccard **0.933**,
    PrepSCTIntegration correlation **0.988**, anchor count ratio **0.989**.
  * Python SCT 1.83× faster (175s vs 320s); RPCA slower due to Numba weight construction.
  * See :doc:`benchmark/rpca_evaluation` for the full tables and
    ``<benchmark-root>/SCTRPCAbench/report.json`` for the raw results.

Version 0.3.43
--------------

* **SCT + RPCA integration** (Seurat v4 parity first):

  * Split API: ``run_sct_integration``, ``integrate_rpca``, ``sct_prep_matrix``;
    ``integrate_sct_rpca`` remains a thin wrapper.
  * Hybrid R step-1 / ``inject_r_scale_data`` path for prep and anchor parity
    (benchmark integrate corr **~0.96** vs Seurat).
  * Explicit opt-in memory release: ``release_sct_integration_cache`` and
    documented manual ``uns`` one-liners (no automatic cache trimming).

* **Performance** (RPCA path):

  * Vectorised / chunked IntegrateData weight construction (Numba).
  * ``n_trees_weight`` (default 10; use 50 for Seurat parity).
  * ``integration_dtype='float32'`` default; ``float64`` for parity.
  * ``sct_batch_n_jobs`` for parallel per-batch SCT.

* **Docs & benchmark**:

  * ``docs/source/usage/batch_integration.rst`` — SCT+RPCA / Harmony / BBKNN.
  * ``docs/source/benchmark/rpca_evaluation.rst`` + ``scripts/run_rpca_benchmark.py``.
  * ``trackcell/benchmark/rpca/`` stepwise and E2E RPCA evaluation suite.

Version 0.3.42
--------------

* Fix ``pyproject.toml`` Poetry extras: declare ``pyglmGamPoi`` / ``xgboost`` as
  optional dependencies and reference bare package names in ``[tool.poetry.extras]``
  (fixes PyPI release workflow for v0.3.41).

Version 0.3.41
--------------

* **SCT v2 parity fixes** (native pyglmGamPoi path):

  * Fix ``vst()`` step-1 cell slicing: use ``cells.get_indexer(cells_step1)`` so
    count columns align with ``cell_attr`` offset (was ``cells.isin``, wrong order).
  * ``mark_suspicious_theta`` uses full-matrix row means/variances (matches R
    ``get_model_pars``).
  * Seurat-parity benchmark loader assigns R rownames by position for duplicate
    genes (``.N`` vs scanpy ``-N``).

* **SCT benchmark** (docs → Benchmark → SCT evaluation):

  * ``trackcell/benchmark/sct/`` — stepwise R→Python injection + end-to-end
    comparison on GSE288946 (3 samples, 6 seeds).
  * ``scripts/run_sct_benchmark.py`` — R export via ``ata`` + Python metrics.
  * QC: ≥100 cells/gene, ≥50 genes/cell. Native HVG Jaccard ~0.98 vs R.

* **Optional dependency** ``trackcell[sct]`` → ``pyglmGamPoi>=0.2.1``.


Version 0.3.40
--------------

* **Fix background image cropping in ROI views**: ``_process_background_image``
  now uses ``data_coords_range`` to crop the background image to the visible
  data region before passing it to ``imshow`` (matching scanpy's behaviour).
  Previously, the full image was always rendered and then zoomed, which could
  cause alignment drift between the background image and cell geometries on
  ROI-subset data due to matplotlib downsampling of large images.


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

