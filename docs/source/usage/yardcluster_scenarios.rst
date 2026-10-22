YardCluster scenario guide
==========================

Practical recipes for **multi-slice Visium HD**, **million-cell STOmics**, and
**Xenium TMA** (including CCHD-style projects).

Concept map (vs DECIPHER / Space Ranger)
----------------------------------------

.. list-table::
   :header-rows: 1

   * - YardCluster
     - DECIPHER (CCHD)
     - Emphasis
   * - ``yardcluster_celltype`` (λ≈0.2)
     - ``leiden_center``
     - Cell-intrinsic / cell typing
   * - ``yardcluster_domain`` (λ≈0.8)
     - ``leiden_nbr``
     - Neighborhood / tissue niche
   * - ``integrate='bbknn'`` + ``expression_source='sct_prep'``
     - ``register_data`` + multi-sample latent
     - Cross-slice expression integration

1. Multi-slice Visium HD
------------------------

**Typical scale:** 1–10 captures per sample, 10k–200k cells per slice.

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   adata = tcl.io.read_hd_cellseg(datapath="segmented_outputs", sample="sample1")
   # Merge multiple captures into one AnnData with batch_key column, e.g. obs['capture_area']

   sc.pp.filter_cells(adata, min_genes=10)
   sc.pp.filter_genes(adata, min_cells=3)

   tcl.tl.run_sct_integration(
       adata,
       batch_key="capture_area",   # or sample_id across merged objects
       layer="counts",
       sct_batch_n_jobs=-1,
   )
   tcl.tl.release_sct_integration_cache(adata, what="scale_data")

   tcl.tl.spatial_cluster(
       adata,
       platform="visium_hd",
       batch_key="capture_area",
       expression_source="sct_prep",
       integrate="bbknn",           # or "joint" for Harmony on embeddings
       graph_mode="hybrid",         # domains: spatial + expression graph
       hybrid_alpha=0.35,
       mode="auto",
   )

   tcl.pl.spatial_cell(adata, color="yardcluster_domain")

**Tips**

* Use **SCT + BBKNN** (not RPCA) when you have many slices (>10).
* ``graph_mode='hybrid'`` often improves domain coherence on HD cell polygons.
* Compare with Space Ranger ``classification`` via ``adjusted_rand_score``.

2. Million-cell STOmics
-----------------------

**Typical scale:** 500k–5M cells per tissue, single or few captures.

.. code-block:: python

   import trackcell as tcl

   adata = tcl.io.read_sto_cellbin("path/to/sample.cellbin.gef", sample="brain1")
   sc.pp.filter_cells(adata, min_genes=50)

   # Multi-sample STOmics (optional SCT front-end)
   # tcl.tl.run_sct_integration(adata, batch_key="sample_id", layer="counts", sct_batch_n_jobs=-1)

   tcl.tl.spatial_cluster(
       adata,
       platform="sto",
       mode="auto",
       sketch=True,                 # or auto when n_obs > 500k
       sketch_n=500_000,
       sketch_stratified=True,
       reuse_neighbors=True,        # fast re-tune resolution
       chunk_size=500,
   )

**Tips**

* **Sketch + stratified sampling** preserves spatial coverage on huge tissues.
* Run on **subset ROI** first to tune ``resolution`` / ``hybrid_alpha``, then full data.
* Single-sample STOmics usually needs ``integrate='none'``; multi-sample → SCT + ``bbknn``.

3. Xenium TMA (CCHD-style)
--------------------------

FLN + SLN 32-slice workflow (SCT + RPCA → YardCluster)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the CCHD cervical TMA **FLN (16) + SLN (16)** layout (``85863_F``,
``85811_S``, ``85863_S``), run integration and spatial clustering in two steps:

.. code-block:: bash

   python scripts/run_cchd_fln_sln_yardcluster.py --step all

Or split steps (recommended on large machines):

.. code-block:: bash

   python scripts/run_cchd_fln_sln_yardcluster.py --step sct_rpca
   python scripts/run_cchd_fln_sln_yardcluster.py --step yardcluster

Step 1 uses ``tcl.tl.integrate_sct_rpca(batch_key='slice_id')`` across 32 slices.
Step 2 uses ``expression_source='sct_integrated'`` (RPCA identity + SCT prep for
spatial M/G) with ``integrate='none'``.

.. code-block:: python

   import trackcell as tcl

   tcl.tl.integrate_sct_rpca(adata, batch_key="slice_id", layer="counts", dims=30)
   tcl.tl.release_sct_integration_cache(adata, what="scale_data")

   tcl.tl.spatial_cluster(
       adata,
       platform="xenium_tma",
       batch_key="slice_id",
       expression_source="sct_integrated",
       integrate="none",
       graph_mode="hybrid",
       hybrid_alpha=0.25,
       mode="auto",
       sketch=True,
       sketch_n=500_000,
   )

**Example project:** cervical cancer Xenium TMA with DBSCAN slice separation,
``obs['slice_id']`` (S001…), ``obs['group']`` (P, LN, FLN, …).

**Typical scale per sample h5ad:** 0.8–2.3M cells, 10–21 slices (see
``xenium_slice_master_summary.csv``).

Preprocessing (geometry)
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   adata = tcl.io.read_xenium_cellseg(
       "/path/to/xenium/region",
       sample="85811_S",
       slice_separate=True,
       slice_eps=80,                # CCHD default; see xenium_slice_separation.rst
   )
   adata = adata[adata.obs["slice_id"] != "debris"].copy()

   # obs['spatial'] from x_centroid / y_centroid (µm) — required for YardCluster

Per-sample YardCluster (recommended first pass)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cluster **within one sample h5ad** (e.g. ``85811_S.h5ad``, 11 LN cores, ~1.87M cells):

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   adata = sc.read_h5ad("/path/to/85811_S.h5ad")
   adata = adata[adata.obs["slice_id"] != "debris"].copy()
   adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

   sc.pp.filter_cells(adata, min_counts=30)
   sc.pp.filter_genes(adata, min_cells=10)

   tcl.tl.run_sct_integration(
       adata,
       batch_key="slice_id",
       layer=None,                  # use adata.X if already counts
       sct_batch_n_jobs=-1,
   )
   tcl.tl.release_sct_integration_cache(adata, what="scale_data")

   tcl.tl.spatial_cluster(
       adata,
       platform="xenium_tma",
       batch_key="slice_id",
       expression_source="sct_prep",
       integrate="bbknn",
       graph_mode="hybrid",
       hybrid_alpha=0.25,
       mode="auto",
       sketch=True,
       sketch_n=400_000,
       resolution_celltype=1.0,
       resolution_domain=0.8,
       key_added="yardcluster",
   )

   # Map to DECIPHER-like column names (optional)
   adata.obs["leiden_center_yard"] = adata.obs["yardcluster_celltype"]
   adata.obs["leiden_nbr_yard"] = adata.obs["yardcluster_domain"]

Downstream (match existing CCHD scripts)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Marker naming: mirror ``02_decipher_primary/step1_markers_*`` on ``yardcluster_*``
* Spatial grids: mirror ``plot_level1_spatial.py`` using ``yardcluster_domain``
* Colony / GC nests: ``yardcluster_domain`` filter → ``tcl.tl.spatial_colony_cluster``

Multi-sample merge (FLN + SLN, etc.)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For combined objects like ``decipher_LN_combined_center_nbr_annotated.h5ad``:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       platform="xenium_tma",
       batch_key="slice_id",        # never cross-slice spatial edges
       expression_source="sct_prep",
       integrate="bbknn",
       mode="auto",
       sketch=True,
       sketch_n=500_000,
   )

Use ``obs['group']`` only for **stratified plotting**, not as ``batch_key`` unless
you intentionally want cross-group expression correction.

4. Topology & interaction analysis (CellCharter-style)
-----------------------------------------------------

Post-clustering analysis for any YardCluster result: resolve physically
disconnected domains, quantify cluster–cluster spatial interactions, and compare
them between conditions (e.g. CCHD Case vs Ctr TMA cores).

.. code-block:: python

   import trackcell as tcl

   # 1. Automatic resolution (stability-based) instead of hand-tuning
   tcl.tl.spatial_cluster(
       adata,
       platform="xenium_tma",
       mode="domain",
       resolution="auto",          # sweep + pick stability peak
       n_resolution_runs=5,
   )
   # sweep table: uns['yardcluster_resolution_sweep'] (resolution/n_clusters/stability)

   # 2. Multi-hop context (larger tissue structures) for domain mode
   tcl.tl.spatial_cluster(adata, mode="domain", n_layers=3)

   # 3. GMM backend with soft assignments / BIC-based k (alternative to Leiden)
   tcl.tl.yard_cluster(
       adata, use_rep="X_yard_mixed", method="gmm", n_clusters=None
   )

   # 4. Split each domain into connected components (per TMA core / tissue piece)
   tcl.tl.connected_components(adata, cluster_key="yardcluster_domain", min_cells=250)
   tcl.pl.spatial_cell(adata, color="component")

   # 5. Neighborhood enrichment (observed vs configuration-model expectation)
   tcl.tl.nhood_enrichment(adata, "yardcluster_domain", log_fold_change=True)
   tcl.pl.nhood_enrichment(adata, "yardcluster_domain")

   # 6. Differential enrichment Case vs Ctr (sample-level permutation p-values)
   tcl.tl.diff_nhood_enrichment(
       adata,
       "yardcluster_domain",
       "group",                    # e.g. Case / Ctr
       ("Case", "Ctr"),
       library_key="slice_id",    # permute at sample level, not cell level
       pvalues=True,
       n_perms=1000,
   )
   tcl.pl.diff_nhood_enrichment(adata, "yardcluster_domain", "group", ("Case", "Ctr"))

Parameter cheat sheet
---------------------

.. list-table::
   :header-rows: 1

   * - Scenario
     - ``platform``
     - ``batch_key``
     - ``integrate``
     - ``graph_mode``
     - Large-data knobs
     - New knobs
   * - Visium HD multi-slice
     - ``visium_hd``
     - ``capture_area`` / ``sample_id``
     - ``bbknn`` or ``joint``
     - ``hybrid``
     - —
     - ``n_layers`` / ``resolution='auto'``
   * - STOmics ≥500k cells
     - ``sto``
     - ``sample_id`` (if multi)
     - ``none`` or ``bbknn``
     - ``expr`` or ``hybrid``
     - ``sketch=True``, ``sketch_stratified=True``
     - ``resolution='auto'``
   * - Xenium TMA (CCHD)
     - ``xenium_tma``
     - ``slice_id``
     - ``bbknn``
     - ``hybrid`` (α≈0.25–0.35)
     - ``sketch_n=400k–500k``, SCT prep
     - ``resolution='auto'``, ``method='gmm'``, ``n_layers``

See also :doc:`spatial_clustering`, :doc:`batch_integration`, :doc:`xenium_slice_separation`.
