Spatial Clustering (YardCluster)
=================================

TrackCell provides **YardCluster**, a lightweight CPU spatial clustering method for
Visium HD cell segmentation and other spatial omics data in ``AnnData`` format.

YardCluster combines ideas from three widely used pipelines:

* **BANKSY** — neighborhood mean (and optional gradient) feature augmentation with a
  mixing parameter ``lambda``
* **DECIPHER** — separate **identity** (cell-intrinsic) and **context**
  (neighborhood) embedding channels before clustering
* **Space Ranger** — optional DE-guided cluster merging after Leiden

Unlike Space Ranger graphclust (expression-only PCA + Louvain), YardCluster uses
spatial coordinates when building neighborhood features, which typically yields
more spatially coherent tissue domains on single-cell HD data.

Prerequisites
-------------

Install TrackCell and clustering dependencies:

.. code-block:: bash

   pip install trackcell
   pip install leidenalg python-igraph

Optional (multi-sample Harmony integration):

.. code-block:: bash

   pip install harmonypy

End-to-End Workflow
-------------------

The typical workflow loads Space Ranger cell segmentation data, optionally runs
scanpy QC, clusters with YardCluster, and visualizes with ``tcl.pl.spatial_cell``.

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   # 1. Load Visium HD cell segmentation (Space Ranger v4)
   adata = tcl.io.read_hd_cellseg(
       datapath="path/to/segmented_outputs",
       sample="sample1",
   )

   # 2. Basic QC (scanpy)
   sc.pp.filter_cells(adata, min_genes=10)
   sc.pp.filter_genes(adata, min_cells=3)

   # 3. YardCluster — includes HVG / normalize / log / scale by default
   tcl.tl.spatial_cluster(adata, mode="auto")

   # 4. Visualize domains and cell-type-oriented clusters
   tcl.pl.spatial_cell(adata, color="yardcluster_domain", figsize=(10, 10))
   tcl.pl.spatial_cell(adata, color="yardcluster_celltype", figsize=(10, 10))

Quick Start (One Function)
--------------------------

``tcl.tl.spatial_cluster()`` runs the full YardCluster pipeline:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       mode="auto",              # celltype + domain in one call
       k_spatial=15,             # spatial neighbors for M/G features
       resolution=1.0,           # Leiden resolution
       key_added="yardcluster",
   )

Clustering modes
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - ``mode``
     - ``lambda``
     - Output column(s)
     - Use case
   * - ``"celltype"``
     - 0.2 (default)
     - ``yardcluster``
     - Cell typing; emphasizes own expression
   * - ``"domain"``
     - 0.8 (default)
     - ``yardcluster``
     - Tissue domains / niches; emphasizes neighborhood
   * - ``"auto"``
     - both
     - ``yardcluster_celltype``, ``yardcluster_domain``
     - Explore both views in one run

Step-by-Step (with scanpy)
--------------------------

For teaching or fine-grained control, call each stage separately. This mirrors
how you would combine scanpy preprocessing with trackcell spatial tools:

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   # --- scanpy preprocessing ---
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)
   sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
   adata = adata[:, adata.var.highly_variable].copy()
   sc.pp.scale(adata, max_value=10)

   # --- YardCluster stages ---
   tcl.tl.spatial_neighbors(adata, k=15, batch_key=None)
   tcl.tl.neighborhood_features(
       adata,
       k_spatial=15,
       use_gradient=False,       # set True for gradient channel G
       gradient_mode="cosphi",   # or "regression"
       key_added="yard",
   )
   tcl.tl.yard_embed(adata, lam=0.2, n_pcs=20, key_added="yard")
   tcl.tl.yard_cluster(
       adata,
       use_rep="X_yard_mixed",
       resolution=1.0,
       key_added="yardcluster",
   )

   # --- optional scanpy visualization on the mixed embedding ---
   sc.pp.neighbors(adata, use_rep="X_yard_mixed")
   sc.tl.umap(adata)
   sc.pl.umap(adata, color="yardcluster")

Function Reference (Summary)
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Function
     - Description
   * - ``tcl.tl.spatial_cluster``
     - One-shot pipeline; recommended entry point
   * - ``tcl.tl.spatial_neighbors``
     - Build within-batch spatial kNN graph (``obsp['spatial_*']``)
   * - ``tcl.tl.neighborhood_features``
     - Compute neighborhood mean ``M`` and optional gradient ``G``
   * - ``tcl.tl.yard_embed``
     - Dual-channel PCA + ``lambda`` mixing → ``obsm['X_yard_mixed']``
   * - ``tcl.tl.yard_cluster``
     - Leiden clustering on a chosen embedding

See the :doc:`../api` page for full parameter lists.

Outputs Stored in ``AnnData``
-----------------------------

After ``spatial_cluster(..., mode='auto', key_added='yardcluster')``:

**Cluster labels (``obs``)**

* ``yardcluster_celltype`` — low-``lambda`` (cell typing)
* ``yardcluster_domain`` — high-``lambda`` (domains)

**Embeddings (``obsm``)** — internal prefix ``yardcluster_features`` when using one-shot:

* ``X_yardcluster_features_identity`` — identity channel PCA
* ``X_yardcluster_features_context`` — context channel PCA
* ``X_yardcluster_features_mixed`` — sqrt-weighted mix used for Leiden

**Graphs (``obsp``)**

* ``spatial_connectivities`` / ``spatial_distances`` — physical neighborhood graph
* ``yardcluster_celltype_connectivities`` — expression kNN for celltype run
* ``yardcluster_domain_connectivities`` — expression kNN for domain run

**Parameters (``uns``)**

* ``yardcluster_params`` — run configuration recap

Visualizing Results
-------------------

Polygon visualization (Visium HD cell segmentation):

.. code-block:: python

   tcl.pl.spatial_cell(
       adata,
       color="yardcluster_domain",
       figsize=(12, 12),
       legend_loc="right margin",
   )

Compare with Space Ranger graphclust (pre-computed at load time):

.. code-block:: python

   # SR labels are in obs['classification'] after read_hd_cellseg()
   tcl.pl.spatial_cell(adata, color="classification", figsize=(12, 12))
   tcl.pl.spatial_cell(adata, color="yardcluster_domain", figsize=(12, 12))

Side-by-side with scanpy UMAP on the YardCluster embedding:

.. code-block:: python

   sc.pp.neighbors(adata, use_rep="X_yardcluster_features_mixed")
   sc.tl.umap(adata)
   sc.pl.umap(adata, color=["yardcluster_domain", "classification"])

Multi-Sample and Batch Handling
---------------------------------

Spatial neighbors are **never drawn across batches** (different tissue slices have
incompatible coordinate systems). Expression integration is controlled separately.

.. code-block:: python

   # Cluster each sample independently (labels prefixed: "sampleA_0", …)
   tcl.tl.spatial_cluster(
       adata,
       batch_key="sample_id",
       integrate="separate",
       mode="domain",
   )

.. code-block:: python

   # Joint integration via Harmony (requires harmonypy)
   tcl.tl.spatial_cluster(
       adata,
       batch_key="sample_id",
       integrate="joint",
       mode="celltype",
   )

.. list-table::
   :header-rows: 1

   * - ``integrate``
     - Behavior
   * - ``"none"``
     - Single sample, or pooled without Harmony (default)
   * - ``"separate"``
     - Per-batch clustering; requires ``batch_key``
   * - ``"joint"``
     - Harmony on identity/context embeddings; requires ``batch_key`` and ``harmonypy``

Advanced Options
----------------

DE-Guided Cluster Merging
~~~~~~~~~~~~~~~~~~~~~~~~~

Space Ranger merges sibling clusters when no genes are differentially expressed.
YardCluster offers a similar optional post-processing step:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       mode="domain",
       merge_clusters=True,
       merge_adj_p=0.05,
   )

Large Datasets (Sketch Mode)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For very large cell counts (>500k by default), YardCluster can cluster a random
subsample and assign remaining cells to the nearest cluster medoid in embedding
space:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       sketch=True,
       sketch_n=500_000,
       sketch_threshold=500_000,
   )

Gradient Features
~~~~~~~~~~~~~~~~~

Neighborhood **mean** ``M`` is always computed. Optional **gradient** ``G``
 augments the context channel:

.. code-block:: python

   tcl.tl.neighborhood_features(
       adata,
       use_gradient=True,
       gradient_mode="cosphi",      # fast BANKSY-like approximation (default)
       # gradient_mode="regression",  # local linear regression magnitude
       k_gradient=30,               # default: 2 * k_spatial
   )

Preprocessing Control
~~~~~~~~~~~~~~~~~~~~~

By default (``preprocess=True``), ``spatial_cluster`` runs NormalizeTotal →
Log1p → HVG → Scale before clustering.  Fine-tune or skip:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       mode="auto",
       preprocess=True,          # default: run standard scanpy pipeline
       n_top_genes=2000,         # HVG count (default 2000)
       hvg_by_batch=True,        # run HVG per batch when batch_key is set
       scale_by_batch=False,     # scale each batch independently
   )

   # If data is already log-normalised and scaled
   tcl.tl.spatial_cluster(adata, preprocess=False, mode="auto")

   # Use a pre-computed layer (e.g. scVI denoised) instead of adata.X
   tcl.tl.spatial_cluster(adata, layer="denoised", preprocess=False, mode="domain")

   # Use adata.raw.X instead of adata.X
   tcl.tl.spatial_cluster(adata, use_raw=True, mode="celltype")

Harmony Integration (per-channel)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``integrate='joint'``, Harmony is run on **both** identity and context
embeddings before mixing.  You can also enable Harmony independently of the
``integrate`` parameter:

.. code-block:: python

   tcl.tl.spatial_cluster(
       adata,
       batch_key="sample_id",
       harmony_integrate=True,   # Harmony on each embedding channel
       integrate="none",         # no separate/joint mode
       mode="domain",
   )

This differs from ``integrate='joint'`` in that ``harmony_integrate`` only
corrects the embedding channels, while ``integrate='joint'`` also triggers
Harmony during the embedding stage.

Return-Copy Mode
~~~~~~~~~~~~~~~~

Set ``copy=True`` to leave the original ``adata`` untouched and return a
modified copy:

.. code-block:: python

   result = tcl.tl.spatial_cluster(adata, mode="auto", copy=True)

When to Use YardCluster vs Space Ranger graphclust
--------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 30 35

   * - Aspect
     - Space Ranger graphclust
     - YardCluster
   * - Uses spatial coordinates
     - No
     - Yes (neighborhood M/G)
   * - Pre-computed in SR output
     - Yes (``obs['classification']``)
     - No — run on demand
   * - Clustering backend
     - Louvain + DE merge
     - Leiden (+ optional DE merge)
   * - GPU required
     - No
     - No
   * - Best for
     - Quick look at SR defaults
     - HD cell-level domains, re-clustering, custom ``lambda``

If you already have satisfactory ``classification`` labels from Space Ranger, you
can use them directly. Run YardCluster when you need spatially informed domains,
different ``lambda`` trade-offs, or to re-cluster after QC/subsetting.

Compare SR vs YardCluster manually when both are available:

.. code-block:: python

   from sklearn.metrics import adjusted_rand_score

   ari = adjusted_rand_score(
       adata.obs["classification"],
       adata.obs["yardcluster_domain"],
   )
   print(f"ARI (SR graphclust vs YardCluster domain): {ari:.3f}")

Troubleshooting
---------------

**``ImportError`` for Leiden**

Install ``leidenalg`` and ``python-igraph`` (see Prerequisites).

**Too many / too few clusters**

Adjust ``resolution`` (e.g. 0.5–1.5) or mode-specific
``resolution_celltype`` / ``resolution_domain``.

**Out of memory on large gene panels**

Neighborhood aggregation uses gene chunks (default 500 genes). Reduce
``n_top_genes`` during preprocessing or increase ``chunk_size`` if you have RAM.

**``integrate='joint'`` fails**

Install ``harmonypy`` and ensure ``batch_key`` has ≥2 batches with sufficient
cells per batch.

**Already preprocessed data**

If ``adata`` is already log-normalized and scaled, disable built-in preprocessing:

.. code-block:: python

   tcl.tl.spatial_cluster(adata, preprocess=False, mode="auto")

See Also
--------

* :doc:`reading` — loading Visium HD and Xenium data
* :doc:`visualization` — ``tcl.pl.spatial_cell`` plotting options
* :doc:`computing_distances` — distance to label after clustering
* :doc:`../api` — autogenerated API reference
