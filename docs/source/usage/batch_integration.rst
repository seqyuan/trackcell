Batch integration (SCT + RPCA / Harmony / BBKNN)
=================================================

trackcell separates **SCT normalization for integration** from the **downstream
correction method**. All paths share the same Seurat v4 front end:

1. Per-batch ``SCTransform`` (``sctransform``, ``vst.flavor='v2'``)
2. ``SelectIntegrationFeatures`` (default 2000 genes)
3. ``PrepSCTIntegration`` → anchor-gene Pearson residuals in ``obsm['X_sct_prep']``

After step 3 you can run **RPCA** (Seurat default), **Harmony**, **BBKNN**, or any
method that consumes a cells × genes embedding matrix.

API overview
------------

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Function
     - Role
     - Output
   * - ``tcl.tl.sctransform``
     - Per-batch SCT only
     - ``uns['sct']``, optional ``obsm['X_sct']``
   * - ``tcl.tl.run_sct_integration``
     - Steps 1–3 (SCT + features + prep)
     - ``obsm['X_sct_prep']``, ``uns['sct_prep']``
   * - ``tcl.tl.sct_prep_matrix``
     - Helper: matrix + gene names
     - ``(ndarray, list[str])``
   * - ``tcl.tl.integrate_rpca``
     - Seurat RPCA anchors + IntegrateData
     - ``obsm['X_sct_integrated']``
   * - ``tcl.tl.release_sct_integration_cache``
     - **Explicit** post-prep memory release (opt-in)
     - Drops ``uns`` SCT caches; prep unchanged
   * - ``tcl.tl.integrate_sct_rpca``
     - Convenience: ``run_sct_integration`` + ``integrate_rpca``
     - Both prep and integrated matrices

Method combinations
-------------------

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Workflow
     - Best for
     - trackcell calls
   * - **SCT + RPCA**
     - Seurat v4 parity; reference mapping; ≤10 batches
     - ``integrate_sct_rpca`` or split API below
   * - **SCT + Harmony**
     - Many slices / batches; fast linear correction in PC space
     - ``run_sct_integration`` → PCA → ``harmonypy``
   * - **SCT + BBKNN**
     - Batch-aware neighbour graph for clustering / UMAP
     - ``run_sct_integration`` → ``bbknn``
   * - **SCT only**
     - Custom integrators (scVI, Scanorama, …)
     - ``run_sct_integration`` → your method on ``X_sct_prep``

When to choose which method
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SCT + RPCA** — closest to Seurat ``FindIntegrationAnchors(reduction='rpca')`` +
``IntegrateData``. Use when you need anchor-based correction or must match a Seurat
reference workflow. Anchor finding scales as **O(B²)** pairwise comparisons; for
**dozens of spatial slices** prefer Harmony or BBKNN on the prep matrix.

**SCT + Harmony** — recommended default for **multi-slice spatial** data: one PCA
on shared anchor residuals, then Harmony removes batch effects in embedding space.
Does not materialise a full integrated expression matrix; downstream clustering /
UMAP uses the Harmony embedding.

**SCT + BBKNN** — builds a batch-balanced kNN graph directly on ``X_sct_prep`` (or
PCA thereof). Good when the goal is **clustering / neighbourhood** rather than a
corrected expression matrix.

   * - ``tcl.tl.integrate_sct_rpca``
     - Convenience: ``run_sct_integration`` + ``integrate_rpca``
     - Both prep and integrated matrices

Seurat parity first
-------------------

**RPCA results matching Seurat take priority.** trackcell uses the same pipeline as
Seurat v4:

* Sample-tree merging (not fixed-reference shortcuts)
* ``n_trees=50`` for anchor finding
* For benchmark-grade parity on IntegrateData, also set ``n_trees_weight=50`` and
  ``integration_dtype='float64'``

Nothing in ``uns`` is deleted automatically. Memory trade-offs are **always**
user-initiated (see below).

SCT + RPCA (Seurat v4)
----------------------

One-shot
~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   tcl.tl.integrate_sct_rpca(
       adata,
       batch_key="orig.ident",
       layer="counts",
       n_features=2000,
       dims=30,
       run_sct=True,
       copy=False,
   )
   # adata.obsm["X_sct_prep"]      — anchor residuals (input to other methods)
   # adata.obsm["X_sct_integrated"] — RPCA-corrected residuals

Split API (reuse prep for Harmony / RPCA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   tcl.tl.run_sct_integration(
       adata,
       batch_key="orig.ident",
       layer="counts",
       n_features=2000,
       sct_batch_n_jobs=4,   # parallel SCT per batch / slice
   )

   # --- explicit memory release (optional; does NOT change RPCA output) ---
   # Pick ONE level after prep is final:

   # (A) Drop per-batch HVG residual caches only (~largest win):
   for _batch, entry in adata.uns["sct"]["batch_models"].items():
       entry.pop("scale_data", None)

   # (B) Drop all batch_models (keep uns['sct']['params']):
   # del adata.uns["sct"]["batch_models"]

   # (C) Drop entire SCT uns block (prep in obsm is enough for RPCA):
   # del adata.uns["sct"]

   # Or use the helper (same as A/B/C):
   # tcl.tl.release_sct_integration_cache(adata, what="scale_data")

   # Option A: RPCA correction (Seurat sample tree)
   tcl.tl.integrate_rpca(
       adata,
       batch_key="orig.ident",
       dims=30,
       n_trees=50,              # anchors — keep 50 for Seurat parity
       n_trees_weight=50,       # integrate — 50 for parity; 10 for speed
       integration_dtype="float64",  # float64 for parity; float32 saves RAM
   )

   # Option B: export prep for external tools
   matrix, genes = tcl.tl.sct_prep_matrix(adata)

SCT + Harmony
-------------

Requires ``harmonypy`` (``pip install harmonypy``).

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   tcl.tl.run_sct_integration(
       adata,
       batch_key="slice_id",
       layer="counts",
       sct_batch_n_jobs=-1,
   )

   sc.pp.pca(adata, use_rep="X_sct_prep", n_comps=30)
   sc.external.pp.harmony_integrate(
       adata,
       key="X_pca",
       basis="X_pca_harmony",
       groupby="slice_id",
   )

   sc.pp.neighbors(adata, use_rep="X_pca_harmony")
   sc.tl.umap(adata)

For spatial clustering that already supports Harmony (``tcl.tl.spatial_cluster``),
see :doc:`spatial_clustering` — that path applies Harmony to **spatial embeddings**,
not SCT prep residuals. Use the recipe above when integrating **expression** across
slices before standard scRNA-seq-style analysis.

SCT + BBKNN
-----------

Requires ``bbknn`` (``pip install bbknn``).

.. code-block:: python

   import scanpy as sc
   import trackcell as tcl

   tcl.tl.run_sct_integration(adata, batch_key="slice_id", layer="counts")

   # On prep residuals directly, or PCA first for speed:
   sc.pp.pca(adata, use_rep="X_sct_prep", n_comps=30)
   sc.external.pp.bbknn(adata, batch_key="slice_id", use_rep="X_pca")

   sc.tl.umap(adata)

Spatial multi-slice tips
------------------------

* Set ``batch_key`` to your slice / capture-area column (e.g. ``slice_id``,
  ``orig.ident``).
* Use ``sct_batch_n_jobs`` to parallelise per-slice SCT (``-1`` = all CPUs).
* For **>10 slices** where Seurat parity is not required, consider **SCT + Harmony**
  instead of RPCA.
* Prep matrix size ≈ ``n_cells × n_features × 4`` bytes (float32 in ``obsm``).
  30 slices × 10k cells × 2000 genes ≈ **2.4 GB** for ``X_sct_prep`` alone.

Explicit memory release (opt-in)
--------------------------------

After ``run_sct_integration``, ``adata.uns['sct']['batch_models']`` may hold large
per-slice ``scale_data`` caches (full HVG residual matrices). **RPCA and Harmony
only need** ``obsm['X_sct_prep']`` — you may delete caches **manually** before
merge to save RAM.

.. important::

   trackcell **never** removes ``uns`` entries automatically. You must run one of
   the commands below (or ``release_sct_integration_cache``) yourself.

.. list-table::
   :header-rows: 1
   :widths: 14 46 40

   * - Level
     - Manual command
     - Changes RPCA vs Seurat?
   * - A
     - ``entry.pop("scale_data", None)`` per batch
     - **No** (prep already in ``obsm``)
   * - B
     - ``del adata.uns["sct"]["batch_models"]``
     - **No**
   * - C
     - ``del adata.uns["sct"]``
     - **No**

Helper (equivalent to A/B/C):

.. code-block:: python

   tcl.tl.release_sct_integration_cache(adata, what="scale_data")   # level A
   tcl.tl.release_sct_integration_cache(adata, what="batch_models") # level B
   tcl.tl.release_sct_integration_cache(adata, what="sct_uns")    # level C

**When you cannot release yet**

* You still need to change ``anchor_features`` and re-run ``prep_sct_integration``
  (missing genes are recomputed from counts + ``model_pars`` — slower without
  ``scale_data``, same numerics).
* You have not finished ``run_sct_integration`` / prep (no ``X_sct_prep`` yet).

**What we do not do (Seurat parity)**

* No automatic cache trimming inside ``integrate_sct_rpca`` or ``integrate_rpca``.
* No ``reference_batch`` / fixed-reference merge shortcuts that alter Seurat's
  sample-tree correction path.

Performance knobs (RPCA path)
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Parameter
     - Default
     - Effect
   * - ``sct_batch_n_jobs``
     - ``1``
     - Parallel per-batch SCT when ``batch_key`` is set
   * - ``n_trees``
     - ``50``
     - Annoy trees for **anchor finding** (keep high for parity)
   * - ``n_trees_weight``
     - ``10``
     - Annoy trees for **IntegrateData** weight PCA. Use **``50``** for Seurat
       parity; ``10`` is faster with small metric impact
   * - ``integration_dtype``
     - ``float32``
     - Halves PCA / correction workspace. Use **``float64``** for Seurat parity
   * - ``weight_query_chunk_size``
     - ``8192``
     - Chunk size for weight matrix construction (lower = less RAM)

Benchmark / Seurat parity mode — pass ``n_trees_weight=50``,
``integration_dtype='float64'``. See :doc:`../benchmark/rpca_evaluation`.

R export bridge (optional)
--------------------------

When R ``SCTransform`` exports exist (e.g. from
``scripts/run_rpca_benchmark.py --export-r``), hybrid parity mode injects R
step-1 subsampling and cached residuals:

.. code-block:: python

   tcl.tl.integrate_sct_rpca(
       adata,
       batch_key="orig.ident",
       r_vst_export_root="/path/to/tcl_test",
       r_vst_export_seed=1448145,
       inject_r_scale_data=True,
   )

See :doc:`../benchmark/rpca_evaluation` for measured parity with and without
R exports.

Further reading
---------------

* :doc:`../benchmark/rpca_evaluation` — Python vs Seurat RPCA benchmark
* :doc:`../benchmark/sct_evaluation` — SCTransform v2 single-sample benchmark
* :doc:`spatial_clustering` — YardCluster (spatial domains; separate from SCT integration)
