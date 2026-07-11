SCT + RPCA integration evaluation
===================================

Abstract
--------

**trackcell** implements Seurat v4 **SCT + reciprocal PCA (RPCA)** integration on a
single ``AnnData``: ``SelectIntegrationFeatures``, ``PrepSCTIntegration``,
``FindIntegrationAnchors(reduction='rpca')``, and ``IntegrateData`` with Annoy
MNN anchors and sample-tree merging.

We evaluate parity on three GSE288946 10x samples (GSM8779707–709) under the same
QC as the SCT benchmark (**≥100 cells per gene**, **≥50 genes per cell**).
After QC, **35,948 cells** remain; integration uses **2,000** anchor features and
ref seed **1448145**.

**Headline (ref seed 1448145, regenerated 2026-07-11):** with hybrid SCT exports
(R step-1 subsampling + R ``residuals_hvg`` cache), end-to-end **IntegrateData**
gene correlation vs Seurat is **~0.959**; stepwise with R prep + R anchors reaches
**~0.959** as well. Anchor local-pair Jaccard remains the main identity gap
(**~0.72** E2E, **~0.87** stepwise).

Motivation
----------

SCT + RPCA is the default Seurat v4 integration workflow. Users need confidence
that a Python ``AnnData`` pipeline reproduces Seurat integration features, anchor
residuals, and corrected embeddings without an R runtime — or, when R exports are
available, that each Python step is numerically sound.

This benchmark complements :doc:`sct_evaluation`: the SCT suite compares against
``sctransform::vst()`` on single samples; the RPCA suite compares against Seurat
``SCTransform()`` inside the full integration pipeline (reference exports under
``steps/rpca_benchmark/r/seed_*/``).

Benchmark design
----------------

Three tiers of comparison:

1. **Integration SCT stepwise** — per-batch native ``vst()`` vs R ``SCTransform``
   (``02_sct``); R ``model_pars_fit`` → Python residuals; native SCT +
   **forced R integration feature list** → prep (isolates prep formula).
2. **RPCA stepwise injection** — R export at step *N* feeds Python step *N+1*
   (features → prep → anchors → integrate).
3. **End-to-end native Python** — ``run_sct_integration`` + ``integrate_rpca``
   (hybrid SCT: Python + R step-1 subsampling + R scale.data overlay) vs full
   Seurat pipeline.

Anchor metrics use **local (dataset, cell) pairs**: Seurat exports batch-local
1-based indices; Python uses global 0-based indices after ``map_anchors_to_global``.

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Handoff
     - Input
     - Compare
   * - SCT per batch
     - Same filtered UMI as R integration
     - Native ``vst()`` HVG vs R ``SCTransform``; R fit → Py residuals
   * - Prep (merged)
     - Native SCT + **R integration feature list**
     - ``PrepSCTIntegration`` gene correlation
   * - 03 features
     - R ``02_sct`` batch models
     - ``SelectIntegrationFeatures`` Jaccard
   * - 04 prep
     - R SCT + R features
     - Prep residual matrix correlation
   * - 05 anchors
     - R prep matrix
     - Local anchor-pair Jaccard
   * - 06 integrate
     - R prep + R anchors
     - Integrated residual correlation

Python vs Seurat summary (ref seed)
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 14 14 14 32

   * - Component
     - Stepwise
     - E2E
     - Target
     - Notes
   * - ``SelectIntegrationFeatures``
     - **1.000**
     - **0.929**
     - 1.0
     - Stepwise injects R SCT; E2E uses hybrid native SCT
   * - ``PrepSCTIntegration``
     - **1.000**
     - **1.000**
     - ~1.0
     - E2E uses R ``residuals_hvg`` scale.data overlay (2026-07)
   * - Anchor local-pair Jaccard
     - **0.868**
     - **0.720**
     - ≥0.95
     - Counts match (~39k); identity differs (Annoy / prep drift)
   * - ``IntegrateData`` gene corr (median)
     - **0.959**
     - **0.959**
     - ≥0.95
     - Shared features only: **0.977** (E2E)
   * - SCT HVG Jaccard (per batch, median)
     - **0.918** native
     - **0.918** hybrid E2E
     - ≥0.95
     - R fit→Py HVG **~0.95** per batch

Detailed tables
---------------

.. include:: ../_generated/rpca_benchmark_table.rst

The generated table above includes a **Python vs Seurat quick comparison** (stepwise
vs E2E) when produced by ``scripts/run_rpca_benchmark.py --run``.

Conclusions (2026-07-11)
------------------------

What is aligned
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 32 18 40

   * - Component
     - Parity
     - Notes
   * - ``PrepSCTIntegration``
     - **~1.0**
     - Stepwise and E2E (with R scale.data injection)
   * - ``IntegrateData`` correction
     - **~0.959**
     - Stepwise 06 and E2E; meets ≥0.95 target
   * - ``SelectIntegrationFeatures`` (R SCT injected)
     - **1.0**
     - Stepwise 03
   * - R ``model_pars_fit`` → Py residuals / HVG
     - **~0.95** HVG Jaccard
     - Per-batch ``02_sct`` export
   * - Anchor **counts**
     - **~1%** diff E2E
     - 39,124 (R) vs ~38,700 (Py E2E)

Remaining gaps
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 32 18 40

   * - Gap
     - Magnitude
     - Mitigation
   * - Native SCT HVG vs R (no export)
     - **~0.92** Jaccard
     - Use ``r_vst_export_root`` + step-1 injection, or accept ~8% HVG drift
   * - E2E feature list vs R
     - **0.929** Jaccard
     - Follows HVG drift; 1926/2000 genes shared
   * - Anchor local-pair identity
     - **0.72** E2E / **0.87** stepwise
     - Annoy tie order; upstream prep differences in pure-native mode
   * - Pure-native E2E (no R exports)
     - Integrate **~0.89** (pre-P2)
     - Hybrid path raises to **~0.959**

Implementation milestones
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Date
     - Change
   * - 2026-07-10
     - Weight formula + sparse ``W^T @ IM`` aligned to Seurat ``FindWeightsC`` /
       ``IntegrateDataC``; integrate corr **0.88 → 0.96** (stepwise).
   * - 2026-07-11
     - **P0** hybrid E2E: ``run_sct_integration`` with R step-1 subsampling +
       ``inject_r_scale_data``; prep **1.0**, integrate **0.959**.
   * - 2026-07-11
     - Split API: ``run_sct_integration`` / ``integrate_rpca``; SCT reusable for
       Harmony / BBKNN.
   * - 2026-07-11
     - **P0 perf**: vectorised chunked ``_build_weights``; ``n_trees_weight``;
       ``integration_dtype='float32'``.

Key diagnostic findings
~~~~~~~~~~~~~~~~~~~~~~~

**Prep drives anchor identity, not the anchor algorithm.** Decomposition on ref
seed 1448145:

.. code-block:: text

   Input prep matrix          Anchor J vs R
   ─────────────────────────────────────────
   Py prep + Py features      0.61
   Py prep + R features       0.60   ← features not the driver
   R prep  + R features       0.87   ← prep is the driver

**Integrate is downstream of prep.** With R prep, integrate corr **~0.959** whether
anchors are R or Python (**0.9595** vs **0.9594** stepwise).

**Reference choice matters.** Single-sample SCT exports use ``sctransform::vst()``;
integration exports use Seurat ``SCTransform()``. Measure RPCA parity against
``steps/rpca_benchmark/``, not ``steps/sct_benchmark/``.

Deployment modes
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Mode
     - API
     - Expected parity
   * - **Hybrid (recommended)**
     - ``integrate_sct_rpca(..., r_vst_export_root=...)``
     - Prep **1.0**, integrate **~0.96**
   * - **Pure Python**
     - ``integrate_sct_rpca`` without exports
     - Integrate **~0.89–0.95**; faster, no R dependency
   * - **SCT + Harmony / BBKNN**
     - ``run_sct_integration`` only
     - No Seurat RPCA benchmark; see :doc:`../usage/batch_integration`

Memory vs parity
~~~~~~~~~~~~~~~~

Seurat RPCA parity is **never** traded for automatic RAM savings. After prep,
``obsm['X_sct_prep']`` is sufficient for ``integrate_rpca``; dropping
``uns['sct']['batch_models'][*]['scale_data']`` (manually or via
``release_sct_integration_cache``) does **not** change anchor or integrate metrics.
See :doc:`../usage/batch_integration` for explicit one-liners.

Regenerate
----------

Requires GSE288946 data under ``/Volumes/process/tmp/tcl_test/data/`` and R with
Seurat + sctransform (``st`` conda env recommended; needs **glmGamPoi**).

.. code-block:: bash

   # Full: R export + Python comparisons + docs/_generated/*.rst
   python scripts/run_rpca_benchmark.py --export-r --run

   # Skip slow E2E (~20 min):
   python scripts/run_rpca_benchmark.py --run --skip-e2e --seeds 1448145

R-only or Python-only:

.. code-block:: bash

   python scripts/run_rpca_benchmark.py --export-r
   python scripts/run_rpca_benchmark.py --run

Implementation references
-------------------------

- Benchmark runner: ``trackcell/benchmark/rpca/``
- Integration SCT stepwise: ``trackcell/benchmark/rpca/sct_stepwise.py``
- Split integration API: ``trackcell/tl/integration/sct_integration.py``,
  ``trackcell/tl/integration/rpca_integration.py``
- R export: ``trackcell/tl/rscripts/export_rpca_benchmark_r.R``
- CLI: ``scripts/run_rpca_benchmark.py``
- Usage: :doc:`../usage/batch_integration`
