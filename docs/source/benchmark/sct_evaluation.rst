SCTransform v2 evaluation
=========================

Abstract
--------

**trackcell** implements Seurat-compatible **SCTransform v2** (``vst.flavor='v2'``)
with native **pyglmGamPoi** offset fitting, kernel-smoothed dispersion
regularization, and Pearson residual HVG selection. We evaluate parity against
``sctransform::vst`` on three GSE288946 10x samples (GSM8779707–709) under
stringent QC (**≥100 cells per gene**, **≥50 genes per cell**).

The benchmark combines:

1. **Stepwise injection** — R exports at each stage feed the next Python step;
   isolates ksmooth, residuals, and HVG logic from the NB fitter.
2. **End-to-end native Python** — full ``vst()`` vs R on identical filtered matrices.
3. **Randomness control** — six R ``set.seed`` values per sample; Python-vs-R HVG
   Jaccard is compared to the R inter-seed distribution.

Across ref-seed runs on three samples (14724–15297 genes × 15297–15800 cells after QC),
native HVG Jaccard median is **~0.98**, R fit→Python residuals HVG Jaccard **~1.0**,
residual matrix correlation **1.0**, and step1 θ correlation **~0.93–0.96**.

Motivation
----------

SCT v2 is the default normalization for Seurat 4/5 single-sample workflows and
the front end of **SCT + RPCA** integration. Users need confidence that a Python
``AnnData`` pipeline reproduces Seurat HVGs and residuals without an R runtime.

trackcell replaces the R subprocess path with **pyglmGamPoi** while keeping
Seurat-compatible ``uns['sct']`` batch models for integration.

Benchmark design
----------------

Data
  GSE288946 — three GSM replicates, 10x gene symbols (``gene.column = 2``).

Filtering (Python and R identical)
  ``CreateSeuratObject(min.cells = 100, min.features = 50)`` — genes detected
  in fewer than 100 cells and cells with fewer than 50 genes are removed before
  ``vst()``.

Tier-B subsampling
  Fixed ``cells_step1`` / ``genes_step1`` from each R export (``n_cells = n_genes
  = 2000``, ``set.seed(seed)``) so stepwise comparisons are not confounded by
  subsampling RNG when injecting R step-1 parameters.

Seeds (R reference)
  ``1448145`` (ref), ``1448146``, ``42``, ``123``, ``999``, ``2024`` — six
  independent R ``vst()`` runs per sample for HVG stability analysis.

Stepwise protocol
  .. list-table::
     :header-rows: 1
     :widths: 28 36 36

     * - Handoff
       - Input
       - Compared output
     * - R step1 → Py ksmooth
       - ``model_pars_step1.csv``
       - ``model_pars_fit`` (θ, Intercept, log_umi)
     * - R fit → Py residuals
       - ``model_pars_fit.csv``
       - Residual matrix + ``residual_variance``
     * - R step1 golden
       - step1 + R ``genes_log_gmean*``
       - Full downstream through HVG
     * - Native Python
       - Filtered counts only
       - Full ``vst()`` vs R export

Metrics
  Pearson *r*, median |Δ|, Spearman RV correlation, HVG Jaccard (top 3000),
  Inf-θ agreement count.

Implementation notes (fixed during this evaluation)
  .. list-table::
     :header-rows: 1
     :widths: 35 65

     * - Issue
       - Fix
     * - ``umi[:, cells.isin(cells_step1)]`` column order
       - Use ``cells.get_indexer(cells_step1)`` so counts align with
         ``cell_attr.loc[cells_step1]`` offset
     * - Scanpy ``-N`` vs Seurat ``.N`` duplicate gene names
       - Assign R rownames by position from ``ref/{sample}_genes_min_cells.txt``

Results
-------

.. include:: ../_generated/sct_benchmark_table.rst
   :start-line: 0

Interpretation
--------------

**Stepwise injection.**
When R ``model_pars_step1`` is injected, Python ksmooth + residuals achieve HVG
Jaccard ≈ **1.0** — downstream trackcell code matches R given identical NB
parameters.

**Native path.**
With input-alignment fixes, pyglmGamPoi batch θ matches R (raw Inf θ 224/224 on
step-1 subsample; after ``mark_suspicious_theta`` 277/277). HVG Jaccard ≈ **0.99**
vs R ref seed.

**Randomness.**
R inter-seed HVG Jaccard (six seeds, three samples) has median **~0.94** (range
~0.87–0.97). Python-vs-R at the ref seed (median **~0.98**) is at or above the
R inter-seed envelope — remaining HVG differences are within Seurat subsampling
variability, not a systematic Python bias.

Reproduction
------------

.. code-block:: bash

   python scripts/run_sct_benchmark.py --export-r --run --threads 3

Outputs:

- ``/Volumes/process/tmp/tcl_test/steps/sct_benchmark/r/seed_*/GSM*/``
- ``docs/source/_generated/sct_benchmark_results.json``
- ``docs/source/_generated/sct_benchmark_table.rst``

Environment
  Python 3.10–3.11, trackcell with pyglmGamPoi ≥ 0.2.1; R reference in conda
  ``st`` (Seurat 4.4, sctransform, glmGamPoi).

What is not covered here
------------------------

- **RPCA integration** (anchors, prep, integrate) — separate benchmark doc planned.
- **Multi-sample merged SCT** — per-sample ``batch_key`` path tested via unit
  tests; full three-sample merge benchmark is future work.
