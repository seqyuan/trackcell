Benchmark
=========

Numerical evaluation of trackcell against Seurat / ``sctransform`` reference
implementations on real 10x Genomics data.

.. toctree::
   :maxdepth: 2

   sct_evaluation
   rpca_evaluation
   yardcluster_evaluation

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Suite
     - Description
   * - :doc:`sct_evaluation`
     - SCTransform v2 (HVG, stepwise injection, end-to-end)
   * - :doc:`rpca_evaluation`
     - SCT + RPCA integration vs Seurat (stepwise + E2E parity tables,
       Python vs R summary, 2026-07-11)
   * - :doc:`yardcluster_evaluation`
     - YardCluster synthetic patch benchmark (ARI + spatial coherence)

Regenerate SCT artifacts
------------------------

Benchmark data is read from ``TRACKCELL_BENCHMARK_ROOT`` (default:
``~/trackcell_benchmark``). Set it before running, e.g.:

.. code-block:: bash

   export TRACKCELL_BENCHMARK_ROOT=/path/to/tcl_test

Requires GSE288946 data under ``<benchmark-root>/data/`` and conda env
``st`` (Seurat + sctransform + glmGamPoi).

.. code-block:: bash

   # R exports (parallel) + Python comparisons + docs/_generated/*.rst
   python scripts/run_sct_benchmark.py --export-r --run --threads 3

R-only or Python-only:

.. code-block:: bash

   python scripts/run_sct_benchmark.py --export-r --threads 3
   python scripts/run_sct_benchmark.py --run

Regenerate RPCA artifacts
-------------------------

Same data root; one R job per seed (full 3-sample integration pipeline).
See :doc:`rpca_evaluation` for the **Python vs Seurat** summary table and E2E
results (integrate gene corr **~0.959** with hybrid SCT exports).

Quick reference (ref seed 1448145):

GSE288946 dataset (35,948 cells, hybrid SCT exports):

.. list-table::
   :header-rows: 1
   :widths: 32 14 14

   * - Metric
     - Stepwise
     - E2E
   * - Features Jaccard
     - 1.000
     - 0.929
   * - Prep gene corr
     - 1.000
     - 1.000
   * - Anchor local-pair Jaccard
     - 0.868
     - 0.720
   * - Integrate gene corr
     - 0.959
     - 0.959

AA/A04 dataset (27,129 cells, pure Python backend, 2026-07-25):

.. list-table::
   :header-rows: 1
   :widths: 32 14 14

   * - Metric
     - Pure Py
     - Target
   * - SCT HVG Jaccard
     - 0.939
     - ≥0.92
   * - Features Jaccard
     - 0.933
     - ≥0.92
   * - Prep gene corr
     - 0.988
     - ~1.0
   * - Anchor count ratio
     - 0.989
     - ~1.0

Full tables: :doc:`rpca_evaluation`.

.. code-block:: bash

   python scripts/run_rpca_benchmark.py --export-r --run
   python scripts/run_rpca_benchmark.py --run --skip-e2e --seeds 1448145

Tests mirroring the benchmark loader live under ``tests/test_sct_*.py``.
