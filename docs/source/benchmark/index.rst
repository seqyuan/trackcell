Benchmark
=========

Numerical evaluation of trackcell against Seurat / ``sctransform`` reference
implementations on real 10x Genomics data.

.. toctree::
   :maxdepth: 2

   sct_evaluation
   rpca_evaluation

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

Regenerate SCT artifacts
------------------------

Requires GSE288946 data under ``/Volumes/process/tmp/tcl_test/data/`` and conda env
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

Full tables: :doc:`rpca_evaluation`.

.. code-block:: bash

   python scripts/run_rpca_benchmark.py --export-r --run
   python scripts/run_rpca_benchmark.py --run --skip-e2e --seeds 1448145

Tests mirroring the benchmark loader live under ``tests/test_sct_*.py``.
