Benchmark
=========

Numerical evaluation of trackcell against Seurat / ``sctransform`` reference
implementations on real 10x Genomics data.

.. toctree::
   :maxdepth: 2

   sct_evaluation

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Suite
     - Description
   * - :doc:`sct_evaluation`
     - SCTransform v2 (HVG, stepwise injection, end-to-end)
   * - RPCA integration
     - Planned — anchor finding & integration (see ``tcl_test`` WIP)

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

Tests mirroring the benchmark loader live under ``tests/test_sct_*.py``.
