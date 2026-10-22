YardCluster evaluation
========================

Synthetic benchmark comparing YardCluster configurations on spatial patch data
with known region labels.

Motivation
----------

Unlike SCT/RPCA, YardCluster had no automated regression suite. This benchmark
tracks **ARI** (vs ground-truth patch labels) and **spatial coherence**
(fraction of spatial neighbors sharing a cluster) across configuration presets.

Configurations
--------------

.. list-table::
   :header-rows: 1

   * - Config
     - Description
   * - ``baseline_no_gradient``
     - ``mode='domain'``, ``use_gradient=False``, ``merge_clusters=False``
   * - ``defaults_auto``
     - ``mode='auto'`` with mode-aware defaults (gradient on, domain DE-merge)
   * - ``domain_gradient``
     - ``mode='domain'``, ``use_gradient=True``, ``merge_clusters=True``

Run locally
-----------

.. code-block:: bash

   python scripts/run_yardcluster_benchmark.py
   python scripts/run_yardcluster_benchmark.py --json --seed 0

Programmatic API:

.. code-block:: python

   from trackcell.benchmark.yardcluster import run_yardcluster_benchmark

   summary = run_yardcluster_benchmark(seed=0)
   print(summary["best_config"], summary["best_ari"])

Real-data evaluation
--------------------

When Space Ranger ``classification`` or anatomical labels exist:

.. code-block:: python

   from sklearn.metrics import adjusted_rand_score
   from trackcell.benchmark.yardcluster.metrics import spatial_coherence

   ari = adjusted_rand_score(adata.obs["classification"], adata.obs["yardcluster_domain"])
   coh = spatial_coherence(adata, "yardcluster_domain")

See also :doc:`../usage/spatial_clustering`.
