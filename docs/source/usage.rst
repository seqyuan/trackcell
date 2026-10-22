Usage Guide
===========

This guide provides examples of how to use TrackCell for reading and analyzing spatial transcriptomics data.

Topics include loading 10x Visium HD / Xenium / **Atera** data, **STOmics (Stereo-seq)** GEF data,
polygon visualization, distance metrics, **DBSCAN TMA slice separation**,
**YardCluster** spatial clustering (``tcl.tl.spatial_cluster``),
**batch integration** (SCT + RPCA / Harmony / BBKNN), and
**CellCharter-style post-analysis** (automatic resolution selection,
multi-hop aggregation, GMM clustering, connected components, neighborhood
enrichment and differential enrichment between conditions).

.. toctree::
   :maxdepth: 2
   :caption: Usage Topics:

   usage/reading
   usage/atera
   usage/annohdcell_conversion
   usage/xenium_slice_separation
   usage/visualization
   usage/computing_distances
   usage/spatial_clustering
   usage/yardcluster_scenarios
   usage/batch_integration
