TrackCell Documentation
========================

.. image:: _static/images/logo.svg
   :alt: TrackCell Logo
   :align: center
   :width: 200px

Welcome to TrackCell's documentation!

TrackCell is a Python package for processing and visualizing single-cell and spatial transcriptomics data,
particularly designed for 10X Genomics HD spatial transcriptomics data.

Key features
------------

* **YardCluster** — lightweight CPU spatial clustering with BANKSY-style identity/context
  embeddings, Leiden or Gaussian-mixture backends, automatic resolution selection,
  multi-hop neighborhood aggregation, and sketch-mode clustering for million-cell data.
* **SCT + RPCA integration** — Seurat v4-parity batch integration across slices.
* **Post-clustering analysis** — connected components, neighborhood enrichment,
  and differential enrichment between conditions (CellCharter-inspired).
* **IO & visualization** — Visium HD cell segmentation, Xenium, **Atera**
  (official Zarr + streaming transcript point cloud), STOmics GEF, TMA slice
  separation, polygon rendering.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   benchmark/index
   api
   examples
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

