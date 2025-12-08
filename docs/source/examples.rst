Examples
========

This section contains example notebooks and code snippets demonstrating TrackCell usage.

Basic Workflow
--------------

A complete workflow example:

.. code-block:: python

   import trackcell as tcl
   import scanpy as sc

   # 1. Read cell segmentation data
   adata = tcl.io.read_hd_cellseg(
       datapath="path/to/segmented_outputs",
       sample="sample1"
   )

   # 2. Compute distance to a specific cell type
   tcl.tl.hd_labeldist(
       adata,
       groupby="classification",
       label="Cluster-2",
       inplace=True
   )

   # 3. Visualize with cell polygons
   tcl.pl.spatial_cell(
       adata,
       color=["classification", "Cluster-2_dist"],
       figsize=(12, 6)
   )


Examples
----------

Colon Cancer Demo

.. nbsphinx:: notebooks/Colon_Cancer_demo.ipynb
   :raw:

.. nbgallery::
    :caption: varies Hi-C map
    :name: rst-gallery1
    :glob:

    notebooks/Colon_Cancer_demo.ipynb