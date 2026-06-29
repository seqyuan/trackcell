Xenium TMA Slice Separation (DBSCAN)
====================================

TrackCell provides **DBSCAN-based spatial slice separation** for 10x Xenium
TMA (tissue microarray) data where multiple tissue cores share one coordinate
system. This is a **geometric preprocessing** step: it splits cells into
physical tissue sections (``S001``, ``S002``, …) before expression-based
clustering (YardCluster, DECIPHER, scanpy).

Use case
--------

* Xenium cassette with many FFPE cores (e.g. 3×7 TMA grid)
* One Xenium **region** directory contains all cores merged in ``cells.parquet``
* Cores are separated by empty space → density clustering on centroids works well

This differs from :doc:`spatial_clustering` (YardCluster), which clusters by
**gene expression + neighborhood** to find cell types or tissue domains *within*
a section.

Quick start
-----------

.. code-block:: python

   import trackcell as tcl

   # Option A: slice separation during load
   adata = tcl.io.read_xenium_cellseg(
       "/path/to/xenium/region",
       sample="85811_S",
       slice_separate=True,
       slice_eps=80,
   )

   # Option B: separate step (works on any AnnData with spatial coords)
   adata = tcl.io.read_xenium_cellseg("/path/to/xenium/region")
   tcl.tl.spatial_slice_cluster(adata, eps=80, min_samples=10)

   # Summarize and split
   summary = tcl.tl.slice_cluster_summary(adata)
   slices = tcl.tl.split_by_slice(adata)  # dict: "S001" -> AnnData

Algorithm
---------

1. Read cell centroids ``(x_centroid, y_centroid)`` or ``obsm['spatial']``.
2. Run ``sklearn.cluster.DBSCAN(eps, min_samples)`` in 2D.
3. Clusters with fewer than ``min_cells`` (default 1000) → ``debris``.
4. Sort remaining clusters by mean Y, then mean X.
5. Assign global slice IDs ``S001``–``SNNN``.

Default parameters match the CCHD cervical cancer Xenium project (FFPE TMA,
pixel size 0.2125 µm, core spacing >100 µm):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Notes
   * - ``eps``
     - 80 µm
     - Increase for fused/close cores; 90 µm for scattered layouts
   * - ``min_samples``
     - 10
     - DBSCAN density threshold
   * - ``min_cells``
     - 1000
     - Filter small fragments / noise

API reference
-------------

``tcl.tl.spatial_slice_cluster(adata, ...)``
    Main entry point. Adds ``obs['slice_id']``, ``obs['slice_cluster']``,
    and ``uns['slice_id_params']`` / ``uns['slice_id_summary']``.

``tcl.tl.dbscan_slice_labels(coords, ...)``
    Low-level function on an N×2 coordinate array. Returns cluster integers
    and string slice IDs.

``tcl.tl.slice_cluster_summary(adata)``
    DataFrame with per-slice cell counts and spatial bounding boxes.

``tcl.tl.split_by_slice(adata)``
    Split into a dict of per-slice AnnData objects; reindexes
    ``uns['cell_boundaries']`` when present.

``tcl.tl.write_slice_annotation(adata, path)``
    Export CCHD-compatible ``slice_annotation.parquet``.

Multi-region TMA workflow
-------------------------

For a full Xenium run with multiple regions (e.g. ``yuanfa``, ``fensan``,
``ronghe``), process each region separately and use a running ``slice_start``
to continue global numbering across regions:

.. code-block:: python

   import trackcell as tcl

   regions = [
       ("85811", "S", "/path/to/fensan", 11),   # slice_start after prior batch
   ]
   next_start = 1
   for slide_id, region, path, _ in regions:
       adata = tcl.io.read_xenium_cellseg(path, sample=f"{slide_id}_{region}")
       tcl.tl.spatial_slice_cluster(adata, slice_start=next_start)
       next_start += adata.uns["slice_id_summary"]["n_slices"]
       adata.write_h5ad(f"{slide_id}_{region}.h5ad")

   # Per-slice files for downstream DECIPHER / YardCluster
   for sid, sub in tcl.tl.split_by_slice(adata).items():
       sub.write_h5ad(f"slices/{sid}.h5ad")

Visualization
-------------

After slice separation, color cells by ``slice_id``:

.. code-block:: python

   tcl.pl.spatial_cell(adata, color="slice_id", figsize=(12, 12))

Combine with YardCluster on a single slice:

.. code-block:: python

   sub = tcl.tl.split_by_slice(adata)["S011"]
   tcl.tl.spatial_cluster(sub, mode="auto")
   tcl.pl.spatial_cell(sub, color="yardcluster_domain")

See also
--------

* :doc:`reading` — ``read_xenium_cellseg`` format details
* :doc:`spatial_clustering` — YardCluster expression + spatial domain clustering

Micro-region colonies (GC / tumor nests)
----------------------------------------

CCHD uses a **second DBSCAN workflow** on **BANKSY-filtered subsets** within
one tissue section — not on the whole cassette. This separates dispersed
**germinal centers (GC)** or **tumor nests** that share the same cell type
label but are spatially disconnected.

Workflow comparison
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Step
     - TMA slice separation
     - GC / tumor nest colonies
   * - Input
     - All cells in a Xenium region
     - Subset after BANKSY (e.g. ``Germinal Center B Cells``, tumor ROI clusters)
   * - Goal
     - Physical cores → ``S001``–``SNNN``
     - Dispersed colonies → ``spatial_colony`` IDs
   * - ``eps``
     - 80–90 µm (Xenium centroids)
     - 50 µm (Xenium GC); 14 (Visium HD bin coords)
   * - ``min_samples``
     - 10
     - 5 (Xenium GC); 10 (Visium tumor/GC)
   * - Size filter
     - ``min_cells=1000`` → ``debris``
     - ``min_cluster_size=40`` → ``NaN``

CCHD notebook sources:

* ``54_1_GC_dist.ipynb`` — Xenium GC colonies (Python)
* ``17.banksy_GCs_identify.ipynb`` — Visium GC nests (R, eps=14)
* ``16.S_P_abnksy_tumor_nest_identify.ipynb`` — Visium tumor nests (R)
* ``dbscan.ipynb`` — k-distance plot for ``eps`` tuning

GC colony example (Xenium)
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import trackcell as tcl

   adata = tcl.io.read_xenium_cellseg("/path/to/xenium/region")
   tcl.tl.spatial_slice_cluster(adata)
   sub = tcl.tl.split_by_slice(adata)["S011"]

   # BANKSY domain clustering first (see spatial_clustering.rst)
   tcl.tl.spatial_cluster(sub, mode="auto")
   gc = sub[sub.obs["yardcluster_domain"] == "Germinal Center B Cells"].copy()

   tcl.tl.spatial_colony_cluster(gc, eps=50, min_samples=5, min_cluster_size=40)
   tcl.tl.mark_colony_centroids(gc, centroid_label="GCC")

   # Distance on full slice for downstream plots
   tcl.tl.distance_to_nearest_centroids(
       sub, centroid_key="cell_centroid_type", distance_key="distance_to_nearest_gcc",
       centroid_label="GCC",
   )
   tcl.pl.spatial_cell(gc, color="spatial_colony")

Visium tumor / GC nests
~~~~~~~~~~~~~~~~~~~~~~~

On Visium HD, CCHD crops a BANKSY tumor or GC ROI, then runs DBSCAN on bin
coordinates with **smaller** ``eps`` (default **14**, ``min_samples=10``).
Tune ``eps`` with a k-distance elbow plot (15th sorted k-NN distance in
``dbscan.ipynb``). Map results back to full tissue for distance-to-border
analysis.

.. code-block:: python

   # After BANKSY on Visium HD bins
   roi = adata[adata.obs["banksy"].isin(["4", "6"])].copy()
   tcl.tl.spatial_colony_cluster(roi, eps=14, min_samples=10, min_cluster_size=20)
   tcl.tl.mark_colony_centroids(roi, centroid_label="TNC")

Parameter tuning tips
~~~~~~~~~~~~~~~~~~~~~

* **Increase ``eps``** when nearby colonies merge; **decrease** when one GC
  splits into fragments.
* **``min_cluster_size``** removes small DBSCAN blobs; GC notebook uses 40
  cells; adjust for sparser platforms (Visium bins).
* Always **subset by cell type / domain first** — raw DBSCAN on all cells
  finds tissue patches, not biological colonies.
* Units follow ``obsm['spatial']`` (µm for Xenium centroids; Visium HD bin
  positions in project coordinates).
