"""
TrackCell: single-cell and spatial transcriptomics toolkit.

Scope
-----
The **reader + spatial-analysis + spatial-plotting** layer for 10x Visium HD,
10x Xenium, 10x Atera and STOmics (Stereo-seq) data. It sits on top of
scanpy/AnnData and complements it: use scanpy for QC, PCA, UMAP, HVG and
plain Leiden, and reach for ``trackcell`` when the task involves spatial
coordinates, cell polygons, physical distances, or platform-specific formats.

Routing table — reach for these instead of hand-rolling
--------------------------------------------------------
============================================  ====================================
Task                                          Use
============================================  ====================================
Read Visium HD bin / cell segmentation        ``io.read_hd_bin`` / ``io.read_hd_cellseg``
Read Xenium                                   ``io.read_xenium_cellseg``
Read Atera (auto-detects bundle flavour)      ``io.read_atera``
Read STOmics GEF                              ``io.read_sto_cellbin`` / ``io.read_sto_bin``
List / stream / merge Atera ROIs              ``io.list_atera_rois`` / ``iter_atera_rois`` / ``concat_atera_rois``
Synchronise polygons after subsetting         ``io.sync_geometries_after_subset``
Distance to nearest cell of a given label     ``tl.hd_labeldist``
Cells inside a circle of a physical radius    ``tl.cells_in_radius``
Microns per coordinate unit                   ``pl.resolve_microns_per_unit``
Split one slide into tissue sections          ``tl.spatial_slice_cluster``
Spatial-domain clustering (YardCluster)       ``tl.spatial_cluster``
Pick a stable Leiden resolution               ``tl.select_resolution``
Neighbourhood enrichment + 2-condition diff   ``tl.nhood_enrichment`` / ``tl.diff_nhood_enrichment``
Multi-gene co-expression blend                ``tl.multigene_blend`` then ``pl.multigene_blend``
Cell-level QC (Atera)                         ``io.atera_cell_qc`` (also ``tl.atera_cell_qc``)
Transcript QC (streams disk, so it is IO)     ``io.atera_transcript_qc``
SCTransform in Python (Seurat 4.4.0 equiv.)   ``tl.sctransform``
SCT + RPCA integration                        ``tl.integrate_sct_rpca``
Plot cell polygons (not scatter)              ``pl.spatial_cell``
Scale bar / ROI marks                         ``pl.add_scale_bar`` / ``pl.mark_circle``
Interactive region selection (Jupyter)        ``pl.select_regions``
============================================  ====================================

Data contract of every ``io.read_*`` result
-------------------------------------------
``.X`` holds raw counts. ``.obs`` holds platform annotation plus a WKT
``geometry`` column on polygon platforms. ``**.obsm['spatial']**`` holds (n, 2)
coordinates — **image pixels** for Visium/Visium HD, **micrometres** for
Xenium/Atera/MERSCOPE. ``.uns['spatial'][sample]`` holds ``images``,
``scalefactors`` (with ``microns_per_pixel``) and ``geometries`` (a
GeoDataFrame). Every ``tl.*``/``pl.*`` function depends only on those keys plus
scanpy conventions, so any AnnData with ``.obsm['spatial']`` can be passed in.

Footguns
--------
* Subsetting invalidates the geometry table — always follow
  ``adata = adata[mask].copy()`` with ``io.sync_geometries_after_subset(adata)``.
  Skipping this silently misplaces every polygon downstream.
* Coordinates are not uniformly metric; pass ``radius_unit='um'`` rather than
  treating a pixel count as a distance.
* Before 0.5.0, ``tl.sctransform`` and ``io.read_atera`` shadowed same-named
  modules, so ``import trackcell.tl.sctransform as m`` handed back a *function*.
  The modules are now ``tl._sct_entry`` and ``io._atera``, and those paths fail
  loudly. Use ``tcl.tl.sctransform(...)`` or ``from trackcell.tl import sctransform``.
* Optional features need extras: ``[clustering]``, ``[integration]``, ``[sct]``,
  ``[ml]``, ``[interactive]``.

Modules
-------
- ``io``: reading, writing, geometry synchronisation
- ``tl``: computation (in-place by default, ``copy=True`` returns a copy)
- ``pl``: plotting (``show=True`` by default, returns ``Axes``)
"""

__version__ = "0.5.1"
__author__ = "Zan Yuan"
__email__ = "yfinddream@gmail.com"

from . import _logging  # noqa: F401  (imported for its NullHandler side effect)
from . import io
from . import pl
from . import tl

__all__ = ["io", "pl", "tl"]
