"""
TL (Tools) module for TrackCell package.

This module provides utility and helper functions for single-cell and spatial transcriptomics data analysis.
"""

from .spatial import cells_in_radius, hd_labeldist, multigene_blend
from .spatial_cluster import (
    neighborhood_features,
    spatial_cluster,
    spatial_neighbors,
    yard_cluster,
    yard_embed,
)
from ._nhood import (
    connected_components,
    diff_nhood_enrichment,
    nhood_enrichment,
)
from ._autores import (
    select_resolution,
    sweep_resolutions,
)
from ._spatial_graph import remove_long_links

# ``atera_cell_qc`` and ``filter_atera_genes`` are pure computation over an
# AnnData — neither touches the filesystem — so by the same rule that puts
# ``nhood_enrichment`` in ``tl`` and only its renderer in ``pl``, they belong
# here too. Definitions stay in ``io`` (the platform reader owns them) and both
# paths keep working: ``tcl.tl.atera_cell_qc is tcl.io.atera_cell_qc``.
#
# ``atera_transcript_qc`` is deliberately *not* re-exported here: its first
# argument is a run directory, and it streams the point cloud from disk. It is
# IO that happens to end in a QC table, not computation.
from ..io import atera_cell_qc, filter_atera_genes
from ._sct_entry import pack_sct_model, sctransform
from .integration import (
    integrate_rpca,
    integrate_sct_rpca,
    prep_sct_integration,
    run_sct_integration,
    release_sct_integration_cache,
    sct_prep_matrix,
    select_integration_features,
)
from .slice_separation import (
    dbscan_slice_labels,
    distance_to_nearest_centroids,
    mark_colony_centroids,
    slice_cluster_summary,
    spatial_colony_cluster,
    spatial_slice_cluster,
    split_by_slice,
    write_slice_annotation,
)

__all__ = [
    # Re-exported from ``io``: pure computation, no filesystem access.
    "atera_cell_qc",
    "filter_atera_genes",
    "sctransform",
    "pack_sct_model",
    "select_integration_features",
    "prep_sct_integration",
    "run_sct_integration",
    "release_sct_integration_cache",
    "sct_prep_matrix",
    "integrate_rpca",
    "integrate_sct_rpca",
    "hd_labeldist",
    "cells_in_radius",
    "multigene_blend",
    "spatial_cluster",
    "spatial_neighbors",
    "neighborhood_features",
    "yard_embed",
    "yard_cluster",
    "connected_components",
    "nhood_enrichment",
    "diff_nhood_enrichment",
    "select_resolution",
    "sweep_resolutions",
    "remove_long_links",
    "dbscan_slice_labels",
    "spatial_slice_cluster",
    "spatial_colony_cluster",
    "mark_colony_centroids",
    "distance_to_nearest_centroids",
    "slice_cluster_summary",
    "split_by_slice",
    "write_slice_annotation",
]
