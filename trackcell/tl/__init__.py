"""
TL (Tools) module for TrackCell package.

This module provides utility and helper functions for single-cell and spatial transcriptomics data analysis.
"""

from .spatial import hd_labeldist, multigene_blend
from .spatial_cluster import (
    neighborhood_features,
    spatial_cluster,
    spatial_neighbors,
    yard_cluster,
    yard_embed,
)
from .sctransform import pack_sct_model, sctransform
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
    "multigene_blend",
    "spatial_cluster",
    "spatial_neighbors",
    "neighborhood_features",
    "yard_embed",
    "yard_cluster",
    "dbscan_slice_labels",
    "spatial_slice_cluster",
    "spatial_colony_cluster",
    "mark_colony_centroids",
    "distance_to_nearest_centroids",
    "slice_cluster_summary",
    "split_by_slice",
    "write_slice_annotation",
]
