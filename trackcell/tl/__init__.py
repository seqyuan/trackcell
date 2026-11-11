"""
TL (Tools) module for TrackCell package.

This module provides utility and helper functions for single-cell and spatial
transcriptomics data analysis.
"""

# Public names are resolved on first access (PEP 562) rather than imported here.
#
# Several of these modules pull in the whole single-cell stack: ``read_data``,
# ``read_sto`` and ``convert_annohdcell`` import ``scanpy`` (2.0 s, and it imports
# ``matplotlib.pyplot``, which locks the backend). Importing them eagerly made
# ``import trackcell.tl`` cost the same as importing every reader, because ``tl``
# re-exports two QC helpers that live here.
#
# ``_LAZY_EXPORTS`` maps each public name to the module that defines it;
# ``scripts/measure_import.py --check`` enforces the result. The mapping is
# generated from the previous import statements and guarded by
# ``tests/test_lazy_imports.py``, which asserts every ``__all__`` name resolves to
# the object the old eager import would have produced.

_LAZY_EXPORTS = {
    "cells_in_radius": ".spatial",
    "hd_labeldist": ".spatial",
    "multigene_blend": ".spatial",
    "neighborhood_features": "._cluster",
    "spatial_cluster": "._cluster",
    "spatial_neighbors": "._cluster",
    "yard_cluster": "._cluster",
    "yard_embed": "._cluster",
    "connected_components": "._nhood",
    "diff_nhood_enrichment": "._nhood",
    "nhood_enrichment": "._nhood",
    "select_resolution": "._autores",
    "sweep_resolutions": "._autores",
    "remove_long_links": "._spatial_graph",
    "atera_cell_qc": "..io",
    "filter_atera_genes": "..io",
    "pack_sct_model": "._sct_entry",
    "sctransform": "._sct_entry",
    "integrate_rpca": ".integration",
    "integrate_sct_rpca": ".integration",
    "prep_sct_integration": ".integration",
    "run_sct_integration": ".integration",
    "release_sct_integration_cache": ".integration",
    "sct_prep_matrix": ".integration",
    "select_integration_features": ".integration",
    "dbscan_slice_labels": ".slice_separation",
    "distance_to_nearest_centroids": ".slice_separation",
    "mark_colony_centroids": ".slice_separation",
    "slice_cluster_summary": ".slice_separation",
    "spatial_colony_cluster": ".slice_separation",
    "spatial_slice_cluster": ".slice_separation",
    "split_by_slice": ".slice_separation",
    "write_slice_annotation": ".slice_separation",
}

__all__ = [
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


def __getattr__(name: str):
    try:
        module_name = _LAZY_EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    import importlib

    value = getattr(importlib.import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY_EXPORTS))
