"""
IO module for TrackCell package.

This module provides functions for reading and writing single-cell and spatial
transcriptomics data.
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
    "read_hd_bin": ".read_data",
    "read_hd_cellseg": ".read_data",
    "read_xenium_cellseg": ".read_data",
    "sync_geometries_after_subset": ".read_data",
    "convert_annohdcell_to_trackcell": ".convert_annohdcell",
    "add_geometries_to_annohdcell_output": ".convert_annohdcell",
    "flatten_geometries": ".serialize",
    "restore_geometries": ".serialize",
    "save_h5ad": ".serialize",
    "read_sto_cellbin": ".read_sto",
    "read_sto_bin": ".read_sto",
    "add_wkt_geometry": "._atera",
    "atera_cell_qc": "._atera",
    "atera_info": "._atera",
    "atera_transcript_qc": "._atera",
    "concat_atera_rois": "._atera",
    "filter_atera_genes": "._atera",
    "iter_atera_rois": "._atera",
    "iter_atera_transcripts": "._atera",
    "list_atera_rois": "._atera",
    "read_atera": "._atera",
    "read_atera_cellseg": "._atera",
    "read_atera_panel": "._atera",
    "read_atera_preview": "._atera",
    "read_atera_transcripts": "._atera",
}

__all__ = [
    "read_hd_bin",
    "read_hd_cellseg",
    "read_xenium_cellseg",
    "sync_geometries_after_subset",
    "convert_annohdcell_to_trackcell",
    "add_geometries_to_annohdcell_output",
    "restore_geometries",
    "flatten_geometries",
    "save_h5ad",
    "read_sto_cellbin",
    "read_sto_bin",
    "atera_info",
    "read_atera",
    "read_atera_cellseg",
    "read_atera_preview",
    "read_atera_panel",
    "list_atera_rois",
    "iter_atera_rois",
    "concat_atera_rois",
    "iter_atera_transcripts",
    "read_atera_transcripts",
    "atera_cell_qc",
    "atera_transcript_qc",
    "add_wkt_geometry",
    "filter_atera_genes",
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
