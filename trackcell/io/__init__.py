"""
IO module for TrackCell package.

This module provides functions for reading and writing single-cell and spatial transcriptomics data.
"""

from .read_data import read_hd_bin, read_hd_cellseg, read_xenium_cellseg, sync_geometries_after_subset
from .convert_annohdcell import convert_annohdcell_to_trackcell, add_geometries_to_annohdcell_output, restore_geometries
from .read_sto import read_sto_cellbin, read_sto_bin
from .read_atera import (
    add_wkt_geometry,
    atera_cell_qc,
    atera_info,
    atera_transcript_qc,
    concat_atera_rois,
    filter_atera_genes,
    iter_atera_rois,
    iter_atera_transcripts,
    list_atera_rois,
    read_atera,
    read_atera_cellseg,
    read_atera_panel,
    read_atera_preview,
    read_atera_transcripts,
)

__all__ = [
    "read_hd_bin",
    "read_hd_cellseg",
    "read_xenium_cellseg",
    "sync_geometries_after_subset",
    "convert_annohdcell_to_trackcell",
    "add_geometries_to_annohdcell_output",
    "restore_geometries",
    "read_sto_cellbin",
    "read_sto_bin",
    # Atera
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
