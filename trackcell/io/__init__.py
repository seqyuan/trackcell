"""
IO module for TrackCell package.

This module provides functions for reading and writing single-cell and spatial transcriptomics data.
"""

from .read_data import read_hd_bin, read_hd_cellseg, read_xenium_cellseg, sync_geometries_after_subset
from .convert_annohdcell import convert_annohdcell_to_trackcell, add_geometries_to_annohdcell_output, restore_geometries
from .read_sto import read_sto_cellbin, read_sto_bin

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
]