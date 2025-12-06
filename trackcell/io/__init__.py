"""
IO module for TrackCell package.

This module provides functions for reading and writing single-cell and spatial transcriptomics data.
"""

from .read_data import read_hd_bin, read_hd_cellseg

__all__ = ["read_hd_bin", "read_hd_cellseg"]