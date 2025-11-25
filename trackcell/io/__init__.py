"""
IO module for TrackCell package.

This module provides functions for reading and writing single-cell and spatial transcriptomics data.
"""

from .spatial import read_hd_cellseg

__all__ = ["read_hd_cellseg"]