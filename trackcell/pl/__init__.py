"""
PL (Plotting) module for TrackCell package.

This module provides functions for plotting and visualizing single-cell and spatial transcriptomics data.
"""

from .plot import spatial_cell, spatial_squarebin, spatial_bin, mark_region

# Lazy import for napari functions (heavy optional dependency)
def __getattr__(name):
    if name in ("napari_view", "select_regions", "napari_extract"):
        from . import napari as _napari
        return getattr(_napari, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "spatial_cell",
    "spatial_squarebin",
    "spatial_bin",
    "mark_region",
    "napari_view",
    "select_regions",
    "napari_extract",
] 