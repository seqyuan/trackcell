"""
TrackCell: A Python package for processing single-cell and spatial transcriptomics data.

This package provides tools for:
- IO: Data input/output operations
- PL: Processing and analysis tools
- TL: Utility and helper functions
"""

__version__ = "0.1.4"
__author__ = "Zan Yuan"
__email__ = "yfinddream@gmail.com"

from . import io
from . import pl
from . import tl

__all__ = ["io", "pl", "tl"] 
