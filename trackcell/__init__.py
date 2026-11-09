"""
TrackCell: A Python package for processing single-cell and spatial transcriptomics data.

This package provides tools for:
- IO: Data input/output operations
- TL: Processing and analysis tools
- PL: Plotting and visualization tools
"""

__version__ = "0.4.4"
__author__ = "Zan Yuan"
__email__ = "yfinddream@gmail.com"

from . import _logging  # noqa: F401  (imported for its NullHandler side effect)
from . import io
from . import pl
from . import tl

__all__ = ["io", "pl", "tl"]
