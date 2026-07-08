"""Tests for annohdcell conversion helpers."""

import numpy as np
from shapely.geometry import Polygon

from trackcell.io.convert_annohdcell import bins_to_cell_polygon


def test_bins_to_cell_polygon_buffers_colinear_bins():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    polygon = bins_to_cell_polygon(coords, microns_per_pixel=None)

    assert isinstance(polygon, Polygon)
    assert polygon.is_valid
    assert not polygon.is_empty
    assert polygon.area > 0
