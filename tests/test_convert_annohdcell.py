"""Tests for annohdcell conversion helpers."""

import anndata as ad
import numpy as np
from shapely.geometry import Polygon

from trackcell.io.convert_annohdcell import (
    add_geometries_to_annohdcell_output,
    bins_to_cell_polygon,
)


def test_bins_to_cell_polygon_buffers_colinear_bins():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    polygon = bins_to_cell_polygon(coords, microns_per_pixel=None)

    assert isinstance(polygon, Polygon)
    assert polygon.is_valid
    assert not polygon.is_empty
    assert polygon.area > 0


def _write_inputs(tmp_path):
    """Two h5ads with the minimum schema ``add_geometries`` needs."""
    # Two bins per cell, one unassigned bin.
    bin_coords = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [10.0, 0.0],
            [12.0, 0.0],
            [20.0, 0.0],
            [22.0, 0.0],
            [40.0, 0.0],
        ]
    )
    labels = [1, 1, 2, 2, 3, 3, 0]
    bins = ad.AnnData(X=np.ones((len(labels), 2), dtype=np.float32))
    bins.obs_names = [f"bin{i}" for i in range(len(labels))]
    bins.var_names = ["g0", "g1"]
    bins.obs["labels_joint"] = labels
    bins.obsm["spatial"] = bin_coords
    bins.uns["spatial"] = {"sample": {"scalefactors": {"microns_per_pixel": 0.5}}}
    bin_path = tmp_path / "b2c_2um.h5ad"
    bins.write_h5ad(bin_path)

    cells = ad.AnnData(X=np.ones((3, 2), dtype=np.float32))
    cells.obs_names = ["cell1", "cell2", "cell3"]
    cells.var_names = ["g0", "g1"]
    cells.obs["object_id"] = [1, 2, 3]
    cell_path = tmp_path / "b2c_cell.h5ad"
    cells.write_h5ad(cell_path)
    return bin_path, cell_path


def test_add_geometries_serial_and_parallel_agree(tmp_path):
    bin_path, cell_path = _write_inputs(tmp_path)

    serial = add_geometries_to_annohdcell_output(
        bin_path, cell_path, sample="sample", n_jobs=1
    )
    parallel = add_geometries_to_annohdcell_output(
        bin_path, cell_path, sample="sample", n_jobs=2
    )

    assert serial.obs["geometry"].notna().all()
    assert parallel.obs["geometry"].notna().all()
    # Every cell keeps its own polygon; the worker state must not mix them.
    assert (
        serial.uns["spatial"]["sample"]["geometries"]
        .geometry.apply(lambda geom: geom.is_valid and not geom.is_empty)
        .all()
    )
    assert (
        serial.obs["geometry"].tolist() == parallel.obs["geometry"].tolist()
    ), "parallel workers returned different polygons than the serial path"


def test_add_geometries_with_all_cells_unassigned_is_empty_not_an_error(tmp_path):
    bin_path, cell_path = _write_inputs(tmp_path)
    bins = ad.read_h5ad(bin_path)
    bins.obs["labels_joint"] = 0
    bins.write_h5ad(bin_path)

    result = add_geometries_to_annohdcell_output(
        bin_path, cell_path, sample="sample", n_jobs=2
    )

    assert result.obs["geometry"].isna().all()
    assert len(result.uns["spatial"]["sample"]["geometries"]) == 0
