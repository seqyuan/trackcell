"""Tests for geometry serialization helpers (h5ad round-trip)."""

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon

import trackcell as tcl


def _square(x: float, y: float) -> Polygon:
    return Polygon([(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)])


def _make_adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A", "B"]}, index=["c0", "c1"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5], [2.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(
                geometry=[_square(0, 0), _square(2, 0)], index=["c0", "c1"]
            ),
            "nucleus_geometries": gpd.GeoDataFrame(
                geometry=[_square(0.2, 0.2), _square(2.2, 0.2)],
                index=["c0", "c1"],
            ),
            "scalefactors": {"spot_diameter_fullres": 1.0},
        }
    }
    return adata


def test_save_h5ad_roundtrip_restores_geometries(tmp_path):
    adata = _make_adata()
    path = tcl.io.save_h5ad(adata, tmp_path / "sample.h5ad")

    # In-memory objects are restored after the write.
    assert isinstance(
        adata.uns["spatial"]["sample"]["geometries"], gpd.GeoDataFrame
    )
    assert isinstance(
        adata.uns["spatial"]["sample"]["nucleus_geometries"], gpd.GeoDataFrame
    )

    restored = tcl.io.restore_geometries(ad.read_h5ad(path))
    expected_bounds = {
        "geometries": (0, 0, 1, 1),
        "nucleus_geometries": (0.2, 0.2, 1.2, 1.2),
    }
    for key, bounds in expected_bounds.items():
        gdf = restored.uns["spatial"]["sample"][key]
        assert isinstance(gdf, gpd.GeoDataFrame)
        assert list(gdf.index) == ["c0", "c1"]
        np.testing.assert_allclose(gdf.loc["c0"].geometry.bounds, bounds)


def test_flatten_geometries_copy_does_not_mutate_input():
    adata = _make_adata()
    flat = tcl.io.flatten_geometries(adata)
    assert isinstance(
        flat.uns["spatial"]["sample"]["geometries"], pd.DataFrame
    )
    assert isinstance(
        adata.uns["spatial"]["sample"]["geometries"], gpd.GeoDataFrame
    )


def test_spatial_cell_reads_restored_geometries(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    adata = _make_adata()
    tcl.io.save_h5ad(adata, tmp_path / "sample.h5ad")
    restored = tcl.io.restore_geometries(ad.read_h5ad(tmp_path / "sample.h5ad"))
    ax = tcl.pl.spatial_cell(
        restored, color="celltype", library_id="sample", show=False
    )
    assert ax is not None
