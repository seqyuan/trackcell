"""Tests for plotting helper behavior."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from anndata import AnnData
from shapely.geometry import Polygon

from trackcell.pl.plot import spatial_cell


def _square(x: float, y: float) -> Polygon:
    return Polygon([(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)])


def test_spatial_cell_syncs_geometries_to_obs_names():
    adata = AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A", "B"]}, index=["c0", "c1"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5], [2.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(
                geometry=[_square(5, 5), _square(0, 0), _square(2, 0)],
                index=["extra", "c0", "c1"],
            )
        }
    }

    ax = spatial_cell(
        adata,
        color="celltype",
        library_id="sample",
        show=False,
        legend=False,
    )

    synced = adata.uns["spatial"]["sample"]["geometries"]
    assert list(synced.index) == ["c0", "c1"]
    plt.close(ax.figure)


def test_spatial_cell_does_not_tight_layout_user_axes(monkeypatch):
    adata = AnnData(
        X=np.ones((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A"]}, index=["c0"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {"geometries": gpd.GeoDataFrame(geometry=[_square(0, 0)], index=["c0"])}
    }
    fig, ax = plt.subplots()

    def fail_tight_layout(*args, **kwargs):
        raise AssertionError("tight_layout should not be called for user axes")

    monkeypatch.setattr(fig, "tight_layout", fail_tight_layout)
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)

    spatial_cell(
        adata,
        color="celltype",
        library_id="sample",
        ax=ax,
        show=True,
        legend=False,
    )
    plt.close(fig)


def test_spatial_cell_flips_colored_background_when_invert_y_false():
    """Colored background image with invert_y=False must be vertically flipped."""
    img = np.arange(12, dtype=np.float32).reshape(3, 4)
    # Polygon that covers the whole 4×3 image area so the crop is a no-op
    full_poly = Polygon([(0, 0), (4, 0), (4, 3), (0, 3)])
    adata = AnnData(
        X=np.ones((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A"]}, index=["c0"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[2.0, 1.5]])
    adata.uns["spatial"] = {
        "sample": {
            "images": {"hires": img},
            "scalefactors": {"tissue_hires_scalef": 1.0},
            "geometries": gpd.GeoDataFrame(geometry=[full_poly], index=["c0"]),
        }
    }

    ax = spatial_cell(
        adata,
        color="celltype",
        library_id="sample",
        invert_y=False,
        show=False,
        legend=False,
    )

    rendered_background = np.asarray(ax.images[0].get_array())
    np.testing.assert_array_equal(rendered_background, img[::-1])
    plt.close(ax.figure)
