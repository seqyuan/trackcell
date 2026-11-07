"""Tests for plotting helper behavior."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from anndata import AnnData
from matplotlib.collections import LineCollection
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
        "sample": {
            "geometries": gpd.GeoDataFrame(geometry=[_square(0, 0)], index=["c0"])
        }
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


def test_spatial_cell_dual_color_draws_category_edges():
    # dual-color mode: color drives fill, edge_color drives boundaries.
    # The categorical edge pass (``_draw_categorical_edges``) must actually
    # render boundaries as LineCollections on top of the fill.
    adata = AnnData(
        X=np.array([[1.0], [2.0], [3.0]], dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["A", "B", "A"]}, index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5], [1.5, 0.5], [2.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(
                geometry=[_square(0, 0), _square(1, 0), _square(2, 0)],
                index=["c0", "c1", "c2"],
            )
        }
    }

    ax = spatial_cell(
        adata,
        color="g0",
        edge_color="celltype",
        edge_palette={"A": "#e41a1c", "B": "#377eb8"},
        library_id="sample",
        show=False,
        legend=False,
    )

    assert any(isinstance(coll, LineCollection) for coll in ax.collections)
    plt.close(ax.figure)


def test_integer_cluster_column_is_categorical():
    """``leiden`` stored as integers must get a legend, not a viridis gradient.

    ``is_categorical = not is_numeric_dtype(...)`` classified every integer
    column as continuous, which dropped the legend and rendered cluster ids as a
    colour ramp.
    """
    from trackcell.pl.plot import _looks_categorical

    clusters = pd.Series([0, 1, 2, 1, 0], dtype="int8")
    assert _looks_categorical(clusters)

    counts = pd.Series(np.arange(1000, dtype="int64"))
    assert not _looks_categorical(counts)

    labels = pd.Series(["A", "B"], dtype="category")
    assert _looks_categorical(labels)


def test_integer_column_can_be_forced_continuous():
    from trackcell.pl.plot import _looks_categorical

    few_levels = pd.Series([0, 1, 2], dtype="int64")
    assert _looks_categorical(few_levels)
    assert not _looks_categorical(few_levels, max_levels=2)


def test_spatial_cell_categorical_override():
    adata = AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"cluster_int": [0, 1]}, index=["c0", "c1"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5], [2.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(
                geometry=[_square(0, 0), _square(2, 0)], index=["c0", "c1"]
            )
        }
    }

    ax = spatial_cell(adata, color="cluster_int", library_id="sample", show=False)
    assert ax is not None
    assert ax.get_legend() is not None  # categorical -> legend

    # Forcing continuous must render a colorbar instead of a legend.
    ax2 = spatial_cell(
        adata, color="cluster_int", library_id="sample", show=False, categorical=False
    )
    assert ax2 is not None
