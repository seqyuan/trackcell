"""Tests for plotting helper behavior."""

import numpy as np
import pytest
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


# ── frameon ──────────────────────────────────────────────────────────────────


def _single_cell_adata(celltype: str = "A") -> AnnData:
    adata = AnnData(
        X=np.ones((1, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": [celltype]}, index=["c0"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.5, 0.5]])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(geometry=[_square(0, 0)], index=["c0"])
        }
    }
    return adata


def test_spatial_cell_frameon_default_keeps_frame():
    adata = _single_cell_adata()
    ax = spatial_cell(adata, color="celltype", library_id="sample", show=False)
    assert ax.axison
    assert ax.spines["left"].get_visible()
    plt.close(ax.figure)


def test_spatial_cell_frameon_false_hides_frame():
    adata = _single_cell_adata()
    ax = spatial_cell(
        adata, color="celltype", library_id="sample", show=False, frameon=False
    )
    assert not ax.axison
    plt.close(ax.figure)


def test_spatial_cell_frameon_false_with_color_none_hides_frame():
    adata = _single_cell_adata()
    ax = spatial_cell(adata, color=None, library_id="sample", show=False, frameon=False)
    assert not ax.axison
    plt.close(ax.figure)


def test_spatial_cell_frameon_none_leaves_user_ax_styling_alone():
    adata = _single_cell_adata()
    fig, ax = plt.subplots()
    ax.axis("off")
    spatial_cell(
        adata, color="celltype", library_id="sample", ax=ax, show=False, legend=False
    )
    assert not ax.axison
    plt.close(fig)


def test_spatial_cell_frameon_is_not_sent_to_geopandas():
    """``frameon`` must not leak into ``**kwargs`` (PatchCollection.set rejects it)."""
    adata = _single_cell_adata()
    ax = spatial_cell(
        adata, color="celltype", library_id="sample", show=False, frameon=False
    )
    assert not ax.axison
    plt.close(ax.figure)


# ── categorical palette: adata.uns[<key>_colors] + legend/fill agreement ─────


def _categorical_adata(labels, colors_key: str = "cell_group") -> AnnData:
    n = len(labels)
    adata = AnnData(
        X=np.ones((n, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {colors_key: pd.Categorical(labels)}, index=[f"c{i}" for i in range(n)]
        ),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[float(i), 0.0] for i in range(n)])
    adata.uns["spatial"] = {
        "sample": {
            "geometries": gpd.GeoDataFrame(
                geometry=[_square(float(i), 0.0) for i in range(n)],
                index=[f"c{i}" for i in range(n)],
            )
        }
    }
    return adata


def _fill_colors(ax, categories):
    """Map each category to the color geopandas actually painted on its polygons."""
    from matplotlib.colors import Normalize, to_hex

    coll = next(
        c
        for c in ax.collections
        if hasattr(c, "get_array") and c.get_array() is not None
    )
    cmap = coll.get_cmap()
    lo, hi = coll.get_clim()
    norm = Normalize(lo, hi)
    return {cat: to_hex(cmap(norm(i))) for i, cat in enumerate(categories)}


def _legend_colors(ax):
    from matplotlib.colors import to_hex

    legend = ax.get_legend()
    handles = getattr(legend, "legend_handles", None)
    if handles is None:  # matplotlib < 3.8
        handles = legend.legendHandles
    return {h.get_label(): to_hex(h.get_facecolor()) for h in handles}


def test_spatial_cell_uses_uns_colors_when_palette_is_none():
    adata = _categorical_adata(["A", "B", "C"])
    adata.uns["cell_group_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    expected = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}
    assert _fill_colors(ax, ["A", "B", "C"]) == expected
    assert _legend_colors(ax) == expected
    plt.close(ax.figure)


def test_spatial_cell_uns_colors_follow_uns_categories():
    """For a plain object column the order has to come from ``uns[..._categories]``."""
    adata = _categorical_adata(["A", "B"])
    adata.obs["cell_group"] = adata.obs["cell_group"].astype(object)  # not categorical
    adata.uns["cell_group_categories"] = ["B", "A"]
    adata.uns["cell_group_colors"] = ["#377eb8", "#e41a1c"]

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    assert _fill_colors(ax, ["A", "B"]) == {"A": "#e41a1c", "B": "#377eb8"}
    plt.close(ax.figure)


def test_spatial_cell_explicit_palette_beats_uns_colors():
    adata = _categorical_adata(["A", "B"])
    adata.uns["cell_group_colors"] = ["#e41a1c", "#377eb8"]

    ax = spatial_cell(
        adata,
        color="cell_group",
        library_id="sample",
        show=False,
        palette={"A": "#000000", "B": "#ffffff"},
    )

    assert _fill_colors(ax, ["A", "B"]) == {"A": "#000000", "B": "#ffffff"}
    plt.close(ax.figure)


def test_spatial_cell_ignores_uns_colors_with_wrong_length():
    adata = _categorical_adata(["A", "B"])
    adata.uns["cell_group_colors"] = ["#e41a1c"]

    with pytest.warns(
        UserWarning, match="Ignoring `adata.uns\\['cell_group_colors'\\]`"
    ):
        ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)
    plt.close(ax.figure)


def test_spatial_cell_legend_matches_fill_without_palette():
    labels = ["A", "B", "C", "D", "E"]
    adata = _categorical_adata(labels)

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    fills = _fill_colors(ax, labels)
    assert _legend_colors(ax) == fills
    plt.close(ax.figure)


def test_spatial_cell_default_palette_does_not_collapse_past_ten_categories():
    labels = [f"g{i}" for i in range(15)]
    adata = _categorical_adata(labels)

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    fills = _fill_colors(ax, labels)
    assert len(set(fills.values())) == 15  # tab10 would have given 5 duplicates
    plt.close(ax.figure)


# ── spatial_squarebin parity ────────────────────────────────────────────────


def test_spatial_squarebin_frameon_and_uns_palette():
    from trackcell.pl.plot import spatial_squarebin

    coords = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0]])
    adata = AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"cell_group": pd.Categorical(["A", "B", "C"])}, index=["b0", "b1", "b2"]
        ),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = coords
    adata.uns["spatial"] = {"sample": {"binsize": 2.0}}
    adata.uns["cell_group_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]

    ax = spatial_squarebin(
        adata,
        color="cell_group",
        library_id="sample",
        show=False,
        frameon=False,
    )
    assert not ax.axison
    facecolors = ax.collections[0].get_facecolors()
    from matplotlib.colors import to_hex

    assert to_hex(facecolors[0]) == "#e41a1c"
    assert to_hex(facecolors[2]) == "#4daf4a"
    plt.close(ax.figure)


def test_spatial_cell_uses_uns_colors_stored_as_dict():
    """``io.read_atera(load_groups=True)`` stores ``{group: hex}``, not a list."""
    adata = _categorical_adata(["A", "B", "C"])
    adata.uns["cell_group_colors"] = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    expected = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}
    assert _fill_colors(ax, ["A", "B", "C"]) == expected
    assert _legend_colors(ax) == expected
    plt.close(ax.figure)


def test_spatial_cell_color_cell_group_picks_up_atera_palette(tmp_path):
    """End-to-end: Atera's own palette + ``color="cell_group"`` needs no palette=."""
    from trackcell.io.read_atera import _attach_cell_groups

    labels = ["Tumour", "Stroma", "Tumour"]
    adata = _categorical_adata(labels)
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs = adata.obs.drop(columns=["cell_group"])
    pd.DataFrame(
        [
            ("c0", "Tumour", "#ff0000"),
            ("c1", "Stroma", "#00ff00"),
            ("c2", "Tumour", "#ff0000"),
        ],
        columns=["cell_id", "group", "color"],
    ).to_csv(tmp_path / "cell_groups.csv", index=False)
    assert _attach_cell_groups(adata, tmp_path) is True

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    expected = {"Tumour": "#ff0000", "Stroma": "#00ff00"}
    assert _fill_colors(ax, ["Stroma", "Tumour"]) == expected
    assert _legend_colors(ax) == expected
    plt.close(ax.figure)


def test_spatial_cell_categorical_with_missing_values():
    """Unassigned cells (NaN) must not break the palette or the legend."""
    labels = ["A", "B", None, "A"]
    adata = _categorical_adata(labels)
    adata.obs["cell_group"] = pd.Categorical(labels)

    ax = spatial_cell(adata, color="cell_group", library_id="sample", show=False)

    assert [h.get_label() for h in ax.get_legend().legend_handles] == ["A", "B"]
    assert len(set(_fill_colors(ax, ["A", "B"]).values())) == 2
    plt.close(ax.figure)


def test_palette_from_uns_survives_subset_and_copy():
    """AnnData packs a ``.uns`` dict into a 0-d object array on subset/copy.

    ``adata[mask].copy().uns['cell_group_colors']`` is then ``array([{...}], dtype=object)``
    — length 1, not a ``Mapping``. Reading it naively made the whole palette fall back
    to the default one after every ``adata[mask]`` (i.e. in every ROI panel).
    """
    from trackcell.pl.plot import _palette_from_uns

    adata = _categorical_adata(["A", "B", "C"])
    adata.uns["cell_group_colors"] = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}
    sub = adata[np.array([True, True, True])].copy()
    assert not isinstance(sub.uns["cell_group_colors"], dict)  # 0-d object array

    assert _palette_from_uns(sub, "cell_group") == {
        "A": "#e41a1c",
        "B": "#377eb8",
        "C": "#4daf4a",
    }

    ax = spatial_cell(sub, color="cell_group", library_id="sample", show=False)
    expected = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a"}
    assert _fill_colors(ax, ["A", "B", "C"]) == expected
    plt.close(ax.figure)


def test_palette_from_uns_survives_anndata_category_pruning():
    """``adata[mask].copy()`` prunes ``uns['<key>_colors']`` to the surviving groups.

    ``annData`` shrinks the color list in lock-step with ``obs[key].cat.categories``
    but leaves ``uns['<key>_categories']`` untouched. Aligning against that stale
    list made every ROI panel (``adata[mask].copy()``) fall back to the default
    palette, so each subset got different colors for the same group.
    """
    from trackcell.pl.plot import _palette_from_uns

    labels = ["A", "B", "C", "A"]
    adata = _categorical_adata(labels)
    adata.uns["cell_group_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]
    adata.uns["cell_group_categories"] = ["A", "B", "C"]

    sub = adata[np.array([True, False, True, True])].copy()  # A and C survive

    assert list(sub.obs["cell_group"].cat.categories) == ["A", "C"]
    assert list(sub.uns["cell_group_colors"]) == ["#e41a1c", "#4daf4a"]  # pruned
    assert _palette_from_uns(sub, "cell_group") == {"A": "#e41a1c", "C": "#4daf4a"}

    ax = spatial_cell(sub, color="cell_group", library_id="sample", show=False)
    assert _fill_colors(ax, ["A", "C"]) == {"A": "#e41a1c", "C": "#4daf4a"}
    plt.close(ax.figure)
