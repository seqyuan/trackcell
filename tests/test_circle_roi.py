"""Tests for radius-based ROI selection and its circle annotation."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.pl.plot import mark_circle, spatial_squarebin
from trackcell.tl.spatial import cells_in_radius


def _make_adata(coords, mpp=None, coord_unit=None):
    coords = np.asarray(coords, dtype=float)
    obs = pd.DataFrame(
        {"grp": pd.Categorical(["a"] * len(coords))},
        index=[f"c{i}" for i in range(len(coords))],
    )
    adata = AnnData(X=np.ones((len(coords), 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = coords
    info = {"metadata": {"sample_id": "s"}, "scalefactors": {}}
    if coord_unit is not None:
        info["metadata"]["coord_unit"] = coord_unit
    if mpp is not None:
        info["scalefactors"]["microns_per_pixel"] = mpp
    adata.uns["spatial"] = {"s": info}
    return adata


# ── cells_in_radius ─────────────────────────────────────────────────────────

# Hand-checkable geometry: distances from (0, 0) are 0, 5, 5, 5, 14.14.
_POINTS = [(0, 0), (3, 4), (5, 0), (-3, -4), (10, 10)]


def test_cells_in_radius_includes_boundary_by_default():
    adata = _make_adata(_POINTS, coord_unit="um")
    mask = cells_in_radius(adata, (0.0, 0.0), radius=5)
    assert mask.dtype == bool
    assert mask.tolist() == [True, True, True, True, False]


def test_cells_in_radius_excludes_boundary_when_asked():
    adata = _make_adata(_POINTS, coord_unit="um")
    mask = cells_in_radius(adata, (0.0, 0.0), radius=5, include_boundary=False)
    assert mask.tolist() == [True, False, False, False, False]


def test_cells_in_radius_mask_subsets_anndata():
    adata = _make_adata(_POINTS, coord_unit="um")
    roi = adata[cells_in_radius(adata, (0.0, 0.0), radius=5)].copy()
    assert roi.n_obs == 4
    assert roi.obs_names.tolist() == ["c0", "c1", "c2", "c3"]


def test_cells_in_radius_data_unit_matches_um_for_micron_coords():
    """Atera/Xenium coordinates are µm, so both units must agree."""
    adata = _make_adata(_POINTS, coord_unit="um")
    as_um = cells_in_radius(adata, (0.0, 0.0), radius=5, radius_unit="um")
    as_data = cells_in_radius(adata, (0.0, 0.0), radius=5, radius_unit="data")
    assert as_um.tolist() == as_data.tolist()


def test_cells_in_radius_converts_physical_radius_for_pixel_coords():
    """150 µm is 300 px at 0.5 µm/px — a pixel-blind radius would be half as big."""
    adata = _make_adata([(0, 0), (8, 0), (11, 0)], mpp=0.5)
    mask = cells_in_radius(adata, (0.0, 0.0), radius=5)  # 5 µm == 10 px
    assert mask.tolist() == [True, True, False]

    raw = cells_in_radius(adata, (0.0, 0.0), radius=5, radius_unit="data")
    assert raw.tolist() == [True, False, False]


def test_cells_in_radius_explicit_microns_per_unit_wins():
    adata = _make_adata([(0, 0), (8, 0), (11, 0)], mpp=0.5)
    mask = cells_in_radius(adata, (0.0, 0.0), radius=5, microns_per_unit=1.0)
    assert mask.tolist() == [True, False, False]


def test_cells_in_radius_returns_micron_distances():
    adata = _make_adata([(0, 0), (8, 0)], mpp=0.5)
    mask, dist = cells_in_radius(adata, (0.0, 0.0), radius=5, return_distances=True)
    assert dist.tolist() == pytest.approx([0.0, 4.0])  # 8 px × 0.5 µm/px
    assert mask.tolist() == [True, True]


def test_cells_in_radius_uses_3d_safe_coordinates():
    adata = _make_adata(_POINTS, coord_unit="um")
    # extra obsm columns beyond x/y must be ignored
    adata.obsm["spatial"] = np.column_stack(
        [np.asarray(_POINTS, dtype=float), np.zeros(len(_POINTS))]
    )
    mask = cells_in_radius(adata, (0.0, 0.0), radius=5)
    assert mask.sum() == 4


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"center": (0.0,), "radius": 5}, "`center` must be two"),
        ({"center": (0.0, 0.0), "radius": 0}, "positive finite"),
        ({"center": (0.0, 0.0), "radius": -1}, "positive finite"),
        ({"center": (0.0, 0.0), "radius": 5, "radius_unit": "mm"}, "radius_unit"),
    ],
)
def test_cells_in_radius_validates_inputs(kwargs, match):
    adata = _make_adata(_POINTS, coord_unit="um")
    with pytest.raises(ValueError, match=match):
        cells_in_radius(adata, **kwargs)


def test_cells_in_radius_requires_spatial_basis():
    adata = _make_adata(_POINTS, coord_unit="um")
    del adata.obsm["spatial"]
    with pytest.raises(ValueError, match="obsm"):
        cells_in_radius(adata, (0.0, 0.0), radius=5)


# ── mark_circle ─────────────────────────────────────────────────────────────


def _axes(lim=20):
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, lim)
    ax.set_ylim(lim, 0)  # inverted, like invert_y=True
    ax.set_aspect("equal")
    return fig, ax


def test_mark_circle_uses_data_units_by_default():
    fig, ax = _axes()
    circle = mark_circle(ax, (10.0, 10.0), radius=3, show=False)
    assert circle.radius == pytest.approx(3.0)
    assert circle.get_zorder() == 100
    assert len(ax.patches) == 1
    plt.close(fig)


def test_mark_circle_converts_microns_to_pixels():
    adata = _make_adata(_POINTS, mpp=0.5)
    fig, ax = _axes(lim=40)
    circle = mark_circle(
        ax, (10.0, 10.0), radius=5, radius_unit="um", adata=adata, show=False
    )
    assert circle.radius == pytest.approx(10.0)  # 5 µm / 0.5 µm-per-px
    plt.close(fig)


def test_mark_circle_for_micron_coords_keeps_radius():
    adata = _make_adata(_POINTS, coord_unit="um")
    fig, ax = _axes(lim=4000)
    circle = mark_circle(
        ax, (2000.0, 2000.0), radius=150, radius_unit="um", adata=adata, show=False
    )
    assert circle.radius == pytest.approx(150.0)
    plt.close(fig)


def test_mark_circle_style_options():
    fig, ax = _axes()
    circle = mark_circle(
        ax,
        (10.0, 10.0),
        radius=3,
        edges_color="blue",
        edges_width=3.0,
        linestyle="-",
        fill_color="red",
        fill_alpha=0.3,
        label="r = 3",
        show=False,
    )
    assert circle.get_edgecolor()[3] > 0
    assert circle.get_linewidth() == pytest.approx(3.0)
    assert circle.get_linestyle() == "-"
    assert circle.get_label() == "r = 3"
    assert circle.get_fill()
    plt.close(fig)


def test_mark_circle_rejects_bad_radius_and_unit():
    fig, ax = _axes()
    with pytest.raises(ValueError, match="positive finite"):
        mark_circle(ax, (10.0, 10.0), radius=0, show=False)
    with pytest.raises(ValueError, match="radius_unit"):
        mark_circle(ax, (10.0, 10.0), radius=1, radius_unit="mm", show=False)
    plt.close(fig)


# ── the notebook workflow, end to end ───────────────────────────────────────


def test_radius_selection_and_scalebar_plot_workflow():
    """The documented workflow: overview -> mask -> ROI plot with a scale bar."""
    rng = np.random.default_rng(0)
    coords = rng.uniform(0, 2000, size=(600, 2))
    adata = _make_adata(coords, coord_unit="um")

    center, radius = (1000.0, 1000.0), 200.0
    mask = cells_in_radius(adata, center, radius)

    expected = (np.linalg.norm(coords - np.array(center), axis=1) <= radius).sum()
    assert mask.sum() == expected
    assert 0 < mask.sum() < adata.n_obs

    roi = adata[mask].copy()
    ax = spatial_squarebin(
        roi,
        color="grp",
        show=False,
        legend=False,
        binsize=10,
        scalebar=True,
    )
    labels = [t.get_text() for t in ax.texts]
    # 400 µm wide view -> 0.18 * 400 = 72 -> 50 µm
    assert labels == ["50 µm"]
    plt.close(ax.figure)

    # the same circle drawn on the overview must match the selection
    ax2 = spatial_squarebin(adata, color="grp", show=False, legend=False, binsize=10)
    circle = mark_circle(ax2, center, radius, radius_unit="um", adata=adata, show=False)
    assert circle.radius == pytest.approx(radius)
    plt.close(ax2.figure)
