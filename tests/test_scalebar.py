"""Tests for the physical scale bar used by spatial plots."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.pl.scalebar import (
    add_scale_bar,
    nice_scale_length,
    resolve_microns_per_unit,
)
from trackcell.pl.plot import spatial_squarebin


def _make_adata(mpp=None, coord_unit=None, span=4000.0):
    """Tiny AnnData with a spatial basis of given extent."""
    n = 50
    grid = np.linspace(0.0, span, n)
    xx, yy = np.meshgrid(grid, grid)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    obs = pd.DataFrame(
        {
            "grp": pd.Categorical(["a", "b"] * ((len(coords) + 1) // 2))[: len(coords)],
            "gene": np.linspace(0, 1, len(coords)),
        },
        index=[f"c{i}" for i in range(len(coords))],
    )
    adata = AnnData(X=np.ones((len(coords), 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = coords
    spatial_info = {"metadata": {"sample_id": "s"}, "scalefactors": {}}
    if coord_unit is not None:
        spatial_info["metadata"]["coord_unit"] = coord_unit
    if mpp is not None:
        spatial_info["scalefactors"]["microns_per_pixel"] = mpp
    adata.uns["spatial"] = {"s": spatial_info}
    adata.uns["spatial"]["s"]["binsize"] = 20.0
    return adata


# ── unit resolution ──────────────────────────────────────────────────────────


def test_resolve_explicit_microns_per_unit_wins():
    adata = _make_adata(mpp=0.25, coord_unit="um")
    assert resolve_microns_per_unit(adata, "s", microns_per_unit=2.0) == (2.0, "um")


def test_resolve_um_metadata_hint_ignores_pixel_size():
    """Atera/Xenium coordinates are µm; pixel size must not be applied."""
    adata = _make_adata(mpp=0.2125, coord_unit="um")
    assert resolve_microns_per_unit(adata, "s") == (1.0, "um")


def test_resolve_atera_uns_hint():
    adata = _make_adata()
    adata.uns["atera"] = {"coord_unit": "um", "pixel_size": 0.21}
    assert resolve_microns_per_unit(adata, "s") == (1.0, "um")


def test_resolve_falls_back_to_microns_per_pixel():
    adata = _make_adata(mpp=0.325)
    factor, family = resolve_microns_per_unit(adata, "s")
    assert factor == pytest.approx(0.325)
    assert family == "px"


def test_resolve_defaults_to_microns_when_nothing_known():
    adata = _make_adata()
    assert resolve_microns_per_unit(adata, "s") == (1.0, "um")


def test_resolve_pixel_unit_without_scalefactors_raises():
    adata = _make_adata(coord_unit="px")
    with pytest.raises(ValueError, match="microns_per_pixel"):
        resolve_microns_per_unit(adata, "s")


def test_resolve_rejects_bad_override():
    with pytest.raises(ValueError, match="positive"):
        resolve_microns_per_unit(None, microns_per_unit=0)


# ── length selection ────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "target,expected",
    [(9.0, 5.0), (12.0, 10.0), (720.0, 500.0), (55.0, 50.0), (2.0, 2.0), (1.0, 1.0)],
)
def test_nice_scale_length(target, expected):
    assert nice_scale_length(target) == pytest.approx(expected)


def test_nice_scale_length_rejects_non_positive():
    with pytest.raises(ValueError):
        nice_scale_length(0)


# ── drawing ─────────────────────────────────────────────────────────────────


def _draw(adata, **kwargs):
    fig, ax = plt.subplots(figsize=(4, 4))
    coords = np.asarray(adata.obsm["spatial"])
    ax.scatter(coords[:, 0], coords[:, 1], s=1)
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    ax.set_aspect("equal")
    info = add_scale_bar(ax, adata, library_id="s", **kwargs)
    return fig, ax, info


def test_add_scale_bar_auto_length_microns():
    adata = _make_adata(coord_unit="um")
    fig, ax, info = _draw(adata, length="auto")
    # 4000 µm wide, frac 0.18 -> 720 µm target -> 500 µm nice value
    assert info["length"] == pytest.approx(500.0)
    assert info["label"] == "500 µm"
    assert info["axes_fraction"] == pytest.approx(500.0 / 4000.0)
    plt.close(fig)


def test_auto_length_shrinks_when_zoomed_in():
    """The whole point: smaller visible extent must yield a shorter physical bar."""
    adata = _make_adata(coord_unit="um")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(1000, 1200)  # 200 µm window
    ax.set_ylim(1200, 1000)
    info = add_scale_bar(ax, adata, library_id="s")
    assert info["length"] == pytest.approx(20.0)  # 200 * 0.18 = 36 -> 20
    plt.close(fig)


def test_auto_length_uses_scalefactors_for_pixel_coords():
    adata = _make_adata(mpp=0.5)  # 4000 px -> 2000 µm
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    info = add_scale_bar(ax, adata, library_id="s")
    assert info["length"] == pytest.approx(200.0)  # 2000*0.18=360 -> 200
    assert info["unit"] == "px"  # coordinates are pixels
    assert info["length_unit"] == "um"  # ... but the bar is labelled in µm
    assert info["label"] == "200 µm"
    # 200 µm = 400 px = 10% of the 4000 px width
    assert info["axes_fraction"] == pytest.approx(0.1)
    plt.close(fig)


def test_fixed_length_and_raw_pixel_units():
    adata = _make_adata(coord_unit="um")
    fig, ax, info = _draw(adata, length="750um")
    assert info["length"] == pytest.approx(750.0)
    plt.close(fig)

    fig, ax, info = _draw(adata, length="200px", unit="px")
    assert info["label"] == "200 px"
    assert info["axes_fraction"] == pytest.approx(0.05)
    plt.close(fig)


def test_length_string_without_space_and_micron_symbol():
    adata = _make_adata(coord_unit="um")
    for value in ("50um", "50 µm", "50μm", " 50 UM "):
        fig, ax, info = _draw(adata, length=value)
        assert info["length"] == pytest.approx(50.0), value
        plt.close(fig)


def test_length_auto_accepts_true_and_none():
    adata = _make_adata(coord_unit="um")
    for value in (True, None, "auto"):
        fig, ax, info = _draw(adata, length=value)
        assert info["length"] == pytest.approx(500.0)
        plt.close(fig)


def test_loc_controls_placement():
    adata = _make_adata(coord_unit="um")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    add_scale_bar(ax, adata, library_id="s", length="100um", loc="lower left")
    line = ax.lines[-1]
    xs, ys = line.get_xdata(), line.get_ydata()
    assert xs[0] < 0.2 and ys[0] < 0.2

    add_scale_bar(ax, adata, library_id="s", length="100um", loc="upper right")
    line = ax.lines[-1]
    xs, ys = line.get_xdata(), line.get_ydata()
    assert xs[0] > 0.6 and ys[0] > 0.6
    plt.close(fig)


def test_invalid_loc_raises():
    adata = _make_adata(coord_unit="um")
    with pytest.raises(ValueError, match="Unknown `loc`"):
        _draw(adata, loc="middle")


def test_bar_is_transform_in_axes_coordinates():
    """Axes-fraction placement survives inverted y-axes and tight bbox saving."""
    adata = _make_adata(coord_unit="um")
    fig, ax, info = _draw(adata, length="100um")
    line = ax.lines[-1]
    assert line.get_transform() is not ax.transData
    ys = line.get_ydata()
    assert 0.0 < ys[0] < 0.5  # lower area of the axes, not affected by invert
    plt.close(fig)


def test_style_options_and_frame():
    adata = _make_adata(coord_unit="um")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    add_scale_bar(
        ax,
        adata,
        library_id="s",
        length="100um",
        color="red",
        linewidth=3.0,
        fontsize=12,
        caps=False,
        label="custom",
        label_position="below",
        frame=True,
        frame_color="white",
    )
    line = ax.lines[-1]
    assert line.get_color() == "red"
    assert line.get_linewidth() == pytest.approx(3.0)
    assert line.get_marker() in (None, "None", "")
    texts = [t.get_text() for t in ax.texts]
    assert "custom" in texts
    # frame rectangle is drawn behind the bar
    assert len(ax.patches) == 1
    plt.close(fig)


def test_label_can_be_hidden():
    adata = _make_adata(coord_unit="um")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    add_scale_bar(ax, adata, library_id="s", length="100um", label_position="none")
    assert not ax.texts
    plt.close(fig)


def test_too_long_bar_raises():
    adata = _make_adata(coord_unit="um")
    with pytest.raises(ValueError, match="spans"):
        _draw(adata, length="5000um")


def test_zero_width_axes_raises(monkeypatch):
    adata = _make_adata(coord_unit="um")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(0, 4000)
    monkeypatch.setattr(ax, "get_xlim", lambda: (5.0, 5.0))
    with pytest.raises(ValueError, match="non-zero axis limits"):
        add_scale_bar(ax, adata, library_id="s")
    plt.close(fig)


# ── integration with the plotting functions ─────────────────────────────────


def test_spatial_squarebin_scalebar_flag():
    adata = _make_adata(coord_unit="um")
    ax = spatial_squarebin(
        adata, color="gene", show=False, legend=False, scalebar=True, binsize=20
    )
    labels = [t.get_text() for t in ax.texts]
    assert any(label.endswith("µm") for label in labels)
    plt.close(ax.figure)


def test_spatial_squarebin_scalebar_false_draws_nothing():
    adata = _make_adata(coord_unit="um")
    ax = spatial_squarebin(adata, color="gene", show=False, legend=False, binsize=20)
    assert not ax.texts
    plt.close(ax.figure)


def test_spatial_squarebin_scalebar_dict_and_crop_zoom():
    adata = _make_adata(coord_unit="um")
    ax = spatial_squarebin(
        adata,
        color="gene",
        show=False,
        legend=False,
        binsize=20,
        crop_coord=(2000, 2200, 2000, 2200),  # 200 µm window
        scalebar=dict(length="auto", loc="upper left", color="white", frame=True),
    )
    labels = [t.get_text() for t in ax.texts]
    assert "20 µm" in labels
    plt.close(ax.figure)


def test_scalebar_dict_options_are_forwarded():
    adata = _make_adata(coord_unit="um")
    ax = spatial_squarebin(
        adata,
        color="gene",
        show=False,
        legend=False,
        binsize=20,
        scalebar=dict(length="100um", color="red", linewidth=4.0, offset=(5, 5)),
    )
    bar_lines = [ln for ln in ax.lines if ln.get_marker() in ("|", "None", None, "")]
    assert any(ln.get_color() == "red" for ln in ax.lines)
    assert bar_lines  # bar exists
    plt.close(ax.figure)
