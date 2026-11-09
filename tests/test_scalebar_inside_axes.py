"""The scale bar must stay inside the axes box, caps and label included.

Regression guard for a real report: with ``loc="lower right"`` the caps and the
label could spill over the spines, landing on top of the colourbar or outside the
figure after ``bbox_inches="tight"``. The bar is drawn with ``clip_on=False`` in
axes coordinates, so nothing in matplotlib stops it.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.pl.scalebar import _LOC_FACTORS, add_scale_bar


def _make_adata(span=4000.0, coord_unit="um"):
    coords = np.array([[0.0, 0.0], [span, span]])
    obs = pd.DataFrame({"g": [0, 1]}, index=["c0", "c1"])
    adata = AnnData(X=np.ones((2, 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = coords
    adata.uns["spatial"] = {
        "s": {"metadata": {"coord_unit": coord_unit}, "scalefactors": {}}
    }
    return adata


def _extents(ax, line, caps_pt, linewidth_pt):
    """Pixel bboxes of the bar (caps included) and of every label.

    ``Line2D.get_window_extent`` pads by a *square* of ``markersize``, which
    overstates a ``|`` cap by an order of magnitude, so the bar box is computed
    from the transformed endpoints and expanded by the real stroke sizes.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_px = fig.dpi / 72.0

    pts = line.get_transform().transform(
        np.column_stack([line.get_xdata(), line.get_ydata()])
    )
    x0, x1 = float(pts[:, 0].min()), float(pts[:, 0].max())
    y0, y1 = float(pts[:, 1].min()), float(pts[:, 1].max())
    bar = (
        x0 - linewidth_pt * to_px / 2.0,
        y0 - caps_pt * to_px / 2.0,
        x1 + linewidth_pt * to_px / 2.0,
        y1 + caps_pt * to_px / 2.0,
    )
    labels = [t.get_window_extent(renderer) for t in ax.texts]
    return [bar] + [(t.x0, t.y0, t.x1, t.y1) for t in labels]


def _assert_inside(ax, extents, tol=1.0):
    box = ax.bbox
    for x0, y0, x1, y1 in extents:
        assert x0 >= box.x0 - tol, f"left overflow: {x0} < {box.x0}"
        assert x1 <= box.x1 + tol, f"right overflow: {x1} > {box.x1}"
        assert y0 >= box.y0 - tol, f"bottom overflow: {y0} < {box.y0}"
        assert y1 <= box.y1 + tol, f"top overflow: {y1} > {box.y1}"


def _draw(adata, figsize=(4, 4), lim=(0, 4000), **kwargs):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(*lim)
    ax.set_ylim(lim[1], lim[0])  # inverted, like invert_y=True
    ax.set_aspect("equal")
    info = add_scale_bar(ax, adata, library_id="s", **kwargs)
    return fig, ax, info


@pytest.mark.parametrize("loc", sorted(_LOC_FACTORS))
def test_bar_and_label_stay_inside_the_axes_for_every_loc(loc):
    adata = _make_adata()
    fig, ax, info = _draw(adata, length="500um", loc=loc)
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


@pytest.mark.parametrize("offset", [(0, 0), (20, -20), (-25, 18), (60, 60), (-60, -60)])
def test_offsets_never_push_the_bar_out(offset):
    adata = _make_adata()
    for loc in ("lower right", "upper left", "center"):
        fig, ax, _ = _draw(adata, length="500um", loc=loc, offset=offset)
        line = ax.lines[-1]
        _assert_inside(
            ax, _extents(ax, line, line.get_markersize(), line.get_linewidth())
        )
        plt.close(fig)


def test_long_label_is_kept_inside_horizontally():
    """A label wider than the bar is recentred, not allowed to hang off the edge."""
    adata = _make_adata()
    fig, ax, _ = _draw(
        adata,
        length="500um",
        loc="lower right",
        label="500 micrometres (measured)",
        fontsize=11,
    )
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


def test_frame_is_inside_too():
    adata = _make_adata()
    fig, ax, _ = _draw(
        adata, length="500um", loc="lower left", frame=True, frame_pad=6, fontsize=12
    )
    frame = ax.patches[-1]
    fb = frame.get_window_extent(ax.figure.canvas.get_renderer())
    _assert_inside(ax, [(fb.x0, fb.y0, fb.x1, fb.y1)])
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


def test_without_caps_the_line_still_stays_inside():
    adata = _make_adata()
    fig, ax, _ = _draw(adata, length="500um", loc="lower right", caps=False)
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, 0.0, line.get_linewidth()))
    plt.close(fig)


def test_label_flips_side_instead_of_escaping_the_top():
    """``upper*`` anchor + positive y offset used to push the label past the top."""
    adata = _make_adata()
    fig, ax, _ = _draw(adata, length="500um", loc="upper right", offset=(0, 8))
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


def test_narrow_axes_keep_bar_inside():
    """Tall/narrow aspect: the bar is a big fraction of the width here."""
    adata = _make_adata()
    fig, ax, info = _draw(adata, figsize=(2.2, 6), length="1000um", loc="lower right")
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


def test_bar_too_long_for_the_axes_warns_instead_of_lying():
    """A truncated bar would report a wrong length; we keep the truth and warn.

    A 3 cm wide axes leaves ~20 % of its width once a framed bar's margins are
    taken out, so a 1000 um bar across a 2000 um span cannot fit.
    """
    adata = _make_adata(span=2000.0)
    fig, ax = plt.subplots(figsize=(1.2, 1.2))
    ax.set_xlim(0, 2000)
    ax.set_ylim(2000, 0)
    ax.set_aspect("equal")
    with pytest.warns(UserWarning, match="does not fit inside the axes"):
        add_scale_bar(
            ax, adata, library_id="s", length="900um", frame=True, frame_pad=20
        )
    plt.close(fig)


def test_bar_stays_inside_on_a_real_spatial_plot():
    """End-to-end: the bar must not overlap the panel border of spatial_cell."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    from trackcell.pl.plot import spatial_cell

    polys = [
        Polygon([(x, y), (x + 10, y), (x + 10, y + 10), (x, y + 10)])
        for x in range(0, 400, 10)
        for y in range(0, 400, 10)
    ]
    obs = pd.DataFrame(
        {"ct": ["a"] * len(polys)}, index=[f"c{i}" for i in range(len(polys))]
    )
    adata = AnnData(X=np.ones((len(polys), 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = np.array([[p.centroid.x, p.centroid.y] for p in polys])
    adata.uns["spatial"] = {
        "s": {
            "metadata": {"coord_unit": "um"},
            "scalefactors": {},
            "geometries": gpd.GeoDataFrame(
                {"geometry": polys}, index=[f"c{i}" for i in range(len(polys))]
            ),
        }
    }
    ax = spatial_cell(
        adata,
        color="ct",
        library_id="s",
        show=False,
        legend=False,
        scalebar=dict(length="50um", loc="lower right", linewidth=3.0),
    )
    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(ax.figure)


def test_geometry_survives_tight_layout_after_placement():
    """The reported bug: layout runs *after* the bar is placed, moving the box.

    ``tight_layout`` (and ``savefig(bbox_inches="tight")``) shrink the axes, so
    geometry clamped against the pre-layout box ends up across the spine. The
    placer artist recomputes at draw time, which is what this checks.
    """
    adata = _make_adata()
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    ax.set_aspect("equal")
    add_scale_bar(ax, adata, library_id="s", length="500um", loc="lower right")
    fig.tight_layout(pad=2.0)  # moves the axes *after* add_scale_bar

    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)


def test_geometry_survives_figure_resize():
    adata = _make_adata()
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0, 4000)
    ax.set_ylim(4000, 0)
    ax.set_aspect("equal")
    add_scale_bar(
        ax, adata, library_id="s", length="800um", loc="upper left", frame=True
    )
    fig.set_size_inches(3.0, 5.0)  # aggressive reshape

    line = ax.lines[-1]
    _assert_inside(ax, _extents(ax, line, line.get_markersize(), line.get_linewidth()))
    plt.close(fig)
