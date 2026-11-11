"""Tests for :func:`trackcell.pl.multigene_blend`.

``tl.multigene_blend`` writes one hex colour per cell into ``obs`` and the
gene-to-colour mapping into ``uns['<key>_params']``. ``spatial_cell`` can draw
that column but cannot label it — a raw-hex column carries no category names, so
the composite previously came out with no legend at all. These tests pin the
legend contract and, just as importantly, pin that nothing is recomputed.

Runs against the committed ``tests/data/hd_cellseg_tiny`` bundle: the renderer
needs real ``uns['spatial']`` metadata and polygons, so a synthetic AnnData
cannot exercise it.
"""

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from matplotlib.patches import Patch

import trackcell as tcl

FIXTURE = Path(__file__).resolve().parent / "data" / "hd_cellseg_tiny"

# One marker per cell-type block, so the three blends land on distinct colours.
GENES = ["EPCAM", "COL1A1", "PTPRC"]

pytestmark = pytest.mark.skipif(
    not FIXTURE.exists(),
    reason="committed fixture missing; run python scripts/make_fixtures.py",
)


def _fixture():
    return tcl.io.read_hd_cellseg(FIXTURE)


def _blended(key: str = "multigene_blend"):
    adata = _fixture()
    tcl.tl.multigene_blend(adata, genes=GENES, mode="blend", key_added=key)
    return adata


def _legend_labels(ax) -> list[str]:
    legend = ax.get_legend()
    assert legend is not None
    return [text.get_text() for text in legend.get_texts()]


# ── the reason this function exists ─────────────────────────────────────────


def test_legend_labels_are_the_input_genes():
    ax = tcl.pl.multigene_blend(_blended(), show=False)

    assert _legend_labels(ax) == GENES
    assert ax.get_legend().get_title().get_text() == "genes"


def test_legend_swatches_use_the_source_colors():
    """The legend must show the *input* colours, not the blended per-cell hex.

    A legend over blended values would list hundreds of near-identical swatches
    and tell a reader nothing.
    """
    adata = _blended()
    ax = tcl.pl.multigene_blend(adata, show=False)

    expected = adata.uns["multigene_blend_params"]["colors"]
    patches = [h for h in ax.get_legend().legend_handles if isinstance(h, Patch)]
    assert len(patches) == len(expected) == len(GENES)


def test_spatial_cell_alone_produces_no_legend():
    """Guard the gap this function fills, so it cannot silently disappear."""
    adata = _blended()
    fig, ax = plt.subplots()
    tcl.pl.spatial_cell(adata, color="multigene_blend", show=False, ax=ax)
    assert ax.get_legend() is None
    plt.close(fig)


# ── nothing is recomputed ───────────────────────────────────────────────────


def test_does_not_recompute_the_blend():
    adata = _blended()
    before = adata.obs["multigene_blend"].copy()
    params_before = dict(adata.uns["multigene_blend_params"])

    tcl.pl.multigene_blend(adata, show=False)

    assert adata.obs["multigene_blend"].equals(before)
    assert adata.uns["multigene_blend_params"] == params_before


def test_expression_matrix_is_not_needed():
    """Zeroing ``X`` must not change what is drawn: the colours come from obs.

    This is the difference between a renderer and a reimplementation of
    ``tl.multigene_blend``.
    """
    adata = _blended()
    expected = adata.obs["multigene_blend"].to_list()

    adata.X.data[:] = 0  # zero stored values without rewiring the sparse indices
    ax = tcl.pl.multigene_blend(adata, show=False)

    assert adata.obs["multigene_blend"].to_list() == expected
    assert _legend_labels(ax) == GENES


def test_series_round_trip_from_inplace_false():
    """``inplace=False`` returns a Series; only ``inplace=True`` leaves something
    for the renderer to read, and the error says so."""
    adata = _fixture()
    series = tcl.tl.multigene_blend(adata, genes=GENES, inplace=False)
    assert "multigene_blend" not in adata.obs
    assert len(series) == adata.n_obs

    with pytest.raises(KeyError, match="tl.multigene_blend"):
        tcl.pl.multigene_blend(adata, show=False)


# ── error paths ─────────────────────────────────────────────────────────────


def test_missing_column_raises_keyerror():
    with pytest.raises(KeyError, match="tl.multigene_blend"):
        tcl.pl.multigene_blend(_fixture(), show=False)


def test_non_hex_column_raises_valueerror():
    adata = _fixture()
    adata.obs["label"] = adata.obs["classification"].astype(str)
    with pytest.raises(ValueError, match="not a column of hex colour strings"):
        tcl.pl.multigene_blend(adata, key="label", show=False)


# ── options ─────────────────────────────────────────────────────────────────


def test_legend_can_be_disabled():
    ax = tcl.pl.multigene_blend(_blended(), show=False, legend=False)
    assert ax.get_legend() is None


def test_passes_through_to_spatial_cell():
    ax = tcl.pl.multigene_blend(_blended(), show=False, xlabel="custom x")
    assert ax.get_xlabel() == "custom x"


def test_custom_key_and_legend_title():
    adata = _blended(key="blend2")
    ax = tcl.pl.multigene_blend(adata, key="blend2", show=False, legend_title=None)
    assert _legend_labels(ax) == GENES
    assert ax.get_legend().get_title().get_text() == ""


def test_returns_the_axes_it_was_given():
    fig, ax = plt.subplots()
    returned = tcl.pl.multigene_blend(_blended(), show=False, ax=ax)
    assert returned is ax
    plt.close(fig)
