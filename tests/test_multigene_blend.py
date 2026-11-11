"""Tests for multigene blend/facet helpers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from anndata import AnnData

import trackcell as tcl


def _make_adata(n_genes: int = 8, n_cells: int = 6) -> AnnData:
    rng = np.random.default_rng(0)
    adata = AnnData(X=rng.random((n_cells, n_genes)).astype(np.float32))
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.layers["scaled"] = adata.X * 2
    return adata


def test_blend_supports_more_than_seven_genes_without_colors():
    adata = _make_adata()
    with pytest.warns(UserWarning, match="muddy"):
        colors = tcl.tl.multigene_blend(
            adata, genes=list(adata.var_names), inplace=False
        )
    assert len(colors) == adata.n_obs
    assert colors.nunique() > 1
    assert "multigene_blend_params" not in adata.uns


def test_blend_inplace_stores_params():
    adata = _make_adata(n_genes=3)
    tcl.tl.multigene_blend(adata, genes=["g0", "g1", "g2"])
    assert "multigene_blend" in adata.obs.columns
    assert adata.uns["multigene_blend_params"]["genes"] == ["g0", "g1", "g2"]


def test_blend_accepts_integer_expression_matrix():
    """Integer counts must be normalised, not truncated to the background.

    ``adata[:, genes].X`` is a dense ``ArrayView`` that also exposes
    ``toarray()``; the old `hasattr` sparse check skipped the float cast, and
    ``np.zeros_like(int_array)`` truncated every normalised value to zero.
    """
    from scipy.sparse import csr_matrix

    values = np.array([[1, 10, 0], [2, 20, 0], [3, 30, 0], [0, 0, 0]], dtype=np.int64)
    expected_background = "#e0e0e0"

    for X in (values, csr_matrix(values)):
        adata = AnnData(X=X.copy())
        adata.var_names = ["g0", "g1", "g2"]
        colors = tcl.tl.multigene_blend(
            adata,
            genes=["g0", "g1"],
            colors=["#ff0000", "#0000ff"],
            inplace=False,
        )
        assert colors.iloc[:3].ne(expected_background).all()
        assert colors.iloc[3] == expected_background
        # The expressed rows are blends, not one source colour or the background.
        assert not set(colors.iloc[:3]) & {"#ff0000", "#0000ff"}


def test_blend_keeps_genes_expressed_below_the_percentile_tail():
    """A gene in 0.3% of cells must not be silently dropped.

    With the default 99th percentile, a gene present in fewer than 1% of cells
    has ``vmax == 0``. The old code zeroed the whole gene, so those cells came
    out as background — exactly the rare markers a co-expression map is used to
    find.
    """
    n = 1000
    X = np.zeros((n, 2), dtype=np.float32)
    X[3:, 0] = np.linspace(1, 100, n - 3)  # common gene, absent in first 3 cells
    X[:3, 1] = 10.0  # rare gene
    adata = AnnData(X=X)
    adata.var_names = ["common", "rare"]
    adata.obs_names = [f"c{i}" for i in range(n)]

    colors = tcl.tl.multigene_blend(
        adata,
        genes=["rare", "common"],
        colors=["#ff0000", "#0000ff"],
        inplace=False,
    )

    assert colors.iloc[:3].tolist() == ["#ff0000"] * 3
    assert colors.iloc[-1] == "#0000ff"
    assert colors.iloc[3:].ne("#e0e0e0").any()


def test_blend_rejects_duplicate_genes():
    adata = _make_adata(n_genes=3)
    with pytest.raises(ValueError, match="duplicates"):
        tcl.tl.multigene_blend(adata, genes=["g0", "g0"])


def test_blend_rejects_negative_expression():
    adata = _make_adata(n_genes=3)
    adata.X[0, 0] = -1.0
    with pytest.raises(ValueError, match="non-negative"):
        tcl.tl.multigene_blend(adata, genes=["g0", "g1"])


def test_blend_normalises_colour_names_and_validates_gamma():
    adata = _make_adata(n_genes=3)
    colors = tcl.tl.multigene_blend(
        adata, genes=["g0", "g1"], colors=["red", "blue"], inplace=False
    )
    assert colors.attrs["colors"] == ["#ff0000", "#0000ff"]

    with pytest.raises(ValueError, match="gamma"):
        tcl.tl.multigene_blend(adata, genes=["g0", "g1"], gamma=0)


def _geometries(adata):
    import geopandas as gpd
    from shapely.geometry import box

    return gpd.GeoDataFrame(
        geometry=[box(i, -0.5, i + 1, 0.5) for i in range(adata.n_obs)],
        index=adata.obs_names,
    )


def test_facet_respects_layer_and_percentiles():
    adata = _make_adata(n_genes=3)
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(adata.n_obs), np.zeros(adata.n_obs)]
    )
    adata.uns["spatial"] = {"sample": {"geometries": _geometries(adata)}}

    axes = tcl.tl.multigene_blend(
        adata,
        genes=["g0", "g1", "g2"],
        mode="facet",
        layer="scaled",
        ncols=3,
        show=False,
    )
    assert np.asarray(axes).size == 3
    # Temporary layer columns must be cleaned up.
    assert not any(col.startswith("__facet_") for col in adata.obs.columns)
    # The colourbar for each panel should exist (vmin/vmax were applied).
    fig = np.asarray(axes).ravel()[0].figure
    assert len(fig.axes) >= 6  # 3 panels + 3 colorbars
    plt.close(fig)
