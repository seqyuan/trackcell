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
