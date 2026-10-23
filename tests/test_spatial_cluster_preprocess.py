"""Tests for YardCluster preprocessing fallbacks and parameter validation."""

import numpy as np
import pytest
from anndata import AnnData

import trackcell as tcl


def _make_counts(n_cells: int = 80, n_genes: int = 50, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    adata = AnnData(
        X=rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    )
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obsm["spatial"] = rng.normal(size=(n_cells, 2))
    return adata


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (dict(mode="celltypes"), "Unknown mode"),
        (dict(method="louvain"), "Unknown method"),
        (dict(integrate="foo"), "Unknown integrate"),
        (dict(graph_mode="foo"), "Unknown graph_mode"),
    ],
)
def test_spatial_cluster_rejects_invalid_parameters(kwargs, match):
    adata = _make_counts()
    with pytest.raises(ValueError, match=match):
        tcl.tl.spatial_cluster(adata, preprocess=False, **kwargs)


def test_preprocess_hvg_falls_back_without_scikit_misc():
    """HVG selection must not require scikit-misc on raw counts."""
    adata = _make_counts()
    tcl.tl.spatial_cluster(
        adata,
        mode="celltype",
        preprocess=True,
        n_top_genes=20,
        k_spatial=5,
        n_pcs=5,
        random_state=0,
    )
    assert adata.n_vars == 20
    assert "yardcluster" in adata.obs.columns


def test_neighborhood_feature_cache_invalidates_on_layer_change():
    rng = np.random.default_rng(1)
    adata = AnnData(X=rng.random((40, 8), dtype=np.float32))
    adata.var_names = [f"g{i}" for i in range(8)]
    adata.obsm["spatial"] = rng.normal(size=(40, 2))
    adata.layers["a"] = adata.X.copy()
    adata.layers["b"] = adata.X + 1.0

    tcl.tl.neighborhood_features(
        adata, k_spatial=3, use_gradient=False, layer="a", key_added="k"
    )
    first = adata.obsm["k_M"].copy()
    assert adata.uns["k_features_params"]["layer"] == "a"

    tcl.tl.neighborhood_features(
        adata, k_spatial=3, use_gradient=False, layer="b", key_added="k"
    )
    second = adata.obsm["k_M"].copy()
    assert adata.uns["k_features_params"]["layer"] == "b"
    assert not np.allclose(first, second)
