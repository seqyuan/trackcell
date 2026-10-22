"""Tests for YardCluster P2 utilities and options."""

import importlib

import numpy as np
import pytest
from anndata import AnnData

from trackcell.tl._yardcluster_utils import (
    PLATFORM_PRESETS,
    mix_spatial_expr_connectivities,
    resolve_platform_params,
    spatial_stratified_sketch_indices,
)
from trackcell.tl.spatial_cluster import neighborhood_features, spatial_cluster


def _patch_adata(n=60, seed=0):
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, 10, size=(n, 2))
    X = rng.normal(size=(n, 30)).astype(np.float32)
    adata = AnnData(X=X)
    adata.obsm["spatial"] = coords
    adata.obs_names = [f"c{i}" for i in range(n)]
    adata.var_names = [f"g{i}" for i in range(30)]
    return adata


def test_platform_preset_sto():
    p = resolve_platform_params("sto", k_spatial=None, k_expr=None, n_pcs=None, resolution=None, n_top_genes=None)
    assert p["k_spatial"] == PLATFORM_PRESETS["sto"]["k_spatial"]
    assert p["resolution"] == PLATFORM_PRESETS["sto"]["resolution"]


def test_platform_preset_xenium_tma():
    p = resolve_platform_params(
        "xenium_tma", k_spatial=None, k_expr=None, n_pcs=None, resolution=None, n_top_genes=None
    )
    assert p == PLATFORM_PRESETS["xenium_tma"]
    assert p["k_spatial"] == 10
    assert p["n_top_genes"] == 3000


def test_platform_preset_xenium_tma_override():
    p = resolve_platform_params(
        "xenium_tma", k_spatial=15, k_expr=None, n_pcs=None, resolution=1.1, n_top_genes=None
    )
    assert p["k_spatial"] == 15
    assert p["resolution"] == 1.1
    assert p["k_expr"] == PLATFORM_PRESETS["xenium_tma"]["k_expr"]

def test_platform_explicit_override():
    p = resolve_platform_params("sto", k_spatial=7, k_expr=None, n_pcs=None, resolution=None, n_top_genes=None)
    assert p["k_spatial"] == 7


def test_neighborhood_features_cache_reuse():
    adata = _patch_adata()
    neighborhood_features(adata, k_spatial=8, key_added="yard", use_gradient=False)
    assert "yard_M" in adata.obsm
    M1 = adata.obsm["yard_M"].copy()
    neighborhood_features(adata, k_spatial=8, key_added="yard", use_gradient=False, reuse_neighbors=True)
    np.testing.assert_array_equal(adata.obsm["yard_M"], M1)


def test_spatial_stratified_sketch_covers_bins():
    coords = np.vstack(
        [
            np.column_stack([np.linspace(0, 1, 20), np.zeros(20)]),
            np.column_stack([np.linspace(0, 1, 20), np.ones(20) * 10]),
        ]
    )
    idx = spatial_stratified_sketch_indices(coords, n_sketch=10, n_bins=5, seed=0)
    assert len(idx) == 10
    assert len(np.unique(idx)) == 10


def test_hybrid_graph_mix():
    pytest.importorskip("leidenalg")
    adata = _patch_adata(n=40)
    spatial_cluster(
        adata,
        mode="domain",
        graph_mode="hybrid",
        hybrid_alpha=0.5,
        use_gradient=False,
        merge_clusters=False,
        preprocess=False,
        k_spatial=6,
        n_pcs=8,
        resolution=0.5,
    )
    cp = adata.uns["yardcluster_cluster_params"]
    assert cp["neighbors_method"] == "hybrid"
    assert cp["hybrid_alpha"] == 0.5


def test_mix_spatial_expr_shapes():
    n = 5
    s = __import__("scipy.sparse", fromlist=["csr_matrix"]).csr_matrix(
        np.eye(n, dtype=np.float64)
    )
    e = __import__("scipy.sparse", fromlist=["csr_matrix"]).csr_matrix(
        np.ones((n, n), dtype=np.float64)
    )
    mixed = mix_spatial_expr_connectivities(s, e, alpha=0.5)
    assert mixed.shape == (n, n)
    row_sums = np.asarray(mixed.sum(axis=1)).ravel()
    np.testing.assert_allclose(row_sums, 1.0, rtol=1e-5)


def test_integrate_bbknn_smoke():
    pytest.importorskip("leidenalg")
    try:
        import bbknn  # noqa: F401
    except ImportError:
        pytest.skip("bbknn not installed")

    adata = _patch_adata(n=50)
    adata.obs["batch"] = ["A"] * 25 + ["B"] * 25
    spatial_cluster(
        adata,
        mode="celltype",
        batch_key="batch",
        integrate="bbknn",
        use_gradient=False,
        merge_clusters=False,
        preprocess=False,
        k_spatial=5,
        n_pcs=6,
        resolution=0.5,
    )
    assert adata.uns["yardcluster_params"]["neighbors_method"] == "bbknn"
