"""Tests for YardCluster spatial clustering."""

import numpy as np
import pytest
import scipy.sparse as sps
from anndata import AnnData

from trackcell.tl._spatial_graph import (
    build_spatial_weights,
    compute_neighborhood_gradient,
    compute_neighborhood_mean,
    mix_embeddings,
)
from trackcell.tl._cluster_merge import merge_clusters_de
from trackcell.tl.spatial_cluster import (
    neighborhood_features,
    spatial_cluster,
    spatial_neighbors,
    yard_cluster,
    yard_embed,
)


def _make_ring_adata(n_per_ring: int = 30, n_genes: int = 50, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    angles_inner = np.linspace(0, 2 * np.pi, n_per_ring, endpoint=False)
    angles_outer = np.linspace(0, 2 * np.pi, n_per_ring, endpoint=False)
    inner = np.column_stack([np.cos(angles_inner), np.sin(angles_inner)]) * 1.0
    outer = np.column_stack([np.cos(angles_outer), np.sin(angles_outer)]) * 2.5
    coords = np.vstack([inner, outer])
    X = rng.normal(size=(2 * n_per_ring, n_genes)).astype(np.float32)
    X[:n_per_ring, :5] += 3.0
    X[n_per_ring:, 5:10] += 3.0
    adata = AnnData(X=X)
    adata.obsm["spatial"] = coords
    adata.obs["ring"] = np.array(["inner"] * n_per_ring + ["outer"] * n_per_ring)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs_names = [f"c{i}" for i in range(2 * n_per_ring)]
    return adata


def test_build_spatial_weights_row_normalized():
    coords = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    w = build_spatial_weights(coords, k=2)
    row_sums = np.array(w.sum(axis=1)).ravel()
    np.testing.assert_allclose(row_sums, 1.0, rtol=1e-5)


def test_gradient_regression():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=12, n_genes=20)
    spatial_neighbors(adata, k=8)
    X = adata.X
    weights = adata.obsp["spatial_connectivities"]
    coords = adata.obsm["spatial"]
    G_reg = compute_neighborhood_gradient(
        X, coords, weights, mode="regression", k_gradient=8, chunk_size=10
    )
    G_cos = compute_neighborhood_gradient(
        X, coords, weights, mode="cosphi", k_gradient=8, chunk_size=10
    )
    assert G_reg.shape == G_cos.shape == adata.shape


def test_yard_pipeline_produces_clusters():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=25, n_genes=40)
    neighborhood_features(adata, k_spatial=8)
    yard_embed(adata, lam=0.2, n_pcs=10)
    yard_cluster(adata, resolution=0.5, key_added="yardcluster")
    assert "yardcluster" in adata.obs
    assert adata.obs["yardcluster"].nunique() >= 2


def test_spatial_cluster_auto_mode():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=25, n_genes=40)
    spatial_cluster(adata, mode="auto", preprocess=False, k_spatial=8, n_pcs=10, resolution=0.5)
    assert "yardcluster_celltype" in adata.obs
    assert "yardcluster_domain" in adata.obs


def test_spatial_cluster_separate_batches():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=15, n_genes=30)
    adata.obs["sample"] = ["A"] * 15 + ["B"] * 15
    spatial_cluster(
        adata,
        mode="celltype",
        batch_key="sample",
        integrate="separate",
        preprocess=False,
        k_spatial=5,
        n_pcs=8,
        resolution=0.5,
    )
    labels = adata.obs["yardcluster"].astype(str)
    assert labels.str.startswith("A_").any() or labels.str.startswith("B_").any()


def test_sketch_clustering():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=40, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.2, n_pcs=8)
    yard_cluster(
        adata,
        resolution=0.5,
        key_added="yardcluster",
        sketch=True,
        sketch_threshold=10,
        sketch_n=20,
    )
    assert "yardcluster" in adata.obs


def test_merge_clusters_de():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=20, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.2, n_pcs=8)
    yard_cluster(adata, resolution=1.5, key_added="yardcluster")
    n_before = adata.obs["yardcluster"].nunique()
    merge_clusters_de(adata, "yardcluster", "X_yard_mixed", adj_p_threshold=0.05)
    n_after = adata.obs["yardcluster_merged"].nunique()
    assert n_after <= n_before


def test_sparse_expression_input():
    pytest.importorskip("leidenalg")
    adata = _make_ring_adata(n_per_ring=12, n_genes=25)
    adata.X = sps.csr_matrix(adata.X)
    spatial_cluster(adata, mode="celltype", preprocess=False, k_spatial=5, n_pcs=5)
