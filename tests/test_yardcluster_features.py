"""Tests for CellCharter-inspired YardCluster features."""

import numpy as np
import pytest
import scipy.sparse as sps
from anndata import AnnData

from trackcell.tl._spatial_graph import (
    build_spatial_weights,
    compute_neighborhood_multi_hop,
    remove_long_links,
)
from trackcell.tl._nhood import (
    connected_components,
    diff_nhood_enrichment,
    nhood_enrichment,
)
from trackcell.tl._autores import select_resolution, sweep_resolutions
from trackcell.tl.spatial_cluster import (
    neighborhood_features,
    spatial_cluster,
    spatial_neighbors,
    yard_cluster,
    yard_embed,
)


def _make_two_blob_adata(n_per_blob: int = 25, n_genes: int = 40, seed: int = 0) -> AnnData:
    """Two well-separated blobs -> two spatial components."""
    rng = np.random.default_rng(seed)
    coords = np.vstack(
        [
            rng.normal(loc=[0, 0], scale=0.3, size=(n_per_blob, 2)),
            rng.normal(loc=[10, 10], scale=0.3, size=(n_per_blob, 2)),
        ]
    )
    X = rng.normal(size=(2 * n_per_blob, n_genes)).astype(np.float32)
    X[:n_per_blob, :5] += 3.0
    X[n_per_blob:, 5:10] += 3.0
    adata = AnnData(X=X)
    adata.obsm["spatial"] = coords
    adata.obs["blob"] = np.array(["A"] * n_per_blob + ["B"] * n_per_blob)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_multi_hop_shapes_and_dimension():
    adata = _make_two_blob_adata(n_per_blob=15, n_genes=20)
    coords = adata.obsm["spatial"]
    w = build_spatial_weights(coords, k=5)
    X = np.asarray(adata.X, dtype=np.float64)
    m1 = compute_neighborhood_multi_hop(X, w, n_layers=1)
    m3 = compute_neighborhood_multi_hop(X, w, n_layers=3)
    assert m1.shape == (adata.n_obs, 20)
    assert m3.shape == (adata.n_obs, 60)
    np.testing.assert_allclose(m3[:, :20], m1, rtol=1e-5, atol=1e-5)


def test_neighborhood_features_n_layers():
    adata = _make_two_blob_adata(n_per_blob=12, n_genes=25)
    neighborhood_features(adata, k_spatial=5, n_layers=2, key_added="yard")
    assert adata.obsm["yard_M"].shape == (adata.n_obs, 50)
    # params recorded for cache invalidation
    assert adata.uns["yard_features_params"]["n_layers"] == 2


def test_remove_long_links_trims_edges():
    coords = np.array([[0, 0], [0.1, 0], [0, 0.1], [50, 50]], dtype=float)
    w = build_spatial_weights(coords, k=3)
    n_edges_before = (w > 0).nnz
    trimmed = remove_long_links(coords, w, distance_percentile=60.0)
    assert (trimmed > 0).nnz < n_edges_before
    # rows still sum to 1 (where any edge remains)
    row_sums = np.asarray(trimmed.sum(axis=1)).ravel()
    assert np.allclose(row_sums[row_sums > 0], 1.0, rtol=1e-5)


def test_connected_components_splits_cluster():
    adata = _make_two_blob_adata(n_per_blob=20, n_genes=20)
    spatial_neighbors(adata, k=5)
    # assign both blobs the same cluster label -> components should split them
    adata.obs["cluster"] = "0"
    connected_components(adata, cluster_key="cluster", min_cells=3, out_key="component")
    comp = adata.obs["component"].astype(str)
    assert comp.nunique() == 2
    assert set(comp.unique()) == {"0", "1"}


def test_connected_components_whole_graph():
    adata = _make_two_blob_adata(n_per_blob=20, n_genes=20)
    spatial_neighbors(adata, k=5)
    connected_components(adata, cluster_key=None, min_cells=3, out_key="component")
    assert adata.obs["component"].astype(str).nunique() == 2


def test_nhood_enrichment_symmetric_and_inter():
    adata = _make_two_blob_adata(n_per_blob=30, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    res = nhood_enrichment(
        adata, "yardcluster", connectivity_key="spatial_connectivities", copy=True
    )
    assert set(res.keys()) == {"enrichment", "observed", "expected"}
    # diagonal should be NaN when only_inter=True (default)
    assert np.isnan(res["enrichment"].values.diagonal()).all()
    # matrix is symmetric
    np.testing.assert_allclose(
        res["enrichment"].values,
        res["enrichment"].values.T,
        rtol=1e-9,
        equal_nan=True,
    )


def test_nhood_enrichment_writes_to_uns():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=20)
    spatial_neighbors(adata, k=6)
    adata.obs["cluster"] = ["0", "1"] * 25
    nhood_enrichment(adata, "cluster")
    assert "cluster_nhood_enrichment" in adata.uns
    assert "params" in adata.uns["cluster_nhood_enrichment"]


def test_diff_nhood_enrichment_with_pvalues():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=40, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    # half cells Case, half Ctr (split each blob)
    adata.obs["condition"] = ["Case"] * 40 + ["Ctr"] * 40
    adata.obs["library"] = ["s1"] * 20 + ["s2"] * 20 + ["s3"] * 20 + ["s4"] * 20
    res = diff_nhood_enrichment(
        adata,
        "yardcluster",
        "condition",
        ("Case", "Ctr"),
        library_key="library",
        pvalues=True,
        n_perms=20,
        copy=True,
    )
    assert "Case_Ctr" in res
    assert "enrichment" in res["Case_Ctr"]
    assert "pvalue" in res["Case_Ctr"]
    pv = res["Case_Ctr"]["pvalue"].values
    assert ((pv >= 0) & (pv <= 1) | np.isnan(pv)).all()


def test_sweep_and_select_resolution():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    sweep = sweep_resolutions(
        adata, use_rep="X_yard_mixed", resolutions=[0.3, 0.5, 0.8], n_runs=3
    )
    assert list(sweep.columns) == ["resolution", "n_clusters", "stability"]
    assert len(sweep) == 3
    assert (sweep["stability"] >= 0).all()

    best, sweep2 = select_resolution(
        adata, use_rep="X_yard_mixed", resolutions=[0.3, 0.5, 0.8], n_runs=3
    )
    assert best in (0.3, 0.5, 0.8)


def test_yard_cluster_auto_resolution():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    yard_cluster(adata, resolution="auto", n_resolution_runs=3, key_added="yardcluster")
    assert "yardcluster" in adata.obs
    assert adata.obs["yardcluster"].nunique() >= 2
    assert "yardcluster_resolution_sweep" in adata.uns
    assert adata.uns["yardcluster_cluster_params"]["resolution"] != "auto"


def test_yard_cluster_gmm():
    pytest.importorskip("sklearn")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    yard_cluster(adata, method="gmm", n_clusters=2, key_added="yardcluster")
    assert "yardcluster" in adata.obs
    assert adata.uns["yardcluster_cluster_params"]["method"] == "gmm"
    assert adata.obs["yardcluster"].nunique() == 2


def test_spatial_cluster_auto_resolution_smoke():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        preprocess=False,
        k_spatial=6,
        n_pcs=10,
        resolution="auto",
        n_resolution_runs=2,
        use_gradient=False,
    )
    assert "yardcluster" in adata.obs


def test_spatial_cluster_n_layers_smoke():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        preprocess=False,
        k_spatial=6,
        n_pcs=10,
        n_layers=2,
        use_gradient=False,
        resolution=0.5,
    )
    assert "yardcluster" in adata.obs


def test_spatial_cluster_gmm_single_mode():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        method="gmm",
        n_clusters=2,
        preprocess=False,
        k_spatial=6,
        n_pcs=8,
        use_gradient=False,
    )
    assert "yardcluster" in adata.obs
    assert adata.uns["yardcluster_cluster_params"]["method"] == "gmm"
    assert adata.obs["yardcluster"].nunique() == 2


def test_spatial_cluster_gmm_auto_mode():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="auto",
        method="gmm",
        n_clusters=None,
        merge_clusters=False,
        preprocess=False,
        k_spatial=6,
        n_pcs=8,
        use_gradient=False,
    )
    assert "yardcluster_celltype" in adata.obs
    assert "yardcluster_domain" in adata.obs
    assert adata.uns["yardcluster_params"]["method"] == "gmm"


def test_pl_nhood_enrichment_heatmap():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from trackcell.pl import diff_nhood_enrichment as pl_diff
    from trackcell.pl import nhood_enrichment as pl_nhood

    adata = _make_two_blob_adata(n_per_blob=25, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    nhood_enrichment(adata, "yardcluster")
    ax = pl_nhood(adata, "yardcluster", show=False)
    assert ax is not None
    assert len(ax.images) == 1
    plt.close("all")

    adata.obs["condition"] = ["Case"] * 25 + ["Ctr"] * 25
    adata.obs["library"] = ["s1"] * 13 + ["s2"] * 12 + ["s3"] * 13 + ["s4"] * 12
    diff_nhood_enrichment(
        adata,
        "yardcluster",
        "condition",
        ("Case", "Ctr"),
        library_key="library",
        pvalues=True,
        n_perms=10,
    )
    ax2 = pl_diff(adata, "yardcluster", "condition", ("Case", "Ctr"), show=False)
    assert ax2 is not None
    plt.close("all")
