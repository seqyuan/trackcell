"""Tests for SCT RPCA integration."""

import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sps
from anndata import AnnData

from trackcell.tl.integration import integrate_sct_rpca, select_integration_features
from trackcell.tl.integration.select_features import _select_integration_features_seurat
from trackcell.tl.integration._annoy_nn import annoy_available
from trackcell.tl.integration._xgboost_integration import xgboost_available
from trackcell.tl.integration.sample_tree import build_sample_tree, count_anchor_similarity
from trackcell.tl.sctransform import sctransform


def _make_batch_adata(n_cells=(40, 35), n_genes=80, seed=0) -> AnnData:
    rng = np.random.default_rng(seed)
    adatas = []
    for i, n in enumerate(n_cells):
        lam = rng.uniform(1.0, 6.0, size=(n, n_genes))
        counts = rng.poisson(lam).astype(np.float64)
        ad = AnnData(X=sps.csr_matrix(counts))
        ad.obs_names = [f"b{i}_c{j}" for j in range(n)]
        ad.var_names = [f"g{j}" for j in range(n_genes)]
        ad.obs["orig.ident"] = f"B{i}"
        adatas.append(ad)
    merged = adatas[0].concatenate(adatas[1], batch_key="orig.ident")
    merged.obs["orig.ident"] = merged.obs["orig.ident"].astype(str).str.replace("-.*", "", regex=True)
    merged.layers["counts"] = merged.X.copy()
    return merged


def _make_three_batch_adata(n_cells=(40, 35, 30), n_genes=80, seed=1) -> AnnData:
    rng = np.random.default_rng(seed)
    adatas = []
    for i, n in enumerate(n_cells):
        lam = rng.uniform(1.0, 6.0, size=(n, n_genes))
        counts = rng.poisson(lam).astype(np.float64)
        ad = AnnData(X=sps.csr_matrix(counts))
        ad.obs_names = [f"b{i}_c{j}" for j in range(n)]
        ad.var_names = [f"g{j}" for j in range(n_genes)]
        ad.obs["orig.ident"] = f"B{i}"
        adatas.append(ad)
    merged = sc.concat(adatas, label="orig.ident", keys=[f"B{i}" for i in range(3)])
    merged.layers["counts"] = merged.X.copy()
    return merged


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_integrate_sct_rpca_auto_run_sct():
    adata = _make_batch_adata()
    out = integrate_sct_rpca(
        adata,
        batch_key="orig.ident",
        n_features=15,
        dims=10,
        k_anchor=3,
        k_weight=20,
        run_sct=True,
        sct_method="python",
        n_sct_hvg=20,
        sct_seed=1,
        copy=True,
    )
    assert "batch_models" in out.uns["sct"]
    assert "X_sct_integrated" in out.obsm


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_integrate_sct_rpca_runs():
    adata = _make_batch_adata()
    sctransform(
        adata,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        seed=1,
    )
    out = integrate_sct_rpca(
        adata,
        batch_key="orig.ident",
        n_features=15,
        dims=10,
        k_anchor=3,
        k_weight=20,
        copy=True,
    )
    assert "X_sct_integrated" in out.obsm
    assert out.obsm["X_sct_integrated"].shape == (adata.n_obs, 15)
    assert np.isfinite(out.obsm["X_sct_integrated"]).all()
    assert out.uns["sct_integrated"]["n_anchors"] > 0
    assert out.uns["sct_integrated"]["params"]["nn_method"] == "annoy"
    assert out.uns["sct_integrated"]["params"]["merge_method"] == "sample_tree"


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_integrate_sct_rpca_three_batches_sample_tree():
    adata = _make_three_batch_adata()
    sctransform(
        adata,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=25,
        n_genes=40,
        seed=2,
    )
    out = integrate_sct_rpca(
        adata,
        batch_key="orig.ident",
        n_features=12,
        dims=8,
        k_anchor=3,
        k_weight=15,
        copy=True,
    )
    assert out.obsm["X_sct_integrated"].shape == (adata.n_obs, 12)
    assert len(out.uns["sct_integrated"]["sample_tree"]) == 2


def test_build_sample_tree_from_anchors():
    import pandas as pd

    anchors = pd.DataFrame(
        {
            "dataset1": [0, 0, 1],
            "dataset2": [1, 2, 2],
            "cell1": [0, 1, 2],
            "cell2": [3, 4, 5],
            "score": [1.0, 0.9, 0.8],
        }
    )
    sim = count_anchor_similarity(anchors, 3, [100, 100, 100])
    tree = build_sample_tree(sim)
    assert tree.shape == (2, 2)


def test_select_integration_features_seurat_tiebreak():
    vf_a = ["g1", "g2", "g3", "g4"]
    vf_b = ["g2", "g1", "g5", "g3"]
    vf_c = ["g1", "g3", "g2", "g6"]
    selected = _select_integration_features_seurat(
        [vf_a, vf_b, vf_c],
        n_features=3,
    )
    assert selected[0] == "g1"
    assert set(selected) == {"g1", "g2", "g3"}


def test_select_integration_features_requires_batch_models():
    adata = _make_batch_adata()
    with pytest.raises(ValueError, match="batch_models|sctransform"):
        select_integration_features(adata, batch_key="orig.ident")


@pytest.mark.skipif(
    not annoy_available() or not xgboost_available(),
    reason="annoy and xgboost required",
)
def test_integrate_sct_rpca_xgboost_method():
    adata = _make_batch_adata()
    sctransform(
        adata,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        seed=1,
    )
    out = integrate_sct_rpca(
        adata,
        batch_key="orig.ident",
        n_features=15,
        dims=10,
        k_anchor=3,
        method="xgboost",
        copy=True,
    )
    assert "X_sct_integrated" in out.obsm
    assert out.obsm["X_sct_integrated"].shape == (adata.n_obs, 15)
    assert np.isfinite(out.obsm["X_sct_integrated"]).all()
    assert out.uns["sct_integrated"]["params"]["method"] == "xgboost"
    assert out.uns["sct_integrated"]["params"]["feature_method"] == "seurat"


@pytest.mark.skipif(not xgboost_available(), reason="xgboost not installed")
def test_select_integration_features_xgboost():
    adata = _make_batch_adata()
    sctransform(
        adata,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        seed=1,
    )
    genes = select_integration_features(
        adata,
        batch_key="orig.ident",
        n_features=10,
        method="xgboost",
    )
    assert len(genes) == 10
    assert len(set(genes)) == 10
