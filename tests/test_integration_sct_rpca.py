"""Tests for SCT RPCA integration."""

import logging

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sps
from anndata import AnnData

from trackcell.tl.integration import integrate_sct_rpca, select_integration_features
from trackcell.tl.integration.select_features import _select_integration_features_seurat
from trackcell.tl.integration._annoy_nn import annoy_available
from trackcell.tl.integration._xgboost_integration import xgboost_available
from trackcell.tl.integration.sample_tree import (
    build_sample_tree,
    count_anchor_similarity,
)
from trackcell.tl._sct_entry import sctransform


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
    merged.obs["orig.ident"] = (
        merged.obs["orig.ident"].astype(str).str.replace("-.*", "", regex=True)
    )
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


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_split_sct_and_rpca_integration():
    """run_sct_integration + integrate_rpca matches integrate_sct_rpca."""
    adata_combined = _make_batch_adata()
    adata_split = _make_batch_adata(seed=0)

    sctransform(
        adata_combined,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        seed=1,
    )
    out_combined = integrate_sct_rpca(
        adata_combined,
        batch_key="orig.ident",
        n_features=15,
        dims=10,
        k_anchor=3,
        k_weight=20,
        run_sct=False,
        copy=True,
    )

    sctransform(
        adata_split,
        layer="counts",
        batch_key="orig.ident",
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        seed=1,
    )
    from trackcell.tl.integration import integrate_rpca, run_sct_integration

    run_sct_integration(
        adata_split,
        batch_key="orig.ident",
        n_features=15,
        run_sct=False,
    )
    out_split = integrate_rpca(
        adata_split,
        batch_key="orig.ident",
        dims=10,
        k_anchor=3,
        k_weight=20,
    )

    assert (
        out_split.obsm["X_sct_integrated"].shape
        == out_combined.obsm["X_sct_integrated"].shape
    )
    np.testing.assert_allclose(
        out_split.obsm["X_sct_integrated"],
        out_combined.obsm["X_sct_integrated"],
        rtol=1e-5,
        atol=1e-5,
    )
    assert "X_sct_prep" in out_split.obsm
    assert out_split.uns["sct_prep"].get("integration_ready") is True


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_run_sct_integration_accepts_explicit_anchor_features():
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
    from trackcell.tl.integration import run_sct_integration

    features = list(adata.uns["sct"]["variable_features"][:15])
    run_sct_integration(
        adata,
        batch_key="orig.ident",
        anchor_features=features,
        run_sct=False,
    )
    assert adata.uns["sct_prep"]["params"]["feature_method"] == "seurat"
    assert adata.uns["sct_prep"]["anchor_features"] == features


@pytest.mark.skipif(not annoy_available(), reason="annoy not installed")
def test_release_sct_integration_cache_then_integrate_rpca():
    """Explicit cache release must not change RPCA output vs Seurat-parity path."""
    from trackcell.tl.integration import (
        integrate_rpca,
        release_sct_integration_cache,
        run_sct_integration,
    )

    adata_ref = _make_batch_adata()
    adata_release = _make_batch_adata(seed=0)

    for ad in (adata_ref, adata_release):
        sctransform(
            ad,
            layer="counts",
            batch_key="orig.ident",
            n_top_genes=20,
            n_cells=30,
            n_genes=40,
            seed=1,
        )

    run_sct_integration(
        adata_ref,
        batch_key="orig.ident",
        n_features=15,
        run_sct=False,
    )
    out_ref = integrate_rpca(
        adata_ref,
        batch_key="orig.ident",
        dims=10,
        k_anchor=3,
        k_weight=20,
    )

    run_sct_integration(
        adata_release,
        batch_key="orig.ident",
        n_features=15,
        run_sct=False,
    )
    assert all(
        "scale_data" in entry
        for entry in adata_release.uns["sct"]["batch_models"].values()
    )
    summary = release_sct_integration_cache(adata_release, what="scale_data")
    assert summary["what"] == "scale_data"
    assert summary["n_batches"] == 2
    assert summary["bytes_estimate"] > 0
    assert all(
        "scale_data" not in entry
        for entry in adata_release.uns["sct"]["batch_models"].values()
    )

    out_release = integrate_rpca(
        adata_release,
        batch_key="orig.ident",
        dims=10,
        k_anchor=3,
        k_weight=20,
    )

    np.testing.assert_allclose(
        out_release.obsm["X_sct_integrated"],
        out_ref.obsm["X_sct_integrated"],
        rtol=1e-5,
        atol=1e-5,
    )


def test_integrate_pair_xgb_reconstructs_full_shift(monkeypatch):
    """The xgb correction must restore ``pca.mean_``, not just the residual part.

    ``delta = y_pred @ pca.components_`` is the projection onto the score space
    without the mean, and the mean *is* the bulk of a batch shift. With a
    predictor that outputs zeros, ``inverse_transform`` must return exactly the
    mean of the anchor deltas — that is the term the old code dropped.
    """
    pytest.importorskip("xgboost")

    from trackcell.tl.integration import _xgboost_integration as xgb_mod

    class _ZeroRegressor:
        def fit(self, x, y):  # noqa: D401 - test double
            self.n_out = np.asarray(y).shape[1]
            return self

        def predict(self, x):
            return np.zeros((len(x), self.n_out))

    class _ZeroMulti:
        def __init__(self, *a, **k):
            pass

        def fit(self, x, y):
            self.n_out = np.asarray(y).shape[1]
            return self

        def predict(self, x):
            return np.zeros((len(x), self.n_out))

    monkeypatch.setattr(xgb_mod, "_require_xgboost", lambda: None)
    monkeypatch.setattr(xgb_mod, "MultiOutputRegressor", _ZeroMulti)

    rng = np.random.default_rng(3)
    n_ref, n_query, n_feat = 60, 80, 12
    data_ref = rng.normal(0, 1, (n_ref, n_feat))
    data_query = rng.normal(0, 1, (n_query, n_feat))

    cell1 = rng.integers(0, n_ref, 40)
    cell2 = rng.integers(0, n_query, 40)
    anchors = pd.DataFrame({"cell1": cell1, "cell2": cell2, "score": np.ones(40)})

    out = xgb_mod.integrate_pair_xgb(data_ref, data_query, anchors, n_delta_pcs=5)
    expected = data_query + (data_query[cell2] - data_ref[cell1]).mean(axis=0)

    np.testing.assert_allclose(out, expected, atol=1e-10)
    # Guard against the regression: the dropped-mean version equals the query
    # matrix untouched under a zero predictor.
    assert not np.allclose(out, data_query)


# ---------------------------------------------------------------------------
# B-1: the two `method` values rank in different directions, and say so
# ---------------------------------------------------------------------------


def _make_batch_driven_adata(seed: int = 7):
    """120 cells, 2 batches, 40 genes, two purpose-built markers.

    ``BATCHONLY`` differs *only* between batches. ``BIO`` varies within a batch,
    which is what ``method='seurat'`` looks for. They pull in different
    directions, which is the whole point of the test below.
    """
    rng = np.random.default_rng(seed)
    n_per, n_genes = 60, 40
    parts = []
    for b, shift in enumerate([0.0, 30.0]):
        lam = rng.uniform(1.0, 5.0, size=(n_per, n_genes))
        lam[:, 0] = 1.0  # BATCHONLY
        lam[:, 1] = rng.uniform(0.5, 25.0, size=n_per)  # BIO
        if b == 1:
            lam[:, 0] += shift
        ad = AnnData(X=sps.csr_matrix(rng.poisson(lam).astype(float)))
        ad.obs_names = [f"b{b}_c{i}" for i in range(n_per)]
        ad.var_names = ["BATCHONLY", "BIO"] + [f"g{i}" for i in range(n_genes - 2)]
        ad.obs["batch"] = f"B{b}"
        parts.append(ad)

    merged = parts[0].concatenate(parts[1], batch_key="batch")
    merged.obs["batch"] = (
        merged.obs["batch"].astype(str).str.replace("-.*", "", regex=True)
    )
    merged.layers["counts"] = merged.X.copy()
    sctransform(
        merged,
        layer="counts",
        batch_key="batch",
        n_top_genes=30,
        n_cells=40,
        n_genes=40,
        seed=1,
    )
    return merged


@pytest.mark.skipif(not xgboost_available(), reason="xgboost not installed")
def test_xgboost_ranks_the_batch_driven_gene_and_seurat_does_not():
    """Pin the documented direction instead of asserting it in prose.

    ``method='xgboost'`` selects genes that *separate* batches; ``'seurat'``
    selects genes that vary *reproducibly within* them. A gene whose only
    variance is a batch shift must therefore be picked by the first and missed by
    the second. If this ever passes with the reverse outcome, the two methods have
    silently become interchangeable and the docstring is wrong.
    """
    adata = _make_batch_driven_adata()

    xgb_genes, table = select_integration_features(
        adata,
        batch_key="batch",
        n_features=10,
        method="xgboost",
        return_importance=True,
    )
    seurat_genes = select_integration_features(
        adata, batch_key="batch", n_features=10, method="seurat"
    )

    assert table["gene"].iloc[0] == "BATCHONLY", table.head().to_string()
    assert "BATCHONLY" in xgb_genes
    assert "BATCHONLY" not in seurat_genes
    # "Not comparable" is a measurable claim, not a disclaimer.
    overlap = set(xgb_genes) & set(seurat_genes)
    assert len(overlap) < len(xgb_genes) // 2, f"the sets agree too much: {overlap}"


@pytest.mark.skipif(not xgboost_available(), reason="xgboost not installed")
def test_return_importance_gives_the_reviewable_table():
    adata = _make_batch_driven_adata()

    genes, table = select_integration_features(
        adata, batch_key="batch", n_features=5, method="xgboost", return_importance=True
    )

    assert list(table.columns) == ["gene", "importance"]
    assert table["importance"].is_monotonic_decreasing
    assert table["gene"].head(len(genes)).tolist() == genes
    assert np.isfinite(table["importance"]).all()


def test_return_importance_is_ignored_for_seurat_with_a_warning(caplog):
    adata = _make_batch_driven_adata()

    with caplog.at_level(logging.WARNING, logger="trackcell"):
        result = select_integration_features(
            adata,
            batch_key="batch",
            n_features=5,
            method="seurat",
            return_importance=True,
        )

    assert isinstance(result, list), "seurat has no single importance score"
    assert any("return_importance" in r.getMessage() for r in caplog.records)


@pytest.mark.skipif(not xgboost_available(), reason="xgboost not installed")
def test_single_batch_yields_nan_importance():
    """No classifier is fitted, so the table must be explicitly empty of scores.

    Built as a genuinely single-batch SCT run: relabelling the batches afterwards
    leaves ``batch_models`` keyed by the old labels, which is a different (and
    invalid) state.
    """
    rng = np.random.default_rng(11)
    counts = rng.poisson(rng.uniform(1.0, 5.0, size=(40, 30))).astype(float)
    adata = AnnData(X=sps.csr_matrix(counts))
    adata.obs_names = [f"c{i}" for i in range(40)]
    adata.var_names = [f"g{i}" for i in range(30)]
    adata.obs["batch"] = "only"
    adata.layers["counts"] = adata.X.copy()
    sctransform(
        adata,
        layer="counts",
        batch_key="batch",
        n_top_genes=25,
        n_cells=30,
        n_genes=30,
        seed=1,
    )

    genes, table = select_integration_features(
        adata, batch_key="batch", n_features=5, method="xgboost", return_importance=True
    )

    assert len(genes) == 5
    assert table["importance"].isna().all()
