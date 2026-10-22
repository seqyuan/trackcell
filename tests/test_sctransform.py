"""Tests for SCTransform normalization."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sps
from anndata import AnnData

from trackcell.tl.sctransform import sctransform


def _make_count_adata(
    n_cells: int = 80,
    n_genes: int = 120,
    seed: int = 0,
    sparse: bool = True,
) -> AnnData:
    rng = np.random.default_rng(seed)
    lam = rng.uniform(0.5, 8.0, size=(n_cells, n_genes))
    counts = rng.poisson(lam).astype(np.float64)
    if sparse:
        x = sps.csr_matrix(counts)
    else:
        x = counts
    adata = AnnData(X=x)
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_sctransform_runs_and_stores_outputs():
    adata = _make_count_adata()
    sctransform(
        adata,
        n_top_genes=20,
        n_cells=40,
        n_genes=50,
        seed=42,
    )
    assert "sct_counts" in adata.layers
    assert "sct_log1p" in adata.layers
    assert "X_sct" in adata.obsm
    assert "sct_highly_variable" in adata.var
    assert "sct" in adata.uns
    assert adata.obsm["X_sct"].shape == (adata.n_obs, 20)
    assert adata.var["sct_highly_variable"].sum() == 20


def test_sctransform_residuals_centered():
    adata = _make_count_adata(n_cells=100, n_genes=150, seed=1)
    sctransform(
        adata,
        n_top_genes=30,
        n_cells=60,
        n_genes=80,
        seed=7,
        clip_range=(-np.inf, np.inf),
    )
    residuals = adata.obsm["X_sct"]
    col_means = residuals.mean(axis=0)
    np.testing.assert_allclose(col_means, 0.0, atol=1e-4)


def test_sctransform_copy_returns_new_object():
    adata = _make_count_adata()
    out = sctransform(
        adata,
        n_top_genes=15,
        n_cells=30,
        n_genes=40,
        copy=True,
        inplace=False,
    )
    assert out is not None
    assert "X_sct" in out.obsm
    assert "X_sct" not in adata.obsm


def test_sctransform_dense_input():
    adata = _make_count_adata(sparse=False)
    sctransform(adata, n_top_genes=10, n_cells=30, n_genes=30)
    assert adata.obsm["X_sct"].shape[1] == 10


def test_sctransform_conserve_memory():
    adata = _make_count_adata(n_cells=80, n_genes=120, seed=2)
    sctransform(
        adata,
        n_top_genes=15,
        n_cells=40,
        n_genes=50,
        conserve_memory=True,
        do_correct_umi=False,
        seed=3,
    )
    assert adata.uns["sct"]["params"]["sct_method"] == "conserve.memory"
    assert adata.obsm["X_sct"].shape == (adata.n_obs, 15)


def test_sctransform_residual_features():
    adata = _make_count_adata(n_cells=60, n_genes=80, seed=4)
    genes = adata.var_names[:10].tolist()
    sctransform(
        adata,
        residual_features=genes,
        n_cells=30,
        n_genes=40,
        seed=5,
    )
    assert adata.uns["sct"]["params"]["sct_method"] == "residual.features"
    assert set(adata.uns["sct"]["variable_features"]) == set(genes)


def test_sctransform_vars_to_regress():
    adata = _make_count_adata(n_cells=70, n_genes=90, seed=6)
    adata.obs["batch"] = np.random.default_rng(0).integers(0, 2, adata.n_obs)
    sctransform(
        adata,
        n_top_genes=12,
        n_cells=35,
        n_genes=45,
        vars_to_regress=["batch"],
        clip_range=(-np.inf, np.inf),
        seed=8,
    )
    assert np.isfinite(adata.obsm["X_sct"]).all()


def test_sctransform_reference_model():
    adata_ref = _make_count_adata(n_cells=60, n_genes=80, seed=9)
    sctransform(adata_ref, n_top_genes=10, n_cells=30, n_genes=40, seed=10)
    ref_model = adata_ref.uns["sct"]["reference_model"]

    adata_query = _make_count_adata(n_cells=50, n_genes=80, seed=11)
    sctransform(
        adata_query,
        reference_sct_model=ref_model,
        n_top_genes=10,
        seed=12,
        save_reference_model=False,
    )
    assert adata_query.uns["sct"]["params"]["sct_method"] == "reference.model"
    assert adata_query.obsm["X_sct"].shape[1] == 10


def test_sctransform_store_in_layer():
    adata = _make_count_adata()
    sctransform(
        adata,
        n_top_genes=12,
        n_cells=30,
        n_genes=30,
        store_residuals_in="layer",
    )
    assert "sct_residuals" in adata.layers
    assert adata.layers["sct_residuals"].shape == (adata.n_obs, adata.n_vars)
    hvg_idx = np.where(adata.var["sct_highly_variable"])[0]
    assert adata.layers["sct_residuals"][:, hvg_idx].nnz > 0


def test_sctransform_batch_var():
    adata = _make_count_adata(n_cells=80, n_genes=100, seed=13)
    adata.obs["batch"] = np.random.default_rng(13).choice(["A", "B"], adata.n_obs)
    sctransform(
        adata,
        n_top_genes=15,
        n_cells=40,
        n_genes=50,
        batch_var="batch",
        seed=14,
    )
    assert "batch" in adata.uns["sct"]["model_str"]
    assert adata.obsm["X_sct"].shape == (adata.n_obs, 15)
    assert np.isfinite(adata.obsm["X_sct"]).all()


def test_sctransform_latent_var():
    adata = _make_count_adata(n_cells=70, n_genes=90, seed=15)
    adata.obs["pct_mito"] = np.random.default_rng(15).uniform(0.01, 0.2, adata.n_obs)
    sctransform(
        adata,
        n_top_genes=12,
        n_cells=35,
        n_genes=45,
        latent_var=["log_umi", "pct_mito"],
        seed=16,
    )
    assert "pct_mito" in adata.uns["sct"]["model_str"]
    assert adata.obsm["X_sct"].shape[1] == 12


def test_sctransform_deviance_residuals():
    adata = _make_count_adata(n_cells=60, n_genes=80, seed=17)
    sctransform(
        adata,
        n_top_genes=10,
        n_cells=30,
        n_genes=40,
        residual_type="deviance",
        seed=18,
    )
    assert adata.uns["sct"]["params"]["residual_type"] == "deviance"
    assert np.isfinite(adata.obsm["X_sct"]).all()


def test_sctransform_vst_flavor_v2_uses_python_by_default(monkeypatch):
    """Default v2 remains pure Python even when pyglmGamPoi is installed."""
    import trackcell.tl._sctransform_v2 as sct_v2

    monkeypatch.setattr(sct_v2, "pyglmGamPoi_available", lambda: True)
    adata = _make_count_adata(n_cells=80, n_genes=100, seed=21)
    sctransform(
        adata,
        n_top_genes=15,
        n_cells=40,
        n_genes=50,
        vst_flavor="v2",
        seed=22,
    )
    assert adata.uns["sct"]["params"]["vst_flavor"] == "v2"
    assert adata.uns["sct"]["params"]["v2_backend"] == "python"
    assert adata.uns["sct"]["params"]["vst_method"] == "nb_offset"
    assert np.isfinite(adata.obsm["X_sct"]).all()


def test_python_v2_offset_handles_multiple_genes():
    from trackcell.tl._sctransform_v2 import _fit_nb_offset_py

    rng = np.random.default_rng(7)
    n = 400
    umi = np.round(np.exp(rng.normal(np.log(3000), 0.4, n)))
    mu = np.exp(np.log(2e-4) + np.log(umi))
    y = np.stack(
        [
            rng.negative_binomial(1.2, 1.2 / (1.2 + mu)),
            rng.negative_binomial(2.0, 2.0 / (2.0 + mu)),
        ]
    ).astype(float)
    intercept, theta = _fit_nb_offset_py(y, np.log10(umi))
    assert intercept.shape == theta.shape == (2,)
    assert np.isfinite(intercept).all()
    assert np.all(theta > 0)
    # R MASS::glm.nb reference for the first generated gene.
    assert intercept[0] == pytest.approx(-8.707863138517, rel=1e-4)
    assert theta[0] == pytest.approx(2.215043455562, rel=1e-3)


def test_sctransform_rejects_unknown_v2_backend():
    adata = _make_count_adata(n_cells=20, n_genes=30, seed=22)
    with pytest.raises(ValueError, match="v2_backend"):
        sctransform(adata, vst_flavor="v2", v2_backend="auto")


def test_explicit_pyglm_backend_requires_dependency(monkeypatch):
    import trackcell.tl._sctransform_v2 as sct_v2

    monkeypatch.setattr(sct_v2, "pyglmGamPoi_available", lambda: False)
    adata = _make_count_adata(n_cells=20, n_genes=30, seed=23)
    with pytest.raises(RuntimeError, match="requires the optional pyglmGamPoi"):
        sctransform(adata, vst_flavor="v2", v2_backend="pyglmGamPoi")


def test_hvg_residual_variance_is_independent_of_scale_data_clip():
    tight = _make_count_adata(n_cells=50, n_genes=60, seed=24)
    wide = tight.copy()
    common = dict(
        n_top_genes=20,
        n_cells=30,
        n_genes=40,
        vst_flavor="v2",
        v2_backend="python",
        seed=25,
    )
    sctransform(tight, clip_range=(-1.0, 1.0), **common)
    sctransform(wide, clip_range=(-40.0, 40.0), **common)
    tight_attr = pd.DataFrame(tight.uns["sct"]["gene_attr"]).T
    wide_attr = pd.DataFrame(wide.uns["sct"]["gene_attr"]).T
    np.testing.assert_allclose(
        tight_attr["residual_variance"],
        wide_attr["residual_variance"],
        rtol=1e-12,
        atol=1e-12,
    )
    assert np.abs(tight.obsm["X_sct"]).max() <= 1.0


def test_sctransform_batch_key():
    adata = _make_count_adata(n_cells=80, n_genes=100, seed=23)
    adata.obs["orig.ident"] = np.random.default_rng(23).choice(["A", "B"], adata.n_obs)
    sctransform(
        adata,
        batch_key="orig.ident",
        n_top_genes=12,
        n_cells=30,
        n_genes=40,
        seed=24,
    )
    assert "batch_models" in adata.uns["sct"]
    assert set(adata.uns["sct"]["batch_models"]) == {"A", "B"}
    assert adata.obsm["X_sct"].shape == (adata.n_obs, 12)
    assert adata.uns["sct"]["params"]["sct_method"] == "batch_key"


def test_sctransform_batch_key_parallel_matches_serial():
    adata_serial = _make_count_adata(n_cells=80, n_genes=100, seed=31)
    adata_parallel = _make_count_adata(n_cells=80, n_genes=100, seed=31)
    for ad in (adata_serial, adata_parallel):
        ad.obs["orig.ident"] = np.random.default_rng(31).choice(["A", "B", "C"], ad.n_obs)
    sctransform(
        adata_serial,
        batch_key="orig.ident",
        n_top_genes=12,
        n_cells=30,
        n_genes=40,
        seed=32,
        batch_n_jobs=1,
    )
    sctransform(
        adata_parallel,
        batch_key="orig.ident",
        n_top_genes=12,
        n_cells=30,
        n_genes=40,
        seed=32,
        batch_n_jobs=3,
    )
    assert set(adata_serial.uns["sct"]["batch_models"]) == set(
        adata_parallel.uns["sct"]["batch_models"]
    )
    assert adata_serial.uns["sct"]["variable_features"] == adata_parallel.uns["sct"][
        "variable_features"
    ]
    np.testing.assert_allclose(
        adata_serial.obsm["X_sct"],
        adata_parallel.obsm["X_sct"],
        rtol=1e-5,
        atol=1e-5,
    )
    assert adata_parallel.uns["sct"]["params"]["batch_n_jobs"] == 3


def test_sctransform_batch_n_jobs_invalid():
    adata = _make_count_adata(n_cells=40, n_genes=50, seed=33)
    adata.obs["orig.ident"] = "A"
    with pytest.raises(ValueError, match="batch_n_jobs"):
        sctransform(
            adata,
            batch_key="orig.ident",
            batch_n_jobs=0,
            n_top_genes=5,
            n_cells=20,
            n_genes=20,
        )


def test_select_integration_features_batch():
    adata = _make_count_adata(n_cells=80, n_genes=100, seed=25)
    adata.obs["orig.ident"] = np.random.default_rng(25).choice(["A", "B"], adata.n_obs)
    sctransform(
        adata,
        batch_key="orig.ident",
        n_top_genes=10,
        n_cells=30,
        n_genes=40,
        seed=26,
    )
    from trackcell.tl.integration import select_integration_features

    features = select_integration_features(
        adata,
        batch_key="orig.ident",
        n_features=8,
    )
    assert len(features) == 8
    assert all(g in adata.var_names for g in features)


def test_glmGamPoi_r_available_removed():
    """R bridge removed — shanuz-inspired pure-Python v2 is always available."""
    import importlib

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("trackcell.tl._r_sctransform")


def test_pyglmGamPoi_available():
    from trackcell.tl._sctransform_v2 import pyglmGamPoi_available

    assert isinstance(pyglmGamPoi_available(), bool)


@pytest.mark.skipif(
    not __import__(
        "trackcell.tl._sctransform_v2", fromlist=["pyglmGamPoi_available"]
    ).pyglmGamPoi_available(),
    reason="pyglmGamPoi not installed",
)
def test_sctransform_vst_flavor_v2_pyglmGamPoi():
    adata = _make_count_adata(n_cells=80, n_genes=100, seed=29)
    sctransform(
        adata,
        n_top_genes=15,
        n_cells=40,
        n_genes=50,
        vst_flavor="v2",
        v2_backend="pyglmGamPoi",
        seed=30,
    )
    assert adata.uns["sct"]["params"]["vst_flavor"] == "v2"
    assert np.isfinite(adata.obsm["X_sct"]).all()



def test_sctransform_model_use_poisson():
    adata = _make_count_adata(n_cells=60, n_genes=80, seed=19)
    adata.obs["batch"] = np.random.default_rng(19).integers(0, 2, adata.n_obs)
    sctransform(
        adata,
        n_top_genes=10,
        n_cells=30,
        n_genes=40,
        vars_to_regress=["batch"],
        model_use="poisson",
        clip_range=(-np.inf, np.inf),
        seed=20,
    )
    assert adata.uns["sct"]["params"]["model_use"] == "poisson"
    assert np.isfinite(adata.obsm["X_sct"]).all()


def test_sctransform_fixed_subsample_indices():
    adata = _make_count_adata(n_cells=80, n_genes=120, seed=31)
    cells_step1 = adata.obs_names[:40].tolist()
    genes_step1 = adata.var_names[:50].tolist()
    sctransform(
        adata,
        n_top_genes=20,
        n_cells=40,
        n_genes=50,
        seed=99,
        cells_step1=cells_step1,
        genes_step1=genes_step1,
    )
    hvg_first = adata.uns["sct"]["variable_features"].copy()
    sctransform(
        adata,
        n_top_genes=20,
        n_cells=40,
        n_genes=50,
        seed=12345,
        cells_step1=cells_step1,
        genes_step1=genes_step1,
    )
    assert hvg_first == adata.uns["sct"]["variable_features"]
