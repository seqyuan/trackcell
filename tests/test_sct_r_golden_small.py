"""Fast R-parity check for SCT, using a small committed golden matrix.

The full Seurat/sctransform parity suite needs ``TRACKCELL_BENCHMARK_ROOT`` and
is therefore skipped in CI. This module ships a 40-gene x 120-cell deterministic
matrix plus the R ``sctransform::vst`` outputs (``tests/data/sct_golden``), so
every PR exercises the step-1 fit, the regularised fit, gene means and Pearson
residuals against R without any external data.

Regenerate the golden files with::

    Rscript tests/data/sct_golden/export.R tests/data/sct_golden

``n_cells`` / ``n_genes`` cover the whole matrix, so R's subsampling (and its
RNG) does not participate; the same matrix is passed to the Python implementation
in full.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from trackcell.tl._sctransform import vst

DATA = Path(__file__).resolve().parent / "data" / "sct_golden"

pytestmark = pytest.mark.skipif(
    not (DATA / "umi.csv").exists(),
    reason="SCT golden vectors missing; run tests/data/sct_golden/export.R",
)


def _read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / name, index_col=0)


def _align_pars(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.rename(columns={"(Intercept)": "Intercept"})
    return frame[["theta", "Intercept", "log_umi"]]


@pytest.fixture(scope="module")
def inputs():
    umi_df = _read_csv("umi.csv")
    umi = csr_matrix(umi_df.to_numpy(dtype=float))
    genes = pd.Index(umi_df.index.astype(str))
    cells = pd.Index(umi_df.columns.astype(str))
    cell_attr = pd.DataFrame(
        {"log_umi": np.log10(np.asarray(umi.sum(axis=0)).ravel())}, index=cells
    )
    return umi, genes, cells, cell_attr


@pytest.fixture(scope="module")
def result(inputs):
    umi, genes, cells, cell_attr = inputs
    return vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        n_cells=len(cells),
        n_genes=len(genes),
        vst_flavor="v2",
        return_corrected_umi=False,
    )


def test_step1_model_pars_match_r(inputs, result):
    """The Poisson offset fit (theta, intercept, log_umi) is R-exact."""
    _umi, genes, _cells, _cell_attr = inputs
    expected = _align_pars(_read_csv("model_pars.csv")).loc[genes]
    observed = _align_pars(result["model_pars"]).loc[genes]

    np.testing.assert_allclose(observed["theta"], expected["theta"], atol=1e-6)
    np.testing.assert_allclose(observed["Intercept"], expected["Intercept"], atol=1e-6)
    np.testing.assert_allclose(observed["log_umi"], expected["log_umi"], atol=1e-12)


def test_regularized_model_pars_match_r(inputs, result):
    _umi, genes, _cells, _cell_attr = inputs
    expected = _align_pars(_read_csv("model_pars_fit.csv")).loc[genes]
    observed = _align_pars(result["model_pars_fit"]).loc[genes]

    np.testing.assert_allclose(observed["theta"], expected["theta"], atol=5e-4)
    np.testing.assert_allclose(observed["Intercept"], expected["Intercept"], atol=5e-4)


def test_genes_log_gmean_matches_r(inputs, result):
    _umi, genes, _cells, _cell_attr = inputs
    expected = _read_csv("genes_log_gmean.csv")["gmean"].loc[genes]
    observed = result["genes_log_gmean_step1"].loc[genes]
    np.testing.assert_allclose(observed, expected, atol=1e-9)


def test_pearson_residuals_match_r(inputs, result):
    _umi, genes, cells, _cell_attr = inputs
    expected = _read_csv("residuals.csv").loc[genes, cells]
    observed = pd.DataFrame(result["y"], index=genes, columns=cells)
    diff = np.abs(observed.to_numpy() - expected.to_numpy())
    assert float(diff.max()) < 2e-3
    assert float(np.median(diff)) < 5e-5


def test_residual_variance_matches_r(inputs, result):
    """``gene_attr`` is what HVG selection reads; keep it within R tolerance."""
    _umi, genes, _cells, _cell_attr = inputs
    expected = _read_csv("gene_attr.csv").loc[genes]
    observed = result["gene_attr"].loc[genes]

    np.testing.assert_allclose(
        observed["residual_variance"],
        expected["residual_variance"],
        rtol=5e-3,
        atol=0,
    )
    np.testing.assert_allclose(observed["gmean"], expected["gmean"], rtol=1e-9)


def test_golden_metadata_matches():
    meta = pd.read_json(DATA / "meta.json", typ="series")
    umi = _read_csv("umi.csv")
    assert int(meta["n_genes"]) == umi.shape[0]
    assert int(meta["n_cells"]) == umi.shape[1]
    assert int(meta["gmean_eps"]) == 1
