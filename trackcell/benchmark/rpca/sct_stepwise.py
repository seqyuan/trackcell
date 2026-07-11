"""Per-batch SCT stepwise parity vs R integration ``SCTransform`` exports."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from trackcell.benchmark.rpca.data import (
    BATCH_KEY,
    N_SCT_HVG,
    RPCA_BENCHMARK_ROOT,
    RPCA_SAMPLES,
    load_integration_batch_umi,
    load_merged_adata_native,
    load_r_features,
    load_r_step1_data,
    r_export_dir,
)
from trackcell.benchmark.rpca.metrics import jaccard, matrix_gene_corr
from trackcell.tl._r_sctransform import normalize_r_model_pars
from trackcell.tl._sctransform import get_residuals, vst
from trackcell.tl._sctransform_v2 import resolve_min_variance
from trackcell.tl.integration.prep_sct import prep_sct_integration
from trackcell.tl.integration.sct_parity import seurat_sct_clip_range
from trackcell.tl.integration.select_features import select_integration_features
from trackcell.tl.sctransform import sctransform


def _load_r_batch_sct(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path,
) -> dict[str, Any]:
    sct_dir = r_export_dir(seed, benchmark_root=benchmark_root) / "02_sct"
    hvg = set((sct_dir / f"{sample}_hvg.txt").read_text().strip().split("\n")) - {""}
    model_pars = normalize_r_model_pars(
        pd.read_csv(sct_dir / f"{sample}_model_pars.csv", index_col=0)
    )
    model_pars.index = model_pars.index.astype(str)
    residuals = pd.read_csv(sct_dir / f"{sample}_residuals_hvg.csv", index_col=0)
    residuals.index = residuals.index.astype(str)
    residuals.columns = residuals.columns.astype(str)
    return {"hvg": hvg, "model_pars_fit": model_pars, "residuals_hvg": residuals}


def _native_hvg(gene_attr: pd.DataFrame, *, n_top: int = N_SCT_HVG) -> set[str]:
    return set(
        gene_attr["residual_variance"]
        .sort_values(ascending=False)
        .index[:n_top]
        .astype(str)
    )


def run_integration_sct_stepwise(
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
) -> dict[str, Any]:
    """
    Compare native Python SCT vs R integration ``02_sct`` per batch.

    Handoffs
    --------
    * native ``vst()`` HVG vs R ``SCTransform`` HVG
    * R ``model_pars_fit`` → Python ``get_residuals`` (HVG + residual matrix)
    * native SCT + R integration features → ``PrepSCTIntegration`` (prep only)
    * native SCT → ``SelectIntegrationFeatures`` vs R feature list
    """
    per_sample: list[dict[str, Any]] = []
    r_features = load_r_features(seed, benchmark_root=benchmark_root)

    for sample in RPCA_SAMPLES:
        umi, genes, cells, cell_attr = load_integration_batch_umi(
            sample, seed, benchmark_root=benchmark_root
        )
        r_batch = _load_r_batch_sct(sample, seed, benchmark_root=benchmark_root)
        r_hvg = r_batch["hvg"]
        r_fit = r_batch["model_pars_fit"]
        r_res = r_batch["residuals_hvg"]

        # --- Native vst with R step-1 injection when available ---
        # Uses R cells_step1 + genes_step1 + model_pars_step1 to get
        # identical step-1 subsampling, then Python ksmooth + residuals.
        r_step1 = load_r_step1_data(sample, seed, benchmark_root=benchmark_root)
        out_rstep1 = None
        rstep1_hvg: set[str] = set()
        rstep1_theta_corr: float | None = None
        if r_step1:
            extra: dict[str, object] = {}
            if "genes_log_gmean" in r_step1:
                extra["genes_log_gmean"] = r_step1["genes_log_gmean"]
            if "genes_log_gmean_step1" in r_step1:
                extra["genes_log_gmean_step1"] = r_step1["genes_log_gmean_step1"]
            out_rstep1 = vst(
                umi,
                genes,
                cells,
                cell_attr=cell_attr,
                vst_flavor="v2",
                return_corrected_umi=False,
                seed=seed,
                model_pars_step1=r_step1.get("model_pars_step1"),
                cells_step1=r_step1.get("cells_step1"),
                genes_step1=r_step1.get("genes_step1"),
                **extra,
            )
            rstep1_hvg = _native_hvg(out_rstep1["gene_attr"])
            py_fit = out_rstep1["model_pars_fit"]
            common = py_fit.index.intersection(r_fit.index)
            py_theta = py_fit.loc[common, "theta"]
            r_theta = r_fit.loc[common, "theta"]
            fin = np.isfinite(py_theta) & np.isfinite(r_theta)
            if fin.sum() > 10:
                rstep1_theta_corr = float(
                    np.corrcoef(
                        np.log10(py_theta[fin].clip(1e-7)),
                        np.log10(r_theta[fin].clip(1e-7)),
                    )[0, 1]
                )

        out_native = vst(
            umi,
            genes,
            cells,
            cell_attr=cell_attr,
            vst_flavor="v2",
            return_corrected_umi=False,
            seed=seed,
        )
        native_hvg = _native_hvg(out_native["gene_attr"])

        min_var = resolve_min_variance("umi_median", umi)
        vst_stub = {
            "model_str": "y ~ log_umi",
            "model_pars_fit": r_fit.reindex(genes),
            "model_pars_nonreg": None,
            "model_str_nonreg": "",
            "cell_attr": cell_attr,
        }
        py_res = get_residuals(
            vst_stub,
            umi,
            genes=genes,
            umi_genes=genes,
            min_variance=float(min_var),
            res_clip_range=seurat_sct_clip_range(umi.shape[1]),
        )
        shared_genes = py_res.index.intersection(r_res.index)
        shared_cells = py_res.columns.intersection(r_res.columns)
        res_corr = float("nan")
        if len(shared_genes) and len(shared_cells):
            pm = py_res.loc[shared_genes, shared_cells].to_numpy(dtype=np.float64)
            rm = r_res.loc[shared_genes, shared_cells].to_numpy(dtype=np.float64)
            mask = np.isfinite(pm) & np.isfinite(rm)
            if mask.sum() >= 3:
                res_corr = float(np.corrcoef(pm[mask].ravel(), rm[mask].ravel())[0, 1])

        py_rfit_hvg = _native_hvg(
            pd.DataFrame({"residual_variance": py_res.var(axis=1)}, index=py_res.index)
        )

        per_sample.append(
            {
                "sample": sample,
                "n_cells": int(umi.shape[1]),
                "n_genes": int(umi.shape[0]),
                "native_hvg_jaccard": jaccard(native_hvg, r_hvg),
                "r_fit_to_py_residual_corr": res_corr,
                "r_fit_to_py_hvg_jaccard": jaccard(py_rfit_hvg, r_hvg),
                "r_step1_to_py_hvg_jaccard": jaccard(rstep1_hvg, r_hvg) if rstep1_hvg else float("nan"),
                "r_step1_theta_corr": rstep1_theta_corr,
            }
        )

    # Merged: native sct + prep/features using R integration feature list
    adata = load_merged_adata_native(benchmark_root=benchmark_root)
    sctransform(
        adata,
        layer="counts",
        batch_key=BATCH_KEY,
        vst_flavor="v2",
        n_top_genes=N_SCT_HVG,
        seed=seed,
        key_added="sct",
    )
    py_features = select_integration_features(
        adata, batch_key=BATCH_KEY, n_features=len(r_features), method="seurat"
    )
    prep_sct_integration(
        adata,
        batch_key=BATCH_KEY,
        anchor_features=r_features,
        sct_key="sct",
        layer="counts",
        key_added="sct_prep",
    )
    from trackcell.benchmark.rpca.data import load_merged_adata_from_r_filter, load_r_prep, inject_r_sct_models

    adata_r = load_merged_adata_from_r_filter(seed, benchmark_root=benchmark_root)
    inject_r_sct_models(adata_r, seed, benchmark_root=benchmark_root)
    prep_sct_integration(
        adata_r,
        batch_key=BATCH_KEY,
        anchor_features=r_features,
        sct_key="sct",
        layer="counts",
        key_added="sct_prep",
    )
    prep_features = list(adata.uns["sct_prep"]["anchor_features"])
    py_prep = pd.DataFrame(adata.obsm["X_sct_prep"], index=adata.obs_names, columns=prep_features)
    r_prep = load_r_prep(adata_r, seed, benchmark_root=benchmark_root)

    merged = {
        "features_jaccard_native_vs_r": jaccard(set(py_features), set(r_features)),
        "prep_native_sct_r_features": matrix_gene_corr(r_prep, py_prep, max_genes=500),
    }

    native_hvgs = [s["native_hvg_jaccard"] for s in per_sample]
    rfit_hvgs = [s["r_fit_to_py_hvg_jaccard"] for s in per_sample]

    return {
        "seed": seed,
        "reference": "R integration SCTransform (steps/rpca_benchmark/r/seed_*/02_sct)",
        "per_sample": per_sample,
        "merged": merged,
        "summary": {
            "native_hvg_jaccard_median": float(np.median(native_hvgs)),
            "r_fit_hvg_jaccard_median": float(np.median(rfit_hvgs)),
            "prep_corr_median_native_sct_r_features": merged["prep_native_sct_r_features"].get(
                "gene_corr_median"
            ),
        },
    }
