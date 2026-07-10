"""Run SCT benchmark: stepwise R→Python injection and end-to-end parity."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from trackcell.benchmark.sct.data import (
    BENCHMARK_ROOT,
    BENCHMARK_SAMPLES,
    BENCHMARK_SEEDS,
    REF_SEED,
    load_benchmark_umi,
    load_r_export,
    r_export_dir,
)
from trackcell.benchmark.sct.metrics import (
    compare_gene_attr,
    compare_model_pars,
    compare_residual_matrices,
    hvg_jaccard,
    jaccard,
    summarize_distribution,
)
from trackcell.tl._r_sctransform import (
    align_r_gmean_series,
    align_r_model_pars,
    normalize_r_model_pars,
)
from trackcell.tl._sctransform import (
    _reg_model_pars,
    _row_gmean_sparse,
    get_residuals,
    vst,
)
from trackcell.tl._sctransform_v2 import resolve_min_variance


def _genes_log_gmean(umi, genes: pd.Index, *, gmean_eps: float = 1.0) -> pd.Series:
    return pd.Series(
        np.log10(np.maximum(_row_gmean_sparse(umi, gmean_eps=gmean_eps), 1e-9)),
        index=genes,
    )


def _run_stepwise_injection(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
) -> dict[str, Any]:
    """R step outputs → Python next step; compare each handoff."""
    umi, genes, cells, cell_attr = load_benchmark_umi(sample, benchmark_root=benchmark_root)
    r = load_r_export(sample, seed, benchmark_root=benchmark_root)
    meta = r.get("meta") or {}
    cells_step1 = pd.Index(r["cells_step1"]).intersection(cells)
    genes_step1 = pd.Index(r["genes_step1"]).intersection(genes)
    gmean_eps = float(meta.get("gmean_eps", 1.0))
    bw_adjust = float(meta.get("bw_adjust", 3.0))

    genes_log_gmean = _genes_log_gmean(umi, genes, gmean_eps=gmean_eps)
    genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]

    r_step1 = normalize_r_model_pars(r["model_pars_step1"]).reindex(genes_step1)
    r_fit = normalize_r_model_pars(r["model_pars_fit"]).reindex(genes)
    r_gene = r["gene_attr"]
    r_res = pd.read_csv(r_export_dir(sample, seed, benchmark_root=benchmark_root) / "residuals.csv", index_col=0)
    r_hvg = set(
        (r_export_dir(sample, seed, benchmark_root=benchmark_root) / "hvg_top3000.txt")
        .read_text()
        .strip()
        .split("\n")
    )

    # Step 1: R step1 → Python ksmooth (model_pars_fit)
    py_fit, _ = _reg_model_pars(
        r_step1.copy(),
        genes_log_gmean_step1,
        genes_log_gmean,
        bw_adjust=bw_adjust,
        theta_regularization="od_factor",
        umi=umi,
        gmean_eps=gmean_eps,
        exclude_poisson=True,
        fix_slope=True,
        genes=genes,
        genes_step1=genes_step1,
        cells_step1=cells_step1,
    )
    step_r1_py_ksmooth = compare_model_pars(py_fit, r_fit)

    # Step 2: R model_pars_fit → Python residuals
    min_var = meta.get("min_variance")
    if min_var is None:
        min_var = resolve_min_variance("umi_median", umi)
    vst_stub = {
        "model_str": "y ~ log_umi",
        "model_pars_fit": r_fit,
        "model_pars_nonreg": None,
        "model_str_nonreg": "",
        "cell_attr": cell_attr,
    }
    py_res = get_residuals(vst_stub, umi, genes=genes, umi_genes=genes, min_variance=float(min_var))
    step_rfit_py_res = compare_residual_matrices(py_res, r_res)
    py_rv = py_res.var(axis=1)
    step_rfit_py_res["residual_variance"] = compare_series_dict(py_rv, r_gene["residual_variance"])
    py_hvg = set(py_rv.sort_values(ascending=False).index[:3000].astype(str))
    step_rfit_py_res["hvg_jaccard"] = jaccard(py_hvg, r_hvg)

    # Step 3: R step1 + R gmean → Python ksmooth (golden r_step1 path)
    out_r1 = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        model_pars_step1=align_r_model_pars(r["model_pars_step1"], genes),
        genes_log_gmean=align_r_gmean_series(r["genes_log_gmean"], genes),
        genes_log_gmean_step1=align_r_gmean_series(r["genes_log_gmean_step1"], genes),
        min_variance=meta.get("min_variance", -np.inf),
        return_corrected_umi=False,
        seed=seed,
    )
    step_r1_full = {
        "model_pars_fit": compare_model_pars(out_r1["model_pars_fit"], r_fit),
        "gene_attr": compare_gene_attr(out_r1["gene_attr"], r_gene),
        "hvg_jaccard": step_rfit_py_res["hvg_jaccard"],
        "note": "HVG from R fit → Py residuals (golden vst returns step1-sized gene_attr only)",
    }

    # Step 4: Native Python full vst
    out_native = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=seed,
    )
    native = {
        "model_pars_step1": compare_model_pars(out_native["model_pars"], r_step1),
        "model_pars_fit": compare_model_pars(out_native["model_pars_fit"], r_fit),
        "gene_attr": compare_gene_attr(out_native["gene_attr"], r_gene),
        "hvg_jaccard": hvg_jaccard(out_native["gene_attr"], r_hvg),
        "residuals": compare_residual_matrices(out_native["y"], r_res),
    }

    return {
        "sample": sample,
        "seed": seed,
        "n_cells": int(umi.shape[1]),
        "n_genes": int(umi.shape[0]),
        "r_step1_to_py_ksmooth": step_r1_py_ksmooth,
        "r_fit_to_py_residuals": step_rfit_py_res,
        "r_step1_golden_path": step_r1_full,
        "native_full_vst": native,
    }


def compare_series_dict(a: pd.Series, b: pd.Series) -> dict[str, float]:
    from trackcell.benchmark.sct.metrics import compare_series

    return compare_series(a, b)


def _run_randomness_analysis(
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    seeds: tuple[int, ...] = BENCHMARK_SEEDS,
    ref_seed: int = REF_SEED,
) -> dict[str, Any]:
    """R inter-seed HVG Jaccard vs Python (ref seed)."""
    r_inter: list[float] = []
    py_vs_ref: list[float] = []
    per_sample: dict[str, Any] = {}

    for sample in BENCHMARK_SAMPLES:
        hvg_by_seed: dict[int, set[str]] = {}
        for seed in seeds:
            path = r_export_dir(sample, seed, benchmark_root=benchmark_root) / "hvg_top3000.txt"
            if path.exists():
                hvg_by_seed[seed] = set(path.read_text().strip().split("\n")) - {""}

        for s1, s2 in combinations(seeds, 2):
            if s1 in hvg_by_seed and s2 in hvg_by_seed:
                r_inter.append(jaccard(hvg_by_seed[s1], hvg_by_seed[s2]))

        ref_hvg = hvg_by_seed.get(ref_seed, set())
        umi, genes, cells, cell_attr = load_benchmark_umi(sample, benchmark_root=benchmark_root)
        r_ref = load_r_export(sample, ref_seed, benchmark_root=benchmark_root)
        cells_step1 = pd.Index(r_ref["cells_step1"]).intersection(cells)
        genes_step1 = pd.Index(r_ref["genes_step1"]).intersection(genes)
        out = vst(
            umi,
            genes,
            cells,
            cell_attr=cell_attr,
            vst_flavor="v2",
            cells_step1=cells_step1,
            genes_step1=genes_step1,
            return_corrected_umi=False,
            seed=ref_seed,
        )
        j = hvg_jaccard(out["gene_attr"], ref_hvg)
        py_vs_ref.append(j)
        per_sample[sample] = {
            "py_vs_r_ref_hvg_jaccard": j,
            "r_ref_n_hvg": len(ref_hvg),
        }

    inter_stats = summarize_distribution(r_inter)
    py_stats = summarize_distribution(py_vs_ref)
    py_median = py_stats.get("median")
    within_range = None
    within_iqr = None
    if py_median is not None and inter_stats:
        within_range = bool(inter_stats["min"] <= py_median <= inter_stats["max"])
        within_iqr = bool(inter_stats["q25"] <= py_median <= inter_stats["q75"])

    return {
        "seeds": list(seeds),
        "ref_seed": ref_seed,
        "r_inter_seed_hvg_jaccard": inter_stats,
        "py_vs_r_ref_hvg_jaccard": py_stats,
        "py_median_within_r_seed_range": within_range,
        "py_median_within_r_inter_seed_iqr": within_iqr,
        "per_sample": per_sample,
    }


def run_sct_benchmark_suite(
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    seeds: tuple[int, ...] = BENCHMARK_SEEDS,
    samples: tuple[str, ...] = BENCHMARK_SAMPLES,
    ref_seed: int = REF_SEED,
) -> dict[str, Any]:
    """Full benchmark payload (stepwise + randomness)."""
    stepwise: list[dict[str, Any]] = []
    missing: list[str] = []

    for seed in seeds:
        for sample in samples:
            export_path = r_export_dir(sample, seed, benchmark_root=benchmark_root) / "meta.json"
            if not export_path.exists():
                missing.append(f"{sample}/seed_{seed}")
                continue
            stepwise.append(_run_stepwise_injection(sample, seed, benchmark_root=benchmark_root))

    randomness = _run_randomness_analysis(
        benchmark_root=benchmark_root, seeds=seeds, ref_seed=ref_seed
    )

    ref_runs = [s for s in stepwise if s["seed"] == ref_seed]
    pooled_native_hvg = [s["native_full_vst"]["hvg_jaccard"] for s in ref_runs]
    pooled_r1_hvg = [s["r_step1_golden_path"]["hvg_jaccard"] for s in ref_runs]

    return {
        "status": "ok" if not missing else "partial",
        "filter": {"min_cells_gene": 100, "min_features": 50},
        "benchmark_root": str(benchmark_root),
        "samples": list(samples),
        "seeds": list(seeds),
        "ref_seed": ref_seed,
        "missing_r_exports": missing,
        "stepwise": stepwise,
        "randomness": randomness,
        "summary": {
            "native_hvg_jaccard": summarize_distribution(pooled_native_hvg),
            "r_step1_golden_hvg_jaccard": summarize_distribution(pooled_r1_hvg),
            "n_stepwise_runs": len(stepwise),
        },
    }


def write_benchmark_json(payload: dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path
