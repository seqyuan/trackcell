"""Load GSE288946 RPCA integration benchmark artifacts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from trackcell.benchmark.sct.data import (
    BENCHMARK_FILTER,
    BENCHMARK_ROOT,
    BENCHMARK_SAMPLES,
    BENCHMARK_SEEDS,
    REF_SEED,
    load_benchmark_umi,
)

RPCA_BENCHMARK_ROOT = BENCHMARK_ROOT
RPCA_SAMPLES = BENCHMARK_SAMPLES
RPCA_SEEDS = BENCHMARK_SEEDS
RPCA_REF_SEED = REF_SEED
BATCH_KEY = "orig.ident"

N_INTEGRATION_FEATURES = 2000
N_SCT_HVG = 3000
DIMS = 30
K_ANCHOR = 5
K_SCORE = 30
K_WEIGHT = 100
N_TREES = 50
N_TREES_WEIGHT = 50
INTEGRATION_DTYPE = "float64"
WEIGHT_QUERY_CHUNK_SIZE = 8192


def r_export_dir(
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
) -> Path:
    return benchmark_root / "steps" / "rpca_benchmark" / "r" / f"seed_{seed}"


def load_merged_adata_from_r_filter(
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
    samples: tuple[str, ...] = RPCA_SAMPLES,
) -> sc.AnnData:
    """Load merged AnnData using R-exported cell lists (post QC)."""
    root = r_export_dir(seed, benchmark_root=benchmark_root)
    filter_dir = root / "01_filter"
    adatas: list[sc.AnnData] = []
    for sample in samples:
        keep = set((filter_dir / f"{sample}_cells.txt").read_text().strip().split("\n"))
        keep.discard("")
        umi, genes, cells, _ = load_benchmark_umi(sample, benchmark_root=benchmark_root)
        from scipy.sparse import csr_matrix

        ad = sc.AnnData(X=csr_matrix(umi.T))
        ad.obs_names = cells.astype(str)
        ad.var_names = genes.astype(str)
        ad = ad[ad.obs_names.isin(keep)].copy()
        ad.obs[BATCH_KEY] = sample
        ad.layers["counts"] = ad.X.copy()
        adatas.append(ad)
    merged = sc.concat(adatas, join="outer", label="batch", keys=list(samples))
    merged.obs[BATCH_KEY] = merged.obs[BATCH_KEY].astype(str)
    merged.layers["counts"] = merged.X.copy()
    return merged


def load_merged_adata_native(
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
    samples: tuple[str, ...] = RPCA_SAMPLES,
) -> sc.AnnData:
    """Load merged AnnData from native benchmark QC (should match R filter)."""
    from scipy.sparse import csr_matrix

    adatas: list[sc.AnnData] = []
    for sample in samples:
        umi, genes, cells, _ = load_benchmark_umi(sample, benchmark_root=benchmark_root)
        ad = sc.AnnData(X=csr_matrix(umi.T))
        ad.obs_names = cells.astype(str)
        ad.var_names = genes.astype(str)
        ad.obs[BATCH_KEY] = sample
        ad.layers["counts"] = ad.X.copy()
        adatas.append(ad)
    merged = sc.concat(adatas, join="outer", label="batch", keys=list(samples))
    merged.obs[BATCH_KEY] = merged.obs[BATCH_KEY].astype(str)
    merged.layers["counts"] = merged.X.copy()
    return merged


def batch_sizes_from_adata(adata: sc.AnnData, batch_key: str = BATCH_KEY) -> dict[str, int]:
    batches = adata.obs[batch_key].astype(str)
    return {label: int((batches == label).sum()) for label in batches.unique()}


def inject_r_sct_models(adata: sc.AnnData, seed: int, *, benchmark_root: Path = RPCA_BENCHMARK_ROOT) -> None:
    """Populate ``adata.uns['sct']`` from R Step 02 CSV exports."""
    sct_dir = r_export_dir(seed, benchmark_root=benchmark_root) / "02_sct"
    batch_models: dict[str, dict] = {}
    union_hvg: set[str] = set()

    for sample in RPCA_SAMPLES:
        hvg = [g for g in (sct_dir / f"{sample}_hvg.txt").read_text().strip().split("\n") if g]
        model_pars = pd.read_csv(sct_dir / f"{sample}_model_pars.csv", index_col=0)
        model_pars.index = model_pars.index.astype(str)
        scale_df = pd.read_csv(sct_dir / f"{sample}_residuals_hvg.csv", index_col=0)
        scale_df.index = scale_df.index.astype(str)
        scale_df.columns = scale_df.columns.astype(str)

        mask = adata.obs[BATCH_KEY] == sample
        sub = adata[mask]
        counts = sub.layers["counts"]
        if hasattr(counts, "toarray"):
            umi = np.asarray(counts.sum(axis=1)).ravel()
        else:
            umi = counts.sum(axis=1)
        cell_attr = pd.DataFrame(
            {"log_umi": np.log10(np.maximum(umi, 1.0))},
            index=sub.obs_names,
        )
        n_cells = int(mask.sum())
        clip = float(np.sqrt(n_cells / 30))

        gene_attr = pd.DataFrame(
            index=model_pars.index,
            data={
                "residual_mean": 0.0,
                "residual_variance": model_pars.index.map(
                    lambda g: float(scale_df.loc[g].var()) if g in scale_df.index else np.nan
                ),
                "gmean": np.nan,
            },
        )

        batch_models[sample] = {
            "model_str": "y ~ log_umi",
            "variable_features": hvg,
            "model_pars_fit": model_pars.to_dict(orient="split"),
            "gene_attr": gene_attr.to_dict(orient="index"),
            "cell_attr": cell_attr.to_dict(orient="index"),
            "scale_data": scale_df.to_dict(orient="split"),
            "arguments": {"sct.clip.range": [-clip, clip], "sct.method": "default"},
        }
        union_hvg.update(hvg)

    adata.uns["sct"] = {
        "params": {"batch_key": BATCH_KEY, "method": "r_export"},
        "batch_models": batch_models,
        "variable_features": sorted(union_hvg),
    }


def load_r_features(seed: int, *, benchmark_root: Path = RPCA_BENCHMARK_ROOT) -> list[str]:
    path = r_export_dir(seed, benchmark_root=benchmark_root) / "03_features" / "integration_features.txt"
    return [g for g in path.read_text().strip().split("\n") if g]


def load_r_prep(
    adata: sc.AnnData,
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
) -> pd.DataFrame:
    prep = pd.read_csv(
        r_export_dir(seed, benchmark_root=benchmark_root) / "04_prep" / "prep_residuals.csv",
        index_col=0,
    )
    prep.index = prep.index.astype(str)
    return prep.loc[adata.obs_names]


def load_r_anchors(seed: int, *, benchmark_root: Path = RPCA_BENCHMARK_ROOT) -> pd.DataFrame:
    return pd.read_csv(r_export_dir(seed, benchmark_root=benchmark_root) / "05_anchors" / "anchors.csv")


def load_r_integrated(seed: int, *, benchmark_root: Path = RPCA_BENCHMARK_ROOT) -> pd.DataFrame:
    df = pd.read_csv(
        r_export_dir(seed, benchmark_root=benchmark_root) / "06_integrated" / "integrated_residuals.csv",
        index_col=0,
    )
    df.index = df.index.astype(str)
    return df


def load_integration_batch_umi(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
):
    """Load genes × cells UMI for one batch using R integration filter cell list."""
    root = r_export_dir(seed, benchmark_root=benchmark_root)
    keep = set((root / "01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")
    umi, genes, cells, cell_attr = load_benchmark_umi(sample, benchmark_root=benchmark_root)
    keep_idx = cells.isin(keep)
    umi = umi[:, keep_idx]
    cells = cells[keep_idx]
    cell_attr = cell_attr.loc[cells]
    return umi, genes, cells, cell_attr


def matrix_from_obsm(adata: sc.AnnData, key: str) -> pd.DataFrame:
    """Build cells × genes DataFrame using feature names stored in ``adata.uns[key]``."""
    features = list(adata.uns[key]["anchor_features"])
    return pd.DataFrame(adata.obsm[f"X_{key}"], index=adata.obs_names, columns=features)


def load_r_step1_data(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
) -> dict[str, object]:
    """Load R step-1 details (cells_step1, genes_step1, model_pars, gmean) for parity."""
    sct_dir = r_export_dir(seed, benchmark_root=benchmark_root) / "02_sct"
    cells_path = sct_dir / f"{sample}_cells_step1.txt"
    genes_path = sct_dir / f"{sample}_genes_step1.txt"
    if not cells_path.exists() or not genes_path.exists():
        return {}
    cells_step1 = pd.Index(cells_path.read_text().strip().split("\n"))
    genes_step1 = pd.Index(genes_path.read_text().strip().split("\n"))

    result: dict[str, object] = {
        "cells_step1": cells_step1,
        "genes_step1": genes_step1,
    }

    gmean_all_path = sct_dir / f"{sample}_genes_log_gmean_all.csv"
    if gmean_all_path.exists():
        gmean_all = pd.read_csv(gmean_all_path)
        result["genes_log_gmean"] = pd.Series(
            gmean_all.set_index("gene")["log10_gmean"].to_dict()
        )

    gmean_s1_path = sct_dir / f"{sample}_genes_log_gmean_step1.csv"
    if gmean_s1_path.exists():
        gmean_s1 = pd.read_csv(gmean_s1_path)
        result["genes_log_gmean_step1"] = pd.Series(
            gmean_s1.set_index("gene")["log10_gmean"].to_dict()
        )

    pars_s1_path = sct_dir / f"{sample}_model_pars_step1.csv"
    if pars_s1_path.exists():
        model_pars = pd.read_csv(pars_s1_path, index_col=0)
        model_pars.index = model_pars.index.astype(str)
        result["model_pars_step1"] = model_pars

    return result
