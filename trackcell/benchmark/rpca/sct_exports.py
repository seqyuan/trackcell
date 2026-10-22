"""Offline helpers for previously exported R SCT/RPCA benchmark artifacts."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Literal, Optional, Sequence

import pandas as pd


def _rpca_data():
    from trackcell.benchmark.rpca.data import load_r_step1_data, r_export_dir

    return load_r_step1_data, r_export_dir


def rpca_step1_available(benchmark_root: Path | str, seed: int, sample: str) -> bool:
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=Path(benchmark_root)) / "02_sct"
    required = (
        sct_dir / f"{sample}_cells_step1.txt",
        sct_dir / f"{sample}_genes_step1.txt",
        sct_dir / f"{sample}_model_pars_step1.csv",
        sct_dir / f"{sample}_genes_log_gmean_step1.csv",
        sct_dir / f"{sample}_genes_log_gmean_all.csv",
    )
    return all(path.exists() for path in required)


def load_rpca_subsample_indices(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> dict[str, dict[str, pd.Index]]:
    out: dict[str, dict[str, pd.Index]] = {}
    load_r_step1_data, _ = _rpca_data()
    for batch in batches:
        step1 = load_r_step1_data(batch, seed, benchmark_root=Path(benchmark_root))
        if step1:
            out[str(batch)] = {
                "cells_step1": pd.Index(step1["cells_step1"]).astype(str),
                "genes_step1": pd.Index(step1["genes_step1"]).astype(str),
            }
    return out


def materialize_rpca_r_step1_export_dir(
    benchmark_root: Path | str,
    seed: int,
    sample: str,
    *,
    force_refresh: bool = False,
) -> Path:
    root = Path(benchmark_root)
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=root) / "02_sct"
    cache = root / ".cache" / "rpca_vst_exports" / f"seed_{seed}" / sample
    if force_refresh and cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True, exist_ok=True)

    def link_or_copy(src: Path, dst: Path) -> None:
        if dst.exists():
            return
        if not src.exists():
            raise FileNotFoundError(f"Missing RPCA SCT export: {src}")
        try:
            dst.symlink_to(src.resolve())
        except OSError:
            shutil.copy2(src, dst)

    link_or_copy(sct_dir / f"{sample}_cells_step1.txt", cache / "cells_step1.txt")
    link_or_copy(sct_dir / f"{sample}_genes_step1.txt", cache / "genes_step1.txt")
    link_or_copy(
        sct_dir / f"{sample}_model_pars_step1.csv",
        cache / "model_pars_step1.csv",
    )
    for source_name, target_name in (
        (f"{sample}_genes_log_gmean_step1.csv", "genes_log_gmean_step1.csv"),
        (f"{sample}_genes_log_gmean_all.csv", "genes_log_gmean.csv"),
    ):
        target = cache / target_name
        if target.exists():
            continue
        frame = pd.read_csv(sct_dir / source_name)
        if "gene" in frame.columns:
            frame = frame.set_index("gene")
        col = "log10_gmean" if "log10_gmean" in frame.columns else frame.columns[0]
        pd.DataFrame(
            {"gmean": frame[col].to_numpy(dtype=float)},
            index=frame.index.astype(str),
        ).to_csv(target)
    meta = cache / "meta.json"
    if not meta.exists():
        meta.write_text('{"min_variance": 0.04, "bw_adjust": 3, "gmean_eps": 1}')
    return cache


def resolve_r_vst_exports(
    adata_batches: Sequence[str],
    *,
    export_root: Path | str,
    export_seed: Optional[int] = None,
    layout: Literal["auto", "per_batch", "rpca_benchmark"] = "auto",
) -> tuple[dict[str, str], Literal["per_batch", "rpca_benchmark"]]:
    root = Path(export_root)
    resolved: Literal["per_batch", "rpca_benchmark"]
    if layout == "auto":
        resolved = "rpca_benchmark" if export_seed is not None else "per_batch"
    else:
        resolved = layout
    exports: dict[str, str] = {}
    if resolved == "rpca_benchmark":
        if export_seed is None:
            raise ValueError("export_seed is required for rpca_benchmark layout.")
        for batch in adata_batches:
            if not rpca_step1_available(root, export_seed, batch):
                raise FileNotFoundError(f"RPCA step-1 export missing for batch '{batch}'.")
            exports[str(batch)] = str(
                materialize_rpca_r_step1_export_dir(root, export_seed, batch)
            )
    else:
        for batch in adata_batches:
            export_dir = root / str(batch)
            if not export_dir.is_dir():
                raise FileNotFoundError(f"Missing R VST export for batch '{batch}'.")
            exports[str(batch)] = str(export_dir)
    return exports, resolved


def rpca_exports_support_r_step1(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> bool:
    return all(rpca_step1_available(benchmark_root, seed, batch) for batch in batches)


def rpca_scale_data_available(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> bool:
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=Path(benchmark_root)) / "02_sct"
    return all((sct_dir / f"{batch}_residuals_hvg.csv").exists() for batch in batches)


def patch_r_sct_scale_data(
    adata,
    *,
    sct_key: str,
    batch_key: str,
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> bool:
    batch_models = adata.uns.get(sct_key, {}).get("batch_models")
    if not batch_models:
        return False
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=Path(benchmark_root)) / "02_sct"
    obs_batches = adata.obs[batch_key].astype(str)
    patched = False
    for batch in batches:
        path = sct_dir / f"{batch}_residuals_hvg.csv"
        if not path.exists() or batch not in batch_models:
            continue
        residuals = pd.read_csv(path, index_col=0)
        residuals.index = residuals.index.astype(str)
        residuals.columns = residuals.columns.astype(str)
        obs_names = pd.Index(adata.obs_names[obs_batches == batch].astype(str))
        if len(obs_names.intersection(residuals.columns)) != len(obs_names):
            continue
        batch_models[batch]["scale_data"] = residuals.loc[:, obs_names].to_dict(
            orient="split"
        )
        patched = True
    return patched
