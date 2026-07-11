"""Seurat SCT parity helpers for integration (clip range, RPCA step-1 exports)."""

from __future__ import annotations

import math
import shutil
from pathlib import Path
from typing import Literal, Optional, Sequence

import pandas as pd


def seurat_sct_clip_range(n_cells: int) -> tuple[float, float]:
    """Seurat ``SCTransform`` / ``sct.clip.range``: ``±sqrt(ncol / 30)``."""
    bound = float(math.sqrt(max(n_cells, 1) / 30.0))
    return (-bound, bound)


def _rpca_data():
    from ...benchmark.rpca.data import load_r_step1_data, r_export_dir

    return load_r_step1_data, r_export_dir


def rpca_step1_available(
    benchmark_root: Path | str,
    seed: int,
    sample: str,
) -> bool:
    """Return True when RPCA ``02_sct`` has step-1 artifacts for ``sample``."""
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=Path(benchmark_root)) / "02_sct"
    required = (
        sct_dir / f"{sample}_cells_step1.txt",
        sct_dir / f"{sample}_genes_step1.txt",
        sct_dir / f"{sample}_model_pars_step1.csv",
        sct_dir / f"{sample}_genes_log_gmean_step1.csv",
    )
    gmean_all = sct_dir / f"{sample}_genes_log_gmean_all.csv"
    return all(path.exists() for path in required) and gmean_all.exists()


def load_rpca_subsample_indices(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> dict[str, dict[str, pd.Index]]:
    """Per-batch ``cells_step1`` / ``genes_step1`` from RPCA ``02_sct`` exports."""
    out: dict[str, dict[str, pd.Index]] = {}
    load_r_step1_data, _ = _rpca_data()
    for batch in batches:
        step1 = load_r_step1_data(batch, seed, benchmark_root=Path(benchmark_root))
        if not step1:
            continue
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
    """
    Stage RPCA ``02_sct`` files into a per-batch directory compatible with
    ``load_r_vst_export`` / ``validate_r_vst_export_artifacts(backend='r_step1')``.
    """
    root = Path(benchmark_root)
    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=root) / "02_sct"
    cache = root / ".cache" / "rpca_vst_exports" / f"seed_{seed}" / sample
    if force_refresh and cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True, exist_ok=True)

    def _link_or_copy(src: Path, dst: Path) -> None:
        if dst.exists():
            return
        if not src.exists():
            raise FileNotFoundError(f"Missing RPCA SCT export: {src}")
        try:
            dst.symlink_to(src.resolve())
        except OSError:
            shutil.copy2(src, dst)

    _link_or_copy(sct_dir / f"{sample}_cells_step1.txt", cache / "cells_step1.txt")
    _link_or_copy(sct_dir / f"{sample}_genes_step1.txt", cache / "genes_step1.txt")
    _link_or_copy(sct_dir / f"{sample}_model_pars_step1.csv", cache / "model_pars_step1.csv")

    gmean_step1_dst = cache / "genes_log_gmean_step1.csv"
    if not gmean_step1_dst.exists():
        src = sct_dir / f"{sample}_genes_log_gmean_step1.csv"
        df = pd.read_csv(src)
        if "gene" in df.columns:
            df = df.set_index("gene")
        col = "log10_gmean" if "log10_gmean" in df.columns else df.columns[0]
        out = pd.DataFrame({"gmean": df[col].to_numpy(dtype=float)}, index=df.index.astype(str))
        out.to_csv(gmean_step1_dst)

    gmean_dst = cache / "genes_log_gmean.csv"
    if not gmean_dst.exists():
        src = sct_dir / f"{sample}_genes_log_gmean_all.csv"
        df = pd.read_csv(src)
        if "gene" in df.columns:
            df = df.set_index("gene")
        col = "log10_gmean" if "log10_gmean" in df.columns else df.columns[0]
        out = pd.DataFrame({"gmean": df[col].to_numpy(dtype=float)}, index=df.index.astype(str))
        out.to_csv(gmean_dst)

    meta_dst = cache / "meta.json"
    if not meta_dst.exists():
        meta_dst.write_text(
            '{"min_variance": 0.04, "bw_adjust": 3, "gmean_eps": 1, "source": "rpca_benchmark"}'
        )

    return cache


def resolve_r_vst_exports(
    adata_batches: Sequence[str],
    *,
    export_root: Path | str,
    export_seed: Optional[int] = None,
    layout: Literal["auto", "per_batch", "rpca_benchmark"] = "auto",
) -> tuple[dict[str, str], Literal["per_batch", "rpca_benchmark"]]:
    """
    Build ``r_vst_exports`` mapping for ``sctransform(..., batch_key=...)``.

    When ``export_seed`` is set (or layout is ``rpca_benchmark``), ``export_root`` is
    treated as the benchmark root and per-batch dirs are staged from ``02_sct``.
    """
    root = Path(export_root)
    resolved_layout = layout
    if layout == "auto":
        if export_seed is not None:
            resolved_layout = "rpca_benchmark"
        elif all((root / str(batch) / "model_pars_step1.csv").exists() for batch in adata_batches):
            resolved_layout = "per_batch"
        else:
            raise ValueError(
                "Could not infer r_vst_export layout. Pass r_vst_export_seed for RPCA "
                "benchmark exports, or use per-batch subdirectories under r_vst_export_root."
            )

    exports: dict[str, str] = {}
    if resolved_layout == "rpca_benchmark":
        if export_seed is None:
            raise ValueError("r_vst_export_seed is required for rpca_benchmark layout.")
        for batch in adata_batches:
            if not rpca_step1_available(root, export_seed, batch):
                raise FileNotFoundError(
                    f"RPCA step-1 export missing for batch '{batch}' "
                    f"(seed={export_seed}). Re-run export_rpca_benchmark_r.R."
                )
            exports[str(batch)] = str(
                materialize_rpca_r_step1_export_dir(root, export_seed, batch)
            )
    else:
        for batch in adata_batches:
            export_dir = root / batch
            if not export_dir.is_dir():
                raise FileNotFoundError(f"Missing R VST export for batch '{batch}': {export_dir}")
            exports[str(batch)] = str(export_dir)

    return exports, resolved_layout


def rpca_exports_support_r_step1(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> bool:
    return all(rpca_step1_available(benchmark_root, seed, b) for b in batches)


def rpca_scale_data_available(
    benchmark_root: Path | str,
    seed: int,
    batches: Sequence[str],
) -> bool:
    """Return True when RPCA ``02_sct`` has per-batch ``residuals_hvg.csv``."""
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
    """
    Overlay R integration ``residuals_hvg`` onto per-batch ``scale_data``.

    Matches Seurat ``PrepSCTIntegration`` cache path when Python SCT residuals
    drift from R ``SCTransform``.
    """
    uns = adata.uns.get(sct_key, {})
    batch_models = uns.get("batch_models")
    if not batch_models:
        return False

    _, r_export_dir = _rpca_data()
    sct_dir = r_export_dir(seed, benchmark_root=Path(benchmark_root)) / "02_sct"
    patched = False
    obs_batches = adata.obs[batch_key].astype(str)

    for batch in batches:
        path = sct_dir / f"{batch}_residuals_hvg.csv"
        if not path.exists() or batch not in batch_models:
            continue
        r_res = pd.read_csv(path, index_col=0)
        r_res.index = r_res.index.astype(str)
        r_res.columns = r_res.columns.astype(str)

        mask = obs_batches == batch
        obs_names = pd.Index(adata.obs_names[mask].astype(str))
        common_cells = obs_names.intersection(r_res.columns)
        if len(common_cells) != len(obs_names):
            continue
        entry = batch_models[batch]
        entry["scale_data"] = r_res.loc[:, obs_names].to_dict(orient="split")
        patched = True

    return patched
