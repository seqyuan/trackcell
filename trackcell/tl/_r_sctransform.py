"""Optional R / glmGamPoi bridge for SCT v2."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, issparse

_RSCRIPT = Path(__file__).resolve().parent / "rscripts" / "fit_glmGamPoi_offset.R"
_RUN_VST = Path(__file__).resolve().parent / "rscripts" / "run_vst.R"


def _candidate_rscript_cmds() -> list[list[str]]:
    cmds: list[list[str]] = []
    override = os.environ.get("TRACKCELL_RSCRIPT")
    if override:
        cmds.append(override.split())
    if shutil.which("conda"):
        cmds.append(["conda", "run", "-n", "st", "Rscript"])
    if shutil.which("Rscript"):
        cmds.append(["Rscript"])
    return cmds


def _rscript_has_glmGamPoi(cmd: list[str]) -> bool:
    try:
        proc = subprocess.run(
            [*cmd, "-e", 'cat(requireNamespace("glmGamPoi", quietly=TRUE))'],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0 and proc.stdout.strip() == "TRUE"


@lru_cache(maxsize=1)
def glmGamPoi_r_available() -> bool:
    """Return True when an Rscript with glmGamPoi can run."""
    return _rscript_cmd() is not None


@lru_cache(maxsize=1)
def _rscript_cmd() -> Optional[list[str]]:
    for cmd in _candidate_rscript_cmds():
        if _rscript_has_glmGamPoi(cmd):
            return cmd
    return None


def fit_glmGamPoi_offset_r(
    umi: csr_matrix,
    regressor_data: pd.DataFrame,
    gene_index: pd.Index,
    *,
    allow_inf_theta: bool = True,
    timeout: int = 3600,
) -> pd.DataFrame:
    """
    Call sctransform::fit_glmGamPoi_offset via Rscript for numerical parity with Seurat v2.
    """
    if "log_umi" not in regressor_data.columns:
        raise ValueError("glmGamPoi_offset requires log_umi in regressor_data.")

    cmd = _rscript_cmd()
    if cmd is None:
        raise RuntimeError(
            "Rscript with glmGamPoi not found. Set TRACKCELL_RSCRIPT, e.g. "
            "'conda run -n st Rscript'."
        )

    if not _RSCRIPT.exists():
        raise FileNotFoundError(f"Missing R helper script: {_RSCRIPT}")

    if issparse(umi):
        umi = umi.tocsr()
    else:
        umi = csr_matrix(umi)

    with tempfile.TemporaryDirectory(prefix="trackcell_ggp_") as tmp:
        tmp_path = Path(tmp)
        umi_path = tmp_path / "umi.mtx"
        genes_path = tmp_path / "genes.txt"
        cells_path = tmp_path / "cells.txt"
        attr_path = tmp_path / "cell_attr.csv"
        out_path = tmp_path / "model_pars.csv"

        mmwrite(str(umi_path), umi)
        genes_path.write_text("\n".join(map(str, gene_index)) + "\n")
        cells_path.write_text("\n".join(map(str, regressor_data.index)) + "\n")
        regressor_data.to_csv(attr_path)

        proc = subprocess.run(
            [
                *cmd,
                str(_RSCRIPT),
                str(umi_path),
                str(genes_path),
                str(cells_path),
                str(attr_path),
                str(out_path),
                "TRUE" if allow_inf_theta else "FALSE",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "glmGamPoi_offset R call failed.\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )

        model_pars = pd.read_csv(out_path, index_col=0)
        model_pars.index = model_pars.index.astype(str)
        return model_pars


def vst_r_available() -> bool:
    """Return True when full sctransform::vst can run via Rscript."""
    cmd = _rscript_cmd()
    if cmd is None or not _RUN_VST.exists():
        return False
    try:
        proc = subprocess.run(
            [
                *cmd,
                "-e",
                'cat(requireNamespace("sctransform", quietly=TRUE))',
            ],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0 and proc.stdout.strip() == "TRUE"


def run_vst_r(
    umi: csr_matrix,
    genes: pd.Index,
    cells: pd.Index,
    *,
    cell_attr: pd.DataFrame,
    n_cells: Optional[int] = 5000,
    n_genes: int = 2000,
    variable_features_n: Optional[int] = 3000,
    variable_features_rv_th: float = 1.3,
    do_correct_umi: bool = True,
    return_only_var_genes: bool = True,
    clip_range: Optional[tuple[float, float]] = None,
    residual_type: str = "pearson",
    vst_flavor: Optional[str] = "v2",
    seed: int = 1448145,
    timeout: int = 7200,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Run full sctransform::vst in R and return assay/vst dicts for run_sctransform.
    """
    if not vst_r_available():
        raise RuntimeError(
            "R sctransform backend unavailable. Install sctransform/glmGamPoi in R "
            "or set TRACKCELL_RSCRIPT, e.g. 'conda run -n st Rscript'."
        )

    cmd = _rscript_cmd()
    assert cmd is not None

    if issparse(umi):
        umi = umi.tocsr()
    else:
        umi = csr_matrix(umi)

    if clip_range is None:
        clip_range = (-np.sqrt(umi.shape[1] / 30), np.sqrt(umi.shape[1] / 30))

    config = {
        "seed": int(seed),
        "n_cells": int(n_cells) if n_cells is not None else None,
        "n_genes": int(n_genes),
        "vst_flavor": vst_flavor,
        "variable_features_n": variable_features_n,
        "variable_features_rv_th": float(variable_features_rv_th),
        "do_correct_umi": bool(do_correct_umi),
        "return_only_var_genes": bool(return_only_var_genes),
        "residual_type": residual_type,
        "clip_range": [float(clip_range[0]), float(clip_range[1])],
    }

    with tempfile.TemporaryDirectory(prefix="trackcell_vst_r_") as tmp:
        tmp_path = Path(tmp)
        umi_path = tmp_path / "umi.mtx"
        genes_path = tmp_path / "genes.txt"
        cells_path = tmp_path / "cells.txt"
        attr_path = tmp_path / "cell_attr.csv"
        config_path = tmp_path / "config.json"
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        mmwrite(str(umi_path), umi)
        genes_path.write_text("\n".join(map(str, genes)) + "\n")
        cells_path.write_text("\n".join(map(str, cells)) + "\n")
        cell_attr.loc[cells].to_csv(attr_path)
        config_path.write_text(json.dumps(config))

        proc = subprocess.run(
            [
                *cmd,
                str(_RUN_VST),
                str(umi_path),
                str(genes_path),
                str(cells_path),
                str(attr_path),
                str(config_path),
                str(out_dir),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "sctransform::vst R call failed.\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )

        meta = json.loads((out_dir / "meta.json").read_text())
        variable_features = pd.Index(
            (out_dir / "variable_features.txt").read_text().strip().split("\n")
        )
        variable_features = variable_features[variable_features != ""]

        scale_data = mmread(out_dir / "scale_data.mtx").tocsr()
        scale_genes = pd.Index(
            (out_dir / "scale_genes.txt").read_text().strip().split("\n")
        )
        scale_cells = pd.Index(
            (out_dir / "scale_cells.txt").read_text().strip().split("\n")
        )
        scale_df = pd.DataFrame(
            scale_data.toarray(),
            index=scale_genes,
            columns=scale_cells,
        )

        gene_attr = pd.read_csv(out_dir / "gene_attr.csv", index_col=0)
        gene_attr.index = gene_attr.index.astype(str)
        model_pars_fit = pd.read_csv(out_dir / "model_pars_fit.csv", index_col=0)
        model_pars_fit.index = model_pars_fit.index.astype(str)
        if "(Intercept)" in model_pars_fit.columns:
            model_pars_fit = model_pars_fit.rename(columns={"(Intercept)": "Intercept"})

        counts_path = out_dir / "counts.mtx"
        if meta.get("has_corrected_counts") and counts_path.exists():
            counts = mmread(counts_path).tocsr().astype(np.float64)
        else:
            counts = umi.astype(np.float64)

        data = counts.copy()
        data.data = np.log1p(data.data)

        vst_out: dict[str, Any] = {
            "y": scale_df,
            "model_str": meta.get("model_str", "y ~ log_umi"),
            "model_pars_fit": model_pars_fit,
            "gene_attr": gene_attr,
            "cell_attr": cell_attr.loc[cells],
            "arguments": {
                "method": meta.get("vst_method"),
                "min_variance": None,
                "vst_flavor": vst_flavor,
            },
        }
        if meta.get("has_corrected_counts"):
            vst_out["umi_corrected"] = counts

        assay_out = {
            "counts": counts,
            "data": data,
            "scale.data": scale_df,
            "variable_features": variable_features.to_numpy(),
            "sct_method": "r.vst",
            "vst_method": meta.get("vst_method"),
        }
        return assay_out, vst_out


def _gene_suffix_aliases(name: str) -> list[str]:
    """Return a gene name plus Seurat ``.N`` / ``-N`` duplicate suffix variants."""
    name = str(name)
    aliases = [name]
    if "." in name:
        base, suffix = name.rsplit(".", 1)
        if suffix.isdigit():
            aliases.append(f"{base}-{suffix}")
    if "-" in name:
        base, suffix = name.rsplit("-", 1)
        if suffix.isdigit():
            aliases.append(f"{base}.{suffix}")
    return aliases


def build_gene_alias_lookup(target_genes: pd.Index) -> dict[str, str]:
    """Map gene names and suffix aliases to canonical names in ``target_genes``."""
    lookup: dict[str, str] = {}
    for gene in map(str, target_genes):
        for alias in _gene_suffix_aliases(gene):
            lookup.setdefault(alias, gene)
    return lookup


def _gene_name_alias(name: str, candidates: set[str] | dict[str, str]) -> Optional[str]:
    """Map Seurat/scanpy duplicate suffix variants to a name in ``candidates``."""
    if isinstance(candidates, dict):
        for alias in _gene_suffix_aliases(name):
            hit = candidates.get(alias)
            if hit is not None:
                return hit
        return None
    if name in candidates:
        return name
    for alias in _gene_suffix_aliases(name)[1:]:
        if alias in candidates:
            return alias
    return None


def resolve_gene_alias(name: str, target_genes: pd.Index) -> Optional[str]:
    """Resolve ``name`` to a canonical entry in ``target_genes`` if possible."""
    return _gene_name_alias(str(name), build_gene_alias_lookup(target_genes))


def align_r_model_pars(model_pars: pd.DataFrame, target_genes: pd.Index) -> pd.DataFrame:
    """Reindex R-exported model parameters onto ``target_genes`` (scanpy var_names)."""
    model_pars = normalize_r_model_pars(model_pars)
    lookup = build_gene_alias_lookup(target_genes)
    rows: list[pd.Series] = []
    index: list[str] = []
    for gene, row in model_pars.iterrows():
        matched = _gene_name_alias(str(gene), lookup)
        if matched is not None and matched not in index:
            rows.append(row)
            index.append(matched)
    if not rows:
        return model_pars.iloc[0:0]
    return pd.DataFrame(rows, index=pd.Index(index))


def align_r_gene_list(genes: Sequence[str], target_genes: pd.Index) -> list[str]:
    lookup = build_gene_alias_lookup(target_genes)
    out: list[str] = []
    for gene in genes:
        matched = _gene_name_alias(str(gene), lookup)
        if matched is not None and matched not in out:
            out.append(matched)
    return out


def validate_r_vst_export_artifacts(
    export_dir: Path | str,
    *,
    backend: str,
) -> None:
    """Ensure R VST export directories contain files required by ``backend``."""
    root = Path(export_dir)
    if backend == "r_fit":
        required = ["model_pars_fit.csv", "model_pars_step1.csv"]
    elif backend == "r_step1":
        required = [
            "model_pars_step1.csv",
            "genes_log_gmean.csv",
            "genes_log_gmean_step1.csv",
        ]
    else:
        raise ValueError(f"Unsupported backend for R export validation: {backend}")
    missing = [name for name in required if not (root / name).exists()]
    if missing:
        raise FileNotFoundError(
            f"R VST export at {root} is missing required files for backend='{backend}': "
            f"{missing}. Re-export with export_sct_stepwise_r.R."
        )


def normalize_r_model_pars(model_pars: pd.DataFrame) -> pd.DataFrame:
    """Normalize R-exported VST parameter tables for Python downstream steps."""
    out = model_pars.copy()
    out.index = out.index.astype(str)
    rename = {col: col.strip("()") for col in out.columns if col.startswith("(") and col.endswith(")")}
    if rename:
        out = out.rename(columns=rename)
    if "Intercept" not in out.columns and "(Intercept)" in out.columns:
        out = out.rename(columns={"(Intercept)": "Intercept"})
    return out


def read_r_gmean_series(path: Path) -> pd.Series:
    """Load R-exported log10 geometric mean series (index = R gene names)."""
    df = pd.read_csv(path, index_col=0)
    col = "gmean" if "gmean" in df.columns else df.columns[0]
    return pd.Series(df[col].to_numpy(dtype=np.float64), index=df.index.astype(str))


def align_r_gmean_series(series: pd.Series, target_genes: pd.Index) -> pd.Series:
    """Align an R gmean series to ``target_genes`` (handles ``.`` / ``-`` aliases)."""
    aligned_genes = align_r_gene_list(series.index.astype(str).tolist(), target_genes)
    if len(aligned_genes) != len(series):
        raise ValueError("Gene alignment length mismatch when aligning R gmean series.")
    return pd.Series(series.to_numpy(dtype=np.float64), index=pd.Index(aligned_genes, dtype=str))


def load_r_vst_export(export_dir: Path | str) -> dict[str, Any]:
    """
    Load R-exported SCT VST intermediates for Python ksmooth/residual continuation.

    Expected layout (per sample or batch)::

        model_pars_step1.csv   # glmGamPoi + shrinkage (step 1)
        model_pars_fit.csv     # optional reference after ksmooth
        cells_step1.txt
        genes_step1.txt
        genes_vst.txt          # optional genes retained after R vst filtering
        genes_log_gmean.csv    # optional full-gene log10 gmean (R gene_attr$gmean)
        genes_log_gmean_step1.csv  # optional step1 log10 gmean
        gene_attr.csv          # optional fallback for gmean / genes_vst
        meta.json              # optional (min_variance, bw_adjust, ...)
    """
    root = Path(export_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"R VST export directory not found: {root}")

    def _read_lines(name: str) -> list[str]:
        path = root / name
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}")
        return [line for line in path.read_text().strip().split("\n") if line]

    def _read_model(name: str) -> Optional[pd.DataFrame]:
        path = root / name
        if not path.exists():
            return None
        return normalize_r_model_pars(pd.read_csv(path, index_col=0))

    meta: dict[str, Any] = {}
    meta_path = root / "meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())

    model_pars_step1 = _read_model("model_pars_step1.csv")
    if model_pars_step1 is None:
        raise FileNotFoundError(f"Missing {root / 'model_pars_step1.csv'}")

    gene_attr_path = root / "gene_attr.csv"
    gene_attr = pd.read_csv(gene_attr_path, index_col=0) if gene_attr_path.exists() else None

    genes_vst: Optional[list[str]] = None
    genes_vst_path = root / "genes_vst.txt"
    if genes_vst_path.exists():
        genes_vst = _read_lines("genes_vst.txt")
    elif gene_attr is not None:
        genes_vst = gene_attr.index.astype(str).tolist()

    genes_log_gmean: Optional[pd.Series] = None
    gmean_path = root / "genes_log_gmean.csv"
    if gmean_path.exists():
        genes_log_gmean = read_r_gmean_series(gmean_path)
    elif gene_attr is not None and "gmean" in gene_attr.columns:
        genes_log_gmean = pd.Series(
            gene_attr["gmean"].to_numpy(dtype=np.float64),
            index=gene_attr.index.astype(str),
        )

    genes_log_gmean_step1: Optional[pd.Series] = None
    gmean_step1_path = root / "genes_log_gmean_step1.csv"
    if gmean_step1_path.exists():
        genes_log_gmean_step1 = read_r_gmean_series(gmean_step1_path)

    return {
        "model_pars_step1": model_pars_step1,
        "model_pars_fit": _read_model("model_pars_fit.csv"),
        "cells_step1": _read_lines("cells_step1.txt"),
        "genes_step1": _read_lines("genes_step1.txt"),
        "genes_vst": genes_vst,
        "genes_log_gmean": genes_log_gmean,
        "genes_log_gmean_step1": genes_log_gmean_step1,
        "gene_attr": gene_attr,
        "meta": meta,
        "export_dir": str(root),
    }
