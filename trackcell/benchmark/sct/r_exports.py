"""Offline readers and alignment helpers for R SCT benchmark exports.

This module never runs R as part of the Trackcell SCT API. It only reads
previously exported benchmark artifacts and optionally detects an R installation
for live parity tests.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

from trackcell.tl._gene_utils import align_gene_list, build_gene_alias_lookup


def normalize_r_model_pars(model_pars: pd.DataFrame) -> pd.DataFrame:
    """Normalize R coefficient names and string gene indexes."""
    out = model_pars.copy()
    out.index = out.index.astype(str)
    rename = {
        col: col.strip("()")
        for col in out.columns
        if col.startswith("(") and col.endswith(")")
    }
    if rename:
        out = out.rename(columns=rename)
    return out


def align_r_gene_list(genes: Sequence[str], target_genes: pd.Index) -> list[str]:
    """Align Seurat duplicate suffixes to names in ``target_genes``."""
    return align_gene_list(genes, target_genes)


def align_r_model_pars(model_pars: pd.DataFrame, target_genes: pd.Index) -> pd.DataFrame:
    """Align an R model parameter table to Trackcell gene names."""
    model_pars = normalize_r_model_pars(model_pars)
    lookup = build_gene_alias_lookup(target_genes)
    rows: list[pd.Series] = []
    index: list[str] = []
    for gene, row in model_pars.iterrows():
        matched = next((lookup[a] for a in _suffix_aliases(str(gene)) if a in lookup), None)
        if matched is not None and matched not in index:
            rows.append(row)
            index.append(matched)
    if not rows:
        return model_pars.iloc[0:0]
    return pd.DataFrame(rows, index=pd.Index(index))


def _suffix_aliases(name: str) -> list[str]:
    aliases = [name]
    for sep, replacement in ((".", "-"), ("-", ".")):
        if sep in name:
            base, suffix = name.rsplit(sep, 1)
            if suffix.isdigit():
                aliases.append(f"{base}{replacement}{suffix}")
    return aliases


def read_r_gmean_series(path: Path) -> pd.Series:
    """Load an R-exported gene geometric-mean series."""
    df = pd.read_csv(path, index_col=0)
    col = "gmean" if "gmean" in df.columns else df.columns[0]
    return pd.Series(df[col].to_numpy(dtype=np.float64), index=df.index.astype(str))


def coerce_genes_log_gmean(series: pd.Series) -> pd.Series:
    """Normalize legacy linear gmean exports to the log10 scale."""
    values = np.asarray(series, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size and np.nanmax(finite) > 5:
        values = np.log10(np.maximum(values, 1e-9))
    return pd.Series(values, index=series.index.astype(str))


def align_r_gmean_series(series: pd.Series, target_genes: pd.Index) -> pd.Series:
    """Align an R gmean series to Trackcell gene names."""
    aligned = align_r_gene_list(series.index.astype(str).tolist(), target_genes)
    if len(aligned) != len(series):
        raise ValueError("Gene alignment length mismatch for R gmean series.")
    return pd.Series(series.to_numpy(dtype=np.float64), index=pd.Index(aligned))


def load_r_vst_export(export_dir: Path | str) -> dict[str, Any]:
    """Read previously exported R SCT intermediates for offline parity tests."""
    root = Path(export_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"R VST export directory not found: {root}")

    def read_lines(name: str, required: bool = True) -> Optional[list[str]]:
        path = root / name
        if not path.exists():
            if required:
                raise FileNotFoundError(f"Missing {path}")
            return None
        return [line for line in path.read_text().splitlines() if line]

    def read_model(name: str) -> Optional[pd.DataFrame]:
        path = root / name
        if not path.exists():
            return None
        return normalize_r_model_pars(pd.read_csv(path, index_col=0))

    meta_path = root / "meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    model_step1 = read_model("model_pars_step1.csv")
    if model_step1 is None:
        raise FileNotFoundError(f"Missing {root / 'model_pars_step1.csv'}")

    gene_attr_path = root / "gene_attr.csv"
    gene_attr = pd.read_csv(gene_attr_path, index_col=0) if gene_attr_path.exists() else None
    genes_vst = read_lines("genes_vst.txt", required=False)
    if genes_vst is None and gene_attr is not None:
        genes_vst = gene_attr.index.astype(str).tolist()

    gmean_path = root / "genes_log_gmean.csv"
    if gmean_path.exists():
        genes_log_gmean = coerce_genes_log_gmean(read_r_gmean_series(gmean_path))
    elif gene_attr is not None and "gmean" in gene_attr.columns:
        genes_log_gmean = coerce_genes_log_gmean(
            pd.Series(gene_attr["gmean"].to_numpy(), index=gene_attr.index.astype(str))
        )
    else:
        genes_log_gmean = None

    gmean_step_path = root / "genes_log_gmean_step1.csv"
    genes_log_gmean_step1 = (
        coerce_genes_log_gmean(read_r_gmean_series(gmean_step_path))
        if gmean_step_path.exists()
        else None
    )
    return {
        "model_pars_step1": model_step1,
        "model_pars_fit": read_model("model_pars_fit.csv"),
        "cells_step1": read_lines("cells_step1.txt"),
        "genes_step1": read_lines("genes_step1.txt"),
        "genes_vst": genes_vst,
        "genes_log_gmean": genes_log_gmean,
        "genes_log_gmean_step1": genes_log_gmean_step1,
        "gene_attr": gene_attr,
        "meta": meta,
        "export_dir": str(root),
    }


@lru_cache(maxsize=1)
def r_glmgampoi_available() -> bool:
    """Return whether an R installation with glmGamPoi exists for live tests."""
    command = shutil.which("Rscript")
    if command is None:
        return False
    try:
        proc = subprocess.run(
            [command, "-e", 'cat(requireNamespace("glmGamPoi", quietly=TRUE))'],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0 and proc.stdout.strip() == "TRUE"
