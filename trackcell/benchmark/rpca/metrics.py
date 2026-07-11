"""Metrics for SCT + RPCA integration benchmark."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def matrix_gene_corr(
    r_df: pd.DataFrame,
    p_df: pd.DataFrame,
    *,
    max_genes: int = 500,
) -> dict[str, Any]:
    """Per-gene Pearson correlation across shared cells (cells × genes)."""
    cells = r_df.index.intersection(p_df.index)
    genes = list(r_df.columns.intersection(p_df.columns))[:max_genes]
    if len(cells) == 0 or len(genes) == 0:
        return {"n_cells": len(cells), "n_genes_compared": len(genes)}
    corrs: list[float] = []
    for g in genes:
        rv = r_df.loc[cells, g].astype(float).to_numpy()
        pv = p_df.loc[cells, g].astype(float).to_numpy()
        if np.std(rv) > 0 and np.std(pv) > 0:
            corrs.append(float(pearsonr(rv, pv).statistic))
    out: dict[str, Any] = {"n_cells": len(cells), "n_genes_compared": len(genes)}
    if corrs:
        out["gene_corr_median"] = float(np.median(corrs))
        out["gene_corr_mean"] = float(np.mean(corrs))
    return out


def batch_offsets(batch_labels: list[str], batch_sizes: dict[str, int]) -> dict[str, int]:
    offsets: dict[str, int] = {}
    off = 0
    for label in batch_labels:
        offsets[label] = off
        off += batch_sizes[label]
    return offsets


def _local_pairs_r_style(anchors: pd.DataFrame) -> set[tuple[tuple[int, int], tuple[int, int]]]:
    pairs: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for row in anchors.itertuples(index=False):
        a = (int(row.dataset1) - 1, int(row.cell1) - 1)
        b = (int(row.dataset2) - 1, int(row.cell2) - 1)
        pairs.add(tuple(sorted((a, b))))
    return pairs


def _local_pairs_global(
    anchors: pd.DataFrame,
    *,
    batch_labels: list[str],
    batch_sizes: dict[str, int],
) -> set[tuple[tuple[int, int], tuple[int, int]]]:
    offsets = batch_offsets(batch_labels, batch_sizes)

    def _to_local(global_idx: int) -> tuple[int, int]:
        for di, label in enumerate(batch_labels):
            start = offsets[label]
            end = start + batch_sizes[label]
            if start <= global_idx < end:
                return di, global_idx - start
        raise ValueError(f"cell index {global_idx} out of range")

    pairs: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for row in anchors.itertuples(index=False):
        pairs.add(tuple(sorted((_to_local(int(row.cell1)), _to_local(int(row.cell2))))))
    return pairs


def compare_anchors(
    r_anchors: pd.DataFrame,
    py_anchors: pd.DataFrame,
    *,
    batch_labels: list[str],
    batch_sizes: dict[str, int],
) -> dict[str, Any]:
    r_pairs = _local_pairs_r_style(r_anchors)
    py_pairs = _local_pairs_global(
        py_anchors, batch_labels=batch_labels, batch_sizes=batch_sizes
    )
    union = r_pairs | py_pairs
    inter = r_pairs & py_pairs
    return {
        "r_n_anchors": int(len(r_anchors)),
        "py_n_anchors": int(len(py_anchors)),
        "r_n_unique_pairs": len(r_pairs),
        "py_n_unique_pairs": len(py_pairs),
        "local_pair_jaccard": len(inter) / len(union) if union else 1.0,
        "local_pair_overlap": len(inter),
        "unique_pair_ratio_py_over_r": len(py_pairs) / max(len(r_pairs), 1),
    }


def r_anchors_to_global(
    r_anchors: pd.DataFrame,
    *,
    batch_labels: list[str],
    batch_offsets_map: dict[str, int],
) -> pd.DataFrame:
    """Convert Seurat local 1-based anchors to Python global 0-based indices."""
    rows: list[dict[str, float | int]] = []
    for row in r_anchors.itertuples(index=False):
        d1, d2 = int(row.dataset1) - 1, int(row.dataset2) - 1
        l1, l2 = batch_labels[d1], batch_labels[d2]
        rows.append(
            {
                "cell1": batch_offsets_map[l1] + int(row.cell1) - 1,
                "cell2": batch_offsets_map[l2] + int(row.cell2) - 1,
                "score": float(row.score),
                "dataset1": d1,
                "dataset2": d2,
            }
        )
    return pd.DataFrame(rows)


def summarize_distribution(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return {}
    return {
        "n": int(arr.size),
        "min": float(arr.min()),
        "q25": float(np.quantile(arr, 0.25)),
        "median": float(np.median(arr)),
        "mean": float(arr.mean()),
        "q75": float(np.quantile(arr, 0.75)),
        "max": float(arr.max()),
    }
