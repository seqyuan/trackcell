"""Integrate batch pairs using anchor vectors (Seurat ``IntegrateData`` SCT path).

Parity with Seurat v4 ``FindWeightsC`` + ``IntegrateDataC``:

* Weight formula: ``1 - exp(-dist * score / (2/sd)^2)``, then column-normalised.
* Correction: ``corrected = expr_query - W^T @ integration_matrix`` (matrix multiply).

Memory notes (multi-slice / spatial)
-------------------------------------
* Weight construction runs in ``weight_query_chunk_size`` cell chunks to cap peak RAM.
* ``integration_dtype='float32'`` halves PCA / correction workspace (default).
* Output remains ``float32`` in ``stack_integrated``; set ``integration_dtype='float64'``
  for maximum Seurat parity on small benchmarks.
"""

from __future__ import annotations

from typing import Literal, Optional

import numba
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, diags
from sklearn.decomposition import PCA

from ._annoy_nn import nn_helper

DEFAULT_WEIGHT_QUERY_CHUNK_SIZE = 8192


def _group_anchors_by_query_cell(
    cell2: np.ndarray,
    n_query: int,
) -> tuple[np.ndarray, np.ndarray]:
    """CSR-like grouping: query cell id -> anchor row indices."""
    cell2 = np.asarray(cell2, dtype=np.int64)
    order = np.argsort(cell2, kind="mergesort")
    anchor_by_cell = order
    counts = np.bincount(cell2, minlength=n_query)
    offsets = np.empty(n_query + 1, dtype=np.int64)
    offsets[0] = 0
    if n_query:
        offsets[1:] = np.cumsum(counts)
    return offsets, anchor_by_cell


@numba.njit(cache=True)
def _build_weights_chunk_numba(
    query_start: int,
    nbr_idx: np.ndarray,
    dists: np.ndarray,
    anchor_query_cells: np.ndarray,
    scores: np.ndarray,
    offsets: np.ndarray,
    anchor_by_cell: np.ndarray,
    k_weight: int,
    denom: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build COO triplets for one query row chunk (Seurat FindWeightsC logic)."""
    n_chunk, k_nn = nbr_idx.shape
    max_nnz = n_chunk * k_weight
    rows = np.empty(max_nnz, dtype=np.int64)
    cols = np.empty(max_nnz, dtype=np.int64)
    data = np.empty(max_nnz, dtype=np.float64)
    nnz = 0

    for qi in range(n_chunk):
        query_j = query_start + qi
        anchor_count = 0
        for ni in range(k_nn):
            if anchor_count >= k_weight:
                break
            anchor_row = nbr_idx[qi, ni]
            if anchor_row < 0:
                continue
            query_cell = anchor_query_cells[anchor_row]
            start = offsets[query_cell]
            end = offsets[query_cell + 1]
            dist = dists[qi, ni]
            for idx in range(start, end):
                if anchor_count >= k_weight:
                    break
                ai = anchor_by_cell[idx]
                w = 1.0 - np.exp(-dist * scores[ai] / denom)
                rows[nnz] = ai
                cols[nnz] = query_j
                data[nnz] = w
                nnz += 1
                anchor_count += 1

    return rows[:nnz], cols[:nnz], data[:nnz]


def _build_weights_seurat_legacy(
    n_anchors: int,
    n_query: int,
    cell2: np.ndarray,
    scores: np.ndarray,
    anchor_query_cells: np.ndarray,
    nbr_idx: np.ndarray,
    dists: np.ndarray,
    *,
    k_weight: int,
    sd_weight: float,
) -> csr_matrix:
    """Original Python-loop implementation (reference / tests)."""
    cell_to_anchors: dict[int, list[int]] = {}
    for ai, qc in enumerate(cell2):
        cell_to_anchors.setdefault(int(qc), []).append(ai)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    denom = (2.0 / max(sd_weight, 1e-9)) ** 2

    for query_j in range(n_query):
        neighbors = nbr_idx[query_j]
        neigh_dists = dists[query_j]
        anchor_count = 0
        for ni in range(len(neighbors)):
            if anchor_count >= k_weight:
                break
            unique_anchor_cell = int(anchor_query_cells[neighbors[ni]])
            anchor_indices = cell_to_anchors.get(unique_anchor_cell)
            if anchor_indices is None:
                continue
            for ai in anchor_indices:
                if anchor_count >= k_weight:
                    break
                w = 1.0 - np.exp(-neigh_dists[ni] * scores[ai] / denom)
                rows.append(ai)
                cols.append(query_j)
                data.append(w)
                anchor_count += 1

    weights_sp = coo_matrix(
        (np.asarray(data, dtype=np.float64), (np.asarray(rows), np.asarray(cols))),
        shape=(n_anchors, n_query),
    ).tocsr()
    col_sums = np.asarray(weights_sp.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    return weights_sp @ diags(1.0 / col_sums, format="csr")


def _build_weights_seurat(
    n_anchors: int,
    n_query: int,
    cell2: np.ndarray,
    scores: np.ndarray,
    anchor_query_cells: np.ndarray,
    nbr_idx: np.ndarray,
    dists: np.ndarray,
    *,
    k_weight: int,
    sd_weight: float,
    query_chunk_size: int = DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
) -> csr_matrix:
    """Build anchor×query weight matrix (chunked + Numba, memory-friendly)."""
    cell2 = np.asarray(cell2, dtype=np.int64)
    scores = np.asarray(scores, dtype=np.float64)
    anchor_query_cells = np.asarray(anchor_query_cells, dtype=np.int64)
    nbr_idx = np.asarray(nbr_idx, dtype=np.int64)
    dists = np.asarray(dists, dtype=np.float64)

    offsets, anchor_by_cell = _group_anchors_by_query_cell(cell2, n_query)
    denom = (2.0 / max(sd_weight, 1e-9)) ** 2
    chunk = max(int(query_chunk_size), 1)

    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []

    for start in range(0, n_query, chunk):
        end = min(start + chunk, n_query)
        rows, cols, data = _build_weights_chunk_numba(
            start,
            nbr_idx[start:end],
            dists[start:end],
            anchor_query_cells,
            scores,
            offsets,
            anchor_by_cell,
            k_weight,
            denom,
        )
        if rows.size:
            row_parts.append(rows)
            col_parts.append(cols)
            data_parts.append(data)

    if not row_parts:
        return csr_matrix((n_anchors, n_query), dtype=np.float64)

    rows = np.concatenate(row_parts)
    cols = np.concatenate(col_parts)
    data = np.concatenate(data_parts)
    weights_sp = coo_matrix((data, (rows, cols)), shape=(n_anchors, n_query)).tocsr()
    col_sums = np.asarray(weights_sp.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    return weights_sp @ diags(1.0 / col_sums, format="csr")


def integrate_pair(
    data_ref: np.ndarray,
    data_query: np.ndarray,
    anchors: pd.DataFrame,
    *,
    k_weight: int = 100,
    sd_weight: float = 1.0,
    dims: int = 30,
    n_trees: int = 50,
    n_trees_weight: Optional[int] = None,
    weight_query_chunk_size: int = DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    integration_dtype: Literal["float32", "float64"] = "float32",
) -> np.ndarray:
    """
    Integrate query batch toward reference using anchors (``RunIntegration`` for one pair).

    ``anchors`` must use local indices: ``cell1`` in ref, ``cell2`` in query.

    Parameters
    ----------
    n_trees_weight
        Annoy trees for weight-PCA neighbor search. Defaults to ``10`` (faster than
        anchor finding). Pass ``50`` for closer Seurat parity on benchmarks.
    weight_query_chunk_size
        Query cells processed per weight-building chunk (limits peak RAM for large
        multi-slice merges).
    integration_dtype
        ``float32`` (default) uses less memory during PCA and correction; ``float64``
        matches Seurat numerics more closely.
    """
    if anchors.empty:
        return data_query.copy()

    work_dtype = np.float32 if integration_dtype == "float32" else np.float64
    trees_weight = 10 if n_trees_weight is None else n_trees_weight

    n_ref = data_ref.shape[0]
    n_query = data_query.shape[0]
    cell1 = anchors["cell1"].astype(int).to_numpy()
    cell2 = anchors["cell2"].astype(int).to_numpy()
    scores = (
        anchors["score"].to_numpy(dtype=np.float64)
        if "score" in anchors.columns
        else np.ones(len(anchors), dtype=np.float64)
    )

    data_ref_w = np.asarray(data_ref, dtype=work_dtype)
    data_query_w = np.asarray(data_query, dtype=work_dtype)

    integration_matrix = data_query_w[cell2] - data_ref_w[cell1]

    merged = np.vstack([data_ref_w, data_query_w])
    merged_centered = merged - merged.mean(axis=0, keepdims=True)
    n_comp = min(dims, merged.shape[0] - 1, merged.shape[1])
    n_comp = max(n_comp, 1)
    dr = PCA(n_components=n_comp, random_state=0).fit_transform(merged_centered)
    dr_query = dr[n_ref:]

    anchor_query_cells = np.unique(cell2)
    if len(anchor_query_cells) == 0:
        return data_query.copy()

    k_nn = min(k_weight, len(anchor_query_cells))
    anchor_dr = dr_query[anchor_query_cells].astype(np.float32, copy=False)
    dr_query_nn = dr_query.astype(np.float32, copy=False)
    search_k = max(trees_weight * k_nn, trees_weight * 2)
    nbr_idx, dists = nn_helper(
        anchor_dr,
        query=dr_query_nn,
        k=k_nn,
        n_trees=trees_weight,
        search_k=search_k,
    )

    weights = _build_weights_seurat(
        n_anchors=len(cell1),
        n_query=n_query,
        cell2=cell2,
        scores=scores,
        anchor_query_cells=anchor_query_cells,
        nbr_idx=nbr_idx,
        dists=dists,
        k_weight=k_weight,
        sd_weight=sd_weight,
        query_chunk_size=weight_query_chunk_size,
    )

    correction = weights.T @ integration_matrix
    corrected = data_query_w - correction
    return np.asarray(corrected, dtype=np.float64)


def stack_integrated(
    integrated: dict[str, np.ndarray],
    batch_labels: list[str],
    cell_offsets: dict[str, int],
) -> np.ndarray:
    n_cells = sum(integrated[b].shape[0] for b in batch_labels)
    n_feat = integrated[batch_labels[0]].shape[1]
    out = np.zeros((n_cells, n_feat), dtype=np.float32)
    for label in batch_labels:
        off = cell_offsets[label]
        block = integrated[label]
        out[off : off + block.shape[0]] = block
    return out
