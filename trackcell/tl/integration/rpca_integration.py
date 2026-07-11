"""Reciprocal PCA integration on SCT prep residuals (Seurat Steps 4–5)."""

from __future__ import annotations

from typing import Literal, Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from .anchors import find_pairwise_anchors, map_anchors_to_global, mirror_anchors
from .integrate_data import DEFAULT_WEIGHT_QUERY_CHUNK_SIZE, stack_integrated
from .sample_tree import integrate_batches_sample_tree
from ._rpca import reciprocal_project, run_batch_pca


def integrate_rpca(
    adata: AnnData,
    *,
    batch_key: str,
    prep_key: str = "sct_prep",
    sct_key: str = "sct",
    anchor_features: Optional[Sequence[str]] = None,
    dims: int = 30,
    k_anchor: int = 5,
    k_score: int = 30,
    k_weight: int = 100,
    sd_weight: float = 1.0,
    l2_norm: bool = True,
    n_trees: int = 50,
    n_trees_weight: Optional[int] = None,
    weight_query_chunk_size: int = DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    integration_dtype: Literal["float32", "float64"] = "float32",
    preserve_order: bool = False,
    method: Literal["seurat", "xgboost"] = "seurat",
    xgb_n_delta_pcs: int = 20,
    xgb_n_estimators: int = 150,
    xgb_max_depth: int = 4,
    seed: int = 0,
    key_added: str = "sct_integrated",
    copy: bool = False,
) -> AnnData:
    """
    RPCA anchor finding + IntegrateData on existing SCT prep residuals.

    Expects ``run_sct_integration`` (or manual ``prep_sct_integration``) to have
    stored anchor residuals in ``adata.obsm[f'X_{prep_key}']``.

    Workflow
    --------
    1. ``FindIntegrationAnchors(reduction='rpca', normalization.method='SCT')``
    2. ``IntegrateData(normalization.method='SCT')``

    Stores integrated residuals in ``adata.obsm[f'X_{key_added}']``.
    """
    if copy:
        adata = adata.copy()
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")
    prep_obs_key = f"X_{prep_key}"
    if prep_obs_key not in adata.obsm:
        raise ValueError(
            f"Missing {prep_obs_key}. Run run_sct_integration(..., prep_key='{prep_key}') "
            "or prep_sct_integration first."
        )

    if anchor_features is None:
        prep_uns = adata.uns.get(prep_key, {})
        anchor_features = prep_uns.get("anchor_features")
        if anchor_features is None:
            raise ValueError(
                f"anchor_features not found in adata.uns['{prep_key}']. "
                "Pass anchor_features explicitly or run run_sct_integration."
            )
    anchor_features = list(anchor_features)

    prep_matrix = np.asarray(adata.obsm[prep_obs_key], dtype=np.float64)
    if not np.isfinite(prep_matrix).all():
        raise ValueError(f"Non-finite values in obsm['{prep_obs_key}'].")

    batches = adata.obs[batch_key].astype(str)
    batch_labels = batches.unique().tolist()
    batch_data: dict[str, np.ndarray] = {}
    batch_pca: dict = {}
    cell_offsets: dict[str, int] = {}
    global_indices: dict[str, np.ndarray] = {}
    label_to_idx = {label: i for i, label in enumerate(batch_labels)}
    offset = 0
    for label in batch_labels:
        mask = batches == label
        batch_data[label] = prep_matrix[mask.to_numpy()]
        batch_pca[label] = run_batch_pca(batch_data[label], n_components=dims, seed=seed)
        cell_offsets[label] = offset
        global_indices[label] = np.where(mask.to_numpy())[0]
        offset += int(mask.sum())

    anchor_frames: list[pd.DataFrame] = []
    for i, label1 in enumerate(batch_labels):
        for j, label2 in enumerate(batch_labels):
            if i >= j:
                continue
            d1 = batch_data[label1]
            d2 = batch_data[label2]
            pca1 = batch_pca[label1]
            pca2 = batch_pca[label2]
            ref_emb, query_emb = reciprocal_project(
                d1, pca1, d2, pca2, dims=dims, l2_norm=l2_norm
            )
            n_ref = d1.shape[0]
            anchors_local = find_pairwise_anchors(
                ref_emb,
                query_emb,
                n_ref=n_ref,
                pca1_embedding=pca1.embeddings,
                pca2_embedding=pca2.embeddings,
                k_anchor=k_anchor,
                k_score=k_score,
                n_trees=n_trees,
            )
            if anchors_local.empty:
                continue
            anchors_global = map_anchors_to_global(
                anchors_local,
                offset1=cell_offsets[label1],
                offset2=cell_offsets[label2],
                n_local1=n_ref,
                dataset1=label_to_idx[label1],
                dataset2=label_to_idx[label2],
            )
            anchors_global["batch1"] = label1
            anchors_global["batch2"] = label2
            anchor_frames.append(anchors_global)

    if not anchor_frames:
        raise RuntimeError("No integration anchors found between batches.")

    anchors = pd.concat(anchor_frames, ignore_index=True)
    anchors = mirror_anchors(anchors)
    integrated_dict, merge_order = integrate_batches_sample_tree(
        batch_data,
        batch_labels,
        anchors,
        global_indices,
        k_weight=k_weight,
        sd_weight=sd_weight,
        dims=dims,
        n_trees=n_trees,
        n_trees_weight=n_trees_weight,
        weight_query_chunk_size=weight_query_chunk_size,
        integration_dtype=integration_dtype,
        preserve_order=preserve_order,
        correction_method=method,
        xgb_n_delta_pcs=xgb_n_delta_pcs,
        xgb_n_estimators=xgb_n_estimators,
        xgb_max_depth=xgb_max_depth,
        xgb_seed=seed,
    )
    integrated = stack_integrated(integrated_dict, batch_labels, cell_offsets)

    adata.obsm[f"X_{key_added}"] = integrated.astype(np.float32)
    adata.uns[key_added] = {
        "batch_key": batch_key,
        "anchor_features": anchor_features,
        "n_anchors": int(len(anchors)),
        "params": {
            "method": method,
            "dims": dims,
            "k_anchor": k_anchor,
            "k_score": k_score,
            "k_weight": k_weight,
            "sd_weight": sd_weight,
            "l2_norm": l2_norm,
            "n_trees": n_trees,
            "n_trees_weight": 10 if n_trees_weight is None else n_trees_weight,
            "weight_query_chunk_size": weight_query_chunk_size,
            "integration_dtype": integration_dtype,
            "nn_method": "annoy",
            "merge_method": "sample_tree",
            "preserve_order": preserve_order,
            "normalization.method": "SCT",
            "reduction": "rpca",
            "sct_method": adata.uns.get(sct_key, {}).get("params", {}).get("method"),
            "xgb_n_delta_pcs": xgb_n_delta_pcs,
            "xgb_n_estimators": xgb_n_estimators,
            "xgb_max_depth": xgb_max_depth,
        },
        "sample_tree": merge_order.tolist(),
        "anchors": anchors.to_dict(orient="list"),
        "sct_key": sct_key,
        "prep_key": prep_key,
    }
    return adata
