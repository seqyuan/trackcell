"""SCT + RPCA integration on a single AnnData (Seurat v4 workflow)."""

from __future__ import annotations

from typing import Literal, Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path

from ..sctransform import sctransform
from .._r_sctransform import validate_r_vst_export_artifacts, vst_r_available
from .._sctransform_v2 import glmGamPoi_offset_available

from .anchors import find_pairwise_anchors, map_anchors_to_global, mirror_anchors
from .integrate_data import stack_integrated
from .sample_tree import integrate_batches_sample_tree
from .prep_sct import prep_sct_integration
from .select_features import select_integration_features
from ._rpca import BatchPCA, reciprocal_project, run_batch_pca


def _resolve_integrate_sct_method(
    sct_method: Literal["python", "r", "r_step1", "r_fit", "auto"],
    *,
    r_vst_export_root: Optional[Union[str, Path]],
) -> Literal["python", "r", "r_step1", "r_fit"]:
    """Resolve ``sct_method='auto'`` for integration (prefer ``r_fit`` when exports exist)."""
    if sct_method != "auto":
        return sct_method
    if r_vst_export_root is not None:
        return "r_fit"
    if vst_r_available():
        return "r"
    if glmGamPoi_offset_available():
        return "python"
    raise RuntimeError(
        "sct_method='auto' could not resolve an SCT backend. Provide r_vst_export_root "
        "for r_fit, install R sctransform, or install pyglmGamPoi / R glmGamPoi."
    )


def integrate_sct_rpca(
    adata: AnnData,
    *,
    batch_key: str,
    anchor_features: Optional[Sequence[str]] = None,
    n_features: int = 2000,
    sct_key: str = "sct",
    prep_key: str = "sct_prep",
    layer: Optional[str] = None,
    dims: int = 30,
    k_anchor: int = 5,
    k_score: int = 30,
    k_weight: int = 100,
    sd_weight: float = 1.0,
    l2_norm: bool = True,
    n_trees: int = 50,
    preserve_order: bool = False,
    method: Literal["seurat", "xgboost"] = "seurat",
    feature_method: Optional[Literal["seurat", "xgboost"]] = None,
    xgb_n_delta_pcs: int = 20,
    xgb_n_estimators: int = 150,
    xgb_max_depth: int = 4,
    seed: int = 0,
    key_added: str = "sct_integrated",
    run_prep: bool = True,
    require_sct_models: bool = True,
    run_sct: bool | None = None,
    sct_method: Literal["python", "r", "r_step1", "r_fit", "auto"] = "auto",
    r_vst_export_root: Optional[Union[str, Path]] = None,
    sct_vst_flavor: str | None = "v2",
    n_sct_hvg: int = 3000,
    sct_seed: int = 1448145,
    copy: bool = False,
) -> AnnData:
    """
    SCT integration with reciprocal PCA anchors (Seurat v4 workflow).

    Workflow
    --------
    1. ``sctransform(..., batch_key=..., vst.flavor='v2')`` — run automatically when
       ``run_sct=True`` or when no per-batch SCT models exist (``run_sct=None``).
       ``sct_method='auto'`` (default) prefers ``method='r_fit'`` when
       ``r_vst_export_root`` is set, otherwise ``method='r'`` when R
       ``sctransform`` is available, otherwise native Python when glmGamPoi is
       available. Use ``sct_method='r_fit'`` with ``r_vst_export_root`` to load
       R ``model_pars_fit`` per batch and run Python residuals only (recommended
       bridge until pyglmGamPoi shrinkage). Use ``sct_method='r_step1'`` to load
       R step-1 ``model_pars`` and run Python ksmooth + residuals (for ksmooth
       development). R exports for ``r_step1`` must include ``genes_log_gmean.csv``
       and ``genes_log_gmean_step1.csv`` (from ``export_sct_stepwise_r.R``).
    2. ``SelectIntegrationFeatures`` → ``anchor_features``
    3. ``PrepSCTIntegration`` → per-batch anchor residuals
    4. ``FindIntegrationAnchors(reduction='rpca', normalization.method='SCT')``
    5. ``IntegrateData(normalization.method='SCT')``

    The default path (``method='seurat'``, ``feature_method='seurat'``) follows
    Seurat v4 SCT + RPCA: Annoy MNN anchors, sample-tree merging, and anchor-weighted
    correction.

    ``method='xgboost'`` keeps the same Seurat RPCA anchor finding and merge tree,
    but replaces the IntegrateData correction step with an optional XGBoost
    low-rank delta model for potentially stronger batch mixing.

    ``feature_method='xgboost'`` optionally replaces ``SelectIntegrationFeatures``
    with an XGBoost importance ranking (Seurat variance prefilter still applied).

    Set ``require_sct_models=False`` when ``obsm['X_{prep_key}']`` was loaded from
    an external source (e.g. R ``PrepSCTIntegration`` export) and SCT models are
    not available in ``adata.uns``.

    Stores integrated anchor residuals in ``adata.obsm[f'X_{key_added}']`` and
    metadata in ``adata.uns[key_added]``.
    """
    if copy:
        adata = adata.copy()
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    has_sct = "batch_models" in adata.uns.get(sct_key, {})
    if run_sct is None:
        run_sct = not has_sct
    if run_sct:
        resolved_sct_method = _resolve_integrate_sct_method(
            sct_method,
            r_vst_export_root=r_vst_export_root,
        )
        r_vst_exports: dict[str, str] | None = None
        if resolved_sct_method in {"r_step1", "r_fit"}:
            if r_vst_export_root is None:
                raise ValueError(
                    f"sct_method='{resolved_sct_method}' requires r_vst_export_root "
                    "(directory containing per-batch export subfolders)."
                )
            root = Path(r_vst_export_root)
            batches = adata.obs[batch_key].astype(str).unique()
            r_vst_exports = {}
            for batch in batches:
                export_dir = root / batch
                if not export_dir.is_dir():
                    raise FileNotFoundError(f"Missing R VST export for batch '{batch}': {export_dir}")
                validate_r_vst_export_artifacts(export_dir, backend=resolved_sct_method)
                r_vst_exports[str(batch)] = str(export_dir)
        sctransform(
            adata,
            layer=layer,
            batch_key=batch_key,
            vst_flavor=sct_vst_flavor,
            n_top_genes=n_sct_hvg,
            seed=sct_seed,
            key_added=sct_key,
            method=resolved_sct_method,
            r_vst_exports=r_vst_exports,
        )
        has_sct = True

    uns = adata.uns.get(sct_key, {})
    if require_sct_models and not has_sct:
        raise ValueError(
            f"Run tcl.tl.sctransform(..., batch_key='{batch_key}') before integration, "
            f"or set run_sct=True."
        )

    if anchor_features is None:
        feat_method = feature_method if feature_method is not None else "seurat"
        anchor_features = select_integration_features(
            adata,
            batch_key=batch_key,
            sct_key=sct_key,
            n_features=n_features,
            method=feat_method,
            layer=layer,
            seed=seed,
        )
    anchor_features = list(anchor_features)

    if run_prep or prep_key not in adata.uns:
        prep_sct_integration(
            adata,
            batch_key=batch_key,
            anchor_features=anchor_features,
            sct_key=sct_key,
            layer=layer,
            key_added=prep_key,
        )

    prep_matrix = adata.obsm[f"X_{prep_key}"]
    if not np.isfinite(prep_matrix).all():
        raise ValueError(f"Non-finite values in obsm['X_{prep_key}'].")

    batches = adata.obs[batch_key].astype(str)
    batch_labels = batches.unique().tolist()
    batch_data: dict[str, np.ndarray] = {}
    batch_pca: dict[str, BatchPCA] = {}
    cell_offsets: dict[str, int] = {}
    global_indices: dict[str, np.ndarray] = {}
    label_to_idx = {label: i for i, label in enumerate(batch_labels)}
    offset = 0
    for label in batch_labels:
        mask = batches == label
        batch_data[label] = np.asarray(prep_matrix[mask], dtype=np.float64)
        batch_pca[label] = run_batch_pca(batch_data[label], n_components=dims, seed=seed)
        cell_offsets[label] = offset
        global_indices[label] = np.where(mask.to_numpy())[0]
        offset += int(mask.sum())

    # Pairwise RPCA + anchors
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
            "feature_method": feature_method if feature_method is not None else "seurat",
            "dims": dims,
            "k_anchor": k_anchor,
            "k_score": k_score,
            "k_weight": k_weight,
            "sd_weight": sd_weight,
            "l2_norm": l2_norm,
            "n_trees": n_trees,
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
