"""SCT + RPCA integration on a single AnnData (Seurat v4 workflow)."""

from __future__ import annotations

from typing import Literal, Optional, Sequence, Union

from anndata import AnnData
from pathlib import Path

from .integrate_data import DEFAULT_WEIGHT_QUERY_CHUNK_SIZE
from .rpca_integration import integrate_rpca
from .sct_integration import resolve_sct_integration_method, run_sct_integration

_resolve_integrate_sct_method = resolve_sct_integration_method


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
    n_trees_weight: Optional[int] = None,
    weight_query_chunk_size: int = DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    integration_dtype: Literal["float32", "float64"] = "float32",
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
    r_vst_export_seed: Optional[int] = None,
    r_vst_export_layout: Literal["auto", "per_batch", "rpca_benchmark"] = "auto",
    inject_r_step1_subsample: bool = True,
    inject_r_scale_data: bool = True,
    prefer_r_step1: bool = False,
    sct_vst_flavor: str | None = "v2",
    n_sct_hvg: int = 3000,
    sct_seed: int = 1448145,
    sct_batch_n_jobs: int = 1,
    clip_range: Optional[tuple[float, float]] = None,
    res_clip_range: Optional[tuple[float, float]] = None,
    copy: bool = False,
) -> AnnData:
    """
    SCT integration with reciprocal PCA anchors (Seurat v4 workflow).

    Convenience wrapper around :func:`run_sct_integration` (Steps 1–3) and
    :func:`integrate_rpca` (Steps 4–5). For Harmony / BBKNN, call
    ``run_sct_integration`` only and use ``sct_prep_matrix`` or
    ``adata.obsm['X_{prep_key}']``.

    See :func:`run_sct_integration` and :func:`integrate_rpca` for parameter
    details.
    """
    if copy:
        adata = adata.copy()

    adata = run_sct_integration(
        adata,
        batch_key=batch_key,
        anchor_features=anchor_features,
        n_features=n_features,
        sct_key=sct_key,
        prep_key=prep_key,
        layer=layer,
        feature_method=feature_method,
        seed=seed,
        run_prep=run_prep,
        require_sct_models=require_sct_models,
        run_sct=run_sct,
        sct_method=sct_method,
        r_vst_export_root=r_vst_export_root,
        r_vst_export_seed=r_vst_export_seed,
        r_vst_export_layout=r_vst_export_layout,
        inject_r_step1_subsample=inject_r_step1_subsample,
        inject_r_scale_data=inject_r_scale_data,
        prefer_r_step1=prefer_r_step1,
        sct_vst_flavor=sct_vst_flavor,
        n_sct_hvg=n_sct_hvg,
        sct_seed=sct_seed,
        sct_batch_n_jobs=sct_batch_n_jobs,
        clip_range=clip_range,
        res_clip_range=res_clip_range,
        copy=False,
    )

    resolved_features = anchor_features
    if resolved_features is None:
        resolved_features = adata.uns[prep_key]["anchor_features"]

    adata = integrate_rpca(
        adata,
        batch_key=batch_key,
        prep_key=prep_key,
        sct_key=sct_key,
        anchor_features=resolved_features,
        dims=dims,
        k_anchor=k_anchor,
        k_score=k_score,
        k_weight=k_weight,
        sd_weight=sd_weight,
        l2_norm=l2_norm,
        n_trees=n_trees,
        n_trees_weight=n_trees_weight,
        weight_query_chunk_size=weight_query_chunk_size,
        integration_dtype=integration_dtype,
        preserve_order=preserve_order,
        method=method,
        xgb_n_delta_pcs=xgb_n_delta_pcs,
        xgb_n_estimators=xgb_n_estimators,
        xgb_max_depth=xgb_max_depth,
        seed=seed,
        key_added=key_added,
        copy=False,
    )

    adata.uns[key_added]["params"]["feature_method"] = (
        feature_method if feature_method is not None else "seurat"
    )
    adata.uns[key_added]["params"]["r_vst_export_seed"] = r_vst_export_seed
    adata.uns[key_added]["params"]["inject_r_step1_subsample"] = inject_r_step1_subsample
    adata.uns[key_added]["params"]["inject_r_scale_data"] = inject_r_scale_data
    adata.uns[key_added]["params"]["prefer_r_step1"] = prefer_r_step1
    adata.uns[key_added]["params"]["sct_batch_n_jobs"] = sct_batch_n_jobs
    return adata
