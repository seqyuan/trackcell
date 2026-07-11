"""SCT normalization pipeline for batch integration (Seurat Steps 1–3).

Produces per-cell anchor-gene residuals in ``obsm`` for RPCA, Harmony, BBKNN,
or other integration backends.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional, Sequence, Union

import numpy as np
from anndata import AnnData

from .._r_sctransform import validate_r_vst_export_artifacts, vst_r_available
from .._sctransform_v2 import glmGamPoi_offset_available
from ..sctransform import sctransform
from .prep_sct import prep_sct_integration
from .select_features import select_integration_features
from .sct_parity import (
    load_rpca_subsample_indices,
    patch_r_sct_scale_data,
    resolve_r_vst_exports,
    rpca_exports_support_r_step1,
    rpca_scale_data_available,
    seurat_sct_clip_range,
)


def resolve_sct_integration_method(
    sct_method: Literal["python", "r", "r_step1", "r_fit", "auto"],
    *,
    r_vst_export_root: Optional[Union[str, Path]],
    r_vst_export_seed: Optional[int],
    batches: Sequence[str],
    prefer_r_step1: bool,
    inject_r_step1_subsample: bool,
) -> Literal["python", "r", "r_step1", "r_fit"]:
    """Resolve ``sct_method='auto'`` for integration."""
    if sct_method != "auto":
        return sct_method
    if r_vst_export_root is not None:
        if (
            prefer_r_step1
            and r_vst_export_seed is not None
            and rpca_exports_support_r_step1(r_vst_export_root, r_vst_export_seed, batches)
        ):
            return "r_step1"
        if (
            inject_r_step1_subsample
            and r_vst_export_seed is not None
            and rpca_exports_support_r_step1(r_vst_export_root, r_vst_export_seed, batches)
            and glmGamPoi_offset_available()
        ):
            return "python"
        return "r_fit"
    if vst_r_available():
        return "r"
    if glmGamPoi_offset_available():
        return "python"
    raise RuntimeError(
        "sct_method='auto' could not resolve an SCT backend. Provide r_vst_export_root "
        "(+ r_vst_export_seed for RPCA step-1 injection) for r_step1/r_fit, install R "
        "sctransform, or install pyglmGamPoi / R glmGamPoi."
    )


def sct_prep_matrix(
    adata: AnnData,
    prep_key: str = "sct_prep",
) -> tuple[np.ndarray, list[str]]:
    """
    Return SCT prep residuals and anchor feature names for external integrators.

    Use with Harmony, BBKNN, scVI, etc.::

        run_sct_integration(adata, batch_key='orig.ident')
        matrix, genes = sct_prep_matrix(adata)
        adata.obsm['X_pca'] = ...  # PCA on matrix, then harmony / bbknn
    """
    if prep_key not in adata.uns:
        raise KeyError(
            f"Prep not found at adata.uns['{prep_key}']. Run run_sct_integration first."
        )
    features = list(adata.uns[prep_key]["anchor_features"])
    matrix = np.asarray(adata.obsm[f"X_{prep_key}"], dtype=np.float64)
    return matrix, features


def release_sct_integration_cache(
    adata: AnnData,
    *,
    sct_key: str = "sct",
    prep_key: str = "sct_prep",
    what: Literal["scale_data", "batch_models", "sct_uns"] = "scale_data",
    require_prep: bool = True,
) -> dict[str, int | str]:
    """
    Explicitly drop SCT caches from ``adata.uns`` after prep is finalised.

    **Opt-in only** — nothing is removed unless you call this (or the manual
    one-liners in :doc:`batch_integration`). ``integrate_rpca`` / Seurat-parity
    RPCA use ``obsm['X_{prep_key}']`` only; releasing caches does **not** change
    integrated output if prep is already stored.

    Parameters
    ----------
    what
        ``'scale_data'`` — remove per-batch ``scale_data`` (largest cache).
        ``'batch_models'`` — remove all ``batch_models`` (keeps ``params``).
        ``'sct_uns'`` — remove entire ``adata.uns[sct_key]``.
    require_prep
        If True (default), require ``prep_key`` and ``obsm['X_{prep_key}']``.

    Returns
    -------
    dict
        Summary: ``what``, ``n_batches``, ``bytes_estimate`` (rough lower bound).
    """
    prep_obs = f"X_{prep_key}"
    if require_prep:
        if prep_key not in adata.uns:
            raise KeyError(
                f"Prep metadata missing at adata.uns['{prep_key}']. "
                "Run run_sct_integration / prep_sct_integration first."
            )
        if prep_obs not in adata.obsm:
            raise KeyError(
                f"Prep matrix missing at adata.obsm['{prep_obs}']. "
                "Run prep before releasing SCT caches."
            )

    uns = adata.uns.get(sct_key, {})
    batch_models = uns.get("batch_models", {})
    n_batches = len(batch_models) if isinstance(batch_models, dict) else 0
    bytes_estimate = 0

    if what == "scale_data":
        if not batch_models:
            raise KeyError(f"No batch_models under adata.uns['{sct_key}'].")
        for entry in batch_models.values():
            sd = entry.pop("scale_data", None)
            if sd is not None and isinstance(sd, dict):
                data = sd.get("data")
                if data is not None:
                    bytes_estimate += len(np.asarray(data, dtype=np.float64)) * 8
    elif what == "batch_models":
        if "batch_models" in uns:
            for entry in batch_models.values():
                sd = entry.get("scale_data")
                if isinstance(sd, dict) and "data" in sd:
                    bytes_estimate += len(np.asarray(sd["data"])) * 8
            del uns["batch_models"]
    elif what == "sct_uns":
        if sct_key in adata.uns:
            bytes_estimate = -1  # unknown
            del adata.uns[sct_key]
    else:
        raise ValueError("what must be 'scale_data', 'batch_models', or 'sct_uns'.")

    return {
        "what": what,
        "sct_key": sct_key,
        "n_batches": n_batches,
        "bytes_estimate": bytes_estimate,
    }


def run_sct_integration(
    adata: AnnData,
    *,
    batch_key: str,
    anchor_features: Optional[Sequence[str]] = None,
    n_features: int = 2000,
    sct_key: str = "sct",
    prep_key: str = "sct_prep",
    layer: Optional[str] = None,
    feature_method: Optional[Literal["seurat", "xgboost"]] = None,
    seed: int = 0,
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
    Run SCT + SelectIntegrationFeatures + PrepSCTIntegration (Seurat Steps 1–3).

  Stores per-batch SCT models in ``adata.uns[sct_key]`` and anchor-gene Pearson
    residuals in ``adata.obsm[f'X_{prep_key}']``. The prep matrix is suitable
    input for :func:`integrate_rpca`, Harmony, BBKNN, or other methods.

    See :func:`sct_prep_matrix` to extract the matrix and feature names.
    """
    if copy:
        adata = adata.copy()
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    batch_labels = adata.obs[batch_key].astype(str).unique().tolist()

    has_sct = "batch_models" in adata.uns.get(sct_key, {})
    if run_sct is None:
        run_sct = not has_sct
    if run_sct:
        resolved_sct_method = resolve_sct_integration_method(
            sct_method,
            r_vst_export_root=r_vst_export_root,
            r_vst_export_seed=r_vst_export_seed,
            batches=batch_labels,
            prefer_r_step1=prefer_r_step1,
            inject_r_step1_subsample=inject_r_step1_subsample,
        )
        r_vst_exports: dict[str, str] | None = None
        subsample_indices: dict[str, dict[str, Sequence[str]]] | None = None

        if (
            inject_r_step1_subsample
            and r_vst_export_root is not None
            and r_vst_export_seed is not None
        ):
            subsample_indices = load_rpca_subsample_indices(
                r_vst_export_root,
                r_vst_export_seed,
                batch_labels,
            )
            if not subsample_indices:
                subsample_indices = None

        if resolved_sct_method in {"r_step1", "r_fit"}:
            if r_vst_export_root is None:
                raise ValueError(
                    f"sct_method='{resolved_sct_method}' requires r_vst_export_root "
                    "(benchmark root or per-batch parent directory)."
                )
            r_vst_exports, _layout = resolve_r_vst_exports(
                batch_labels,
                export_root=r_vst_export_root,
                export_seed=r_vst_export_seed,
                layout=r_vst_export_layout,
            )
            for export_dir in r_vst_exports.values():
                validate_r_vst_export_artifacts(export_dir, backend=resolved_sct_method)

        resolved_clip = clip_range
        resolved_res_clip = res_clip_range
        if resolved_clip is None and resolved_res_clip is None:
            max_cells = max(
                int((adata.obs[batch_key].astype(str) == label).sum()) for label in batch_labels
            )
            resolved_clip = seurat_sct_clip_range(max_cells)
            resolved_res_clip = resolved_clip

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
            subsample_indices=subsample_indices,
            clip_range=resolved_clip,
            res_clip_range=resolved_res_clip,
            batch_n_jobs=sct_batch_n_jobs,
        )
        if (
            inject_r_scale_data
            and r_vst_export_root is not None
            and r_vst_export_seed is not None
            and rpca_scale_data_available(r_vst_export_root, r_vst_export_seed, batch_labels)
        ):
            patch_r_sct_scale_data(
                adata,
                sct_key=sct_key,
                batch_key=batch_key,
                benchmark_root=r_vst_export_root,
                seed=r_vst_export_seed,
                batches=batch_labels,
            )
        has_sct = True

    if require_sct_models and not has_sct:
        raise ValueError(
            f"Run sctransform(..., batch_key='{batch_key}') before integration, "
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

    adata.uns.setdefault(prep_key, {})
    adata.uns[prep_key]["integration_ready"] = True
    adata.uns[prep_key]["sct_key"] = sct_key
    adata.uns[prep_key]["params"] = {
        "batch_key": batch_key,
        "n_features": len(anchor_features),
        "feature_method": feature_method if feature_method is not None else "seurat",
        "sct_method": adata.uns.get(sct_key, {}).get("params", {}).get("method"),
        "r_vst_export_seed": r_vst_export_seed,
        "inject_r_step1_subsample": inject_r_step1_subsample,
        "inject_r_scale_data": inject_r_scale_data,
        "prefer_r_step1": prefer_r_step1,
        "sct_batch_n_jobs": sct_batch_n_jobs,
    }
    return adata
