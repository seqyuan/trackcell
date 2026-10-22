"""SCT normalization pipeline for batch integration (Seurat Steps 1–3).

Produces per-cell anchor-gene residuals in ``obsm`` for RPCA, Harmony, BBKNN,
or other integration backends.
"""

from __future__ import annotations

from typing import Literal, Optional, Sequence

import numpy as np
from anndata import AnnData

from ..sctransform import sctransform
from .prep_sct import prep_sct_integration
from .select_features import select_integration_features
from .sct_parity import seurat_sct_clip_range


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
            f"Prep not found at adata.uns['{prep_key}']. "
            "Run run_sct_integration first."
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

    **Opt-in only** — nothing is removed unless you call this. ``integrate_rpca``
    / Seurat-parity RPCA use ``obsm['X_{prep_key}']`` only; releasing caches does
    **not** change integrated output if prep is already stored.

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
    sct_vst_flavor: str | None = "v2",
    sct_v2_backend: Literal["python", "pyglmGamPoi"] = "python",
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
        resolved_clip = clip_range
        resolved_res_clip = res_clip_range
        if resolved_clip is None and resolved_res_clip is None:
            max_cells = max(
                int((adata.obs[batch_key].astype(str) == label).sum())
                for label in batch_labels
            )
            resolved_clip = seurat_sct_clip_range(max_cells)
            resolved_res_clip = resolved_clip

        sctransform(
            adata,
            layer=layer,
            batch_key=batch_key,
            vst_flavor=sct_vst_flavor,
            v2_backend=sct_v2_backend,
            n_top_genes=n_sct_hvg,
            seed=sct_seed,
            key_added=sct_key,
            clip_range=resolved_clip,
            res_clip_range=resolved_res_clip,
            batch_n_jobs=sct_batch_n_jobs,
        )
        has_sct = True

    if require_sct_models and not has_sct:
        raise ValueError(
            f"Run sctransform(..., batch_key='{batch_key}') before integration, "
            f"or set run_sct=True."
        )

    feat_method = feature_method if feature_method is not None else "seurat"
    if anchor_features is None:
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
        "feature_method": feat_method,
        "sct_v2_backend": sct_v2_backend,
        "sct_batch_n_jobs": sct_batch_n_jobs,
    }
    return adata
