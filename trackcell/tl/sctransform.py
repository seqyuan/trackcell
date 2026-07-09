"""
SCTransform normalization for AnnData objects.

Replicates the Seurat 4.4.0 ``SCTransform`` workflow using regularized
negative binomial regression (Hafemeister & Satija, Genome Biology 2019).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sps
from anndata import AnnData

from ._sctransform import pack_sct_model, run_sctransform


def _serialize_batch_model(
    vst_out: dict[str, Any],
    clip_range: tuple[float, float],
    variable_features: Sequence[str],
    scale_data: Optional[pd.DataFrame] = None,
) -> dict[str, Any]:
    entry = pack_sct_model(vst_out, clip_range=clip_range)
    entry["cell_attr"] = vst_out["cell_attr"].to_dict(orient="index")
    entry["arguments"] = {
        **entry.get("arguments", {}),
        **vst_out.get("arguments", {}),
    }
    if isinstance(vst_out.get("model_pars_nonreg"), pd.DataFrame):
        entry["model_pars_nonreg"] = vst_out["model_pars_nonreg"].to_dict(orient="split")
    entry["model_pars_fit"] = vst_out["model_pars_fit"].to_dict(orient="split")
    entry["gene_attr"] = vst_out["gene_attr"].to_dict(orient="index")
    entry["variable_features"] = list(variable_features)
    sd = scale_data
    if sd is None and isinstance(vst_out.get("y"), pd.DataFrame):
        sd = vst_out["y"]
    if isinstance(sd, pd.DataFrame):
        entry["scale_data"] = sd.to_dict(orient="split")
    return entry


def _write_sctransform_results(
    adata: AnnData,
    assay_out: dict[str, Any],
    vst_out: dict[str, Any],
    *,
    key_added: str,
    store_residuals_in: str,
    n_top_genes: Optional[int],
    n_cells: int,
    n_genes: int,
    do_correct_umi: bool,
    do_center: bool,
    do_scale: bool,
    clip_range: Optional[tuple[float, float]],
    res_clip_range: Optional[tuple[float, float]],
    conserve_memory: bool,
    residual_features: Optional[Sequence[str]],
    vars_to_regress: Optional[Sequence[str]],
    latent_var: Optional[Sequence[str]],
    batch_var: Optional[str],
    latent_var_nonreg: Optional[Sequence[str]],
    residual_type: str,
    model_use: str,
    seed: int,
    layer: Optional[str],
    vst_flavor: Optional[str],
    save_reference_model: bool,
    reference_sct_model: Optional[dict],
    batch_models: Optional[dict[str, dict]] = None,
    batch_key: Optional[str] = None,
    method: str = "python",
) -> None:
    counts_cxg = assay_out["counts"].T.tocsr()
    data_cxg = assay_out["data"].T.tocsr()
    residuals_gxc = assay_out["scale.data"]
    hvg_names = assay_out["variable_features"]
    hvg_mask = np.asarray(adata.var_names.isin(hvg_names))

    adata.layers[f"{key_added}_counts"] = counts_cxg
    adata.layers[f"{key_added}_log1p"] = data_cxg
    adata.var[f"{key_added}_highly_variable"] = hvg_mask

    if store_residuals_in == "obsm":
        adata.obsm[f"X_{key_added}"] = np.asarray(
            residuals_gxc.T.to_numpy(), dtype=np.float32
        )
    elif store_residuals_in == "layer":
        hvg_indices = np.where(hvg_mask)[0]
        residual_layer = sps.lil_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)
        residual_layer[:, hvg_indices] = residuals_gxc.T.to_numpy()
        adata.layers[f"{key_added}_residuals"] = residual_layer.tocsr()
    else:
        raise ValueError("store_residuals_in must be 'obsm' or 'layer'.")

    gene_attr = vst_out["gene_attr"].loc[
        adata.var_names.intersection(vst_out["gene_attr"].index)
    ]
    resolved_clip = clip_range or (
        -np.sqrt(adata.n_obs / 30),
        np.sqrt(adata.n_obs / 30),
    )
    uns_entry: dict[str, Any] = {
        "params": {
            "n_top_genes": n_top_genes,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "do_correct_umi": do_correct_umi,
            "do_center": do_center,
            "do_scale": do_scale,
            "clip_range": clip_range,
            "res_clip_range": res_clip_range,
            "conserve_memory": conserve_memory,
            "residual_features": list(residual_features) if residual_features else None,
            "vars_to_regress": list(vars_to_regress) if vars_to_regress else None,
            "latent_var": list(latent_var) if latent_var else None,
            "batch_var": batch_var,
            "latent_var_nonreg": list(latent_var_nonreg) if latent_var_nonreg else None,
            "residual_type": residual_type,
            "model_use": model_use,
            "seed": seed,
            "layer": layer,
            "sct_method": assay_out["sct_method"],
            "vst_flavor": vst_flavor,
            "vst_method": assay_out.get("vst_method"),
            "method": method,
            "batch_key": batch_key,
        },
        "model_str": vst_out["model_str"],
        "variable_features": hvg_names.tolist(),
        "gene_attr": gene_attr.to_dict(orient="index"),
    }
    if batch_models is not None:
        uns_entry["batch_models"] = batch_models
    if save_reference_model and reference_sct_model is None and batch_key is None:
        uns_entry["reference_model"] = pack_sct_model(vst_out, clip_range=resolved_clip)
    adata.uns[key_added] = uns_entry


def _run_sctransform_on_counts(
    counts,
    adata: AnnData,
    *,
    variable_features_rv_th: float,
    n_top_genes: Optional[int],
    n_cells: int,
    n_genes: int,
    do_correct_umi: bool,
    do_center: bool,
    do_scale: bool,
    clip_range: Optional[tuple[float, float]],
    res_clip_range: Optional[tuple[float, float]],
    return_only_var_genes: bool,
    residual_features: Optional[Sequence[str]],
    conserve_memory: bool,
    reference_sct_model: Optional[dict],
    vars_to_regress: Optional[Sequence[str]],
    latent_var: Optional[Sequence[str]],
    batch_var: Optional[str],
    latent_var_nonreg: Optional[Sequence[str]],
    residual_type: str,
    model_use: str,
    seed: int,
    vst_flavor: Optional[str],
    cell_attr: Optional[pd.DataFrame] = None,
    latent_data: Optional[pd.DataFrame] = None,
    cells_step1: Optional[Sequence[str]] = None,
    genes_step1: Optional[Sequence[str]] = None,
    backend: str = "python",
    r_vst_export_dir: Optional[Union[str, Path]] = None,
):
    if sps.issparse(counts):
        counts = counts.tocsr()
    else:
        counts = np.asarray(counts)

    return run_sctransform(
        counts,
        adata.var_names.to_list(),
        adata.obs_names.to_list(),
        variable_features_n=n_top_genes,
        variable_features_rv_th=variable_features_rv_th,
        n_cells=n_cells,
        n_genes=n_genes,
        do_correct_umi=do_correct_umi,
        do_center=do_center,
        do_scale=do_scale,
        clip_range=clip_range,
        res_clip_range=res_clip_range,
        return_only_var_genes=return_only_var_genes,
        residual_features=residual_features,
        conserve_memory=conserve_memory,
        reference_sct_model=reference_sct_model,
        vars_to_regress=vars_to_regress,
        latent_data=latent_data,
        cell_attr=cell_attr if cell_attr is not None else adata.obs.copy(),
        latent_var=latent_var,
        batch_var=batch_var,
        latent_var_nonreg=latent_var_nonreg,
        residual_type=residual_type,
        model_use=model_use,
        seed=seed,
        vst_flavor=vst_flavor,
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        backend=backend,
        r_vst_export_dir=r_vst_export_dir,
    )


def sctransform(
    adata: AnnData,
    *,
    layer: Optional[str] = None,
    n_top_genes: Optional[int] = 3000,
    variable_features_rv_th: float = 1.3,
    n_cells: int = 5000,
    n_genes: int = 2000,
    do_correct_umi: bool = True,
    do_center: bool = True,
    do_scale: bool = False,
    clip_range: Optional[tuple[float, float]] = None,
    res_clip_range: Optional[tuple[float, float]] = None,
    return_only_var_genes: bool = True,
    residual_features: Optional[Sequence[str]] = None,
    conserve_memory: bool = False,
    reference_sct_model: Optional[dict] = None,
    vars_to_regress: Optional[Sequence[str]] = None,
    latent_var: Optional[Sequence[str]] = None,
    batch_var: Optional[str] = None,
    batch_key: Optional[str] = None,
    latent_var_nonreg: Optional[Sequence[str]] = None,
    residual_type: str = "pearson",
    model_use: str = "linear",
    vst_flavor: Optional[str] = None,
    seed: int = 1448145,
    key_added: str = "sct",
    store_residuals_in: str = "obsm",
    save_reference_model: bool = True,
    inplace: bool = True,
    copy: bool = False,
    cells_step1: Optional[Sequence[str]] = None,
    genes_step1: Optional[Sequence[str]] = None,
    subsample_indices: Optional[dict[str, dict[str, Sequence[str]]]] = None,
    r_vst_export_dir: Optional[Union[str, Path]] = None,
    r_vst_exports: Optional[dict[str, Union[str, Path]]] = None,
    method: str = "python",
) -> Optional[AnnData]:
    """
    Normalize UMI counts with SCTransform (Seurat 4.4.0).

    Parameters
    ----------
    adata
        Annotated data matrix with raw UMI counts in ``adata.X`` or ``layer``.
    layer
        Layer containing raw counts. Uses ``adata.X`` when ``None``.
    batch_key
        Column in ``adata.obs`` identifying biological batches/replicates.
        Runs independent SCT models per batch (Seurat ``SplitObject`` equivalent)
        and stores them in ``adata.uns[key_added]['batch_models']``.
    vst_flavor
        When ``"v2"``, use SCT v2 defaults (``glmGamPoi_offset``, ``exclude_poisson``,
        ``min_variance='umi_median'``, ``n_cells=2000``). Requires pyglmGamPoi or
        R glmGamPoi; use ``method='r_fit'`` / ``method='r'`` when unavailable.
    batch_var
        Column in ``adata.obs`` for batch-aware VST fitting within a single object.
        Do not combine with ``batch_key``.
    cells_step1, genes_step1
        Fixed cell/gene subsample indices for VST step 1 (Seurat ``vst`` parity).
        When set, random subsampling is skipped for the provided dimension(s).
    subsample_indices
        Per-batch fixed subsample indices when using ``batch_key``. Each entry
        maps a batch label to ``{"cells_step1": [...], "genes_step1": [...]}``.
    method
        ``"python"`` (default) native VST; ``"r"`` full R ``sctransform::vst``;
        ``"r_step1"`` loads R step-1 ``model_pars`` then runs Python ksmooth +
        residuals; ``"r_fit"`` loads R ``model_pars_fit`` and runs Python
        residuals only.
    r_vst_export_dir
        Directory with R-exported VST intermediates (single-sample mode).
    r_vst_exports
        Per-batch export directories when using ``batch_key`` with
        ``method='r_step1'`` or ``method='r_fit'``.
    """
    if copy:
        adata = adata.copy()
    elif not inplace:
        raise ValueError("Use copy=True or inplace=True.")

    if batch_key is not None and batch_var is not None:
        raise ValueError("Use batch_key or batch_var, not both.")
    if subsample_indices is not None and batch_key is None:
        raise ValueError("subsample_indices requires batch_key.")
    if (cells_step1 is not None or genes_step1 is not None) and batch_key is not None:
        raise ValueError("Use subsample_indices with batch_key, not cells_step1/genes_step1.")
    if method not in {"python", "r", "r_step1", "r_fit"}:
        raise ValueError("method must be 'python', 'r', 'r_step1', or 'r_fit'.")
    if method in {"r_step1", "r_fit"} and batch_key is None and r_vst_export_dir is None:
        raise ValueError(f"method='{method}' requires r_vst_export_dir.")
    if method in {"r_step1", "r_fit"} and batch_key is not None and not r_vst_exports:
        raise ValueError(f"method='{method}' with batch_key requires r_vst_exports.")
    if method == "r" and (
        subsample_indices is not None or cells_step1 is not None or genes_step1 is not None
    ):
        raise ValueError("method='r' does not support cells_step1/genes_step1/subsample_indices.")
    if method in {"r_step1", "r_fit"} and (
        cells_step1 is not None or genes_step1 is not None or subsample_indices is not None
    ):
        raise ValueError(
            f"method='{method}' uses r_vst_export_dir(s); do not pass cells_step1/genes_step1."
        )

    if layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"layer '{layer}' not in adata.layers.")

    obs_cols = set(adata.obs.columns)
    if batch_key is not None and batch_key not in obs_cols:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")
    if batch_var is not None and batch_var not in obs_cols:
        raise ValueError(f"batch_var '{batch_var}' not found in adata.obs.")
    for name in list(latent_var or ()) + list(latent_var_nonreg or ()):
        if name not in obs_cols and name not in {
            "umi",
            "gene",
            "log_umi",
            "log_gene",
            "umi_per_gene",
            "log_umi_per_gene",
        }:
            raise ValueError(f"latent covariate '{name}' not found in adata.obs.")

    common_kwargs = dict(
        variable_features_rv_th=variable_features_rv_th,
        n_top_genes=n_top_genes,
        n_cells=n_cells,
        n_genes=n_genes,
        do_correct_umi=do_correct_umi,
        do_center=do_center,
        do_scale=do_scale,
        clip_range=clip_range,
        res_clip_range=res_clip_range,
        return_only_var_genes=return_only_var_genes,
        residual_features=residual_features,
        conserve_memory=conserve_memory,
        reference_sct_model=reference_sct_model,
        vars_to_regress=vars_to_regress,
        latent_var=latent_var,
        batch_var=batch_var,
        latent_var_nonreg=latent_var_nonreg,
        residual_type=residual_type,
        model_use=model_use,
        seed=seed,
        vst_flavor=vst_flavor,
    )

    if batch_key is not None:
        if reference_sct_model is not None:
            raise ValueError("reference_sct_model is not supported with batch_key.")
        batch_labels = adata.obs[batch_key].astype(str)
        batch_models: dict[str, dict] = {}
        per_batch: dict[str, tuple[dict, dict]] = {}
        union_hvg: set[str] = set()

        for batch_label in batch_labels.unique():
            mask = batch_labels == batch_label
            sub = adata[mask]
            counts = sub.layers[layer] if layer is not None else sub.X
            latent_data = sub.obs[list(vars_to_regress)].copy() if vars_to_regress else None
            batch_subsample = (subsample_indices or {}).get(str(batch_label), {})
            batch_export = (r_vst_exports or {}).get(str(batch_label)) if r_vst_exports else None
            assay_out, vst_out = _run_sctransform_on_counts(
                counts,
                sub,
                cell_attr=sub.obs.copy(),
                latent_data=latent_data,
                cells_step1=batch_subsample.get("cells_step1"),
                genes_step1=batch_subsample.get("genes_step1"),
                backend=method,
                r_vst_export_dir=batch_export,
                **common_kwargs,
            )
            resolved_clip = clip_range or (
                -np.sqrt(sub.n_obs / 30),
                np.sqrt(sub.n_obs / 30),
            )
            batch_models[batch_label] = _serialize_batch_model(
                vst_out,
                resolved_clip,
                assay_out["variable_features"].tolist(),
                scale_data=assay_out.get("scale.data"),
            )
            per_batch[batch_label] = (assay_out, vst_out)
            union_hvg.update(assay_out["variable_features"].tolist())

        combined_hvg = sorted(union_hvg)
        if n_top_genes is not None and len(combined_hvg) > n_top_genes:
            rv_scores: dict[str, float] = {}
            for entry in batch_models.values():
                gene_attr = pd.DataFrame(entry["gene_attr"]).T
                for gene, row in gene_attr.iterrows():
                    rv = row.get("residual_variance", np.nan)
                    if np.isfinite(rv):
                        rv_scores[gene] = max(rv_scores.get(gene, -np.inf), float(rv))
            combined_hvg = sorted(rv_scores, key=rv_scores.get, reverse=True)[:n_top_genes]

        gene_to_col = {g: i for i, g in enumerate(combined_hvg)}
        combined_residuals = np.zeros((adata.n_obs, len(combined_hvg)), dtype=np.float32)

        for batch_label in batch_labels.unique():
            mask = batch_labels == batch_label
            assay_out, _ = per_batch[batch_label]
            batch_genes = [g for g in combined_hvg if g in assay_out["scale.data"].index]
            if batch_genes:
                block = assay_out["scale.data"].loc[batch_genes].T.to_numpy(dtype=np.float32)
                cols = [gene_to_col[g] for g in batch_genes]
                combined_residuals[np.where(mask)[0][:, None], cols] = block

        if layer is not None:
            combined_counts = adata.layers[layer]
            combined_data = adata.layers[layer].copy()
        else:
            combined_counts = adata.X
            combined_data = adata.X.copy()
        if sps.issparse(combined_data):
            combined_data = combined_data.tocsr(copy=True)
            combined_data.data = np.log1p(combined_data.data)
        else:
            combined_data = np.log1p(np.asarray(combined_data))
        scale_data = pd.DataFrame(
            combined_residuals.T,
            index=combined_hvg,
            columns=adata.obs_names,
        )
        first_vst = per_batch[batch_labels.unique()[0]][1]
        _write_sctransform_results(
            adata,
            {
                "counts": combined_counts.T.tocsr() if sps.issparse(combined_counts) else sps.csr_matrix(combined_counts.T),
                "data": combined_data.T.tocsr() if sps.issparse(combined_data) else sps.csr_matrix(combined_data.T),
                "scale.data": scale_data,
                "variable_features": np.array(combined_hvg),
                "sct_method": "batch_key",
            },
            first_vst,
            key_added=key_added,
            store_residuals_in="obsm",
            layer=layer,
            batch_models=batch_models,
            batch_key=batch_key,
            save_reference_model=False,
            reference_sct_model=None,
            **{k: common_kwargs[k] for k in (
                "n_top_genes", "n_cells", "n_genes", "do_correct_umi", "do_center",
                "do_scale", "clip_range", "res_clip_range", "conserve_memory",
                "residual_features", "vars_to_regress", "latent_var", "batch_var",
                "latent_var_nonreg", "residual_type", "model_use", "seed", "vst_flavor",
            )},
            method=method,
        )
        adata.uns[key_added]["variable_features"] = combined_hvg
        return None if inplace else adata

    counts = adata.layers[layer] if layer is not None else adata.X
    latent_data = None
    if vars_to_regress:
        missing = [v for v in vars_to_regress if v not in adata.obs.columns]
        if missing:
            raise ValueError(f"vars_to_regress not found in adata.obs: {missing}")
        latent_data = adata.obs[list(vars_to_regress)].copy()

    assay_out, vst_out = _run_sctransform_on_counts(
        counts,
        adata,
        cell_attr=adata.obs.copy(),
        latent_data=latent_data,
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        backend=method,
        r_vst_export_dir=r_vst_export_dir,
        **common_kwargs,
    )
    _write_sctransform_results(
        adata,
        assay_out,
        vst_out,
        key_added=key_added,
        store_residuals_in=store_residuals_in,
        layer=layer,
        save_reference_model=save_reference_model,
        reference_sct_model=reference_sct_model,
        batch_key=None,
        **{k: common_kwargs[k] for k in (
            "n_top_genes", "n_cells", "n_genes", "do_correct_umi", "do_center",
            "do_scale", "clip_range", "res_clip_range", "conserve_memory",
            "residual_features", "vars_to_regress", "latent_var", "batch_var",
            "latent_var_nonreg", "residual_type", "model_use", "seed", "vst_flavor",
        )},
        method=method,
    )
    return None if inplace else adata
