"""PrepSCTIntegration — extract anchor-gene residuals per batch."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from .._sctransform import _as_genes_by_cells, get_residuals


def _split_dict_to_dataframe(obj: dict | pd.DataFrame) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, dict) and "data" in obj and "index" in obj:
        return pd.DataFrame(obj["data"], index=obj["index"], columns=obj["columns"])
    return pd.DataFrame.from_dict(obj, orient="index")


def _unpack_vst_out(entry: dict) -> dict:
    """Rebuild a vst_out dict from stored uns entry."""
    model_pars_fit = _split_dict_to_dataframe(entry["model_pars_fit"])
    gene_attr = _split_dict_to_dataframe(entry["gene_attr"])
    cell_attr = _split_dict_to_dataframe(entry["cell_attr"])
    vst_out = {
        "model_str": entry["model_str"],
        "model_pars_fit": model_pars_fit,
        "gene_attr": gene_attr,
        "cell_attr": cell_attr,
        "model_str_nonreg": entry.get("model_str_nonreg", ""),
        "arguments": entry.get("arguments", {}),
    }
    if entry.get("model_pars_nonreg"):
        vst_out["model_pars_nonreg"] = _split_dict_to_dataframe(entry["model_pars_nonreg"])
    scale_data = entry.get("scale_data")
    if scale_data is not None:
        vst_out["y"] = _split_dict_to_dataframe(scale_data)
    return vst_out


def _sct_clip_range(entry: dict, n_cells: int) -> tuple[float, float]:
    args = entry.get("arguments", {})
    clip = args.get("sct.clip.range")
    if clip is not None and len(clip) == 2:
        return float(clip[0]), float(clip[1])
    bound = float(np.sqrt(n_cells / 30))
    return -bound, bound


def _batch_anchor_residuals(
    entry: dict,
    *,
    counts,
    var_names: Sequence[str],
    obs_names: Sequence[str],
    anchor_features: Sequence[str],
) -> pd.DataFrame:
    """Seurat GetResidualSCTModel parity: reuse scale.data, compute missing pearson residuals."""
    vst_out = _unpack_vst_out(entry)
    obs_index = pd.Index(obs_names)
    genes = pd.Index(anchor_features).intersection(vst_out["model_pars_fit"].index)
    if len(genes) == 0:
        raise ValueError("No anchor features overlap the SCT model.")

    scale_data = vst_out.get("y")
    out = pd.DataFrame(index=genes, columns=obs_index, dtype=np.float64)

    cached_genes = pd.Index([])
    if isinstance(scale_data, pd.DataFrame):
        cached_genes = genes.intersection(scale_data.index)
        common_cells = obs_index.intersection(scale_data.columns)
        if len(cached_genes) > 0 and len(common_cells) == len(obs_index):
            out.loc[cached_genes, obs_index] = scale_data.loc[cached_genes, obs_index].to_numpy()

    missing_genes = genes.difference(cached_genes)
    if len(missing_genes) > 0:
        umi = _as_genes_by_cells(counts, var_names, obs_names)
        clip_range = _sct_clip_range(entry, umi.shape[1])
        cell_attr = vst_out["cell_attr"].reindex(obs_index)
        res = get_residuals(
            vst_out,
            umi,
            missing_genes,
            pd.Index(var_names),
            cell_attr=cell_attr,
            res_clip_range=clip_range,
            residual_type="pearson",
        )
        res = res.reindex(index=missing_genes, columns=obs_index)
        res = res.sub(res.mean(axis=1), axis=0)
        out.loc[missing_genes, obs_index] = res.to_numpy()

    return out


def prep_sct_integration(
    adata: AnnData,
    *,
    batch_key: str,
    anchor_features: Sequence[str],
    sct_key: str = "sct",
    layer: Optional[str] = None,
    key_added: str = "sct_integrated",
) -> AnnData:
    """
    Compute Pearson residuals for anchor genes per batch (Seurat ``PrepSCTIntegration``).

    Matches Seurat ``GetResidualSCTModel`` / ``PrepSCTIntegration``: reuse stored
    per-batch ``scale.data`` for anchor genes already present, compute pearson
    residuals (row-centered) only for missing genes.

    Stores a cells x genes residual matrix in ``adata.obsm[f'X_{key_added}']`` and
    per-batch metadata in ``adata.uns[key_added]``.
    """
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    uns = adata.uns.get(sct_key, {})
    batch_models = uns.get("batch_models")
    if not batch_models:
        raise ValueError(
            f"Per-batch SCT models not found at adata.uns['{sct_key}']['batch_models']. "
            f"Run sctransform(..., batch_key='{batch_key}') first."
        )

    anchor_features = list(anchor_features)
    batches = adata.obs[batch_key].astype(str)
    residual_blocks: list[tuple[np.ndarray, np.ndarray]] = []
    batch_info: dict[str, dict] = {}

    for batch_label in batches.unique():
        mask = batches == batch_label
        if batch_label not in batch_models:
            raise ValueError(f"Missing SCT model for batch '{batch_label}'.")
        entry = batch_models[batch_label]
        sub = adata[mask]
        counts = sub.layers[layer] if layer is not None else sub.X
        res = _batch_anchor_residuals(
            entry,
            counts=counts,
            var_names=adata.var_names.to_list(),
            obs_names=sub.obs_names.to_list(),
            anchor_features=anchor_features,
        )
        full = np.full((mask.sum(), len(anchor_features)), np.nan, dtype=np.float64)
        gene_to_col = {g: i for i, g in enumerate(anchor_features)}
        for gene in res.index:
            if gene in gene_to_col:
                full[:, gene_to_col[gene]] = res.loc[gene].to_numpy()
        residual_blocks.append((np.where(mask)[0], full))
        batch_info[batch_label] = {"n_anchor_genes": int(res.shape[0])}

    out = np.full((adata.n_obs, len(anchor_features)), np.nan, dtype=np.float32)
    for rows, mat in residual_blocks:
        out[rows] = mat.astype(np.float32)

    adata.obsm[f"X_{key_added}"] = out
    adata.uns[key_added] = {
        "batch_key": batch_key,
        "anchor_features": anchor_features,
        "batch_info": batch_info,
        "sct_key": sct_key,
    }
    return adata
