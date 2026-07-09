"""SelectIntegrationFeatures for SCT-normalized AnnData."""

from __future__ import annotations

from collections import Counter
from typing import Literal, Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from .._r_sctransform import resolve_gene_alias
from ._xgboost_integration import select_integration_features_xgb


def _median_vf_rank(gene: str, vf_lists: Sequence[Sequence[str]]) -> float:
    ranks = [vf.index(gene) + 1 for vf in vf_lists if gene in vf]
    return float(np.median(ranks)) if ranks else np.inf


def _select_integration_features_seurat(
    vf_lists: Sequence[Sequence[str]],
    *,
    n_features: int,
    gene_universe: Optional[Sequence[str]] = None,
) -> list[str]:
    """
    Seurat v4 ``SelectIntegrationFeatures`` for SCT objects.

    Rank genes by how many objects list them as ``VariableFeatures``, break ties
    by median rank within each object's variable-feature list.
    """
    if not vf_lists:
        return []

    counts = Counter(g for vf in vf_lists for g in vf)
    if not counts:
        return []

    var_features = pd.Series(dict(counts), dtype=int).sort_values(ascending=False)
    if gene_universe is not None:
        universe = set(gene_universe)
        var_features = var_features[var_features.index.isin(universe)]

    if var_features.empty:
        return []

    n_pick = min(n_features, len(var_features))
    tie_val = int(var_features.iloc[n_pick - 1])
    features = var_features[var_features > tie_val].index.tolist()
    if features:
        features = sorted(features, key=lambda g: _median_vf_rank(g, vf_lists))

    features_tie = var_features[var_features == tie_val].index.tolist()
    if len(features) < n_features and features_tie:
        tie_sorted = sorted(features_tie, key=lambda g: _median_vf_rank(g, vf_lists))
        features.extend(tie_sorted[: n_features - len(features)])

    return features[:n_features]


def select_integration_features(
    adata: AnnData,
    *,
    batch_key: str,
    sct_key: str = "sct",
    n_features: int = 2000,
    features: Optional[Sequence[str]] = None,
    method: Literal["seurat", "xgboost"] = "seurat",
    layer: Optional[str] = None,
    seed: int = 0,
) -> list[str]:
    """
    Choose anchor genes for SCT integration (Seurat ``SelectIntegrationFeatures``).

    Default (``method='seurat'``): genes most frequently selected as
    ``VariableFeatures`` across per-batch SCT models, with Seurat tie-breaking.

    Optional enhancement (``method='xgboost'``): same prefilter, then rank
    by XGBoost batch-discrimination importance on per-cell SCT residuals.
    """
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not found in adata.obs.")

    if features is not None:
        return list(features)

    if method == "xgboost":
        return select_integration_features_xgb(
            adata,
            batch_key=batch_key,
            sct_key=sct_key,
            n_features=n_features,
            layer=layer,
            seed=seed,
        )

    uns = adata.uns.get(sct_key, {})
    batch_models = uns.get("batch_models")
    if batch_models:
        batch_labels = adata.obs[batch_key].astype(str).unique()
        vf_lists: list[list[str]] = []
        gene_sets: list[set[str]] = []
        for batch_label in batch_labels:
            entry = batch_models.get(str(batch_label))
            if entry is None:
                raise ValueError(f"Missing SCT batch model for '{batch_label}'.")
            vf = entry.get("variable_features") or []
            vf_lists.append(
                [resolve_gene_alias(g, adata.var_names) or g for g in vf]
            )
            gene_attr = entry.get("gene_attr", {})
            if gene_attr:
                gene_sets.append(
                    {
                        canon
                        for g in gene_attr.keys()
                        if (canon := resolve_gene_alias(g, adata.var_names)) is not None
                    }
                )
        gene_universe = set.intersection(*gene_sets) if gene_sets else set(adata.var_names)
        return _select_integration_features_seurat(
            vf_lists,
            n_features=n_features,
            gene_universe=gene_universe,
        )

    vf = uns.get("variable_features")
    if vf:
        return _select_integration_features_seurat([list(vf)], n_features=n_features)

    raise ValueError(
        f"Run tcl.tl.sctransform with batch_key='{batch_key}' first, "
        "or provide features= explicitly."
    )
