"""SelectIntegrationFeatures for SCT-normalized AnnData."""

from __future__ import annotations

from collections import Counter
from typing import Literal, Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from ..._logging import get_logger
from .._gene_utils import resolve_gene_alias
from ._xgboost_integration import select_integration_features_xgb

logger = get_logger(__name__)


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
    return_importance: bool = False,
):
    """
    Choose anchor genes for SCT integration.

    The two ``method`` values answer **different questions** and their results are
    not comparable. This is the most important thing to know about this function:

    ``method='seurat'`` (default) — Seurat ``SelectIntegrationFeatures``: genes
    most often selected as ``VariableFeatures`` across the per-batch SCT models,
    with Seurat tie-breaking. Selects genes that vary *reproducibly within*
    batches.

    ``method='xgboost'`` — genes that best *separate* batches, ranked by gain
    importance from a classifier trained on per-cell SCT residuals. On the same
    input this ranking points in close to the opposite direction to the default:
    a gene whose variance is purely batch-driven ranks high here and low under
    ``'seurat'``, and vice versa. Running both and comparing the resulting
    integrated embeddings therefore compares two different feature spaces, not two
    implementations of one method.

    ``'xgboost'`` is an experiment, not a Seurat equivalent. Its importances are
    in-sample (no held-out fold, no cross-validation), and anchors found on
    batch-discriminating features can match cells by batch structure rather than
    by cell identity. ``'seurat'`` remains the supported path for reciprocal PCA.
    Use ``return_importance=True`` to inspect the ranking it produces.

    Parameters
    ----------
    return_importance
        Only meaningful for ``method='xgboost'``. Return ``(genes, table)``, where
        ``table`` is a ``DataFrame`` of ``gene`` / ``importance`` sorted
        descending. For ``method='seurat'`` the flag is ignored and a plain list
        is returned, because that method has no single importance score.

    Returns
    -------
    list of str, or (list of str, DataFrame)
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
            return_importance=return_importance,
        )

    if return_importance and method == "seurat":
        logger.warning(
            "return_importance is only produced by method='xgboost'; method=%r has "
            "no single importance score, returning a plain list.",
            method,
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
            vf_lists.append([resolve_gene_alias(g, adata.var_names) or g for g in vf])
            gene_attr = entry.get("gene_attr", {})
            if gene_attr:
                gene_sets.append(
                    {
                        canon
                        for g in gene_attr.keys()
                        if (canon := resolve_gene_alias(g, adata.var_names)) is not None
                    }
                )
        gene_universe = (
            set.intersection(*gene_sets) if gene_sets else set(adata.var_names)
        )
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
