"""
DE-guided cluster merging (Space Ranger-inspired).
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.stats import false_discovery_control


def _cluster_medoids(adata: AnnData, cluster_key: str, use_rep: str) -> tuple[np.ndarray, list[str]]:
    labels = adata.obs[cluster_key].astype(str)
    clusters = sorted(labels.unique(), key=lambda x: (len(x), x))
    emb = np.asarray(adata.obsm[use_rep])
    medoids = []
    for cl in clusters:
        idx = np.where(labels.values == cl)[0]
        sub = emb[idx]
        if len(idx) == 1:
            medoids.append(sub[0])
            continue
        dists = pdist(sub)
        triu = np.triu_indices(len(idx), k=1)
        sum_d = np.zeros(len(idx))
        for k, (i, j) in enumerate(zip(triu[0], triu[1])):
            sum_d[i] += dists[k]
            sum_d[j] += dists[k]
        medoids.append(sub[np.argmin(sum_d)])
    return np.vstack(medoids), clusters


def _has_de_genes(
    adata: AnnData,
    cluster_key: str,
    cl_a: str,
    cl_b: str,
    adj_p_threshold: float,
) -> bool:
    mask = adata.obs[cluster_key].astype(str).isin([cl_a, cl_b])
    sub = adata[mask].copy()
    sub.obs["_pair"] = sub.obs[cluster_key].astype(str)
    sc.tl.rank_genes_groups(
        sub, groupby="_pair", groups=[cl_a], reference=cl_b, method="wilcoxon"
    )
    result = sc.get.rank_genes_groups_df(sub, group=cl_a)
    if result.empty:
        return False
    if "pvals_adj" in result.columns:
        adj = result["pvals_adj"].fillna(1.0).values
    else:
        pvals = result.get("pvals", pd.Series(dtype=float)).fillna(1.0).values
        adj = false_discovery_control(pvals) if len(pvals) else np.array([])
    return bool(len(adj) and np.any(adj < adj_p_threshold))


def _leaves_at(Z, clusters: list[str], idx: int, n_leaves: int) -> list[str]:
    if idx < n_leaves:
        return [clusters[idx]]
    row = Z[int(idx - n_leaves)]
    left = _leaves_at(Z, clusters, int(row[0]), n_leaves)
    right = _leaves_at(Z, clusters, int(row[1]), n_leaves)
    return left + right


def merge_clusters_de(
    adata: AnnData,
    cluster_key: str,
    use_rep: str,
    adj_p_threshold: float = 0.05,
    key_added: Optional[str] = None,
    max_de_tests: Optional[int] = 2000,
) -> str:
    """
    Merge sibling clusters with no significant DE genes (hierarchical medoids + Wilcoxon).

    ``max_de_tests`` caps the number of pairwise DE comparisons. When the
    possible cluster pairs exceed this cap, labels are left unchanged instead of
    launching an unexpectedly expensive merge pass.
    """
    out_key = key_added or f"{cluster_key}_merged"
    labels = adata.obs[cluster_key].astype(str).copy()
    unique = sorted(labels.unique(), key=lambda x: (len(x), x))
    if len(unique) <= 1:
        adata.obs[out_key] = labels.values
        return out_key

    max_possible_tests = len(unique) * (len(unique) - 1) // 2
    if max_de_tests is not None and max_possible_tests > max_de_tests:
        warnings.warn(
            f"Skipping DE-guided cluster merge: {len(unique)} clusters could require "
            f"up to {max_possible_tests} pairwise DE tests, above max_de_tests={max_de_tests}. "
            "Increase max_de_tests to run the merge anyway."
        )
        adata.obs[out_key] = labels.values
        adata.uns[f"{out_key}_merge_params"] = {
            "source": cluster_key,
            "use_rep": use_rep,
            "adj_p_threshold": adj_p_threshold,
            "n_before": len(unique),
            "n_after": len(unique),
            "skipped": True,
            "max_de_tests": max_de_tests,
            "max_possible_tests": max_possible_tests,
        }
        return out_key

    medoids, clusters = _cluster_medoids(adata, cluster_key, use_rep)
    Z = linkage(medoids, method="average")
    n = len(clusters)
    parent = {c: c for c in clusters}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    checked: set[tuple[str, str]] = set()
    for i in range(len(Z)):
        left_clusters = _leaves_at(Z, clusters, int(Z[i, 0]), n)
        right_clusters = _leaves_at(Z, clusters, int(Z[i, 1]), n)
        for ca in left_clusters:
            for cb in right_clusters:
                if ca == cb:
                    continue
                pair = tuple(sorted([ca, cb]))
                if pair in checked:
                    continue
                checked.add(pair)
                if not _has_de_genes(adata, cluster_key, ca, cb, adj_p_threshold):
                    union(ca, cb)

    merged_map = {c: find(c) for c in clusters}
    adata.obs[out_key] = labels.map(merged_map).values
    adata.uns[f"{out_key}_merge_params"] = {
        "source": cluster_key,
        "use_rep": use_rep,
        "adj_p_threshold": adj_p_threshold,
        "n_before": len(unique),
        "n_after": len(set(merged_map.values())),
        "skipped": False,
        "max_de_tests": max_de_tests,
    }
    return out_key
