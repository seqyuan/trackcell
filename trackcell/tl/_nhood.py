"""
Spatial topology and neighborhood-enrichment analysis for YardCluster.

Adapted from the CellCharter design (``cc.gr.connected_components``,
``cc.gr.nhood_enrichment``, ``cc.gr.diff_nhood_enrichment``), reimplemented on
top of scipy/numpy so no extra dependencies are required.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sps
from anndata import AnnData

_DEFAULT_CONN_KEY = "spatial_connectivities"


def _bh_adjust(pvalue: pd.DataFrame) -> pd.DataFrame:
    """Benjamini-Hochberg FDR over the flattened (cluster x cluster) p-values.

    A permutation matrix is dozens-to-hundreds of simultaneous tests; without
    correction every cell of the heatmap invites over-reading. Computed over
    finite entries only, so the shape and NaN layout of the input survive.
    """
    values = pvalue.to_numpy(dtype=np.float64, copy=True)
    finite = np.isfinite(values)
    if not finite.any():
        return pd.DataFrame(values, index=pvalue.index, columns=pvalue.columns)

    flat = values[finite]
    order = np.argsort(flat)
    n = flat.size
    ranked = flat[order] * n / np.arange(1, n + 1)
    # Step-up monotonicity: adjusted p is the running minimum from the largest.
    ranked = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    adjusted = np.empty_like(flat)
    adjusted[order] = ranked
    values[finite] = adjusted
    return pd.DataFrame(values, index=pvalue.index, columns=pvalue.columns)


def _as_binary(adj: sps.spmatrix) -> sps.csr_matrix:
    adj = adj.tocsr().astype(bool).astype(np.float64)
    return adj


def _component_labels(adj: sps.csr_matrix, min_cells: int) -> tuple[np.ndarray, int]:
    """Connected-component labels per node; tiny components -> -1.

    Returns
    -------
    labels
        Component id per node (-1 for components smaller than ``min_cells``).
    n_kept
        Number of components kept.
    """
    if adj.shape[0] == 0:
        return np.array([], dtype=np.int64), 0
    n_components, labels = sps.csgraph.connected_components(
        adj, directed=False, return_labels=True
    )
    counts = np.bincount(labels, minlength=n_components)
    kept = counts >= min_cells
    out = np.full(labels.shape[0], -1, dtype=np.int64)
    next_id = 0
    for lbl in range(n_components):
        if kept[lbl]:
            out[labels == lbl] = next_id
            next_id += 1
    return out, next_id


def connected_components(
    adata: AnnData,
    cluster_key: Optional[str] = None,
    min_cells: int = 250,
    connectivity_key: str = _DEFAULT_CONN_KEY,
    out_key: str = "component",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Split spatial clusters into connected components of the spatial graph.

    A single cluster label can correspond to several physically disconnected
    regions (e.g. the same domain across multiple TMA cores). This function
    relabels each connected region as a distinct component, discarding regions
    smaller than ``min_cells`` (assigned ``NaN``).

    Parameters
    ----------
    adata
        Annotated data object.
    cluster_key
        Key in ``adata.obs`` with cluster labels. If ``None``, components are
        computed on the whole graph.
    min_cells
        Minimum number of cells for a component to be kept.
    connectivity_key
        Key in ``adata.obsp`` with the spatial connectivities.
    out_key
        Key in ``adata.obs`` where component labels are stored.
    copy
        If ``True``, return a modified copy.

    Returns
    -------
    ``None`` or :class:`anndata.AnnData` if ``copy=True``.
    """
    adata = adata.copy() if copy else adata
    if connectivity_key not in adata.obsp:
        raise ValueError(
            f"`adata.obsp['{connectivity_key}']` not found. "
            "Run `tcl.tl.spatial_neighbors` first."
        )
    adj = _as_binary(adata.obsp[connectivity_key])

    comp = pd.Series(np.nan, index=adata.obs.index, dtype=object)

    if cluster_key is not None:
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"cluster_key '{cluster_key}' not in adata.obs.")
        labels = adata.obs[cluster_key].astype(str)
        offset = 0
        for cl in sorted(labels.unique()):
            mask = (labels == cl).to_numpy()
            sub_adj = adj[mask][:, mask]
            cl_labels, n_kept = _component_labels(sub_adj, min_cells)
            vals = np.full(mask.sum(), np.nan, dtype=object)
            vals[cl_labels >= 0] = [
                str(int(x) + offset) for x in cl_labels[cl_labels >= 0]
            ]
            comp[mask] = vals
            offset += n_kept
    else:
        labels, _ = _component_labels(adj, min_cells)
        vals = np.full(len(labels), np.nan, dtype=object)
        vals[labels >= 0] = [str(int(x)) for x in labels[labels >= 0]]
        comp.loc[:] = vals

    adata.obs[out_key] = comp.astype("category")
    adata.uns[f"{out_key}_params"] = {
        "cluster_key": cluster_key,
        "min_cells": min_cells,
        "connectivity_key": connectivity_key,
    }
    return adata if copy else None


def _without_self_loops(adj: sps.csr_matrix) -> sps.csr_matrix:
    """Drop the diagonal of a spatial graph.

    Self loops in ``obsp`` are not spatial links. ``build_spatial_weights``
    writes ``(i, i, 1.0)`` for a cell with no usable neighbour, and
    ``remove_long_links`` does the same for cells whose edges were all trimmed,
    so that such a cell's neighbourhood context is its own expression rather
    than zero. Counting those entries as observed links inflates the diagonal of
    the observed matrix while the configuration-model null only halves it (one
    self loop is counted once, not twice), which can flip the sign of the
    same-cluster enrichment. Removing the diagonal keeps both sides comparable.
    """
    out = adj.tocsr().copy()
    out.setdiag(0.0)
    out.eliminate_zeros()
    return out


def _observed_n_links(
    adj: sps.csr_matrix,
    labels: pd.Series,
    symmetric: bool,
) -> pd.DataFrame:
    cats = list(labels.cat.categories)
    n = len(cats)
    obs = np.zeros((n, n))
    for i, c1 in enumerate(cats):
        row_mask = (labels == c1).to_numpy()
        sub = adj[row_mask, :]
        for j, c2 in enumerate(cats):
            col_mask = (labels == c2).to_numpy()
            total = sub[:, col_mask].sum()
            if not symmetric and i != j:
                denom = sub.sum()
                obs[i, j] = total / denom if denom > 0 else 0.0
            else:
                obs[i, j] = total
    return pd.DataFrame(obs, index=cats, columns=cats)


def _expected_n_links(
    adj: sps.csr_matrix,
    labels: pd.Series,
    symmetric: bool,
) -> pd.DataFrame:
    cats = list(labels.cat.categories)
    n = len(cats)
    degrees = np.asarray(adj.sum(axis=1)).ravel()
    total = degrees.sum()
    if total <= 0:
        return pd.DataFrame(np.zeros((n, n)), index=cats, columns=cats)
    cluster_degree = np.array([degrees[(labels == c).to_numpy()].sum() for c in cats])
    # Configuration model: E[i,j] = k_i * k_j / total (undirected).
    exp = np.outer(cluster_degree, cluster_degree) / total
    np.fill_diagonal(exp, np.diag(exp) / 2.0)
    if not symmetric:
        # Normalize each row by the source cluster degree.
        row_denom = cluster_degree
        exp = exp / np.maximum(row_denom[:, None], 1e-12)
    return pd.DataFrame(exp, index=cats, columns=cats)


def nhood_enrichment(
    adata: AnnData,
    cluster_key: str,
    connectivity_key: str = _DEFAULT_CONN_KEY,
    log_fold_change: bool = False,
    only_inter: bool = True,
    symmetric: bool = True,
    copy: bool = False,
) -> Optional[dict]:
    """
    Neighborhood enrichment between clusters on the spatial graph.

    Compares the observed weighted links between cluster pairs with the expected
    links under a degree-preserving (configuration-model) null. Enrichment is
    ``observed - expected`` by default, or ``log2(observed / expected)`` when
    ``log_fold_change=True``.

    Parameters
    ----------
    adata
        Annotated data object.
    cluster_key
        Key in ``adata.obs`` with cluster labels.
    connectivity_key
        Key in ``adata.obsp`` with the spatial connectivities.
    log_fold_change
        If ``True``, report ``log2(observed / expected)`` instead of the difference.
    only_inter
        If ``True``, ignore the diagonal (same-cluster links).
    symmetric
        If ``True``, the enrichment matrix is symmetric; otherwise rows are
        normalized by the source cluster degree (directed view).
    copy
        If ``True``, return the result dict instead of writing to ``adata.uns``.

    Returns
    -------
    ``None`` or a dict with keys ``'enrichment'``, ``'observed'``, ``'expected'``.
    """
    adata = adata.copy() if copy else adata
    if connectivity_key not in adata.obsp:
        raise ValueError(
            f"`adata.obsp['{connectivity_key}']` not found. "
            "Run `tcl.tl.spatial_neighbors` first."
        )
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"cluster_key '{cluster_key}' not in adata.obs.")

    adj = _without_self_loops(adata.obsp[connectivity_key].tocsr().astype(np.float64))
    labels = adata.obs[cluster_key].astype("category")

    observed = _observed_n_links(adj, labels, symmetric=symmetric)
    expected = _expected_n_links(adj, labels, symmetric=symmetric)

    if log_fold_change:
        with np.errstate(divide="ignore", invalid="ignore"):
            enrichment = np.log2(observed.values / expected.values)
        enrichment = pd.DataFrame(
            enrichment, index=observed.index, columns=observed.columns
        )
        enrichment = enrichment.replace([np.inf, -np.inf], np.nan)
    else:
        enrichment = observed - expected

    if only_inter:
        np.fill_diagonal(observed.values, np.nan)
        np.fill_diagonal(expected.values, np.nan)
        np.fill_diagonal(enrichment.values, np.nan)

    result = {"enrichment": enrichment, "observed": observed, "expected": expected}
    if copy:
        return result
    adata.uns[f"{cluster_key}_nhood_enrichment"] = result
    adata.uns[f"{cluster_key}_nhood_enrichment"]["params"] = {
        "connectivity_key": connectivity_key,
        "log_fold_change": log_fold_change,
        "only_inter": only_inter,
        "symmetric": symmetric,
    }
    return None


def _diff_enrichment(
    adj: sps.csr_matrix,
    labels: pd.Series,
    condition_groups: tuple[str, str],
    conditions: pd.Series,
    log_fold_change: bool,
) -> pd.DataFrame:
    enrichments = {}
    for cond in condition_groups:
        mask = np.asarray(conditions == cond)
        sub_adj = adj[mask][:, mask]
        sub_labels = labels[mask].cat.remove_unused_categories()
        observed = _observed_n_links(sub_adj, sub_labels, symmetric=True)
        expected = _expected_n_links(sub_adj, sub_labels, symmetric=True)
        if log_fold_change:
            with np.errstate(divide="ignore", invalid="ignore"):
                enr = np.log2(observed.values / expected.values)
            enr = pd.DataFrame(
                enr, index=observed.index, columns=observed.columns
            ).replace([np.inf, -np.inf], np.nan)
        else:
            enr = observed - expected
        # Align to the union of categories.
        enrichments[cond] = enr.reindex(
            index=labels.cat.categories, columns=labels.cat.categories
        )
    return enrichments[condition_groups[0]] - enrichments[condition_groups[1]]


def diff_nhood_enrichment(
    adata: AnnData,
    cluster_key: str,
    condition_key: str,
    condition_groups: tuple[str, str],
    connectivity_key: str = _DEFAULT_CONN_KEY,
    library_key: Optional[str] = None,
    pvalues: bool = False,
    n_perms: int = 1000,
    random_state: int = 0,
    log_fold_change: bool = False,
    copy: bool = False,
) -> Optional[dict]:
    """
    Differential neighborhood enrichment between two conditions.

    Computes the per-condition enrichment matrices and subtracts them. When
    ``pvalues=True``, a permutation test shuffles the sample (``library_key``)
    labels between conditions to build the null distribution of each
    cluster-pair difference.

    Parameters
    ----------
    adata
        Annotated data object.
    cluster_key
        Key in ``adata.obs`` with cluster labels.
    condition_key
        Key in ``adata.obs`` with the condition label.
    condition_groups
        Two condition values to compare, e.g. ``("Case", "Ctr")``.
    connectivity_key
        Key in ``adata.obsp`` with the spatial connectivities.
    library_key
        Key in ``adata.obs`` with the sample/library id used for permutation.
        Required when ``pvalues=True``.
    pvalues
        If ``True``, compute empirical p-values via sample-label permutation.
    n_perms
        Number of permutations.
    random_state
        Seed for the permutation RNG.
    log_fold_change
        If ``True``, use log2 ratio enrichment before differencing.
    copy
        If ``True``, return the result dict.

    Returns
    -------
    ``None`` or dict with key ``'{c1}_{c2}'`` containing ``'enrichment'`` and
    optionally ``'pvalue'``.
    """
    adata = adata.copy() if copy else adata
    if connectivity_key not in adata.obsp:
        raise ValueError(
            f"`adata.obsp['{connectivity_key}']` not found. "
            "Run `tcl.tl.spatial_neighbors` first."
        )
    if cluster_key not in adata.obs.columns or condition_key not in adata.obs.columns:
        raise ValueError("cluster_key and condition_key must both be in adata.obs.")
    if pvalues and library_key is None:
        raise ValueError("library_key is required when pvalues=True.")
    if pvalues and library_key not in adata.obs.columns:
        raise ValueError(f"library_key '{library_key}' not in adata.obs.")

    conditions_all = adata.obs[condition_key].astype(str)
    missing_groups = [g for g in condition_groups if g not in set(conditions_all)]
    if missing_groups:
        raise ValueError(
            f"condition_groups not found in adata.obs['{condition_key}']: "
            f"{missing_groups}. Available: {sorted(set(conditions_all))}"
        )

    # Restrict both the observed statistic and the permutation null to the two
    # compared conditions.  Otherwise cells from other conditions silently
    # join the "0" group during permutations and bias the p-values.
    keep = conditions_all.isin(list(condition_groups)).to_numpy()
    adj_full = adata.obsp[connectivity_key].tocsr().astype(np.float64)
    adj = adj_full[keep][:, keep]
    labels = (
        adata.obs.loc[keep, cluster_key]
        .astype("category")
        .cat.remove_unused_categories()
    )
    conditions = pd.Categorical(
        conditions_all.to_numpy()[keep], categories=list(condition_groups)
    )
    libraries = (
        adata.obs.loc[keep, library_key].astype(str)
        if library_key is not None
        else None
    )

    # Hoisted: ``adj`` never changes across permutations, and copying it inside
    # the loop cost ~6.6 s per 200 permutations on a 20k-cell graph.
    adj_no_loops = _without_self_loops(adj)

    obs_diff = _diff_enrichment(
        adj_no_loops,
        labels,
        condition_groups,
        conditions,
        log_fold_change=log_fold_change,
    )

    result_key = f"{condition_groups[0]}_{condition_groups[1]}"
    result: dict = {result_key: {"enrichment": obs_diff}}

    if pvalues:
        rng = np.random.default_rng(random_state)
        libs1 = libraries[conditions == condition_groups[0]].unique()
        libs2 = libraries[conditions == condition_groups[1]].unique()
        all_libs = np.concatenate([libs1, libs2])
        n1 = len(libs1)

        # Two-sided permutation test: p = P(|perm_diff| >= |obs_diff|).
        #
        # ``perm_diff`` / ``obs_diff`` are reindexed onto the full category list,
        # so a category absent from one condition yields NaN. Comparing NaN
        # yields False and never increments the counter, which biased every
        # such entry towards a spuriously small p-value. Missing categories have
        # no links on either side, so their difference is genuinely 0.
        obs_abs = np.nan_to_num(np.abs(obs_diff.values), nan=0.0, posinf=np.inf)
        counts = np.zeros_like(obs_diff.values, dtype=np.float64)
        for _ in range(n_perms):
            perm = rng.permutation(all_libs)
            group1 = set(perm[:n1].tolist())
            permuted_cond = pd.Categorical(
                libraries.isin(group1).astype(int).astype(str)
            )
            perm_diff = _diff_enrichment(
                adj_no_loops,
                labels,
                ("1", "0"),
                permuted_cond,
                log_fold_change=log_fold_change,
            )
            perm_abs = np.nan_to_num(np.abs(perm_diff.values), nan=0.0, posinf=np.inf)
            counts += (perm_abs >= obs_abs).astype(np.float64)

        # +1 correction: with a finite number of permutations, p == 0 is not
        # attainable; (0+1)/(n+1) is the standard estimator and keeps the value
        # usable by multiple-testing correction.
        pvalue = (counts + 1.0) / (n_perms + 1.0)
        pvalue = pd.DataFrame(pvalue, index=obs_diff.index, columns=obs_diff.columns)
        result[result_key]["pvalue"] = pvalue
        result[result_key]["pvalue_adj"] = _bh_adjust(pvalue)

    if copy:
        return result
    adata.uns[f"{cluster_key}_{condition_key}_diff_nhood_enrichment"] = result
    adata.uns[f"{cluster_key}_{condition_key}_diff_nhood_enrichment"]["params"] = {
        "connectivity_key": connectivity_key,
        "condition_groups": condition_groups,
        "library_key": library_key,
        "pvalues": pvalues,
        "n_perms": n_perms,
        "log_fold_change": log_fold_change,
    }
    return None
