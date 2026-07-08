"""
YardCluster: lightweight spatial clustering for trackcell.
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from typing import Literal, Optional, Union

import numpy as np
import scanpy as sc
import scipy.sparse as sps
from anndata import AnnData

from ._cluster_merge import merge_clusters_de
from ._spatial_graph import (
    build_spatial_weights,
    compute_neighborhood_gradient,
    compute_neighborhood_mean,
    expr_k_neighbors,
    mix_embeddings,
    spatial_distances_from_weights,
)

ModeType = Literal["celltype", "domain", "auto"]
IntegrateType = Literal["none", "separate", "joint"]
GradientMode = Literal["cosphi", "regression"]

_SKETCH_THRESHOLD = 500_000
_SKETCH_N = 500_000


def _get_coords(adata: AnnData) -> np.ndarray:
    if "spatial" not in adata.obsm:
        raise ValueError("`adata.obsm['spatial']` is required.")
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("`adata.obsm['spatial']` must be (n_obs, >=2).")
    return coords[:, :2]


def _get_expression(
    adata: AnnData,
    layer: Optional[str] = None,
    use_raw: bool = False,
) -> Union[np.ndarray, sps.spmatrix]:
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        return adata.raw.X
    if layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"layer '{layer}' not in adata.layers.")
        return adata.layers[layer]
    return adata.X


def _batch_labels(adata: AnnData, batch_key: Optional[str]) -> Optional[np.ndarray]:
    if batch_key is None:
        return None
    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs.")
    return adata.obs[batch_key].to_numpy()


def _store_params(adata: AnnData, key: str, params: dict) -> None:
    adata.uns[key] = deepcopy(params)


def _run_harmony(
    adata: AnnData,
    use_rep: str,
    batch_key: str,
    key_out: Optional[str] = None,
) -> str:
    """Harmony batch correction on an obsm embedding (optional dependency)."""
    key_out = key_out or f"{use_rep}_harmony"
    try:
        import scanpy.external as sce
    except ImportError as exc:
        raise ImportError(
            "Harmony integration requires harmonypy. "
            "Install with: pip install harmonypy"
        ) from exc

    adata.obsm[use_rep] = np.asarray(adata.obsm[use_rep], dtype=np.float32)
    sce.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis=use_rep,
        adjusted_basis=key_out,
    )
    return key_out


def spatial_neighbors(
    adata: AnnData,
    k: int = 15,
    batch_key: Optional[str] = None,
    key_added: str = "spatial",
    copy: bool = False,
) -> Optional[AnnData]:
    """Build within-batch spatial kNN graph with Gaussian weights."""
    adata = adata.copy() if copy else adata
    coords = _get_coords(adata)
    batch = _batch_labels(adata, batch_key)
    weights = build_spatial_weights(coords, k=k, batch_labels=batch)
    adata.obsp[f"{key_added}_connectivities"] = weights
    adata.obsp[f"{key_added}_distances"] = spatial_distances_from_weights(weights)
    adata.uns[f"{key_added}_neighbors_params"] = {"k": k, "batch_key": batch_key}
    return adata if copy else None


def neighborhood_features(
    adata: AnnData,
    k_spatial: int = 15,
    batch_key: Optional[str] = None,
    use_gradient: bool = False,
    gradient_mode: GradientMode = "cosphi",
    k_gradient: Optional[int] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    chunk_size: int = 500,
    key_added: str = "yard",
    copy: bool = False,
) -> Optional[AnnData]:
    """Compute neighborhood mean M and optional gradient G features."""
    adata = adata.copy() if copy else adata
    spatial_neighbors(adata, k=k_spatial, batch_key=batch_key, key_added="spatial")

    X = _get_expression(adata, layer=layer, use_raw=use_raw)
    weights = adata.obsp["spatial_connectivities"]
    M = compute_neighborhood_mean(X, weights, chunk_size=chunk_size)
    adata.obsm[f"{key_added}_M"] = M.astype(np.float32)

    if use_gradient:
        if k_gradient is None:
            k_gradient = 2 * k_spatial
        coords = _get_coords(adata)
        G = compute_neighborhood_gradient(
            X,
            coords,
            weights,
            mode=gradient_mode,
            k_gradient=k_gradient,
            chunk_size=chunk_size,
        )
        adata.obsm[f"{key_added}_G"] = G.astype(np.float32)

    adata.uns[f"{key_added}_features_params"] = {
        "k_spatial": k_spatial,
        "k_gradient": k_gradient if use_gradient else None,
        "use_gradient": use_gradient,
        "gradient_mode": gradient_mode if use_gradient else None,
        "batch_key": batch_key,
    }
    return adata if copy else None


def _pca_embed(adata: AnnData, matrix: np.ndarray, n_pcs: int, key: str) -> np.ndarray:
    tmp = AnnData(X=matrix.astype(np.float32))
    n_comps = max(1, min(n_pcs, matrix.shape[0] - 1, matrix.shape[1]))
    sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack")
    adata.uns[f"{key}_pca_variance_ratio"] = tmp.uns["pca"]["variance_ratio"].tolist()
    return tmp.obsm["X_pca"].astype(np.float32)


def yard_embed(
    adata: AnnData,
    lam: float = 0.2,
    n_pcs: int = 20,
    key_added: str = "yard",
    layer: Optional[str] = None,
    use_raw: bool = False,
    batch_key: Optional[str] = None,
    harmony_integrate: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """Dual-channel YardCluster embeddings and lambda-mixed representation."""
    adata = adata.copy() if copy else adata
    m_key = f"{key_added}_M"
    if m_key not in adata.obsm:
        raise ValueError(f"Run `neighborhood_features` first (missing `{m_key}`).")

    X = _get_expression(adata, layer=layer, use_raw=use_raw)
    if sps.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    M = np.asarray(adata.obsm[m_key], dtype=np.float32)
    g_key = f"{key_added}_G"
    context_matrix = np.hstack([M, adata.obsm[g_key]]) if g_key in adata.obsm else M

    id_emb = _pca_embed(adata, X, n_pcs, f"{key_added}_identity")
    ctx_emb = _pca_embed(adata, context_matrix, n_pcs, f"{key_added}_context")

    id_key = f"X_{key_added}_identity"
    ctx_key = f"X_{key_added}_context"
    adata.obsm[id_key] = id_emb
    adata.obsm[ctx_key] = ctx_emb

    if harmony_integrate and batch_key is not None:
        id_key = _run_harmony(adata, id_key, batch_key)
        ctx_key = _run_harmony(adata, ctx_key, batch_key)

    mixed = mix_embeddings(adata.obsm[id_key], adata.obsm[ctx_key], lam)
    adata.obsm[f"X_{key_added}_mixed"] = mixed
    adata.uns[f"{key_added}_embed_params"] = {
        "lam": lam,
        "n_pcs": n_pcs,
        "harmony_integrate": harmony_integrate,
        "identity_rep": id_key,
        "context_rep": ctx_key,
    }
    return adata if copy else None


def _leiden(
    adata: AnnData,
    resolution: float,
    key_added: str,
    neighbors_key: Optional[str],
) -> None:
    leiden_kw = dict(resolution=resolution, key_added=key_added, neighbors_key=neighbors_key)
    try:
        sc.tl.leiden(
            adata,
            flavor="igraph",
            directed=False,
            n_iterations=2,
            **leiden_kw,
        )
    except (ImportError, TypeError, ValueError):
        sc.tl.leiden(adata, **leiden_kw)


def _assign_by_nearest_medoid(
    adata_full: AnnData,
    adata_sketch: AnnData,
    cluster_key: str,
    use_rep: str,
) -> None:
    """Assign non-sketch cells to nearest sketch cluster medoid in embedding space."""
    labels = adata_sketch.obs[cluster_key].astype(str)
    emb_sk = np.asarray(adata_sketch.obsm[use_rep])
    emb_full = np.asarray(adata_full.obsm[use_rep])
    clusters = sorted(labels.unique(), key=lambda x: (len(x), x))
    medoids = []
    for cl in clusters:
        idx = np.where(labels.values == cl)[0]
        sub = emb_sk[idx]
        if len(idx) == 1:
            medoids.append(sub[0])
        else:
            dists = np.linalg.norm(sub[:, None, :] - sub[None, :, :], axis=2)
            medoids.append(sub[np.argmin(dists.sum(axis=1))])
    medoids = np.vstack(medoids)
    d = np.linalg.norm(emb_full[:, None, :] - medoids[None, :, :], axis=2)
    assigned = np.array([clusters[i] for i in np.argmin(d, axis=1)])
    adata_full.obs[cluster_key] = assigned


def yard_cluster(
    adata: AnnData,
    use_rep: str = "X_yard_mixed",
    k_expr: Optional[int] = 50,
    resolution: float = 1.0,
    key_added: str = "yardcluster",
    neighbors_key: Optional[str] = None,
    merge_clusters: bool = False,
    merge_adj_p: float = 0.05,
    merge_max_de_tests: Optional[int] = 2000,
    sketch: bool = False,
    sketch_threshold: int = _SKETCH_THRESHOLD,
    sketch_n: int = _SKETCH_N,
    copy: bool = False,
) -> Optional[AnnData]:
    """Cluster using YardCluster embedding and Leiden; optional DE merge and sketch."""
    adata = adata.copy() if copy else adata
    if use_rep not in adata.obsm:
        raise ValueError(f"`adata.obsm['{use_rep}']` not found.")

    n_neighbors = expr_k_neighbors(adata.n_obs, k_expr)
    do_sketch = sketch or (adata.n_obs > sketch_threshold)

    if do_sketch and adata.n_obs > sketch_n:
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, size=sketch_n, replace=False)
        sketch_mask = np.zeros(adata.n_obs, dtype=bool)
        sketch_mask[idx] = True
        adata_sk = adata[sketch_mask].copy()
        sc.pp.neighbors(adata_sk, n_neighbors=n_neighbors, use_rep=use_rep, key_added=neighbors_key)
        _leiden(adata_sk, resolution, key_added, neighbors_key)
        if merge_clusters:
            merge_clusters_de(
                adata_sk,
                key_added,
                use_rep,
                adj_p_threshold=merge_adj_p,
                max_de_tests=merge_max_de_tests,
            )
            key_added = f"{key_added}_merged"
        _assign_by_nearest_medoid(adata, adata_sk, key_added, use_rep)
        adata.uns[f"{key_added}_cluster_params"] = {
            "use_rep": use_rep,
            "sketch": True,
            "sketch_n": sketch_n,
            "n_obs": int(adata.n_obs),
        }
        return adata if copy else None

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, key_added=neighbors_key)
    _leiden(adata, resolution, key_added, neighbors_key)
    if merge_clusters:
        merged_key = merge_clusters_de(
            adata,
            key_added,
            use_rep,
            adj_p_threshold=merge_adj_p,
            max_de_tests=merge_max_de_tests,
        )
        adata.obs[key_added] = adata.obs[merged_key].astype(str).values
    adata.uns[f"{key_added}_cluster_params"] = {
        "use_rep": use_rep,
        "k_expr": n_neighbors,
        "resolution": resolution,
        "merge_clusters": merge_clusters,
        "merge_max_de_tests": merge_max_de_tests,
        "sketch": False,
    }
    return adata if copy else None


def _run_single_yardcluster(
    adata: AnnData,
    *,
    lam: float,
    mode_label: str,
    k_spatial: int,
    k_gradient: Optional[int],
    batch_key: Optional[str],
    use_gradient: bool,
    gradient_mode: GradientMode,
    layer: Optional[str],
    use_raw: bool,
    chunk_size: int,
    n_pcs: int,
    k_expr: Optional[int],
    resolution: float,
    key_added: str,
    neighbors_key: str,
    merge_clusters: bool,
    merge_adj_p: float,
    merge_max_de_tests: Optional[int],
    harmony_integrate: bool,
    sketch: bool,
    sketch_threshold: int,
    sketch_n: int,
) -> None:
    n_key = f"{key_added}_features"
    neighborhood_features(
        adata,
        k_spatial=k_spatial,
        batch_key=batch_key,
        use_gradient=use_gradient,
        gradient_mode=gradient_mode,
        k_gradient=k_gradient,
        layer=layer,
        use_raw=use_raw,
        chunk_size=chunk_size,
        key_added=n_key,
    )
    yard_embed(
        adata,
        lam=lam,
        n_pcs=n_pcs,
        key_added=n_key,
        layer=layer,
        use_raw=use_raw,
        batch_key=batch_key,
        harmony_integrate=harmony_integrate,
    )
    use_rep = f"X_{n_key}_mixed"
    obs_key = key_added if mode_label == "default" else f"{key_added}_{mode_label}"
    yard_cluster(
        adata,
        use_rep=use_rep,
        k_expr=k_expr,
        resolution=resolution,
        key_added=obs_key,
        neighbors_key=neighbors_key,
        merge_clusters=merge_clusters,
        merge_adj_p=merge_adj_p,
        merge_max_de_tests=merge_max_de_tests,
        sketch=sketch,
        sketch_threshold=sketch_threshold,
        sketch_n=sketch_n,
    )


def spatial_cluster(
    adata: AnnData,
    mode: ModeType = "auto",
    lambda_celltype: float = 0.2,
    lambda_domain: float = 0.8,
    k_spatial: int = 15,
    k_gradient: Optional[int] = None,
    k_expr: Optional[int] = 50,
    use_gradient: bool = False,
    gradient_mode: GradientMode = "cosphi",
    n_pcs: int = 20,
    resolution: float = 1.0,
    resolution_celltype: Optional[float] = None,
    resolution_domain: Optional[float] = None,
    batch_key: Optional[str] = None,
    integrate: IntegrateType = "none",
    harmony_integrate: bool = False,
    hvg_by_batch: bool = True,
    scale_by_batch: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    chunk_size: int = 500,
    preprocess: bool = True,
    n_top_genes: int = 2000,
    merge_clusters: bool = False,
    merge_adj_p: float = 0.05,
    merge_max_de_tests: Optional[int] = 2000,
    sketch: bool = False,
    sketch_threshold: int = _SKETCH_THRESHOLD,
    sketch_n: int = _SKETCH_N,
    key_added: str = "yardcluster",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    YardCluster spatial clustering (one-shot or dual-mode).

    Phase 2 options: ``merge_clusters``, ``integrate='joint'`` + ``harmony_integrate``,
    ``sketch`` for large datasets, ``gradient_mode='regression'``.
    """
    adata = adata.copy() if copy else adata

    if integrate == "joint" and batch_key is None:
        raise ValueError("integrate='joint' requires batch_key.")
    use_harmony = harmony_integrate or integrate == "joint"

    if preprocess:
        _preprocess_for_yardcluster(
            adata,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            hvg_by_batch=hvg_by_batch,
            scale_by_batch=scale_by_batch,
            layer=layer,
        )

    if integrate == "separate":
        if batch_key is None:
            raise ValueError("integrate='separate' requires batch_key.")
        batches = adata.obs[batch_key].astype(str)
        for b in batches.unique():
            sub = adata[batches == b].copy()
            spatial_cluster(
                sub,
                mode=mode,
                lambda_celltype=lambda_celltype,
                lambda_domain=lambda_domain,
                k_spatial=k_spatial,
                k_gradient=k_gradient,
                k_expr=k_expr,
                use_gradient=use_gradient,
                gradient_mode=gradient_mode,
                n_pcs=n_pcs,
                resolution=resolution,
                resolution_celltype=resolution_celltype,
                resolution_domain=resolution_domain,
                batch_key=None,
                integrate="none",
                preprocess=False,
                merge_clusters=merge_clusters,
                merge_adj_p=merge_adj_p,
                merge_max_de_tests=merge_max_de_tests,
                sketch=sketch,
                sketch_threshold=sketch_threshold,
                sketch_n=sketch_n,
                key_added=key_added,
            )
            for col in sub.obs.columns:
                if col.startswith(key_added):
                    adata.obs.loc[sub.obs_names, col] = [
                        f"{b}_{label}" for label in sub.obs[col].astype(str).values
                    ]
        _store_params(adata, f"{key_added}_params", {"integrate": "separate", "batch_key": batch_key})
        return adata if copy else None

    res_ct = resolution if resolution_celltype is None else resolution_celltype
    res_dm = resolution if resolution_domain is None else resolution_domain
    run_kw = dict(
        k_spatial=k_spatial,
        k_gradient=k_gradient,
        batch_key=batch_key,
        use_gradient=use_gradient,
        gradient_mode=gradient_mode,
        layer=layer,
        use_raw=use_raw,
        chunk_size=chunk_size,
        n_pcs=n_pcs,
        k_expr=k_expr,
        key_added=key_added,
        merge_clusters=merge_clusters,
        merge_adj_p=merge_adj_p,
        merge_max_de_tests=merge_max_de_tests,
        harmony_integrate=use_harmony,
        sketch=sketch,
        sketch_threshold=sketch_threshold,
        sketch_n=sketch_n,
    )

    if mode in ("celltype", "auto"):
        _run_single_yardcluster(
            adata,
            lam=lambda_celltype,
            mode_label="celltype" if mode == "auto" else "default",
            resolution=res_ct,
            neighbors_key=f"{key_added}_celltype" if mode == "auto" else key_added,
            **run_kw,
        )
    if mode in ("domain", "auto"):
        _run_single_yardcluster(
            adata,
            lam=lambda_domain,
            mode_label="domain" if mode == "auto" else "default",
            resolution=res_dm,
            neighbors_key=f"{key_added}_domain" if mode == "auto" else key_added,
            **run_kw,
        )

    _store_params(
        adata,
        f"{key_added}_params",
        {
            "mode": mode,
            "integrate": integrate,
            "harmony_integrate": use_harmony,
            "merge_clusters": merge_clusters,
            "merge_max_de_tests": merge_max_de_tests,
            "gradient_mode": gradient_mode,
            "sketch": sketch,
        },
    )
    return adata if copy else None


def _scale_for_yardcluster(adata: AnnData, max_value: float = 10) -> None:
    if sps.issparse(adata.X):
        warnings.warn(
            "Sparse adata.X detected during YardCluster preprocessing; using "
            "zero_center=False in sc.pp.scale() to avoid densifying the matrix.",
            stacklevel=3,
        )
        sc.pp.scale(adata, max_value=max_value, zero_center=False)
    else:
        sc.pp.scale(adata, max_value=max_value)


def _preprocess_for_yardcluster(
    adata: AnnData,
    n_top_genes: int,
    batch_key: Optional[str],
    hvg_by_batch: bool,
    scale_by_batch: bool,
    layer: Optional[str],
) -> None:
    if layer is not None:
        warnings.warn(
            "preprocess=True with layer=... uses adata.X for HVG/scale; "
            "set preprocess=False if data is already prepared.",
            stacklevel=3,
        )
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if adata.n_vars > n_top_genes:
        if "highly_variable" not in adata.var.columns:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor="seurat_v3",
                batch_key=batch_key if hvg_by_batch else None,
            )
        if "highly_variable" in adata.var.columns:
            adata._inplace_subset_var(adata.var["highly_variable"].values)
    if scale_by_batch and batch_key is not None:
        for b in adata.obs[batch_key].unique():
            idx = adata.obs[batch_key] == b
            sub = adata[idx].copy()
            _scale_for_yardcluster(sub, max_value=10)
            adata.X[idx] = sub.X
    else:
        _scale_for_yardcluster(adata, max_value=10)
