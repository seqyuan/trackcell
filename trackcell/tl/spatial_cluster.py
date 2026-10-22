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

from ._autores import select_resolution
from ._cluster_merge import merge_clusters_de
from ._yardcluster_utils import (
    GraphModeType,
    PlatformType,
    mix_spatial_expr_connectivities,
    neighborhood_features_cached,
    resolve_platform_params,
    spatial_neighbors_cached,
    spatial_stratified_sketch_indices,
)
from ._spatial_graph import (
    build_spatial_weights,
    compute_neighborhood_gradient,
    compute_neighborhood_mean,
    compute_neighborhood_multi_hop,
    expr_k_neighbors,
    mix_embeddings,
    spatial_distances_from_weights,
)

ModeType = Literal["celltype", "domain", "auto"]
IntegrateType = Literal["none", "separate", "joint", "bbknn"]
GradientMode = Literal["cosphi", "regression"]
ExpressionSourceType = Literal["x", "sct_prep", "sct_integrated"]
ResolutionType = Union[float, Literal["auto"]]
MethodType = Literal["leiden", "gmm"]

_SKETCH_THRESHOLD = 500_000
_SKETCH_N = 500_000


def _context_expression_source(expression_source: ExpressionSourceType) -> ExpressionSourceType:
    """Gene-level matrix for spatial M/G; identity may use RPCA embedding."""
    if expression_source == "sct_integrated":
        return "sct_prep"
    return expression_source


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


def _get_yard_expression(
    adata: AnnData,
    *,
    expression_source: ExpressionSourceType = "x",
    layer: Optional[str] = None,
    use_raw: bool = False,
    sct_prep_key: str = "sct_prep",
    sct_integrated_key: str = "sct_integrated",
) -> Union[np.ndarray, sps.spmatrix]:
    """Return expression matrix for YardCluster (sparse when possible)."""
    if expression_source == "sct_integrated":
        obsm_key = f"X_{sct_integrated_key}"
        if obsm_key not in adata.obsm:
            raise ValueError(
                f"expression_source='sct_integrated' requires adata.obsm['{obsm_key}']. "
                "Run tcl.tl.integrate_sct_rpca first."
            )
        mat = np.asarray(adata.obsm[obsm_key], dtype=np.float32)
        if mat.shape[0] != adata.n_obs:
            raise ValueError(
                f"SCT integrated rows ({mat.shape[0]}) != n_obs ({adata.n_obs})."
            )
        return mat
    if expression_source == "sct_prep":
        obsm_key = f"X_{sct_prep_key}"
        if obsm_key not in adata.obsm:
            raise ValueError(
                f"expression_source='sct_prep' requires adata.obsm['{obsm_key}']. "
                "Run tcl.tl.run_sct_integration first."
            )
        prep_uns = adata.uns.get(sct_prep_key)
        if not prep_uns or "anchor_features" not in prep_uns:
            raise ValueError(
                f"expression_source='sct_prep' requires adata.uns['{sct_prep_key}']['anchor_features']."
            )
        mat = np.asarray(adata.obsm[obsm_key], dtype=np.float32)
        if mat.shape[0] != adata.n_obs:
            raise ValueError(
                f"SCT prep matrix rows ({mat.shape[0]}) != n_obs ({adata.n_obs})."
            )
        n_features = len(prep_uns["anchor_features"])
        if mat.shape[1] != n_features:
            raise ValueError(
                f"SCT prep matrix columns ({mat.shape[1]}) != "
                f"len(anchor_features) ({n_features})."
            )
        return mat
    return _get_expression(adata, layer=layer, use_raw=use_raw)


def _resolve_yard_expression(
    adata: AnnData,
    *,
    expression_source: ExpressionSourceType = "x",
    layer: Optional[str] = None,
    use_raw: bool = False,
    sct_prep_key: str = "sct_prep",
    sct_integrated_key: str = "sct_integrated",
) -> np.ndarray:
    """Return a dense (n_obs, n_genes) matrix when densification is required."""
    X = _get_yard_expression(
        adata,
        expression_source=expression_source,
        layer=layer,
        use_raw=use_raw,
        sct_prep_key=sct_prep_key,
        sct_integrated_key=sct_integrated_key,
    )
    if sps.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def _default_use_gradient(mode: ModeType, use_gradient: Optional[bool]) -> bool:
    if use_gradient is not None:
        return use_gradient
    return mode in ("domain", "auto")


def _resolve_merge_clusters(
    mode_label: str,
    merge_clusters: Optional[bool],
    merge_clusters_celltype: Optional[bool],
    merge_clusters_domain: Optional[bool],
) -> bool:
    """Per-run DE-merge flag with mode-aware defaults (domain merges by default)."""
    if mode_label == "celltype" and merge_clusters_celltype is not None:
        return merge_clusters_celltype
    if mode_label == "domain" and merge_clusters_domain is not None:
        return merge_clusters_domain
    if merge_clusters is not None:
        return merge_clusters
    return mode_label == "domain"


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


def _resolve_neighbors_method(
    integrate: IntegrateType,
    graph_mode: GraphModeType,
) -> GraphModeType | Literal["bbknn"]:
    if integrate == "bbknn":
        return "bbknn"
    return graph_mode


def _ensure_neighbors_uns(adata: AnnData, nk: str) -> None:
    """Ensure scanpy Leiden can resolve a custom neighbor graph key."""
    conn_key = f"{nk}_connectivities"
    dist_key = f"{nk}_distances"
    if conn_key not in adata.obsp:
        raise ValueError(f"Missing adata.obsp['{conn_key}'] for neighbors_key='{nk}'.")
    adata.uns[nk] = {
        "connectivities_key": conn_key,
        "distances_key": dist_key,
        "params": adata.uns.get(nk, {}).get("params", {}),
    }


def _build_cluster_neighbors(
    adata: AnnData,
    *,
    use_rep: str,
    n_neighbors: int,
    neighbors_key: Optional[str],
    batch_key: Optional[str],
    neighbors_method: GraphModeType | Literal["bbknn"],
    hybrid_alpha: float,
    spatial_conn_key: str = "spatial_connectivities",
) -> str:
    """Build clustering graph; returns effective ``neighbors_key`` for Leiden."""
    nk = neighbors_key or "yard_neighbors"
    conn_key = f"{nk}_connectivities"
    dist_key = f"{nk}_distances"

    if neighbors_method == "bbknn":
        if batch_key is None:
            raise ValueError("integrate='bbknn' requires batch_key.")
        try:
            import scanpy.external as sce
        except ImportError as exc:
            raise ImportError(
                "integrate='bbknn' requires bbknn. Install with: pip install bbknn"
            ) from exc
        sce.pp.bbknn(
            adata,
            batch_key=batch_key,
            use_rep=use_rep,
            neighbors_within_batch=max(3, n_neighbors // 3),
        )
        adata.obsp[conn_key] = adata.obsp["connectivities"].copy()
        adata.obsp[dist_key] = adata.obsp["distances"].copy()
        _ensure_neighbors_uns(adata, nk)
        return nk

    if neighbors_method == "spatial":
        if spatial_conn_key not in adata.obsp:
            raise ValueError(
                f"graph_mode='spatial' requires adata.obsp['{spatial_conn_key}']."
            )
        adata.obsp[conn_key] = adata.obsp[spatial_conn_key].copy()
        adata.obsp[dist_key] = spatial_distances_from_weights(adata.obsp[conn_key])
        _ensure_neighbors_uns(adata, nk)
        return nk

    expr_nk = f"{nk}_expr" if neighbors_method == "hybrid" else nk
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=use_rep,
        key_added=expr_nk,
    )
    if neighbors_method == "hybrid":
        if spatial_conn_key not in adata.obsp:
            raise ValueError(
                f"graph_mode='hybrid' requires adata.obsp['{spatial_conn_key}']."
            )
        expr_conn = adata.obsp[f"{expr_nk}_connectivities"]
        spatial_conn = adata.obsp[spatial_conn_key]
        adata.obsp[conn_key] = mix_spatial_expr_connectivities(
            spatial_conn, expr_conn, hybrid_alpha
        )
        adata.obsp[dist_key] = spatial_distances_from_weights(adata.obsp[conn_key])
        _ensure_neighbors_uns(adata, nk)
    else:
        _ensure_neighbors_uns(adata, expr_nk)
    return nk


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
    n_layers: int = 1,
    layer: Optional[str] = None,
    use_raw: bool = False,
    expression_source: ExpressionSourceType = "x",
    sct_prep_key: str = "sct_prep",
    sct_integrated_key: str = "sct_integrated",
    chunk_size: int = 500,
    key_added: str = "yard",
    reuse_neighbors: bool = True,
    copy: bool = False,
) -> Optional[AnnData]:
    """Compute neighborhood mean M and optional gradient G features.

    ``n_layers`` controls multi-hop aggregation (CellCharter-style): when >1,
    ``M`` becomes the concatenation of hop 1..``n_layers`` means instead of the
    1-hop mean.
    """
    adata = adata.copy() if copy else adata

    def _build_spatial() -> None:
        spatial_neighbors(adata, k=k_spatial, batch_key=batch_key, key_added="spatial")

    if reuse_neighbors:
        spatial_neighbors_cached(
            adata,
            k=k_spatial,
            batch_key=batch_key,
            key_added="spatial",
            build_fn=lambda ad, **kw: spatial_neighbors(ad, **kw),
        )
    else:
        _build_spatial()

    resolved_k_gradient = (
        (k_gradient if k_gradient is not None else 2 * k_spatial) if use_gradient else None
    )
    expected_params = {
        "k_spatial": k_spatial,
        "k_gradient": resolved_k_gradient,
        "use_gradient": use_gradient,
        "gradient_mode": gradient_mode if use_gradient else None,
        "n_layers": n_layers,
        "batch_key": batch_key,
        "expression_source": expression_source,
        "sct_prep_key": sct_prep_key if expression_source in ("sct_prep", "sct_integrated") else None,
        "sct_integrated_key": sct_integrated_key if expression_source == "sct_integrated" else None,
    }

    context_source = _context_expression_source(expression_source)

    def _compute_mg() -> None:
        X = _get_yard_expression(
            adata,
            expression_source=context_source,
            layer=layer,
            use_raw=use_raw,
            sct_prep_key=sct_prep_key,
            sct_integrated_key=sct_integrated_key,
        )
        weights = adata.obsp["spatial_connectivities"]
        if n_layers > 1:
            M = compute_neighborhood_multi_hop(
                X, weights, n_layers=n_layers, chunk_size=chunk_size
            )
        else:
            M = compute_neighborhood_mean(X, weights, chunk_size=chunk_size)
        adata.obsm[f"{key_added}_M"] = M.astype(np.float32)

        if use_gradient:
            coords = _get_coords(adata)
            G = compute_neighborhood_gradient(
                X,
                coords,
                weights,
                mode=gradient_mode,
                k_gradient=resolved_k_gradient,
                chunk_size=chunk_size,
            )
            adata.obsm[f"{key_added}_G"] = G.astype(np.float32)
        adata.uns[f"{key_added}_features_params"] = expected_params

    if reuse_neighbors:
        neighborhood_features_cached(
            adata,
            expected_params=expected_params,
            key_added=key_added,
            compute_fn=_compute_mg,
        )
    else:
        _compute_mg()

    return adata if copy else None


def _pca_embed(
    adata: AnnData,
    matrix: Union[np.ndarray, sps.spmatrix],
    n_pcs: int,
    key: str,
) -> np.ndarray:
    if sps.issparse(matrix):
        tmp = AnnData(X=matrix)
        n_genes = matrix.shape[1]
    else:
        matrix = np.asarray(matrix, dtype=np.float32)
        tmp = AnnData(X=matrix)
        n_genes = matrix.shape[1]
    n_comps = max(1, min(n_pcs, tmp.n_obs - 1, n_genes))
    sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack")
    adata.uns[f"{key}_pca_variance_ratio"] = tmp.uns["pca"]["variance_ratio"].tolist()
    return tmp.obsm["X_pca"].astype(np.float32)


def _yard_embed_channels(
    adata: AnnData,
    expression: Union[np.ndarray, sps.spmatrix],
    n_pcs: int,
    key_added: str,
    batch_key: Optional[str] = None,
    harmony_integrate: bool = False,
) -> tuple[str, str]:
    """Compute identity and context PCA channels (no lambda mixing)."""
    m_key = f"{key_added}_M"
    if m_key not in adata.obsm:
        raise ValueError(f"Run `neighborhood_features` first (missing `{m_key}`).")

    M = np.asarray(adata.obsm[m_key], dtype=np.float32)
    g_key = f"{key_added}_G"
    context_matrix = np.hstack([M, adata.obsm[g_key]]) if g_key in adata.obsm else M

    id_emb = _pca_embed(adata, expression, n_pcs, f"{key_added}_identity")
    ctx_emb = _pca_embed(adata, context_matrix, n_pcs, f"{key_added}_context")

    id_key = f"X_{key_added}_identity"
    ctx_key = f"X_{key_added}_context"
    adata.obsm[id_key] = id_emb
    adata.obsm[ctx_key] = ctx_emb

    if harmony_integrate and batch_key is not None:
        id_key = _run_harmony(adata, id_key, batch_key)
        ctx_key = _run_harmony(adata, ctx_key, batch_key)

    return id_key, ctx_key


def yard_embed(
    adata: AnnData,
    lam: float = 0.2,
    n_pcs: int = 20,
    key_added: str = "yard",
    layer: Optional[str] = None,
    use_raw: bool = False,
    expression_source: ExpressionSourceType = "x",
    sct_prep_key: str = "sct_prep",
    sct_integrated_key: str = "sct_integrated",
    batch_key: Optional[str] = None,
    harmony_integrate: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """Dual-channel YardCluster embeddings and lambda-mixed representation."""
    adata = adata.copy() if copy else adata
    expression = _get_yard_expression(
        adata,
        expression_source=expression_source,
        layer=layer,
        use_raw=use_raw,
        sct_prep_key=sct_prep_key,
    )
    id_key, ctx_key = _yard_embed_channels(
        adata,
        expression,
        n_pcs,
        key_added,
        batch_key=batch_key,
        harmony_integrate=harmony_integrate,
    )

    mixed = mix_embeddings(adata.obsm[id_key], adata.obsm[ctx_key], lam)
    adata.obsm[f"X_{key_added}_mixed"] = mixed
    adata.uns[f"{key_added}_embed_params"] = {
        "lam": lam,
        "n_pcs": n_pcs,
        "harmony_integrate": harmony_integrate,
        "identity_rep": id_key,
        "context_rep": ctx_key,
        "expression_source": expression_source,
        "sct_prep_key": sct_prep_key if expression_source in ("sct_prep", "sct_integrated") else None,
        "sct_integrated_key": sct_integrated_key if expression_source == "sct_integrated" else None,
    }
    return adata if copy else None


def _leiden(
    adata: AnnData,
    resolution: float,
    key_added: str,
    neighbors_key: Optional[str],
    random_state: int = 0,
) -> None:
    leiden_kw = dict(
        resolution=resolution,
        key_added=key_added,
        neighbors_key=neighbors_key,
        random_state=random_state,
    )
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


def _leiden_with_optional_auto(
    adata: AnnData,
    *,
    use_rep: str,
    resolution: Union[float, str],
    key_added: str,
    neighbors_key: str,
    n_resolution_runs: int,
    random_state: int,
) -> float:
    """Run Leiden, selecting ``resolution`` by stability when it is ``'auto'``."""
    if isinstance(resolution, str) and resolution == "auto":
        best, sweep = select_resolution(
            adata,
            use_rep=use_rep,
            neighbors_key=neighbors_key,
            n_runs=n_resolution_runs,
            random_state=random_state,
        )
        adata.uns[f"{key_added}_resolution_sweep"] = sweep
        resolution = float(best)
    _leiden(adata, resolution, key_added, neighbors_key, random_state=random_state)
    return resolution


def _gmm_cluster(
    adata: AnnData,
    use_rep: str,
    n_clusters: Optional[int],
    key_added: str,
    random_state: int,
) -> int:
    """Gaussian-mixture clustering on an embedding (sklearn, BIC for auto k)."""
    from sklearn.mixture import GaussianMixture

    X = np.asarray(adata.obsm[use_rep], dtype=np.float64)
    n_obs = X.shape[0]
    if n_clusters is None:
        ks = [k for k in range(2, min(16, n_obs) + 1)]
        if not ks:
            ks = [1]
        best_k = None
        best_bic = float("inf")
        model = None
        for k in ks:
            g = GaussianMixture(
                n_components=k,
                covariance_type="full",
                random_state=random_state,
                n_init=3,
                reg_covar=1e-6,
            )
            g.fit(X)
            bic = float(g.bic(X))
            if bic < best_bic:
                best_bic = bic
                best_k = k
                model = g
    else:
        best_k = max(1, int(n_clusters))
        model = GaussianMixture(
            n_components=best_k,
            covariance_type="full",
            random_state=random_state,
            n_init=5,
            reg_covar=1e-6,
        )
        model.fit(X)
    adata.obs[key_added] = model.predict(X).astype(str)
    return int(best_k)


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
    resolution: Union[float, Literal["auto"]] = 1.0,
    key_added: str = "yardcluster",
    neighbors_key: Optional[str] = None,
    merge_clusters: bool = False,
    merge_adj_p: float = 0.05,
    merge_max_de_tests: Optional[int] = 2000,
    sketch: bool = False,
    sketch_threshold: int = _SKETCH_THRESHOLD,
    sketch_n: int = _SKETCH_N,
    sketch_stratified: bool = True,
    batch_key: Optional[str] = None,
    neighbors_method: GraphModeType | Literal["bbknn"] = "expr",
    hybrid_alpha: float = 0.3,
    spatial_conn_key: str = "spatial_connectivities",
    method: Literal["leiden", "gmm"] = "leiden",
    n_clusters: Optional[int] = None,
    n_resolution_runs: int = 5,
    random_state: int = 0,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Cluster using the YardCluster embedding.

    ``method='leiden'`` (default) runs Leiden on a neighbor graph; set
    ``resolution='auto'`` to select the resolution by clustering stability.
    ``method='gmm'`` fits a Gaussian mixture (sklearn) on ``use_rep``, with
    ``n_clusters=None`` picking k by BIC. Optional DE merge and sketch apply to
    the Leiden path.
    """
    adata = adata.copy() if copy else adata
    if use_rep not in adata.obsm:
        raise ValueError(f"`adata.obsm['{use_rep}']` not found.")

    if method == "gmm":
        chosen_k = _gmm_cluster(adata, use_rep, n_clusters, key_added, random_state)
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
            "method": "gmm",
            "use_rep": use_rep,
            "n_clusters": chosen_k,
            "merge_clusters": merge_clusters,
        }
        return adata if copy else None

    n_neighbors = expr_k_neighbors(adata.n_obs, k_expr)
    do_sketch = sketch or (adata.n_obs > sketch_threshold)

    if do_sketch and adata.n_obs > sketch_n:
        if sketch_stratified and "spatial" in adata.obsm:
            idx = spatial_stratified_sketch_indices(
                adata.obsm["spatial"], sketch_n, seed=0
            )
        else:
            rng = np.random.default_rng(0)
            idx = rng.choice(adata.n_obs, size=sketch_n, replace=False)
        sketch_mask = np.zeros(adata.n_obs, dtype=bool)
        sketch_mask[idx] = True
        adata_sk = adata[sketch_mask].copy()
        leiden_nk = _build_cluster_neighbors(
            adata_sk,
            use_rep=use_rep,
            n_neighbors=n_neighbors,
            neighbors_key=neighbors_key,
            batch_key=batch_key,
            neighbors_method=neighbors_method,
            hybrid_alpha=hybrid_alpha,
            spatial_conn_key=spatial_conn_key,
        )
        resolution = _leiden_with_optional_auto(
            adata_sk,
            use_rep=use_rep,
            resolution=resolution,
            key_added=key_added,
            neighbors_key=leiden_nk,
            n_resolution_runs=n_resolution_runs,
            random_state=random_state,
        )
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
            "sketch_stratified": sketch_stratified,
            "sketch_n": sketch_n,
            "n_obs": int(adata.n_obs),
            "neighbors_method": neighbors_method,
        }
        return adata if copy else None

    leiden_nk = _build_cluster_neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        neighbors_key=neighbors_key,
        batch_key=batch_key,
        neighbors_method=neighbors_method,
        hybrid_alpha=hybrid_alpha,
        spatial_conn_key=spatial_conn_key,
    )
    resolution = _leiden_with_optional_auto(
        adata,
        use_rep=use_rep,
        resolution=resolution,
        key_added=key_added,
        neighbors_key=leiden_nk,
        n_resolution_runs=n_resolution_runs,
        random_state=random_state,
    )
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
        "neighbors_method": neighbors_method,
        "hybrid_alpha": hybrid_alpha if neighbors_method == "hybrid" else None,
    }
    return adata if copy else None


def _run_shared_yard_features(
    adata: AnnData,
    *,
    n_key: str,
    k_spatial: int,
    k_gradient: Optional[int],
    n_layers: int,
    batch_key: Optional[str],
    use_gradient: bool,
    gradient_mode: GradientMode,
    layer: Optional[str],
    use_raw: bool,
    expression_source: ExpressionSourceType,
    sct_prep_key: str,
    sct_integrated_key: str,
    chunk_size: int,
    n_pcs: int,
    harmony_integrate: bool,
    reuse_neighbors: bool,
) -> tuple[str, str]:
    """Compute spatial neighborhood features and dual-channel PCA once."""
    neighborhood_features(
        adata,
        k_spatial=k_spatial,
        batch_key=batch_key,
        use_gradient=use_gradient,
        gradient_mode=gradient_mode,
        k_gradient=k_gradient,
        n_layers=n_layers,
        layer=layer,
        use_raw=use_raw,
        expression_source=expression_source,
        sct_prep_key=sct_prep_key,
        sct_integrated_key=sct_integrated_key,
        chunk_size=chunk_size,
        key_added=n_key,
        reuse_neighbors=reuse_neighbors,
    )
    expression = _get_yard_expression(
        adata,
        expression_source=expression_source,
        layer=layer,
        use_raw=use_raw,
        sct_prep_key=sct_prep_key,
        sct_integrated_key=sct_integrated_key,
    )
    return _yard_embed_channels(
        adata,
        expression,
        n_pcs,
        n_key,
        batch_key=batch_key,
        harmony_integrate=harmony_integrate,
    )


def _run_yardcluster_from_channels(
    adata: AnnData,
    *,
    n_key: str,
    id_key: str,
    ctx_key: str,
    lam: float,
    mode_label: str,
    k_expr: Optional[int],
    resolution: ResolutionType,
    key_added: str,
    neighbors_key: str,
    merge_clusters: bool,
    merge_adj_p: float,
    merge_max_de_tests: Optional[int],
    sketch: bool,
    sketch_threshold: int,
    sketch_n: int,
    sketch_stratified: bool,
    batch_key: Optional[str],
    neighbors_method: GraphModeType | Literal["bbknn"],
    hybrid_alpha: float,
    method: MethodType,
    n_clusters: Optional[int],
    n_resolution_runs: int,
    random_state: int,
) -> None:
    """Mix precomputed channels with ``lam`` and run clustering."""
    mixed_key = f"X_{n_key}_mixed"
    adata.obsm[mixed_key] = mix_embeddings(adata.obsm[id_key], adata.obsm[ctx_key], lam)
    obs_key = key_added if mode_label == "default" else f"{key_added}_{mode_label}"
    yard_cluster(
        adata,
        use_rep=mixed_key,
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
        sketch_stratified=sketch_stratified,
        batch_key=batch_key,
        neighbors_method=neighbors_method,
        hybrid_alpha=hybrid_alpha,
        method=method,
        n_clusters=n_clusters,
        n_resolution_runs=n_resolution_runs,
        random_state=random_state,
    )


def _run_single_yardcluster(
    adata: AnnData,
    *,
    lam: float,
    mode_label: str,
    k_spatial: int,
    k_gradient: Optional[int],
    n_layers: int,
    batch_key: Optional[str],
    use_gradient: bool,
    gradient_mode: GradientMode,
    layer: Optional[str],
    use_raw: bool,
    expression_source: ExpressionSourceType,
    sct_prep_key: str,
    sct_integrated_key: str = "sct_integrated",
    chunk_size: int,
    n_pcs: int,
    k_expr: Optional[int],
    resolution: ResolutionType,
    key_added: str,
    neighbors_key: str,
    merge_clusters: bool,
    merge_adj_p: float,
    merge_max_de_tests: Optional[int],
    harmony_integrate: bool,
    sketch: bool,
    sketch_threshold: int,
    sketch_n: int,
    sketch_stratified: bool,
    neighbors_method: GraphModeType | Literal["bbknn"],
    hybrid_alpha: float,
    reuse_neighbors: bool,
    method: MethodType,
    n_clusters: Optional[int],
    n_resolution_runs: int,
    random_state: int,
) -> None:
    n_key = f"{key_added}_features"
    id_key, ctx_key = _run_shared_yard_features(
        adata,
        n_key=n_key,
        k_spatial=k_spatial,
        k_gradient=k_gradient,
        n_layers=n_layers,
        batch_key=batch_key,
        use_gradient=use_gradient,
        gradient_mode=gradient_mode,
        layer=layer,
        use_raw=use_raw,
        expression_source=expression_source,
        sct_prep_key=sct_prep_key,
        sct_integrated_key=sct_integrated_key,
        chunk_size=chunk_size,
        n_pcs=n_pcs,
        harmony_integrate=harmony_integrate,
        reuse_neighbors=reuse_neighbors,
    )
    _run_yardcluster_from_channels(
        adata,
        n_key=n_key,
        id_key=id_key,
        ctx_key=ctx_key,
        lam=lam,
        mode_label=mode_label,
        k_expr=k_expr,
        resolution=resolution,
        key_added=key_added,
        neighbors_key=neighbors_key,
        merge_clusters=merge_clusters,
        merge_adj_p=merge_adj_p,
        merge_max_de_tests=merge_max_de_tests,
        sketch=sketch,
        sketch_threshold=sketch_threshold,
        sketch_n=sketch_n,
        sketch_stratified=sketch_stratified,
        batch_key=batch_key,
        neighbors_method=neighbors_method,
        hybrid_alpha=hybrid_alpha,
        method=method,
        n_clusters=n_clusters,
        n_resolution_runs=n_resolution_runs,
        random_state=random_state,
    )


def spatial_cluster(
    adata: AnnData,
    mode: ModeType = "auto",
    lambda_celltype: float = 0.2,
    lambda_domain: float = 0.8,
    platform: Optional[PlatformType] = None,
    k_spatial: Optional[int] = None,
    k_gradient: Optional[int] = None,
    k_expr: Optional[int] = None,
    use_gradient: Optional[bool] = None,
    gradient_mode: GradientMode = "cosphi",
    n_pcs: Optional[int] = None,
    n_layers: int = 1,
    resolution: Optional[ResolutionType] = None,
    resolution_celltype: Optional[ResolutionType] = None,
    resolution_domain: Optional[ResolutionType] = None,
    method: MethodType = "leiden",
    n_clusters: Optional[int] = None,
    n_resolution_runs: int = 5,
    random_state: int = 0,
    batch_key: Optional[str] = None,
    integrate: IntegrateType = "none",
    harmony_integrate: bool = False,
    graph_mode: GraphModeType = "expr",
    hybrid_alpha: float = 0.3,
    hvg_by_batch: bool = True,
    scale_by_batch: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    expression_source: ExpressionSourceType = "x",
    sct_prep_key: str = "sct_prep",
    sct_integrated_key: str = "sct_integrated",
    chunk_size: int = 500,
    preprocess: bool = True,
    n_top_genes: Optional[int] = None,
    merge_clusters: Optional[bool] = None,
    merge_clusters_celltype: Optional[bool] = None,
    merge_clusters_domain: Optional[bool] = None,
    merge_adj_p: float = 0.05,
    merge_max_de_tests: Optional[int] = 2000,
    sketch: bool = False,
    sketch_threshold: int = _SKETCH_THRESHOLD,
    sketch_n: int = _SKETCH_N,
    sketch_stratified: bool = True,
    reuse_neighbors: bool = True,
    key_added: str = "yardcluster",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    YardCluster spatial clustering (one-shot or dual-mode).

    Multi-slice: ``expression_source='sct_prep'`` after SCT; ``expression_source='sct_integrated'``
    after SCT+RPCA (M/G from prep, identity from RPCA). Use ``integrate='bbknn'``
    or ``'joint'`` (Harmony) when batch correction is not already in the embedding.

    Platform presets (``platform='visium_hd'|'sto'|'xenium'|'xenium_tma'``) fill ``k_spatial``,
    ``n_pcs``, ``resolution``, ``k_expr``, ``n_top_genes`` when those arguments are
    ``None``.

    ``graph_mode``: ``'expr'`` (default), ``'spatial'`` (Leiden on spatial graph),
    or ``'hybrid'`` (weighted mix; tune with ``hybrid_alpha``).

    ``method``: ``'leiden'`` (default) or ``'gmm'`` (sklearn Gaussian mixture).
    ``resolution='auto'`` (or ``resolution_celltype``/``resolution_domain``)
    selects the resolution by clustering stability. ``n_layers`` enables
    CellCharter-style multi-hop neighborhood aggregation.
    """
    adata = adata.copy() if copy else adata

    resolved = resolve_platform_params(
        platform,
        k_spatial=k_spatial,
        k_expr=k_expr,
        n_pcs=n_pcs,
        resolution=resolution,
        n_top_genes=n_top_genes,
    )
    k_spatial = int(resolved["k_spatial"])
    k_expr = resolved["k_expr"]
    n_pcs = int(resolved["n_pcs"])
    _res = resolved["resolution"]
    resolution = _res if isinstance(_res, str) else float(_res)
    n_top_genes = int(resolved["n_top_genes"])

    resolved_use_gradient = _default_use_gradient(mode, use_gradient)
    neighbors_method = _resolve_neighbors_method(integrate, graph_mode)

    if integrate in ("joint", "bbknn") and batch_key is None:
        raise ValueError(f"integrate='{integrate}' requires batch_key.")
    if integrate == "bbknn" and harmony_integrate:
        warnings.warn(
            "harmony_integrate is ignored when integrate='bbknn'.",
            stacklevel=2,
        )
    use_harmony = harmony_integrate or integrate == "joint"

    if expression_source in ("sct_prep", "sct_integrated") and preprocess:
        warnings.warn(
            f"preprocess=True ignored when expression_source='{expression_source}' "
            "(using SCT matrices from obsm).",
            stacklevel=2,
        )
        preprocess = False

    if expression_source == "sct_integrated" and integrate != "none":
        warnings.warn(
            "integrate != 'none' ignored when expression_source='sct_integrated' "
            "(batch correction already applied via RPCA).",
            stacklevel=2,
        )
        integrate = "none"
        neighbors_method = _resolve_neighbors_method("none", graph_mode)
        use_harmony = False

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
                n_layers=n_layers,
                resolution=resolution,
                resolution_celltype=resolution_celltype,
                resolution_domain=resolution_domain,
                method=method,
                n_clusters=n_clusters,
                n_resolution_runs=n_resolution_runs,
                random_state=random_state,
                batch_key=None,
                integrate="none",
                preprocess=False,
                expression_source=expression_source,
                sct_prep_key=sct_prep_key,
                merge_clusters=merge_clusters,
                merge_clusters_celltype=merge_clusters_celltype,
                merge_clusters_domain=merge_clusters_domain,
                merge_adj_p=merge_adj_p,
                merge_max_de_tests=merge_max_de_tests,
                sketch=sketch,
                sketch_threshold=sketch_threshold,
                sketch_n=sketch_n,
                sketch_stratified=sketch_stratified,
                platform=platform,
                graph_mode=graph_mode,
                hybrid_alpha=hybrid_alpha,
                reuse_neighbors=reuse_neighbors,
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
    feature_kw = dict(
        k_spatial=k_spatial,
        k_gradient=k_gradient,
        n_layers=n_layers,
        batch_key=batch_key,
        use_gradient=resolved_use_gradient,
        gradient_mode=gradient_mode,
        layer=layer,
        use_raw=use_raw,
        expression_source=expression_source,
        sct_prep_key=sct_prep_key,
        sct_integrated_key=sct_integrated_key,
        chunk_size=chunk_size,
        n_pcs=n_pcs,
        harmony_integrate=use_harmony,
        reuse_neighbors=reuse_neighbors,
    )
    cluster_base_kw = dict(
        k_expr=k_expr,
        key_added=key_added,
        merge_adj_p=merge_adj_p,
        merge_max_de_tests=merge_max_de_tests,
        sketch=sketch,
        sketch_threshold=sketch_threshold,
        sketch_n=sketch_n,
        sketch_stratified=sketch_stratified,
        batch_key=batch_key,
        neighbors_method=neighbors_method,
        hybrid_alpha=hybrid_alpha,
        method=method,
        n_clusters=n_clusters,
        n_resolution_runs=n_resolution_runs,
        random_state=random_state,
    )

    def _cluster_kw_for(mode_label: str) -> dict:
        return {
            **cluster_base_kw,
            "merge_clusters": _resolve_merge_clusters(
                mode_label,
                merge_clusters,
                merge_clusters_celltype,
                merge_clusters_domain,
            ),
        }

    if mode == "auto":
        n_key = f"{key_added}_features"
        id_key, ctx_key = _run_shared_yard_features(adata, n_key=n_key, **feature_kw)
        _run_yardcluster_from_channels(
            adata,
            n_key=n_key,
            id_key=id_key,
            ctx_key=ctx_key,
            lam=lambda_celltype,
            mode_label="celltype",
            resolution=res_ct,
            neighbors_key=f"{key_added}_celltype",
            **_cluster_kw_for("celltype"),
        )
        _run_yardcluster_from_channels(
            adata,
            n_key=n_key,
            id_key=id_key,
            ctx_key=ctx_key,
            lam=lambda_domain,
            mode_label="domain",
            resolution=res_dm,
            neighbors_key=f"{key_added}_domain",
            **_cluster_kw_for("domain"),
        )
    else:
        mode_label = "celltype" if mode == "celltype" else "domain"
        run_kw = {
            **feature_kw,
            **_cluster_kw_for(mode_label),
        }
        if mode == "celltype":
            _run_single_yardcluster(
                adata,
                lam=lambda_celltype,
                mode_label="default",
                resolution=res_ct,
                neighbors_key=key_added,
                **run_kw,
            )
        elif mode == "domain":
            _run_single_yardcluster(
                adata,
                lam=lambda_domain,
                mode_label="default",
                resolution=res_dm,
                neighbors_key=key_added,
                **run_kw,
            )

    _store_params(
        adata,
        f"{key_added}_params",
        {
            "mode": mode,
            "integrate": integrate,
            "harmony_integrate": use_harmony,
            "use_gradient": resolved_use_gradient,
            "merge_clusters": merge_clusters,
            "merge_clusters_celltype": merge_clusters_celltype,
            "merge_clusters_domain": merge_clusters_domain,
            "merge_max_de_tests": merge_max_de_tests,
            "gradient_mode": gradient_mode,
            "n_layers": n_layers,
            "method": method,
            "n_clusters": n_clusters,
            "n_resolution_runs": n_resolution_runs,
            "random_state": random_state,
            "sketch": sketch,
            "expression_source": expression_source,
            "sct_prep_key": sct_prep_key if expression_source == "sct_prep" else None,
            "platform": platform,
            "graph_mode": graph_mode,
            "hybrid_alpha": hybrid_alpha,
            "neighbors_method": neighbors_method,
            "reuse_neighbors": reuse_neighbors,
            "sketch_stratified": sketch_stratified,
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
