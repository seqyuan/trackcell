"""
Spatial neighborhood graph and neighborhood feature computation for YardCluster.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sps
from scipy.spatial import cKDTree
from typing import Optional, Tuple, Union

ArrayLike = Union[np.ndarray, sps.spmatrix]


def _as_dense(X: ArrayLike) -> np.ndarray:
    if sps.issparse(X):
        return X.toarray()
    return np.asarray(X, dtype=np.float64)


def _batch_groups(batch_labels: Optional[np.ndarray], n_obs: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return (global_indices, coords_subset) per batch."""
    if batch_labels is None:
        idx = np.arange(n_obs, dtype=np.int64)
        return [(idx, idx)]
    batch_labels = np.asarray(batch_labels)
    batches = np.unique(batch_labels)
    groups = []
    for b in batches:
        idx = np.where(batch_labels == b)[0]
        groups.append((idx, idx))
    return groups


def build_spatial_weights(
    coords: np.ndarray,
    k: int = 15,
    batch_labels: Optional[np.ndarray] = None,
) -> sps.csr_matrix:
    """
    Build a row-normalized sparse spatial weight matrix using Gaussian-kernel kNN.

    Neighbors are restricted to within-batch when ``batch_labels`` is provided.
    Weights follow BANKSY-style scale-invariant Gaussian modulation using the
    distance to the k-th neighbor as bandwidth.

    Parameters
    ----------
    coords
        (n_obs, 2) spatial coordinates.
    k
        Number of spatial neighbors (excluding self).
    batch_labels
        Optional batch/sample id per observation; edges do not cross batches.

    Returns
    -------
    scipy.sparse.csr_matrix
        Row-normalized weights, shape (n_obs, n_obs).
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_obs = coords.shape[0]
    if n_obs == 0:
        raise ValueError("coords must contain at least one point.")
    if coords.shape[1] < 2:
        raise ValueError("coords must have at least 2 dimensions.")
    k = int(k)
    if k < 1:
        raise ValueError("k must be >= 1.")
    k_query = min(k + 1, n_obs)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    if batch_labels is None:
        batch_iter = [(np.arange(n_obs, dtype=np.int64), coords)]
    else:
        batch_labels = np.asarray(batch_labels)
        batch_iter = []
        for b in np.unique(batch_labels):
            idx = np.where(batch_labels == b)[0]
            batch_iter.append((idx, coords[idx]))

    for global_idx, batch_coords in batch_iter:
        n_batch = batch_coords.shape[0]
        if n_batch == 1:
            rows.append(int(global_idx[0]))
            cols.append(int(global_idx[0]))
            data.append(1.0)
            continue

        k_local = min(k_query, n_batch)
        tree = cKDTree(batch_coords)
        dists, local_nbrs = tree.query(batch_coords, k=k_local)

        if k_local == 1:
            dists = dists[:, np.newaxis]
            local_nbrs = local_nbrs[:, np.newaxis]

        for i, gi in enumerate(global_idx):
            neighbor_dists = dists[i]
            neighbor_local = local_nbrs[i]

            # Drop self (zero distance)
            mask = neighbor_dists > 0
            if not np.any(mask):
                rows.append(int(gi))
                cols.append(int(gi))
                data.append(1.0)
                continue

            nd = neighbor_dists[mask]
            nl = neighbor_local[mask]
            r_k = nd[-1] if nd.size >= k else nd.max()
            if r_k <= 0:
                r_k = 1.0

            weights = np.exp(-(nd / r_k) ** 2)
            weights /= weights.sum()

            for lj, w in zip(nl, weights):
                rows.append(int(gi))
                cols.append(int(global_idx[lj]))
                data.append(float(w))

    adj = sps.csr_matrix((data, (rows, cols)), shape=(n_obs, n_obs))
    return adj


def spatial_distances_from_weights(weights: sps.csr_matrix) -> sps.csr_matrix:
    """Convert weight matrix to a distance-like matrix (1 - weight on edges)."""
    dist = weights.copy()
    dist.data = 1.0 - dist.data
    return dist.tocsr()


def compute_neighborhood_mean(
    X: ArrayLike,
    weights: sps.csr_matrix,
    chunk_size: int = 500,
) -> np.ndarray:
    """
    Compute neighborhood mean expression M = weights @ X in gene chunks.

    Parameters
    ----------
    X
        (n_obs, n_genes) expression matrix.
    weights
        Row-normalized spatial weight matrix.
    chunk_size
        Genes per chunk for memory control.

    Returns
    -------
    np.ndarray
        Neighborhood mean matrix, same shape as X.
    """
    if sps.issparse(X):
        n_obs, n_genes = X.shape
        M = np.empty((n_obs, n_genes), dtype=np.float64)
        for start in range(0, n_genes, chunk_size):
            end = min(start + chunk_size, n_genes)
            block = X[:, start:end].toarray()
            M[:, start:end] = weights @ block
        return M

    X = np.asarray(X, dtype=np.float64)
    n_obs, n_genes = X.shape
    M = np.empty((n_obs, n_genes), dtype=np.float64)
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        M[:, start:end] = weights @ X[:, start:end]
    return M


def compute_neighborhood_gradient_cosphi(
    X: ArrayLike,
    coords: np.ndarray,
    weights: sps.csr_matrix,
    k_gradient: Optional[int] = None,
    chunk_size: int = 500,
) -> np.ndarray:
    """
    Approximate BANKSY AGF with a real-valued cos(phi) directional weighting.

    For each cell u and gene q:
        G_uq = sum_v w_uv * x_vq * cos(phi_uv)

    where phi_uv is the azimuth of neighbor v in u's local polar frame.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_obs = coords.shape[0]
    X_dense = _as_dense(X)
    _, n_genes = X_dense.shape

    if k_gradient is not None:
        weights = truncate_weights_per_row(weights, k_gradient)

    G = np.zeros((n_obs, n_genes), dtype=np.float64)
    weights = weights.tocsr()

    for u in range(n_obs):
        start, end = weights.indptr[u], weights.indptr[u + 1]
        if start == end:
            continue
        nbr_cols = weights.indices[start:end]
        nbr_w = weights.data[start:end]

        dx = coords[nbr_cols, 0] - coords[u, 0]
        dy = coords[nbr_cols, 1] - coords[u, 1]
        phi = np.arctan2(dy, dx)
        dir_w = nbr_w * np.cos(phi)
        dir_w /= np.abs(dir_w).sum() if np.abs(dir_w).sum() > 0 else 1.0

        for g_start in range(0, n_genes, chunk_size):
            g_end = min(g_start + chunk_size, n_genes)
            G[u, g_start:g_end] = dir_w @ X_dense[nbr_cols, g_start:g_end]

    return G


def truncate_weights_per_row(weights: sps.csr_matrix, k: int) -> sps.csr_matrix:
    """Keep only top-k neighbors per row (by weight), re-normalize rows."""
    weights = weights.tocsr()
    n_obs = weights.shape[0]
    rows, cols, data = [], [], []
    for u in range(n_obs):
        start, end = weights.indptr[u], weights.indptr[u + 1]
        if start == end:
            rows.append(u)
            cols.append(u)
            data.append(1.0)
            continue
        c = weights.indices[start:end]
        d = weights.data[start:end]
        if len(d) > k:
            top = np.argsort(d)[-k:]
            c = c[top]
            d = d[top]
        d = d / d.sum()
        rows.extend([u] * len(c))
        cols.extend(c.tolist())
        data.extend(d.tolist())
    return sps.csr_matrix((data, (rows, cols)), shape=(n_obs, n_obs))


def mix_embeddings(
    identity: np.ndarray,
    context: np.ndarray,
    lam: float,
) -> np.ndarray:
    """
    BANKSY-style sqrt-weighted concatenation of identity and context PCA spaces.

    Z = sqrt(1 - lam) * identity || sqrt(lam) * context
    """
    lam = float(lam)
    if not 0.0 <= lam <= 1.0:
        raise ValueError("lam must be in [0, 1].")
    return np.hstack(
        [np.sqrt(1.0 - lam) * identity, np.sqrt(lam) * context]
    ).astype(np.float32)


def compute_neighborhood_gradient_regression(
    X: ArrayLike,
    coords: np.ndarray,
    weights: sps.csr_matrix,
    k_gradient: Optional[int] = None,
    chunk_size: int = 500,
    ridge: float = 1e-6,
) -> np.ndarray:
    """
    Local linear regression gradient magnitude: ||(beta_x, beta_y)|| per gene.

    For each cell u, fit expr ~ beta0 + beta1*dx + beta2*dy across spatial neighbors
    in local coordinates centered at u.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_obs = coords.shape[0]
    X_dense = _as_dense(X)
    _, n_genes = X_dense.shape

    if k_gradient is not None:
        weights = truncate_weights_per_row(weights, k_gradient)

    G = np.zeros((n_obs, n_genes), dtype=np.float64)
    weights = weights.tocsr()
    eye3 = np.eye(3) * ridge

    for u in range(n_obs):
        start, end = weights.indptr[u], weights.indptr[u + 1]
        if start == end:
            continue
        nbr_cols = weights.indices[start:end]
        if len(nbr_cols) < 3:
            continue

        dx = coords[nbr_cols, 0] - coords[u, 0]
        dy = coords[nbr_cols, 1] - coords[u, 1]
        A = np.column_stack([np.ones(len(nbr_cols)), dx, dy])
        ata_inv = np.linalg.inv(A.T @ A + eye3)

        for g_start in range(0, n_genes, chunk_size):
            g_end = min(g_start + chunk_size, n_genes)
            y = X_dense[nbr_cols, g_start:g_end]
            beta = ata_inv @ A.T @ y
            # beta shape (3, n_genes_chunk); rows 1,2 are spatial slopes
            G[u, g_start:g_end] = np.linalg.norm(beta[1:3, :], axis=0)

    return G


def compute_neighborhood_gradient(
    X: ArrayLike,
    coords: np.ndarray,
    weights: sps.csr_matrix,
    mode: str = "cosphi",
    k_gradient: Optional[int] = None,
    chunk_size: int = 500,
    ridge: float = 1e-6,
) -> np.ndarray:
    """Dispatch gradient computation by mode ('cosphi' or 'regression')."""
    mode = mode.lower()
    if mode == "cosphi":
        return compute_neighborhood_gradient_cosphi(
            X, coords, weights, k_gradient=k_gradient, chunk_size=chunk_size
        )
    if mode == "regression":
        return compute_neighborhood_gradient_regression(
            X, coords, weights, k_gradient=k_gradient, chunk_size=chunk_size, ridge=ridge
        )
    raise ValueError(f"Unknown gradient mode: {mode}. Use 'cosphi' or 'regression'.")


def expr_k_neighbors(n_obs: int, k_expr: Optional[int] = None) -> int:
    """Expression-space kNN count (Space Ranger-inspired log scaling)."""
    if k_expr is not None:
        return max(2, int(k_expr))
    if n_obs <= 1:
        return 2
    return max(15, int(np.floor(np.log2(n_obs))))
