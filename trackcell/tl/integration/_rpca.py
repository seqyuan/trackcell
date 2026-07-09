"""Reciprocal PCA helpers for Seurat-style integration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.decomposition import PCA


@dataclass
class BatchPCA:
    embeddings: np.ndarray  # cells x dims
    loadings: np.ndarray  # features x dims
    mean: np.ndarray  # features,


def run_batch_pca(
    data: np.ndarray,
    *,
    n_components: int,
    seed: int = 0,
) -> BatchPCA:
    """
    PCA on centered residual matrix (cells x features), Seurat ``RunPCA``-like.
    """
    if data.ndim != 2:
        raise ValueError("data must be 2D (cells x features).")
    n_components = min(n_components, data.shape[0] - 1, data.shape[1])
    if n_components < 1:
        raise ValueError("Cannot compute PCA with fewer than 1 component.")
    centered = data - data.mean(axis=0, keepdims=True)
    pca = PCA(n_components=n_components, random_state=seed, svd_solver="randomized")
    embeddings = pca.fit_transform(centered)
    loadings = pca.components_.T  # features x dims
    return BatchPCA(embeddings=embeddings, loadings=loadings, mean=pca.mean_)


def project_svd(
    loadings: np.ndarray,
    data: np.ndarray,
    *,
    feature_idx: np.ndarray | None = None,
    center: bool = False,
) -> np.ndarray:
    """
    Project cells x features with reference loadings (Seurat ``ProjectSVD``).

    Seurat SCT+RPCA uses ``do.center=FALSE`` and ``do.scale=FALSE``.
    """
    if feature_idx is not None:
        loadings = loadings[feature_idx]
        data = data[:, feature_idx]
    if center:
        data = data - data.mean(axis=0, keepdims=True)
    return data @ loadings


def l2_normalize_rows(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.maximum(norms, eps)


def l2_normalize_cols_by_sd(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    sd = x.std(axis=0, keepdims=True)
    return x / np.maximum(sd, eps)


def reciprocal_project(
    data1: np.ndarray,
    pca1: BatchPCA,
    data2: np.ndarray,
    pca2: BatchPCA,
    *,
    dims: int,
    l2_norm: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reciprocal PCA projections for a batch pair (Seurat ``ReciprocalProject``).

    Returns
    -------
    ref_embedding, query_embedding
        Both are (n1+n2) x dims matrices used for anchor finding.
    """
    dims = min(dims, pca1.embeddings.shape[1], pca2.embeddings.shape[1])
    n_feat = min(data1.shape[1], data2.shape[1], pca1.loadings.shape[0], pca2.loadings.shape[0])
    feature_idx = np.arange(n_feat)
    loadings1 = pca1.loadings[feature_idx, :dims]
    loadings2 = pca2.loadings[feature_idx, :dims]
    emb1 = pca1.embeddings[:, :dims]
    emb2 = pca2.embeddings[:, :dims]
    proj1 = project_svd(loadings2, data1, feature_idx=feature_idx, center=False)[:, :dims]
    proj2 = project_svd(loadings1, data2, feature_idx=feature_idx, center=False)[:, :dims]
    ref_embedding = np.vstack([emb1, proj2])
    query_embedding = np.vstack([proj1, emb2])
    if l2_norm:
        ref_embedding = l2_normalize_cols_by_sd(ref_embedding)
        query_embedding = l2_normalize_cols_by_sd(query_embedding)
        ref_embedding = l2_normalize_rows(ref_embedding)
        query_embedding = l2_normalize_rows(query_embedding)
    return ref_embedding, query_embedding
