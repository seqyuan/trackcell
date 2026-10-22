"""Synthetic spatial datasets for YardCluster evaluation."""

from __future__ import annotations

import numpy as np
from anndata import AnnData


def make_patch_adata(
    *,
    n_per_patch: int = 40,
    n_patches: int = 3,
    patch_spacing: float = 6.0,
    n_genes: int = 80,
    seed: int = 0,
) -> AnnData:
    """
    Build a grid of spatial expression patches with known region labels.

    Each patch is a square cluster of cells with distinct marker genes.
    """
    rng = np.random.default_rng(seed)
    coords_list = []
    labels = []
    for p in range(n_patches):
        row, col = divmod(p, 2)
        cx = col * patch_spacing
        cy = row * patch_spacing
        jitter = rng.uniform(-0.4, 0.4, size=(n_per_patch, 2))
        patch_coords = np.column_stack(
            [
                cx + rng.uniform(-1.2, 1.2, n_per_patch),
                cy + rng.uniform(-1.2, 1.2, n_per_patch),
            ]
        ) + jitter
        coords_list.append(patch_coords)
        labels.extend([f"patch_{p}"] * n_per_patch)

    coords = np.vstack(coords_list)
    n_obs = coords.shape[0]
    X = rng.normal(size=(n_obs, n_genes)).astype(np.float32)
    genes_per_patch = max(3, n_genes // (n_patches + 2))
    for p in range(n_patches):
        start = p * n_per_patch
        end = start + n_per_patch
        g0 = (p * genes_per_patch) % n_genes
        g1 = min(g0 + genes_per_patch, n_genes)
        X[start:end, g0:g1] += 4.0

    adata = AnnData(X=X)
    adata.obsm["spatial"] = coords
    adata.obs["patch"] = np.array(labels, dtype=str)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    return adata
