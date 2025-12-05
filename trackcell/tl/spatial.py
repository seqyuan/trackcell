"""
Tools for spatial analyses.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


def hd_labeldist(adata, group: str, label: str, inplace: bool = True, method: str = "kdtree"):
    """
    Compute the distance from every cell to the nearest cell annotated with a specific label (10x HD data).
    
    Distances are reported both in the pixel coordinate system stored in `.obsm["spatial"]`
    (SpaceRanger target/hires layer) and in microns using the scalefactors embedded in
    `adata.uns["spatial"]`.
    
    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix with spatial coordinates and SpaceRanger scalefactors.
    group : str
        Column name in `adata.obs` that contains the annotation labels.
    label : str
        Target label within `group` for which distances will be computed.
    inplace : bool, default True
        If True, the function adds two columns to `adata.obs`:
        `{label}_px` (pixel distance on the hires/registered image) and
        `{label}_dist` (physical distance in microns).
        If False, the function returns a dataframe with the two columns.
    method : str, default "kdtree"
        Method to use for distance computation:
        - "kdtree": Use KDTree spatial indexing (recommended, O(n log n) time, O(n) memory).
                   Best for large datasets with many cells.
        - "cdist": Use scipy's cdist function (O(n*m) time and memory, where m is number of label cells).
                   Faster for small datasets but memory-intensive for large ones.
    
    Returns
    -------
    pandas.DataFrame | None
        Returns a DataFrame with `{label}_px` and `{label}_dist` when `inplace=False`.
        Otherwise, modifies `adata.obs` in place and returns None.
    """
    if "spatial" not in adata.obsm:
        raise ValueError("`adata.obsm['spatial']` is required but missing.")
    if group not in adata.obs.columns:
        raise ValueError(f"`{group}` not found in `adata.obs`.")
    if label not in adata.obs[group].unique():
        raise ValueError(f"`{label}` not present in `adata.obs['{group}']`.")
    
    coords = adata.obsm["spatial"]
    if coords is None or len(coords) == 0:
        raise ValueError("Spatial coordinates are empty.")
    
    mask_label = (adata.obs[group] == label).to_numpy()
    if not mask_label.any():
        raise ValueError(f"No observations with label `{label}` were found.")
    
    coords_all = np.asarray(coords, dtype=float)
    coords_label = coords_all[mask_label]
    
    # Compute distances using selected method
    if method == "kdtree":
        # Method 1: KDTree (recommended for large datasets)
        # Memory: O(n + m), Time: O(n log m) where n=all cells, m=label cells
        tree = cKDTree(coords_label)
        dist_px, _ = tree.query(coords_all, k=1)
    elif method == "cdist":
        # Method 2: cdist (faster for small datasets, but memory-intensive)
        # Memory: O(n*m), Time: O(n*m)
        dist_matrix = cdist(coords_all, coords_label, metric='euclidean')
        dist_px = dist_matrix.min(axis=1)
    else:
        raise ValueError(f"Unknown method: {method}. Choose 'kdtree' or 'cdist'.")
    
    spatial_meta = adata.uns.get("spatial")
    if not spatial_meta:
        raise ValueError("`adata.uns['spatial']` is missing scalefactor information.")
    sample_key = next(iter(spatial_meta))
    scalefactors = spatial_meta[sample_key].get("scalefactors", {})
    
    microns_per_pixel = scalefactors.get("microns_per_pixel")
    if microns_per_pixel is None:
        raise ValueError("`microns_per_pixel` not found in scalefactors.")
    
    hires_scale = scalefactors.get("tissue_hires_scalef") or scalefactors.get("regist_target_img_scalef") or 1.0
    dist_um = dist_px * (microns_per_pixel / hires_scale)
    
    col_px = f"{label}_px"
    col_um = f"{label}_dist"
    
    if inplace:
        adata.obs[col_px] = dist_px
        adata.obs[col_um] = dist_um
        return None
    
    return pd.DataFrame({col_px: dist_px, col_um: dist_um}, index=adata.obs.index)

