"""
Tools for spatial analyses.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


def hd_labeldist(adata, groupby: str, label: str, inplace: bool = True, method: str = "kdtree"):
    """
    Compute the distance from every cell to the nearest cell annotated with a specific label (10x HD data).
    
    Distances are reported both in the pixel coordinate system stored in `.obsm["spatial"]`
    (SpaceRanger target/hires layer) and in microns using the scalefactors embedded in
    `adata.uns["spatial"]`.
    
    The function automatically detects whether coordinates are in hires or full-res resolution
    by comparing the computed tissue size with the expected chip size (6.5mm for 10X HD).
    This ensures compatibility with both SpaceRanger output and bin2cell-processed data.
    
    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix with spatial coordinates and SpaceRanger scalefactors.
    groupby : str
        Column name in `adata.obs` that contains the annotation labels.
    label : str
        Target label within `groupby` for which distances will be computed.
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
    if groupby not in adata.obs.columns:
        raise ValueError(f"`{groupby}` not found in `adata.obs`.")
    if label not in adata.obs[groupby].unique():
        raise ValueError(f"`{label}` not present in `adata.obs['{groupby}']`.")
    
    coords = adata.obsm["spatial"]
    if coords is None or len(coords) == 0:
        raise ValueError("Spatial coordinates are empty.")
    
    mask_label = (adata.obs[groupby] == label).to_numpy()
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
    
    # Auto-detect coordinate resolution by comparing computed size with expected chip size
    # 10X HD chip is 6.5mm x 6.5mm, VisiumHD is similar
    # Calculate coordinate range in both possible resolutions
    coord_range = np.ptp(coords_all, axis=0)  # [max_x - min_x, max_y - min_y]
    max_range_px = np.max(coord_range)
    
    # Test both hypotheses:
    # 1. Coordinates are hires: need to divide by hires_scale
    size_um_hires = max_range_px * (microns_per_pixel / hires_scale)
    # 2. Coordinates are full-res: use directly
    size_um_fullres = max_range_px * microns_per_pixel
    
    # Expected chip size: 6.5mm = 6500um (with some tolerance for tissue coverage)
    # Typical tissue coverage is 60-90% of chip, so reasonable range is 4000-7000um
    expected_max_size_um = 7000  # Upper bound for reasonable chip size
    
    # Determine which resolution gives more reasonable results
    # If hires calculation gives reasonable size (< expected_max_size_um), coordinates are hires
    # If fullres calculation gives reasonable size, coordinates are full-res
    # If hires_scale is 1.0 or very close, assume full-res
    if hires_scale < 1.01:  # hires_scale close to 1.0 means likely full-res
        is_hires = False
    elif size_um_hires <= expected_max_size_um and size_um_fullres > expected_max_size_um:
        # hires calculation gives reasonable result, fullres doesn't
        is_hires = True
    elif size_um_fullres <= expected_max_size_um and size_um_hires > expected_max_size_um:
        # fullres calculation gives reasonable result, hires doesn't
        is_hires = False
    else:
        # Both or neither give reasonable results, prefer hires if hires_scale is significantly < 1.0
        # This handles edge cases where tissue might be larger than chip
        is_hires = (hires_scale < 0.5)  # If hires_scale is significantly < 1, likely hires coords
    
    # Calculate physical distance based on detected resolution
    if is_hires:
        dist_um = dist_px * (microns_per_pixel / hires_scale)
    else:
        dist_um = dist_px * microns_per_pixel
    
    col_px = f"{label}_px"
    col_um = f"{label}_dist"
    
    if inplace:
        adata.obs[col_px] = dist_px
        adata.obs[col_um] = dist_um
        return None
    
    return pd.DataFrame({col_px: dist_px, col_um: dist_um}, index=adata.obs.index)