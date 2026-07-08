"""
Tools for spatial analyses.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from typing import Optional, Union, List
import warnings


# Perceptually distinct colors for up to 7 genes (ColorBrewer-inspired)
DEFAULT_GENE_COLORS = [
    '#e41a1c',  # Red
    '#377eb8',  # Blue
    '#4daf4a',  # Green
    '#ff7f00',  # Orange
    '#984ea3',  # Purple
    '#a65628',  # Brown
    '#f781bf',  # Pink
]


def _hex_to_rgb(hex_color: str) -> tuple:
    """Convert hex color string to (R, G, B) tuple in [0, 1]."""
    h = hex_color.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def _rgb_to_hex(r: float, g: float, b: float) -> str:
    """Convert (R, G, B) in [0, 1] to hex string."""
    def clamp(v):
        return max(0, min(255, int(round(v * 255))))
    return f'#{clamp(r):02x}{clamp(g):02x}{clamp(b):02x}'


def _multigene_facet(
    adata,
    genes: List[str],
    colors: List[str],
    layer: Optional[str] = None,
    vmin_percentile: float = 1.0,
    vmax_percentile: float = 99.0,
    ncols: int = 3,
    figsize: Optional[tuple] = None,
    edges_width: float = 0.3,
    edges_color: str = '#cccccc30',
    alpha: float = 0.8,
    show: bool = True,
    ax=None,
    library_id: Optional[str] = None,
):
    """Faceted multi-gene subplots, one gene per panel with single-hue colormap."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    n_genes = len(genes)
    nrows = int(np.ceil(n_genes / ncols))
    
    if figsize is None:
        figsize = (5 * ncols, 5 * nrows)
    
    if ax is None:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    else:
        axes = np.atleast_2d(np.asarray(ax))
        fig = axes.flat[0].figure
    
    axes_flat = axes.flatten()
    
    # Temporarily suppress auto-show in spatial_cell
    for i, (gene, color) in enumerate(zip(genes, colors)):
        ax_i = axes_flat[i]
        # Single-hue colormap: white → gene color
        gene_cmap = LinearSegmentedColormap.from_list(
            f'_{gene}', ['#fafafa', color]
        )
        
        # Call spatial_cell with gene expression
        from trackcell.pl.plot import spatial_cell
        spatial_cell(
            adata,
            color=gene,
            cmap=gene_cmap,
            ax=ax_i,
            show=False,
            legend=True,
            edges_width=edges_width,
            edges_color=edges_color,
            alpha=alpha,
            library_id=library_id,
        )
        ax_i.set_title(gene, fontsize=12, fontweight='bold')
    
    # Hide unused axes
    for j in range(n_genes, len(axes_flat)):
        axes_flat[j].set_visible(False)
    
    fig.tight_layout(rect=[0, 0, 1, 1])
    
    if show:
        plt.show()
    
    return axes


def multigene_blend(
    adata,
    genes: List[str],
    mode: str = 'blend',
    colors: Optional[List[str]] = None,
    layer: Optional[str] = None,
    vmin_percentile: float = 1.0,
    vmax_percentile: float = 99.0,
    gamma: float = 1.0,
    background: str = '#e0e0e0',
    key_added: str = 'multigene_blend',
    # Facet-specific
    ncols: int = 3,
    figsize: Optional[tuple] = None,
    edges_width: float = 0.3,
    edges_color: str = '#cccccc30',
    alpha_facet: float = 0.8,
    show: bool = True,
    ax=None,
    library_id: Optional[str] = None,
    inplace: bool = True,
):
    """
    Multi-gene co-expression visualization with two modes.
    
    **Mode 'blend'** (default): Computes blended hex colors by weighted RGB averaging
    of per-gene base colors. The result is a hex color per cell stored in
    ``adata.obs[key_added]``, usable with ``spatial_cell(color=key_added)``.
    
    **Mode 'facet'**: Draws faceted subplots, one gene per panel, each using
    a single-hue colormap (white → gene color). Returns matplotlib axes.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data with gene expression.
    genes : list of str
        Gene names. Up to 7 recommended for blend mode.
    mode : str
        'blend' (default) or 'facet'.
    colors : list of str, optional
        Base colors for each gene (hex). Default: perceptually distinct palette.
    layer : str, optional
        Use adata.layers[layer] instead of adata.X.
    vmin_percentile : float
        Lower percentile for clipping (blend & facet, default 1).
    vmax_percentile : float
        Upper percentile for clipping (blend & facet, default 99).
    gamma : float
        Gamma correction for blend mode (>1 darkens, <1 brightens).
    background : str
        Background hex color for blend mode (zero-expression cells).
    key_added : str
        obs column name for blend mode output.
    inplace : bool
        If True (blend), stores in adata.obs. If False, returns Series.
    
    **Facet-specific parameters:**
    ncols : int
        Number of columns in facet grid (default 3).
    figsize : tuple, optional
        Figure size. Auto-computed from ncols if None.
    edges_width : float
        Cell edge width for facet plots (default 0.3).
    edges_color : str
        Cell edge color for facet plots (default '#cccccc30').
    alpha : float
        Cell transparency for facet plots (default 0.8).
    show : bool
        Whether to call plt.show() (default True).
    ax : matplotlib.Axes or array, optional
        Pre-existing axes to plot on.
    library_id : str, optional
        Library id passed to spatial_cell.
    
    Returns
    -------
    blend mode: pd.Series or None
    facet mode: np.ndarray of matplotlib.Axes
    
    Examples
    --------
    >>> import trackcell as tcl
    >>>
    >>> # Mode 1: Blend → single composite image
    >>> tcl.tl.multigene_blend(adata, genes=['EPCAM', 'PECAM1', 'VWF'])
    >>> tcl.pl.spatial_cell(adata, color='multigene_blend')
    >>>
    >>> # Mode 2: Facet → subplots per gene
    >>> tcl.tl.multigene_blend(
    ...     adata, genes=['EPCAM', 'PECAM1', 'VWF'],
    ...     mode='facet', ncols=3,
    ... )
    """
    n_genes = len(genes)
    if n_genes < 2:
        raise ValueError(f'Need at least 2 genes, got {n_genes}')
    
    # Validate genes
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise ValueError(f'Genes not found in adata.var_names: {missing}')
    
    # Default colors
    if colors is None:
        colors = DEFAULT_GENE_COLORS[:n_genes]
    elif len(colors) != n_genes:
        raise ValueError(
            f'colors length ({len(colors)}) must match genes length ({n_genes})'
        )
    
    if mode == 'facet':
        return _multigene_facet(
            adata=adata, genes=genes, colors=colors,
            layer=layer, vmin_percentile=vmin_percentile,
            vmax_percentile=vmax_percentile,
            ncols=ncols, figsize=figsize,
            edges_width=edges_width, edges_color=edges_color,
            alpha=alpha_facet, show=show, ax=ax, library_id=library_id,
        )
    
    if mode != 'blend':
        raise ValueError(
            f"Unknown mode '{mode}'. Valid modes: 'blend', 'facet'."
        )
    
    # ── Blend mode ──
    if n_genes > 7:
        warnings.warn(
            f'{n_genes} genes may produce muddy blends. '
            f'Consider using <= 7 genes or mode="facet".'
        )
    
    # ── Extract expression matrix ──
    X = adata[:, genes].X if layer is None else adata[:, genes].layers[layer]
    if hasattr(X, 'toarray'):
        X = X.toarray()
    else:
        X = np.asarray(X, dtype=float)
    
    n_cells = X.shape[0]
    
    # ── Convert base colors to RGB ──
    base_rgb = np.array([_hex_to_rgb(c) for c in colors])  # (n_genes, 3)
    
    # ── Per-gene normalization (percentile clipping) ──
    X_norm = np.zeros_like(X)
    for j in range(n_genes):
        col = X[:, j]
        vmin = np.percentile(col, vmin_percentile)
        vmax = np.percentile(col, vmax_percentile)
        if vmax > vmin:
            col_clipped = np.clip(col, vmin, vmax)
            X_norm[:, j] = (col_clipped - vmin) / (vmax - vmin)
        else:
            X_norm[:, j] = 0.0
    
    # ── Weighted RGB blending ──
    # weight_i = norm_expr_i / (sum(norm_expr) + eps)
    # RGB = Σ(weight_i * base_rgb_i)
    total = X_norm.sum(axis=1, keepdims=True)  # (n_cells, 1)
    weights = np.divide(X_norm, total, where=total > 0)  # (n_cells, n_genes)
    rgb = weights @ base_rgb  # (n_cells, 3)
    
    # ── Gamma correction ──
    if gamma != 1.0:
        rgb = np.power(rgb, gamma)
    
    # ── Convert to hex ──
    bg_rgb = np.array(_hex_to_rgb(background))
    hex_colors = []
    for i in range(n_cells):
        if total[i, 0] > 0:
            hex_colors.append(_rgb_to_hex(rgb[i, 0], rgb[i, 1], rgb[i, 2]))
        else:
            hex_colors.append(background)
    
    result = pd.Series(hex_colors, index=adata.obs_names, name=key_added)
    
    # ── Store metadata ──
    adata.uns[f'{key_added}_params'] = {
        'genes': genes,
        'colors': colors,
        'vmin_percentile': vmin_percentile,
        'vmax_percentile': vmax_percentile,
        'gamma': gamma,
        'background': background,
        'n_unique_colors': result.nunique(),
    }
    
    if inplace:
        adata.obs[key_added] = result
        print(
            f'[multigene_blend] {n_genes} genes → '
            f'{result.nunique():,} unique colors in adata.obs["{key_added}"]'
        )
        return None
    
    return result


def _resolve_labeldist_scalefactors(
    adata,
    library_id: Optional[str],
    microns_per_pixel: Optional[float],
    hires_scale: Optional[float],
) -> tuple[str, float, Optional[float]]:
    spatial_meta = adata.uns.get("spatial")
    if not spatial_meta:
        raise ValueError("`adata.uns['spatial']` is missing scalefactor information.")

    if library_id is None:
        available = list(spatial_meta.keys())
        if len(available) == 0:
            raise ValueError("No library_id found in `adata.uns['spatial']`.")
        library_id = available[0]
        if len(available) > 1:
            warnings.warn(
                f"Multiple library_ids found: {available}. Using '{library_id}'. "
                "Specify `library_id` explicitly to use a different one."
            )
    elif library_id not in spatial_meta:
        raise ValueError(
            f"`library_id` '{library_id}' not found in `adata.uns['spatial']`. "
            f"Available library_ids: {list(spatial_meta.keys())}"
        )

    scalefactors = spatial_meta[library_id].get("scalefactors", {})
    resolved_mpp = microns_per_pixel
    if resolved_mpp is None:
        resolved_mpp = scalefactors.get("microns_per_pixel")
    if resolved_mpp is None or not np.isfinite(resolved_mpp) or resolved_mpp <= 0:
        raise ValueError(
            "A positive `microns_per_pixel` is required. Provide it explicitly "
            "or include `microns_per_pixel` in scalefactors."
        )

    resolved_hires_scale = hires_scale
    if resolved_hires_scale is None:
        resolved_hires_scale = (
            scalefactors.get("tissue_hires_scalef")
            or scalefactors.get("regist_target_img_scalef")
        )
    if resolved_hires_scale is not None and (
        not np.isfinite(resolved_hires_scale) or resolved_hires_scale <= 0
    ):
        raise ValueError("`hires_scale` must be positive when provided.")

    return library_id, float(resolved_mpp), (
        float(resolved_hires_scale) if resolved_hires_scale is not None else None
    )


def _infer_labeldist_coordinate_system(
    coords_all: np.ndarray,
    microns_per_pixel: float,
    hires_scale: Optional[float],
) -> str:
    if hires_scale is None or hires_scale >= 0.99:
        return "fullres"

    coord_range = np.ptp(coords_all[:, :2], axis=0)
    max_range_px = float(np.nanmax(coord_range))
    if not np.isfinite(max_range_px) or max_range_px <= 0:
        warnings.warn(
            "Could not infer spatial coordinate resolution from coordinate range; "
            "using coordinate_system='fullres'."
        )
        return "fullres"

    size_um_fullres = max_range_px * microns_per_pixel
    size_um_hires = max_range_px * (microns_per_pixel / hires_scale)
    expected_max_size_um = 7000

    fullres_reasonable = size_um_fullres <= expected_max_size_um
    hires_reasonable = size_um_hires <= expected_max_size_um

    if hires_reasonable and not fullres_reasonable:
        return "hires"
    if fullres_reasonable and not hires_reasonable:
        return "fullres"

    resolved = "hires" if hires_scale < 0.5 else "fullres"
    warnings.warn(
        "Could not confidently infer whether `adata.obsm['spatial']` is in fullres "
        f"or hires pixels (fullres span={size_um_fullres:.1f} um, "
        f"hires span={size_um_hires:.1f} um). Using coordinate_system='{resolved}'. "
        "Pass `coordinate_system='fullres'` or `coordinate_system='hires'` explicitly "
        "to avoid this heuristic."
    )
    return resolved


def hd_labeldist(
    adata,
    groupby: str,
    label: str,
    inplace: bool = True,
    method: str = "kdtree",
    coordinate_system: str = "auto",
    library_id: Optional[str] = None,
    microns_per_pixel: Optional[float] = None,
    hires_scale: Optional[float] = None,
):
    """
    Compute the distance from every cell to the nearest cell annotated with a specific label (10x HD data).
    
    Distances are reported both in the pixel coordinate system stored in `.obsm["spatial"]`
    (SpaceRanger target/hires layer) and in microns using the scalefactors embedded in
    `adata.uns["spatial"]`.
    
    By default, the function preserves the previous auto-detection heuristic for
    hires vs full-res coordinates. For reproducible distances, pass
    ``coordinate_system='fullres'`` or ``coordinate_system='hires'`` explicitly.
    
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
    coordinate_system : {'auto', 'fullres', 'hires'}, default 'auto'
        Pixel coordinate system used by `adata.obsm['spatial']`. `fullres` uses
        `microns_per_pixel` directly; `hires` divides by `hires_scale`; `auto`
        uses a heuristic and warns when ambiguous.
    library_id : str, optional
        Spatial library id in `adata.uns['spatial']`. Defaults to the first one.
    microns_per_pixel : float, optional
        Override the value stored in scalefactors.
    hires_scale : float, optional
        Override `tissue_hires_scalef` / `regist_target_img_scalef`.
    
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

    coordinate_system = coordinate_system.lower()
    if coordinate_system not in {"auto", "fullres", "hires"}:
        raise ValueError("`coordinate_system` must be 'auto', 'fullres', or 'hires'.")
    
    coords = adata.obsm["spatial"]
    if coords is None or len(coords) == 0:
        raise ValueError("Spatial coordinates are empty.")
    
    mask_label = (adata.obs[groupby] == label).to_numpy()
    if not mask_label.any():
        raise ValueError(f"No observations with label `{label}` were found.")
    
    coords_all = np.asarray(coords, dtype=float)
    if coords_all.ndim != 2 or coords_all.shape[1] < 2:
        raise ValueError("`adata.obsm['spatial']` must be an n_obs x >=2 array.")
    coords_all = coords_all[:, :2]
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
    
    resolved_library_id, resolved_mpp, resolved_hires_scale = _resolve_labeldist_scalefactors(
        adata,
        library_id=library_id,
        microns_per_pixel=microns_per_pixel,
        hires_scale=hires_scale,
    )

    resolved_coordinate_system = coordinate_system
    if coordinate_system == "auto":
        resolved_coordinate_system = _infer_labeldist_coordinate_system(
            coords_all,
            microns_per_pixel=resolved_mpp,
            hires_scale=resolved_hires_scale,
        )
    elif coordinate_system == "hires" and resolved_hires_scale is None:
        raise ValueError(
            "coordinate_system='hires' requires `hires_scale` or a hires scale "
            "in scalefactors."
        )

    if resolved_coordinate_system == "hires":
        dist_um = dist_px * (resolved_mpp / resolved_hires_scale)
    else:
        dist_um = dist_px * resolved_mpp
    
    col_px = f"{label}_px"
    col_um = f"{label}_dist"
    params = {
        "groupby": groupby,
        "label": label,
        "method": method,
        "coordinate_system": resolved_coordinate_system,
        "requested_coordinate_system": coordinate_system,
        "library_id": resolved_library_id,
        "microns_per_pixel": resolved_mpp,
        "hires_scale": resolved_hires_scale,
    }
    
    if inplace:
        adata.obs[col_px] = dist_px
        adata.obs[col_um] = dist_um
        adata.uns[f"{label}_dist_params"] = params
        return None
    
    result = pd.DataFrame({col_px: dist_px, col_um: dist_um}, index=adata.obs.index)
    result.attrs["params"] = params
    return result