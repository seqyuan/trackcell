"""
Automatic Leiden resolution selection via clustering stability.

Adapted from CellCharter's ``ClusterAutoK`` idea (repeated fits + stability
curve), translated to Leiden's resolution parameter: sweep a resolution grid,
run Leiden several times with different seeds, and score each resolution by the
mean pairwise ARI of its repeated runs.
"""

from __future__ import annotations

from itertools import combinations
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score


def _leiden_with_seed(
    adata: AnnData,
    resolution: float,
    neighbors_key: str,
    key_added: str,
    random_state: int,
) -> None:
    kw = dict(
        resolution=float(resolution),
        key_added=key_added,
        neighbors_key=neighbors_key,
        random_state=int(random_state),
    )
    try:
        sc.tl.leiden(
            adata,
            flavor="igraph",
            directed=False,
            n_iterations=2,
            **kw,
        )
    except (ImportError, TypeError, ValueError):
        sc.tl.leiden(adata, flavor="leidenalg", **kw)


def sweep_resolutions(
    adata: AnnData,
    *,
    use_rep: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    resolutions: Optional[Sequence[float]] = None,
    n_runs: int = 5,
    n_neighbors: int = 15,
    random_state: int = 0,
    copy: bool = False,
) -> pd.DataFrame:
    """
    Sweep Leiden resolutions and score each by clustering stability.

    For each resolution the pipeline runs Leiden ``n_runs`` times with
    different seeds and reports the mean pairwise ARI (stability) plus the
    median number of clusters.

    Parameters
    ----------
    adata
        Annotated data object.
    use_rep
        Key in ``adata.obsm`` used to build the neighbor graph when
        ``neighbors_key`` is not provided.
    neighbors_key
        Key of an existing neighbor graph (``adata.uns[neighbors_key]``). If
        ``None``, ``sc.pp.neighbors`` is run on ``use_rep`` first.
    resolutions
        Resolution grid. Default ``np.arange(0.2, 2.01, 0.2)``.
    n_runs
        Number of repeated Leiden runs per resolution.
    n_neighbors
        kNN size used when building a new neighbor graph.
    random_state
        Base seed; run ``i`` uses ``random_state + i``.
    copy
        If ``True``, run on a copy of ``adata``.

    Returns
    -------
    :class:`pandas.DataFrame` with columns ``resolution``, ``n_clusters``,
    ``stability``.
    """
    adata = adata.copy() if copy else adata

    if resolutions is None:
        resolutions = np.round(np.arange(0.2, 2.01, 0.2), 2).tolist()
    resolutions = [float(r) for r in resolutions]

    if neighbors_key is None:
        if use_rep is None:
            raise ValueError("use_rep is required when neighbors_key is None.")
        # Use a non-reserved key: scanpy's NeighborsView treats "neighbors"
        # specially and would otherwise look for unprefixed "connectivities".
        neighbors_key = "yard_autores"
        sc.pp.neighbors(
            adata, use_rep=use_rep, n_neighbors=n_neighbors, key_added=neighbors_key
        )
    conn_key = f"{neighbors_key}_connectivities"
    if conn_key not in adata.obsp:
        raise ValueError(
            f"`adata.obsp['{conn_key}']` not found for neighbors_key='{neighbors_key}'."
        )

    rows: list[dict] = []
    temp_keys: list[str] = []
    for res in resolutions:
        labels: list[np.ndarray] = []
        for run in range(n_runs):
            key = f"__yard_autores_{res}_{run}"
            temp_keys.append(key)
            _leiden_with_seed(adata, res, neighbors_key, key, random_state + run)
            labels.append(adata.obs[key].astype(str).to_numpy())

        n_clusters = int(np.median([len(np.unique(l)) for l in labels]))
        if n_runs > 1:
            stability = float(
                np.mean(
                    [
                        adjusted_rand_score(labels[i], labels[j])
                        for i, j in combinations(range(n_runs), 2)
                    ]
                )
            )
        else:
            stability = float("nan")
        rows.append(
            {"resolution": res, "n_clusters": n_clusters, "stability": stability}
        )

    for key in temp_keys:
        if key in adata.obs.columns:
            del adata.obs[key]

    return pd.DataFrame(rows)


def select_resolution(
    adata: AnnData,
    *,
    use_rep: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    resolutions: Optional[Sequence[float]] = None,
    n_runs: int = 5,
    n_neighbors: int = 15,
    random_state: int = 0,
    copy: bool = False,
) -> tuple[float, pd.DataFrame]:
    """
    Pick the Leiden resolution with the highest stability.

    Returns
    -------
    best_resolution
        Resolution maximizing stability among resolutions producing >1 cluster.
    sweep
        The full sweep :class:`pandas.DataFrame`.
    """
    sweep = sweep_resolutions(
        adata,
        use_rep=use_rep,
        neighbors_key=neighbors_key,
        resolutions=resolutions,
        n_runs=n_runs,
        n_neighbors=n_neighbors,
        random_state=random_state,
        copy=copy,
    )
    candidates = sweep[sweep["n_clusters"] > 1]
    if candidates.empty or candidates["stability"].isna().all():
        best = float(sweep["resolution"].iloc[-1])
    else:
        best = float(candidates.loc[candidates["stability"].idxmax(), "resolution"])
    return best, sweep
