"""Tests for CellCharter-inspired YardCluster features."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scipy import sparse as sps

from trackcell.tl._spatial_graph import (
    _row_normalize,
    build_spatial_weights,
    compute_neighborhood_mean,
    compute_neighborhood_multi_hop,
    remove_long_links,
)
from trackcell.tl._nhood import (
    connected_components,
    diff_nhood_enrichment,
    nhood_enrichment,
)
from trackcell.tl._autores import select_resolution, sweep_resolutions
from trackcell.tl._cluster import (
    neighborhood_features,
    spatial_cluster,
    spatial_neighbors,
    yard_cluster,
    yard_embed,
)


def _make_two_blob_adata(
    n_per_blob: int = 25, n_genes: int = 40, seed: int = 0
) -> AnnData:
    """Two well-separated blobs -> two spatial components."""
    rng = np.random.default_rng(seed)
    coords = np.vstack(
        [
            rng.normal(loc=[0, 0], scale=0.3, size=(n_per_blob, 2)),
            rng.normal(loc=[10, 10], scale=0.3, size=(n_per_blob, 2)),
        ]
    )
    X = rng.normal(size=(2 * n_per_blob, n_genes)).astype(np.float32)
    X[:n_per_blob, :5] += 3.0
    X[n_per_blob:, 5:10] += 3.0
    adata = AnnData(X=X)
    adata.obsm["spatial"] = coords
    adata.obs["blob"] = np.array(["A"] * n_per_blob + ["B"] * n_per_blob)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_multi_hop_shapes_and_dimension():
    adata = _make_two_blob_adata(n_per_blob=15, n_genes=20)
    coords = adata.obsm["spatial"]
    w = build_spatial_weights(coords, k=5)
    X = np.asarray(adata.X, dtype=np.float64)
    m1 = compute_neighborhood_multi_hop(X, w, n_layers=1)
    m3 = compute_neighborhood_multi_hop(X, w, n_layers=3)
    assert m1.shape == (adata.n_obs, 20)
    assert m3.shape == (adata.n_obs, 60)
    np.testing.assert_allclose(m3[:, :20], m1, rtol=1e-5, atol=1e-5)


def test_neighborhood_features_n_layers():
    adata = _make_two_blob_adata(n_per_blob=12, n_genes=25)
    neighborhood_features(adata, k_spatial=5, n_layers=2, key_added="yard")
    assert adata.obsm["yard_M"].shape == (adata.n_obs, 50)
    # params recorded for cache invalidation
    assert adata.uns["yard_features_params"]["n_layers"] == 2


def test_remove_long_links_trims_edges():
    coords = np.array([[0, 0], [0.1, 0], [0, 0.1], [50, 50]], dtype=float)
    w = build_spatial_weights(coords, k=3)
    n_edges_before = (w > 0).nnz
    trimmed = remove_long_links(coords, w, distance_percentile=60.0)
    assert (trimmed > 0).nnz < n_edges_before
    # rows still sum to 1 (where any edge remains)
    row_sums = np.asarray(trimmed.sum(axis=1)).ravel()
    assert np.allclose(row_sums[row_sums > 0], 1.0, rtol=1e-5)


def test_connected_components_splits_cluster():
    adata = _make_two_blob_adata(n_per_blob=20, n_genes=20)
    spatial_neighbors(adata, k=5)
    # assign both blobs the same cluster label -> components should split them
    adata.obs["cluster"] = "0"
    connected_components(adata, cluster_key="cluster", min_cells=3, out_key="component")
    comp = adata.obs["component"].astype(str)
    assert comp.nunique() == 2
    assert set(comp.unique()) == {"0", "1"}


def test_connected_components_whole_graph():
    adata = _make_two_blob_adata(n_per_blob=20, n_genes=20)
    spatial_neighbors(adata, k=5)
    connected_components(adata, cluster_key=None, min_cells=3, out_key="component")
    assert adata.obs["component"].astype(str).nunique() == 2


def test_nhood_enrichment_symmetric_and_inter():
    adata = _make_two_blob_adata(n_per_blob=30, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    res = nhood_enrichment(
        adata, "yardcluster", connectivity_key="spatial_connectivities", copy=True
    )
    assert set(res.keys()) == {"enrichment", "observed", "expected"}
    # diagonal should be NaN when only_inter=True (default)
    assert np.isnan(res["enrichment"].values.diagonal()).all()
    # matrix is symmetric
    np.testing.assert_allclose(
        res["enrichment"].values,
        res["enrichment"].values.T,
        rtol=1e-9,
        equal_nan=True,
    )


def test_nhood_enrichment_writes_to_uns():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=20)
    spatial_neighbors(adata, k=6)
    adata.obs["cluster"] = ["0", "1"] * 25
    nhood_enrichment(adata, "cluster")
    assert "cluster_nhood_enrichment" in adata.uns
    assert "params" in adata.uns["cluster_nhood_enrichment"]


def test_diff_nhood_enrichment_with_pvalues():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=40, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    # half cells Case, half Ctr (split each blob)
    adata.obs["condition"] = ["Case"] * 40 + ["Ctr"] * 40
    adata.obs["library"] = ["s1"] * 20 + ["s2"] * 20 + ["s3"] * 20 + ["s4"] * 20
    res = diff_nhood_enrichment(
        adata,
        "yardcluster",
        "condition",
        ("Case", "Ctr"),
        library_key="library",
        pvalues=True,
        n_perms=20,
        copy=True,
    )
    assert "Case_Ctr" in res
    assert "enrichment" in res["Case_Ctr"]
    assert "pvalue" in res["Case_Ctr"]
    pv = res["Case_Ctr"]["pvalue"].values
    assert ((pv >= 0) & (pv <= 1) | np.isnan(pv)).all()
    # Permutation p-values use the (count+1)/(n+1) estimator, so 0 is not
    # attainable — otherwise a downstream "p < 0.001" filter is meaningless.
    assert (pv > 0).all()
    # Multiple-testing correction is provided alongside the raw matrix.
    assert "pvalue_adj" in res["Case_Ctr"]
    adj = res["Case_Ctr"]["pvalue_adj"].values
    assert (adj >= pv - 1e-12).all()
    assert (adj <= 1.0 + 1e-12).all()


def test_bh_adjust_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests

    from trackcell.tl._nhood import _bh_adjust

    p = pd.DataFrame([[0.001, 0.01], [0.02, 0.5]], index=["a", "b"], columns=["a", "b"])
    adj = _bh_adjust(p)
    _, expected, _, _ = multipletests(p.to_numpy().ravel(), method="fdr_bh")
    np.testing.assert_allclose(adj.to_numpy().ravel(), expected)


def test_bh_adjust_keeps_nan_layout():
    from trackcell.tl._nhood import _bh_adjust

    p = pd.DataFrame([[0.01, np.nan], [0.2, 0.4]])
    adj = _bh_adjust(p)
    assert np.isnan(adj.iloc[0, 1])
    assert np.isfinite(adj.iloc[0, 0])


def test_sweep_and_select_resolution():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    sweep = sweep_resolutions(
        adata, use_rep="X_yard_mixed", resolutions=[0.3, 0.5, 0.8], n_runs=3
    )
    assert list(sweep.columns) == ["resolution", "n_clusters", "stability"]
    assert len(sweep) == 3
    assert (sweep["stability"] >= 0).all()

    best, sweep2 = select_resolution(
        adata, use_rep="X_yard_mixed", resolutions=[0.3, 0.5, 0.8], n_runs=3
    )
    assert best in (0.3, 0.5, 0.8)


def test_yard_cluster_auto_resolution():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    yard_cluster(adata, resolution="auto", n_resolution_runs=3, key_added="yardcluster")
    assert "yardcluster" in adata.obs
    assert adata.obs["yardcluster"].nunique() >= 2
    assert "yardcluster_resolution_sweep" in adata.uns
    assert adata.uns["yardcluster_cluster_params"]["resolution"] != "auto"


def test_yard_cluster_gmm():
    pytest.importorskip("sklearn")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=10)
    yard_cluster(adata, method="gmm", n_clusters=2, key_added="yardcluster")
    assert "yardcluster" in adata.obs
    assert adata.uns["yardcluster_cluster_params"]["method"] == "gmm"
    assert adata.obs["yardcluster"].nunique() == 2


def test_spatial_cluster_auto_resolution_smoke():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        preprocess=False,
        k_spatial=6,
        n_pcs=10,
        resolution="auto",
        n_resolution_runs=2,
        use_gradient=False,
    )
    assert "yardcluster" in adata.obs


def test_spatial_cluster_n_layers_smoke():
    pytest.importorskip("leidenalg")
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        preprocess=False,
        k_spatial=6,
        n_pcs=10,
        n_layers=2,
        use_gradient=False,
        resolution=0.5,
    )
    assert "yardcluster" in adata.obs


def test_spatial_cluster_gmm_single_mode():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="domain",
        method="gmm",
        n_clusters=2,
        preprocess=False,
        k_spatial=6,
        n_pcs=8,
        use_gradient=False,
    )
    assert "yardcluster" in adata.obs
    assert adata.uns["yardcluster_cluster_params"]["method"] == "gmm"
    assert adata.obs["yardcluster"].nunique() == 2


def test_spatial_cluster_gmm_auto_mode():
    adata = _make_two_blob_adata(n_per_blob=25, n_genes=30)
    spatial_cluster(
        adata,
        mode="auto",
        method="gmm",
        n_clusters=None,
        merge_clusters=False,
        preprocess=False,
        k_spatial=6,
        n_pcs=8,
        use_gradient=False,
    )
    assert "yardcluster_celltype" in adata.obs
    assert "yardcluster_domain" in adata.obs
    assert adata.uns["yardcluster_params"]["method"] == "gmm"


def test_pl_nhood_enrichment_heatmap():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from trackcell.pl import diff_nhood_enrichment as pl_diff
    from trackcell.pl import nhood_enrichment as pl_nhood

    adata = _make_two_blob_adata(n_per_blob=25, n_genes=20)
    neighborhood_features(adata, k_spatial=6)
    yard_embed(adata, lam=0.8, n_pcs=8)
    yard_cluster(adata, resolution=1.0, key_added="yardcluster")
    nhood_enrichment(adata, "yardcluster")
    ax = pl_nhood(adata, "yardcluster", show=False)
    assert ax is not None
    assert len(ax.images) == 1
    plt.close("all")

    adata.obs["condition"] = ["Case"] * 25 + ["Ctr"] * 25
    adata.obs["library"] = ["s1"] * 13 + ["s2"] * 12 + ["s3"] * 13 + ["s4"] * 12
    diff_nhood_enrichment(
        adata,
        "yardcluster",
        "condition",
        ("Case", "Ctr"),
        library_key="library",
        pvalues=True,
        n_perms=10,
    )
    ax2 = pl_diff(adata, "yardcluster", "condition", ("Case", "Ctr"), show=False)
    assert ax2 is not None
    plt.close("all")


def test_gene_alias_lookup_exact_names_win_over_aliases():
    """A real gene named ``A-1`` must not resolve to the alias of ``A.1``.

    A single ``setdefault`` pass let whichever gene came first claim the alias,
    so ``align_gene_list(['A-1'])`` returned ``A.1`` — a different gene, silently.
    """
    from trackcell.tl._gene_utils import (
        align_gene_list,
        build_gene_alias_lookup,
        resolve_gene_alias,
    )

    targets = pd.Index(["A.1", "A-1", "B", "C.2"])

    lookup = build_gene_alias_lookup(targets)
    assert lookup["A.1"] == "A.1"
    assert lookup["A-1"] == "A-1"

    assert resolve_gene_alias("A-1", targets) == "A-1"
    assert resolve_gene_alias("A.1", targets) == "A.1"
    assert align_gene_list(["A-1"], targets) == ["A-1"]

    # When only one spelling exists, the alias still resolves to it.
    only_dotted = pd.Index(["A.1", "B"])
    assert resolve_gene_alias("A-1", only_dotted) == "A.1"


def test_align_gene_list_deduplicates_matches():
    from trackcell.tl._gene_utils import align_gene_list

    targets = pd.Index(["A.1", "B"])
    assert align_gene_list(["A.1", "A-1", "A.1"], targets) == ["A.1"]


def test_align_gene_frame_matches_list_variant_and_keeps_dtypes():
    """The frame variant must behave like the list variant and stay cheap.

    It deduplicated through ``matched not in index`` on a *list*, an O(n^2) scan
    that cost 2.4 s on a 20k-gene model frame (147x the list version), and it
    rebuilt the frame from ``Series`` rows, which can upcast column dtypes.
    """
    from trackcell.tl._gene_utils import align_gene_frame

    targets = pd.Index(["A.1", "B", "G0", "G1"])
    frame = pd.DataFrame(
        {"theta": [1.0, 2.0, 3.0], "name": ["x", "y", "z"]},
        index=["A.1", "A-1", "ZZZ"],
    )
    out = align_gene_frame(frame, targets)

    # Alias resolves to the canonical gene, duplicate rows collapse to one,
    # unmatched genes are dropped.
    assert list(out.index) == ["A.1"]
    assert out["theta"].tolist() == [1.0]
    # dtypes survive the realignment
    assert out["theta"].dtype == frame["theta"].dtype
    assert out["name"].dtype == frame["name"].dtype


def test_align_gene_frame_empty_and_paren_columns():
    from trackcell.tl._gene_utils import align_gene_frame

    targets = pd.Index(["G0", "G1"])
    # No gene matches: still returns a frame with the original columns.
    empty = align_gene_frame(pd.DataFrame({"v": [1.0]}, index=["ZZZ"]), targets)
    assert empty.shape == (0, 1)
    assert list(empty.columns) == ["v"]

    # ``normalize_model_pars`` still strips the R intercept parentheses first.
    intercept = align_gene_frame(
        pd.DataFrame({"(Intercept)": [1.0]}, index=["G0"]), targets
    )
    assert list(intercept.columns) == ["Intercept"]
    assert list(intercept.index) == ["G0"]


def test_remove_long_links_trims_the_only_edge_of_a_degree_one_cell():
    """A degree-1 cell normalises its sole edge to weight 1.0.

    The old self-loop test was ``data >= 1 - 1e-12``, which mistook that genuine
    edge for a self loop and exempted it from trimming — so the edge most in need
    of removal, the one long enough to bridge two regions, always survived.
    ``build_spatial_weights`` emits no self loops at all, so the weight test was
    never detecting what it claimed.
    """
    coords = np.array([[0.0, 0.0], [100.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    adjacency = np.array(
        [
            [0.0, 1.0, 1.0, 1.0],
            [1.0, 0.0, 0.0, 0.0],  # degree 1, neighbour 100 units away
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0],
        ]
    )
    w = _row_normalize(sps.csr_matrix(adjacency))
    assert w[1, 0] == pytest.approx(1.0)  # the weight that fooled the old test

    trimmed = remove_long_links(coords, w, distance_percentile=50.0)

    assert trimmed[1, 0] == 0.0, "the 100-unit edge from the degree-1 cell survived"


def test_remove_long_links_keeps_explicit_self_loops():
    """Explicit diagonal entries are still exempt — now by position, not weight."""
    coords = np.array([[0.0, 0.0], [0.1, 0.0], [50.0, 0.0]])
    w = sps.csr_matrix(
        np.array(
            [
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 1.0],
            ]
        )
    )

    trimmed = remove_long_links(coords, w, distance_percentile=50.0).toarray()

    # Self loops are never trimmed, and a self loop's zero distance must not
    # enter the percentile sample either.
    assert trimmed[0, 0] > 0.0
    assert trimmed[1, 1] > 0.0
    assert trimmed[2, 2] > 0.0


def _leaf_graph():
    """4 cells: 0/2/3 clustered, cell 1 alone 100 units away (only edge: 1->0)."""
    coords = np.array([[0.0, 0.0], [100.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    adjacency = np.array(
        [
            [0.0, 1.0, 1.0, 1.0],
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0],
        ]
    )
    return coords, _row_normalize(sps.csr_matrix(adjacency))


def test_remove_long_links_leaves_no_zero_row_sums():
    """A dangling row would feed zeros into ``weights @ X``."""
    coords, w = _leaf_graph()
    trimmed = remove_long_links(coords, w, distance_percentile=50.0).toarray()

    assert np.allclose(trimmed.sum(axis=1), 1.0)


def test_trimmed_away_cell_context_becomes_its_own_expression():
    """Isolating a cell must not silently zero its neighbourhood features.

    ``build_spatial_weights`` already gives a neighbour-less cell ``(i, i, 1.0)``
    — context = own expression. Trimming used to break that invariant: the row
    normalised to nothing, so ``compute_neighborhood_mean`` reported zeros and
    downstream PCA read the cell as "expresses nothing".
    """
    coords, w = _leaf_graph()
    trimmed = remove_long_links(coords, w, distance_percentile=50.0)

    rng = np.random.default_rng(0)
    X = rng.random((4, 5))
    M = compute_neighborhood_mean(X, trimmed)

    assert np.allclose(M[1], X[1]), "isolated cell's context is not its own expression"
    assert np.all(np.abs(M) > 0), "some cell still has an all-zero context"


def test_remove_long_links_adds_no_self_loops_when_nothing_is_isolated():
    """The self loop is a repair, not a new global convention."""
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    w = build_spatial_weights(coords, k=3)
    trimmed = remove_long_links(coords, w, distance_percentile=99.9).toarray()

    assert np.allclose(np.diag(trimmed), 0.0)
    assert np.allclose(trimmed.sum(axis=1), 1.0)


def test_isolated_cell_fallback_matches_the_graph_builder_convention():
    """Same fallback as ``build_spatial_weights`` gives a single-cell batch."""
    from trackcell.tl._spatial_graph import build_spatial_weights

    # A batch with one cell: the builder's own "no neighbour" branch.
    solo = build_spatial_weights(np.array([[0.0, 0.0]]), k=5).toarray()
    assert solo[0, 0] == 1.0

    coords, w = _leaf_graph()
    trimmed = remove_long_links(coords, w, distance_percentile=50.0).toarray()
    assert trimmed[1, 1] == solo[0, 0]


def test_isolated_cell_features_are_finite_in_every_mode():
    """A self loop must not leak NaN into any neighbourhood feature.

    The self loop gives the cell a neighbour at offset (0, 0), which is exactly
    the kind of degenerate input a directional weighting can turn into NaN
    (``arctan2(0, 0)`` -> 0, so ``cos(phi)`` -> 1; the regression variant skips
    rows with fewer than three neighbours). All three consumers are checked
    because only checking the mean would leave two untested.
    """
    from trackcell.tl._spatial_graph import (
        compute_neighborhood_gradient_cosphi,
        compute_neighborhood_gradient_regression,
    )

    coords, w = _leaf_graph()
    trimmed = remove_long_links(coords, w, distance_percentile=50.0)

    rng = np.random.default_rng(1)
    X = rng.random((4, 6))

    mean = compute_neighborhood_mean(X, trimmed)
    cosphi = compute_neighborhood_gradient_cosphi(X, coords, trimmed)
    regression = compute_neighborhood_gradient_regression(X, coords, trimmed)

    for name, matrix in (
        ("mean", mean),
        ("cosphi", cosphi),
        ("regression", regression),
    ):
        assert np.all(np.isfinite(matrix)), f"{name} produced non-finite values"

    # Multi-hop diffusion must also stay finite: the self loop persists under W^k.
    multi = compute_neighborhood_multi_hop(X, trimmed, n_layers=3)
    assert np.all(np.isfinite(multi))
    assert np.allclose(multi[:, :6], mean, rtol=1e-5, atol=1e-5)
