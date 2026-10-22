"""Run YardCluster benchmark on synthetic patch data."""

from __future__ import annotations

from typing import Any, Literal

from anndata import AnnData

from trackcell.benchmark.yardcluster.data import make_patch_adata
from trackcell.benchmark.yardcluster.metrics import evaluate_clustering
from trackcell.tl.spatial_cluster import spatial_cluster

ConfigName = Literal["baseline_no_gradient", "defaults_auto", "domain_gradient"]


def _run_config(
    adata: AnnData,
    config: ConfigName,
    *,
    truth_key: str = "patch",
    key_added: str = "yardcluster",
) -> dict[str, Any]:
    ad = adata.copy()
    if config == "baseline_no_gradient":
        spatial_cluster(
            ad,
            mode="domain",
            use_gradient=False,
            merge_clusters=False,
            preprocess=False,
            k_spatial=10,
            n_pcs=15,
            resolution=0.8,
            key_added=key_added,
        )
        pred_key = key_added
    elif config == "defaults_auto":
        spatial_cluster(
            ad,
            mode="auto",
            preprocess=False,
            k_spatial=10,
            n_pcs=15,
            resolution=0.8,
            key_added=key_added,
        )
        pred_key = f"{key_added}_domain"
    elif config == "domain_gradient":
        spatial_cluster(
            ad,
            mode="domain",
            use_gradient=True,
            merge_clusters=True,
            preprocess=False,
            k_spatial=10,
            n_pcs=15,
            resolution=0.8,
            key_added=key_added,
        )
        pred_key = key_added
    else:
        raise ValueError(f"Unknown config: {config}")

    metrics = evaluate_clustering(ad, pred_key, truth_key)
    metrics["config"] = config
    return metrics


def run_yardcluster_benchmark(
    *,
    seed: int = 0,
    n_per_patch: int = 40,
    n_patches: int = 3,
    configs: tuple[ConfigName, ...] = (
        "baseline_no_gradient",
        "defaults_auto",
        "domain_gradient",
    ),
) -> dict[str, Any]:
    """
    Compare YardCluster configurations on synthetic spatial patch data.

    Returns per-config ARI and spatial coherence vs known patch labels.
    """
    adata = make_patch_adata(
        n_per_patch=n_per_patch,
        n_patches=n_patches,
        seed=seed,
    )
    results = [_run_config(adata, cfg) for cfg in configs]
    best = max(results, key=lambda r: (r["ari"], r["spatial_coherence"]))
    return {
        "seed": seed,
        "n_obs": int(adata.n_obs),
        "n_patches": n_patches,
        "configs": results,
        "best_config": best["config"],
        "best_ari": best["ari"],
        "best_spatial_coherence": best["spatial_coherence"],
    }


def format_benchmark_table(summary: dict[str, Any]) -> str:
    """Render a markdown table from :func:`run_yardcluster_benchmark` output."""
    lines = [
        "| Config | ARI | Spatial coherence | Clusters |",
        "|--------|-----|-------------------|----------|",
    ]
    for row in summary["configs"]:
        lines.append(
            f"| {row['config']} | {row['ari']:.3f} | "
            f"{row['spatial_coherence']:.3f} | {row['n_clusters']} |"
        )
    return "\n".join(lines)
