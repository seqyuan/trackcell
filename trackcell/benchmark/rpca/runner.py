"""Run RPCA integration benchmark: stepwise R→Python and end-to-end parity."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from trackcell.benchmark.rpca.data import (
    BATCH_KEY,
    DIMS,
    K_ANCHOR,
    K_SCORE,
    K_WEIGHT,
    N_INTEGRATION_FEATURES,
    N_SCT_HVG,
    N_TREES,
    N_TREES_WEIGHT,
    INTEGRATION_DTYPE,
    WEIGHT_QUERY_CHUNK_SIZE,
    RPCA_BENCHMARK_ROOT,
    RPCA_REF_SEED,
    RPCA_SAMPLES,
    RPCA_SEEDS,
    batch_sizes_from_adata,
    inject_r_sct_models,
    load_merged_adata_from_r_filter,
    load_merged_adata_native,
    load_r_anchors,
    load_r_features,
    load_r_integrated,
    load_r_prep,
    matrix_from_obsm,
    r_export_dir,
)
from trackcell.benchmark.rpca.metrics import (
    batch_offsets,
    compare_anchors,
    jaccard,
    matrix_gene_corr,
    r_anchors_to_global,
    summarize_distribution,
)
from trackcell.benchmark.rpca.sct_stepwise import run_integration_sct_stepwise
from trackcell.tl.integration.anchors import (
    find_pairwise_anchors,
    map_anchors_to_global,
    mirror_anchors,
)
from trackcell.tl.integration.integrate_data import stack_integrated
from trackcell.tl.integration.prep_sct import prep_sct_integration
from trackcell.tl.integration.sample_tree import integrate_batches_sample_tree
from trackcell.tl.integration.select_features import select_integration_features
from trackcell.tl.integration._rpca import reciprocal_project, run_batch_pca

def _run_py_anchors(
    prep: pd.DataFrame,
    adata,
    features: list[str],
    *,
    seed: int,
) -> pd.DataFrame:
    batches = adata.obs[BATCH_KEY].astype(str)
    batch_labels = batches.unique().tolist()
    batch_data: dict[str, np.ndarray] = {}
    batch_pca: dict = {}
    cell_offsets: dict[str, int] = {}
    label_to_idx = {label: i for i, label in enumerate(batch_labels)}
    offset = 0

    for label in batch_labels:
        mask = batches == label
        batch_data[label] = prep.loc[mask, features].to_numpy(dtype=np.float64)
        batch_pca[label] = run_batch_pca(batch_data[label], n_components=DIMS, seed=seed)
        cell_offsets[label] = offset
        offset += int(mask.sum())

    anchor_frames: list[pd.DataFrame] = []
    for i, label1 in enumerate(batch_labels):
        for j, label2 in enumerate(batch_labels):
            if i >= j:
                continue
            ref_emb, query_emb = reciprocal_project(
                batch_data[label1],
                batch_pca[label1],
                batch_data[label2],
                batch_pca[label2],
                dims=DIMS,
                l2_norm=True,
            )
            n_ref = batch_data[label1].shape[0]
            anchors_local = find_pairwise_anchors(
                ref_emb,
                query_emb,
                n_ref=n_ref,
                pca1_embedding=batch_pca[label1].embeddings,
                pca2_embedding=batch_pca[label2].embeddings,
                k_anchor=K_ANCHOR,
                k_score=K_SCORE,
                n_trees=N_TREES,
            )
            if anchors_local.empty:
                continue
            anchors_global = map_anchors_to_global(
                anchors_local,
                offset1=cell_offsets[label1],
                offset2=cell_offsets[label2],
                n_local1=n_ref,
                dataset1=label_to_idx[label1],
                dataset2=label_to_idx[label2],
            )
            anchor_frames.append(anchors_global)

    if not anchor_frames:
        return pd.DataFrame(columns=["cell1", "cell2", "score", "dataset1", "dataset2"])
    return mirror_anchors(pd.concat(anchor_frames, ignore_index=True))


def _run_py_integrate(
    prep: pd.DataFrame,
    anchors: pd.DataFrame,
    adata,
    features: list[str],
) -> pd.DataFrame:
    batches = adata.obs[BATCH_KEY].astype(str)
    batch_labels = batches.unique().tolist()
    batch_data: dict[str, np.ndarray] = {}
    global_indices: dict[str, np.ndarray] = {}
    cell_offsets: dict[str, int] = {}
    offset = 0
    for label in batch_labels:
        mask = batches == label
        batch_data[label] = prep.loc[mask, features].to_numpy(dtype=np.float64)
        global_indices[label] = np.where(mask.to_numpy())[0]
        cell_offsets[label] = offset
        offset += int(mask.sum())

    integrated_dict, _ = integrate_batches_sample_tree(
        batch_data,
        batch_labels,
        anchors,
        global_indices,
        k_weight=K_WEIGHT,
        sd_weight=1.0,
        dims=DIMS,
        n_trees=N_TREES,
        n_trees_weight=N_TREES_WEIGHT,
        weight_query_chunk_size=WEIGHT_QUERY_CHUNK_SIZE,
        integration_dtype=INTEGRATION_DTYPE,
    )
    stacked = stack_integrated(integrated_dict, batch_labels, cell_offsets)
    return pd.DataFrame(stacked, index=adata.obs_names, columns=features)


def _run_stepwise(
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
) -> dict[str, Any]:
    adata = load_merged_adata_from_r_filter(seed, benchmark_root=benchmark_root)
    inject_r_sct_models(adata, seed, benchmark_root=benchmark_root)

    r_features = load_r_features(seed, benchmark_root=benchmark_root)
    r_feat_set = set(r_features)
    r_prep = load_r_prep(adata, seed, benchmark_root=benchmark_root)
    r_anchors_raw = load_r_anchors(seed, benchmark_root=benchmark_root)
    r_int = load_r_integrated(seed, benchmark_root=benchmark_root)

    batch_labels = adata.obs[BATCH_KEY].astype(str).unique().tolist()
    batch_sizes = batch_sizes_from_adata(adata)
    offsets = batch_offsets(batch_labels, batch_sizes)

    steps: dict[str, Any] = {}

    py_features = select_integration_features(
        adata, batch_key=BATCH_KEY, n_features=N_INTEGRATION_FEATURES, method="seurat"
    )
    steps["03_features"] = {
        "input": "R Step 02 SCT batch_models",
        "jaccard": jaccard(set(py_features), r_feat_set),
        "r_n": len(r_feat_set),
        "py_n": len(py_features),
    }

    prep_sct_integration(
        adata,
        batch_key=BATCH_KEY,
        anchor_features=r_features,
        sct_key="sct",
        layer="counts",
        key_added="sct_prep",
    )
    py_prep = pd.DataFrame(adata.obsm["X_sct_prep"], index=adata.obs_names, columns=r_features)
    steps["04_prep"] = {
        "input": "R Step 02 SCT + R Step 03 features",
        **matrix_gene_corr(r_prep, py_prep),
    }

    py_anchors = _run_py_anchors(r_prep, adata, r_features, seed=seed)
    steps["05_anchors"] = {
        "input": "R Step 04 prep_residuals",
        **compare_anchors(
            r_anchors_raw,
            py_anchors,
            batch_labels=batch_labels,
            batch_sizes=batch_sizes,
        ),
    }

    r_anchors_global = r_anchors_to_global(
        r_anchors_raw, batch_labels=batch_labels, batch_offsets_map=offsets
    )
    py_int_r = _run_py_integrate(r_prep, r_anchors_global, adata, r_features)
    py_int_py = _run_py_integrate(r_prep, py_anchors, adata, r_features)
    steps["06_integrated_r_anchors"] = {
        "input": "R Step 04 prep + R Step 05 anchors",
        **matrix_gene_corr(r_int, py_int_r),
    }
    steps["06_integrated_py_anchors"] = {
        "input": "R Step 04 prep + Python Step 05 anchors",
        **matrix_gene_corr(r_int, py_int_py),
    }

    return {
        "seed": seed,
        "n_cells": int(adata.n_obs),
        "batch_sizes": batch_sizes,
        "steps": steps,
    }


def _run_e2e_native(
    seed: int,
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
    sct_method: str = "auto",
) -> dict[str, Any]:
    adata = load_merged_adata_native(benchmark_root=benchmark_root)

    from trackcell.tl.integration import integrate_sct_rpca

    out = integrate_sct_rpca(
        adata,
        batch_key=BATCH_KEY,
        layer="counts",
        n_features=N_INTEGRATION_FEATURES,
        dims=DIMS,
        k_anchor=K_ANCHOR,
        k_score=K_SCORE,
        k_weight=K_WEIGHT,
        n_trees=N_TREES,
        n_trees_weight=N_TREES_WEIGHT,
        weight_query_chunk_size=WEIGHT_QUERY_CHUNK_SIZE,
        integration_dtype=INTEGRATION_DTYPE,
        seed=seed,
        sct_seed=seed,
        sct_method=sct_method,
        r_vst_export_root=benchmark_root,
        r_vst_export_seed=seed,
        run_sct=True,
        copy=True,
    )

    r_features = load_r_features(seed, benchmark_root=benchmark_root)
    r_prep = load_r_prep(out, seed, benchmark_root=benchmark_root)
    r_anchors_raw = load_r_anchors(seed, benchmark_root=benchmark_root)
    r_int = load_r_integrated(seed, benchmark_root=benchmark_root)

    py_features = out.uns["sct_integrated"]["anchor_features"]
    py_prep = matrix_from_obsm(out, "sct_prep")
    py_int = matrix_from_obsm(out, "sct_integrated")
    py_anchors = pd.DataFrame(out.uns["sct_integrated"]["anchors"])

    batch_labels = out.obs[BATCH_KEY].astype(str).unique().tolist()
    batch_sizes = batch_sizes_from_adata(out)

    hvg_js: dict[str, float] = {}
    sct_dir = r_export_dir(seed, benchmark_root=benchmark_root) / "02_sct"
    for sample in RPCA_SAMPLES:
        r_hvg = set((sct_dir / f"{sample}_hvg.txt").read_text().strip().split("\n")) - {""}
        py_hvg = set(out.uns["sct"]["batch_models"][sample]["variable_features"])
        hvg_js[sample] = jaccard(py_hvg, r_hvg)

    shared_features = sorted(set(r_features) & set(py_features))

    return {
        "seed": seed,
        "n_cells": int(out.n_obs),
        "sct_hvg_jaccard": hvg_js,
        "features_jaccard": jaccard(set(py_features), set(r_features)),
        "n_shared_features": len(shared_features),
        "prep": matrix_gene_corr(r_prep, py_prep, max_genes=500),
        "prep_shared_features_only": matrix_gene_corr(
            r_prep[shared_features], py_prep[shared_features], max_genes=500
        ),
        "anchors": compare_anchors(
            r_anchors_raw,
            py_anchors,
            batch_labels=batch_labels,
            batch_sizes=batch_sizes,
        ),
        "integrated": matrix_gene_corr(r_int, py_int, max_genes=500),
        "integrated_shared_features_only": matrix_gene_corr(
            r_int[shared_features], py_int[shared_features], max_genes=500
        ),
    }


def run_rpca_benchmark_suite(
    *,
    benchmark_root: Path = RPCA_BENCHMARK_ROOT,
    seeds: tuple[int, ...] = RPCA_SEEDS,
    ref_seed: int = RPCA_REF_SEED,
    skip_e2e: bool = False,
) -> dict[str, Any]:
    stepwise: list[dict[str, Any]] = []
    sct_stepwise: list[dict[str, Any]] = []
    e2e: list[dict[str, Any]] = []
    missing: list[str] = []

    for seed in seeds:
        export_marker = r_export_dir(seed, benchmark_root=benchmark_root) / "pipeline_summary.json"
        if not export_marker.exists():
            missing.append(f"seed_{seed}")
            continue
        stepwise.append(_run_stepwise(seed, benchmark_root=benchmark_root))
        sct_stepwise.append(run_integration_sct_stepwise(seed, benchmark_root=benchmark_root))
        if not skip_e2e:
            e2e.append(_run_e2e_native(seed, benchmark_root=benchmark_root))

    ref_step = next((s for s in stepwise if s["seed"] == ref_seed), None)
    ref_sct = next((s for s in sct_stepwise if s["seed"] == ref_seed), None)
    ref_e2e = next((s for s in e2e if s["seed"] == ref_seed), None)

    summary: dict[str, Any] = {}
    if ref_sct:
        summary["sct_stepwise_ref"] = ref_sct["summary"]
    if ref_step:
        summary["stepwise_ref"] = {
            "features_jaccard": ref_step["steps"]["03_features"]["jaccard"],
            "prep_corr_median": ref_step["steps"]["04_prep"].get("gene_corr_median"),
            "anchor_local_jaccard": ref_step["steps"]["05_anchors"]["local_pair_jaccard"],
            "integrate_corr_median_r_anchors": ref_step["steps"]["06_integrated_r_anchors"].get(
                "gene_corr_median"
            ),
            "integrate_corr_median_py_anchors": ref_step["steps"]["06_integrated_py_anchors"].get(
                "gene_corr_median"
            ),
        }
    if ref_e2e:
        summary["e2e_ref"] = {
            "features_jaccard": ref_e2e["features_jaccard"],
            "n_shared_features": ref_e2e["n_shared_features"],
            "prep_corr_median": ref_e2e["prep"].get("gene_corr_median"),
            "prep_shared_corr_median": ref_e2e["prep_shared_features_only"].get(
                "gene_corr_median"
            ),
            "anchor_local_jaccard": ref_e2e["anchors"]["local_pair_jaccard"],
            "integrate_corr_median": ref_e2e["integrated"].get("gene_corr_median"),
            "integrate_shared_corr_median": ref_e2e["integrated_shared_features_only"].get(
                "gene_corr_median"
            ),
            "sct_hvg_jaccard_median": summarize_distribution(
                list(ref_e2e["sct_hvg_jaccard"].values())
            ).get("median"),
        }

    return {
        "status": "ok" if not missing else "partial",
        "filter": {"min_cells_gene": 100, "min_features": 50},
        "benchmark_root": str(benchmark_root),
        "samples": list(RPCA_SAMPLES),
        "seeds": list(seeds),
        "ref_seed": ref_seed,
        "missing_r_exports": missing,
        "sct_stepwise": sct_stepwise,
        "stepwise": stepwise,
        "e2e": e2e,
        "summary": summary,
    }


def write_benchmark_json(payload: dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path
