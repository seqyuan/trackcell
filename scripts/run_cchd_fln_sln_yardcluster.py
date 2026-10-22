#!/usr/bin/env python3
"""FLN + SLN (32 Xenium TMA slices): SCT+RPCA then YardCluster.

Mirrors the sample layout in ``02_scvi_FLN_SLN/run_scvi.py``:

* 85863_F  — 16 FLN slices (S022–S037)
* 85811_S  — 11 SLN slices (S011–S021)
* 85863_S  —  5 SLN slices (S038–S042)

Step 1 — ``integrate_sct_rpca`` with ``batch_key='slice_id'`` (32 batches).
Step 2 — ``spatial_cluster`` with ``expression_source='sct_integrated'``.

Example::

    # Full pipeline (needs large RAM; ~4.8M cells)
    python scripts/run_cchd_fln_sln_yardcluster.py --step all

    # Resume YardCluster from saved SCT+RPCA checkpoint
    python scripts/run_cchd_fln_sln_yardcluster.py --step yardcluster \\
        --input-h5ad /path/to/fln_sln_sct_rpca.h5ad
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import trackcell as tcl  # noqa: E402

DEFAULT_DATA = Path("/Volumes/data/co_project/CCHD/xenium_analysis/01_data")
DEFAULT_WORK = Path("/Volumes/data/co_project/CCHD/xenium_analysis/03_yardcluster_fln_sln")

SAMPLE_CONFIG = [
    ("85863_F", "FLN", "F"),
    ("85811_S", "SLN1", "LN"),
    ("85863_S", "SLN2", "LN"),
]

MIN_COUNTS = 30
MIN_CELLS = 10


def _load_fln_sln(data_dir: Path) -> ad.AnnData:
    """Load 3 samples, QC, intersect genes, concat (~4.8M cells)."""
    adatas: dict[str, ad.AnnData] = {}
    common_genes: set[str] | None = None

    for fname, sample_label, tissue in SAMPLE_CONFIG:
        path = data_dir / f"{fname}.h5ad"
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"  Loading {path.name} ...", end=" ", flush=True)
        t0 = time.time()
        adata = sc.read_h5ad(path)
        n_pre = adata.n_obs
        if "slice_id" in adata.obs:
            adata = adata[adata.obs["slice_id"].astype(str) != "debris"].copy()
            print(f"{n_pre:,} → {adata.n_obs:,} cells (debris removed), {time.time()-t0:.0f}s")
        else:
            print(f"{adata.n_obs:,} cells, {time.time()-t0:.0f}s")

        adata.obs["sample"] = sample_label
        adata.obs["tissue_type"] = tissue
        adata.obs["sample_id"] = fname
        adata.obs["group"] = "FLN" if tissue == "F" else "SLN"

        sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
        sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

        if "x_centroid" in adata.obs and "y_centroid" in adata.obs:
            adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

        adatas[sample_label] = adata
        genes = set(adata.var_names)
        common_genes = genes if common_genes is None else (common_genes & genes)

    assert common_genes is not None
    common_genes_sorted = sorted(common_genes)
    print(f"  Common genes: {len(common_genes_sorted):,}")

    for label in adatas:
        adatas[label] = adatas[label][:, common_genes_sorted].copy()

    keys = list(adatas.keys())
    merged = ad.concat(
        [adatas[k] for k in keys],
        join="inner",
        label="sample_batch",
        keys=keys,
    )
    merged.obs_names_make_unique()
    for label in adatas:
        mask = merged.obs["sample_batch"] == label
        for col in ("sample", "tissue_type", "sample_id", "group"):
            merged.obs.loc[mask, col] = adatas[label].obs[col].values

    if "counts" not in merged.layers:
        merged.layers["counts"] = merged.X.copy()

    slices = sorted(merged.obs["slice_id"].astype(str).unique())
    print(f"  Combined: {merged.n_obs:,} cells × {merged.n_vars:,} genes")
    print(f"  Slices (batch_key): {len(slices)} — {slices[0]} … {slices[-1]}")
    print(f"  FLN: {(merged.obs['group'] == 'FLN').sum():,}  SLN: {(merged.obs['group'] == 'SLN').sum():,}")
    del adatas
    gc.collect()
    return merged


def step_sct_rpca(
    adata: ad.AnnData,
    *,
    batch_key: str,
    n_features: int,
    dims: int,
    sct_batch_n_jobs: int,
    rpca_method: str,
    release_cache: bool,
) -> ad.AnnData:
    print("\n=== Step 1: SCT + RPCA ===")
    print(f"  batch_key={batch_key!r}, n_slices={adata.obs[batch_key].nunique()}")
    t0 = time.time()
    tcl.tl.integrate_sct_rpca(
        adata,
        batch_key=batch_key,
        layer="counts",
        n_features=n_features,
        dims=dims,
        sct_batch_n_jobs=sct_batch_n_jobs,
        method=rpca_method,
        n_sct_hvg=3000,
    )
    print(f"  integrate_sct_rpca done in {time.time() - t0:.0f}s")
    print(f"  X_sct_prep: {adata.obsm['X_sct_prep'].shape}")
    print(f"  X_sct_integrated: {adata.obsm['X_sct_integrated'].shape}")
    print(f"  anchors: {adata.uns['sct_integrated']['n_anchors']}")

    if release_cache:
        summary = tcl.tl.release_sct_integration_cache(adata, what="scale_data")
        print(f"  released SCT scale_data cache (~{summary.get('bytes_estimate', 0) / 1e9:.1f} GB est.)")
        gc.collect()
    return adata


def step_yardcluster(
    adata: ad.AnnData,
    *,
    batch_key: str,
    sketch_n: int,
    hybrid_alpha: float,
    resolution_celltype: float,
    resolution_domain: float,
) -> ad.AnnData:
    print("\n=== Step 2: YardCluster ===")
    if "X_sct_integrated" not in adata.obsm:
        raise ValueError("Missing X_sct_integrated — run step sct_rpca first.")

    sketch = adata.n_obs > sketch_n
    print(
        f"  platform=xenium_tma, expression_source=sct_integrated, "
        f"sketch={sketch}, sketch_n={sketch_n if sketch else adata.n_obs}"
    )
    t0 = time.time()
    tcl.tl.spatial_cluster(
        adata,
        platform="xenium_tma",
        batch_key=batch_key,
        expression_source="sct_integrated",
        integrate="none",
        graph_mode="hybrid",
        hybrid_alpha=hybrid_alpha,
        mode="auto",
        sketch=sketch,
        sketch_n=sketch_n,
        resolution_celltype=resolution_celltype,
        resolution_domain=resolution_domain,
        key_added="yardcluster",
    )
    print(f"  spatial_cluster done in {time.time() - t0:.0f}s")

    adata.obs["leiden_center_yard"] = adata.obs["yardcluster_celltype"]
    adata.obs["leiden_nbr_yard"] = adata.obs["yardcluster_domain"]
    print(
        f"  yardcluster_celltype: {adata.obs['yardcluster_celltype'].nunique()} clusters"
    )
    print(f"  yardcluster_domain: {adata.obs['yardcluster_domain'].nunique()} clusters")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description="CCHD FLN+SLN SCT+RPCA → YardCluster")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    parser.add_argument(
        "--step",
        choices=("load", "sct_rpca", "yardcluster", "all"),
        default="all",
        help="Pipeline step (default: all)",
    )
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=None,
        help="Input for yardcluster step (default: work-dir/fln_sln_sct_rpca.h5ad)",
    )
    parser.add_argument("--batch-key", default="slice_id")
    parser.add_argument("--n-features", type=int, default=2000)
    parser.add_argument("--rpca-dims", type=int, default=30)
    parser.add_argument("--sct-jobs", type=int, default=-1)
    parser.add_argument("--rpca-method", choices=("seurat", "xgboost"), default="seurat")
    parser.add_argument("--no-release-cache", action="store_true")
    parser.add_argument("--sketch-n", type=int, default=500_000)
    parser.add_argument("--hybrid-alpha", type=float, default=0.25)
    parser.add_argument("--resolution-celltype", type=float, default=1.0)
    parser.add_argument("--resolution-domain", type=float, default=0.8)
    args = parser.parse_args()

    work = args.work_dir.expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    rpca_path = work / "fln_sln_sct_rpca.h5ad"
    final_path = work / "fln_sln_yardcluster.h5ad"

    if args.step in ("load", "all"):
        print("Loading FLN + SLN samples ...")
        adata = _load_fln_sln(args.data_dir.expanduser().resolve())
        if args.step == "load":
            out = work / "fln_sln_loaded.h5ad"
            print(f"Writing {out} ...")
            adata.write_h5ad(out)
            return

    if args.step in ("sct_rpca", "all"):
        if args.step == "all":
            pass  # adata from load
        else:
            inp = args.input_h5ad or (work / "fln_sln_loaded.h5ad")
            print(f"Loading {inp} ...")
            adata = sc.read_h5ad(inp.expanduser().resolve())
            if "spatial" not in adata.obsm and "x_centroid" in adata.obs:
                adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
            if "counts" not in adata.layers:
                adata.layers["counts"] = adata.X.copy()

        adata = step_sct_rpca(
            adata,
            batch_key=args.batch_key,
            n_features=args.n_features,
            dims=args.rpca_dims,
            sct_batch_n_jobs=args.sct_jobs,
            rpca_method=args.rpca_method,
            release_cache=not args.no_release_cache,
        )
        print(f"Writing checkpoint {rpca_path} ...")
        adata.write_h5ad(rpca_path)
        if args.step == "sct_rpca":
            print("Step 1 complete.")
            return

    if args.step in ("yardcluster", "all"):
        if args.step == "yardcluster":
            inp = args.input_h5ad or rpca_path
            print(f"Loading {inp} ...")
            adata = sc.read_h5ad(inp.expanduser().resolve())
            if "spatial" not in adata.obsm and "x_centroid" in adata.obs:
                adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

        adata = step_yardcluster(
            adata,
            batch_key=args.batch_key,
            sketch_n=args.sketch_n,
            hybrid_alpha=args.hybrid_alpha,
            resolution_celltype=args.resolution_celltype,
            resolution_domain=args.resolution_domain,
        )
        print(f"Writing {final_path} ...")
        adata.write_h5ad(final_path)
        print("Done.")


if __name__ == "__main__":
    main()
