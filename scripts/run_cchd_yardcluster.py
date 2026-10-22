#!/usr/bin/env python3
"""Run YardCluster on CCHD Xenium TMA h5ad (per-sample workflow).

Example::

    python scripts/run_cchd_yardcluster.py \\
        --h5ad /Volumes/data/co_project/CCHD/xenium_analysis/01_data/85811_S.h5ad \\
        --out 85811_S_yardcluster.h5ad

Requires bbknn for cross-slice integration (``pip install bbknn``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import trackcell as tcl  # noqa: E402


DEFAULT_CCHD_DATA = Path("/Volumes/data/co_project/CCHD/xenium_analysis/01_data")


def _prepare_adata(adata, min_counts: int, min_cells: int):
    if "slice_id" in adata.obs.columns:
        adata = adata[adata.obs["slice_id"] != "debris"].copy()
    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    elif "spatial" not in adata.obsm:
        raise ValueError("Need obs x_centroid/y_centroid or obsm['spatial'].")

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description="CCHD Xenium TMA YardCluster")
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=DEFAULT_CCHD_DATA / "85811_S.h5ad",
        help="Input h5ad (default: CCHD 85811_S)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output h5ad path (default: <stem>_yardcluster.h5ad next to input)",
    )
    parser.add_argument("--batch-key", default="slice_id")
    parser.add_argument("--sketch-n", type=int, default=400_000)
    parser.add_argument("--no-sketch", action="store_true")
    parser.add_argument("--hybrid-alpha", type=float, default=0.25)
    parser.add_argument("--resolution-celltype", type=float, default=1.0)
    parser.add_argument("--resolution-domain", type=float, default=0.8)
    parser.add_argument("--min-counts", type=int, default=30)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--dry-run", action="store_true", help="Load + filter only")
    args = parser.parse_args()

    h5ad = args.h5ad.expanduser().resolve()
    if not h5ad.is_file():
        raise SystemExit(f"Input not found: {h5ad}")

    out = args.out
    if out is None:
        out = h5ad.with_name(f"{h5ad.stem}_yardcluster.h5ad")
    else:
        out = out.expanduser().resolve()

    print(f"Loading {h5ad} ...")
    adata = sc.read_h5ad(h5ad)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    adata = _prepare_adata(adata, args.min_counts, args.min_cells)
    print(f"  after filter: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if args.dry_run:
        print("Dry run — skipping SCT and YardCluster.")
        return

    print("SCT integration (batch_key=%s) ..." % args.batch_key)
    tcl.tl.run_sct_integration(
        adata,
        batch_key=args.batch_key,
        layer=None,
        sct_batch_n_jobs=-1,
    )
    tcl.tl.release_sct_integration_cache(adata, what="scale_data")

    sketch = not args.no_sketch and adata.n_obs > args.sketch_n
    print(
        "YardCluster (platform=xenium_tma, sketch=%s, n_sketch=%s) ..."
        % (sketch, args.sketch_n if sketch else adata.n_obs)
    )
    tcl.tl.spatial_cluster(
        adata,
        platform="xenium_tma",
        batch_key=args.batch_key,
        expression_source="sct_prep",
        integrate="bbknn",
        graph_mode="hybrid",
        hybrid_alpha=args.hybrid_alpha,
        mode="auto",
        sketch=sketch,
        sketch_n=args.sketch_n,
        resolution_celltype=args.resolution_celltype,
        resolution_domain=args.resolution_domain,
        key_added="yardcluster",
    )

    adata.obs["leiden_center_yard"] = adata.obs["yardcluster_celltype"]
    adata.obs["leiden_nbr_yard"] = adata.obs["yardcluster_domain"]

    print(f"Writing {out} ...")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print("Done.")
    print("  yardcluster_celltype:", adata.obs["yardcluster_celltype"].nunique(), "clusters")
    print("  yardcluster_domain:", adata.obs["yardcluster_domain"].nunique(), "clusters")


if __name__ == "__main__":
    main()
