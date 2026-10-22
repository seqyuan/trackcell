#!/usr/bin/env python3
"""CLI for YardCluster synthetic benchmark."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from trackcell.benchmark.yardcluster.runner import format_benchmark_table, run_yardcluster_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run YardCluster synthetic benchmark")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-per-patch", type=int, default=40)
    parser.add_argument("--n-patches", type=int, default=3)
    parser.add_argument("--json", action="store_true", help="Print JSON instead of table")
    args = parser.parse_args()

    summary = run_yardcluster_benchmark(
        seed=args.seed,
        n_per_patch=args.n_per_patch,
        n_patches=args.n_patches,
    )
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(format_benchmark_table(summary))
        print()
        print(
            f"Best: {summary['best_config']} "
            f"(ARI={summary['best_ari']:.3f}, "
            f"coherence={summary['best_spatial_coherence']:.3f})"
        )


if __name__ == "__main__":
    main()
