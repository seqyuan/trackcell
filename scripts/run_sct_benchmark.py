#!/usr/bin/env python3
"""
Run trackcell SCT benchmark (GSE288946, 3 samples).

1. Optionally export R reference via ata (``--export-r``).
2. Stepwise: R step N → Python step N+1, compare metrics.
3. End-to-end native Python vs R; randomness across 6 R seeds.

Example::

    python scripts/run_sct_benchmark.py --export-r --run
    python scripts/run_sct_benchmark.py --run --benchmark-root /Volumes/process/tmp/tcl_test
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from trackcell.benchmark.sct import (  # noqa: E402
    BENCHMARK_ROOT,
    BENCHMARK_SAMPLES,
    BENCHMARK_SEEDS,
    run_sct_benchmark_suite,
    write_sct_benchmark_docs,
)

RSCRIPT = REPO / "trackcell" / "tl" / "rscripts" / "export_sct_benchmark_r.R"
ATA = Path("/Volumes/data/GOPATH/bin/ata")


def export_r_commands(benchmark_root: Path, seeds: tuple[int, ...], samples: tuple[str, ...]) -> Path:
    """Write ata command file for R exports."""
    cmd_file = benchmark_root / "steps" / "sct_benchmark" / "export_r_commands.txt"
    cmd_file.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for seed in seeds:
        for sample in samples:
            lines.append(
                "SCT_BENCHMARK_ROOT="
                f"{benchmark_root} SCT_SAMPLE={sample} SCT_SEED={seed} "
                f"conda run -n st Rscript {RSCRIPT}"
            )
    cmd_file.write_text("\n".join(lines) + "\n")
    return cmd_file


def main() -> None:
    parser = argparse.ArgumentParser(description="trackcell SCT benchmark")
    parser.add_argument("--benchmark-root", type=Path, default=BENCHMARK_ROOT)
    parser.add_argument("--export-r", action="store_true", help="Run R exports via ata")
    parser.add_argument("--run", action="store_true", help="Run Python comparisons")
    parser.add_argument("--threads", type=int, default=3, help="ata concurrency for R export")
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=list(BENCHMARK_SEEDS),
        help="R/Python vst seeds",
    )
    args = parser.parse_args()

    seeds = tuple(args.seeds)
    samples = BENCHMARK_SAMPLES
    root = args.benchmark_root

    if args.export_r:
        cmd_file = export_r_commands(root, seeds, samples)
        if not ATA.exists():
            print(f"ata not found at {ATA}; run commands manually from {cmd_file}")
        else:
            n_tasks = len(seeds) * len(samples)
            threads = min(args.threads, n_tasks)
            print(f"Exporting R reference: {n_tasks} jobs, {threads} threads …")
            subprocess.run(
                [str(ATA), "-i", str(cmd_file), "-t", str(threads)],
                check=True,
            )

    if args.run:
        print("Running Python benchmark suite …")
        payload = run_sct_benchmark_suite(
            benchmark_root=root,
            seeds=seeds,
            samples=samples,
        )
        json_path, rst_path = write_sct_benchmark_docs(payload)
        print(f"Wrote {json_path}")
        print(f"Wrote {rst_path}")
        print(f"Status: {payload['status']}")
        if payload.get("missing_r_exports"):
            print("Missing R exports:", payload["missing_r_exports"][:5], "…")
        summary = payload.get("summary", {})
        print("Native HVG Jaccard (ref seed):", summary.get("native_hvg_jaccard"))

    if not args.export_r and not args.run:
        parser.print_help()


if __name__ == "__main__":
    main()
