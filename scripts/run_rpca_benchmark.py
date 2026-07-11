#!/usr/bin/env python3
"""
Run trackcell SCT + RPCA integration benchmark (GSE288946, 3 samples).

1. Optionally export R reference via ata (``--export-r``).
2. Stepwise: R step N → Python step N+1.
3. End-to-end native Python vs R (filter: 100 cells/gene, 50 genes/cell).

Example::

    python scripts/run_rpca_benchmark.py --export-r --run
    python scripts/run_rpca_benchmark.py --run --benchmark-root /Volumes/process/tmp/tcl_test
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from trackcell.benchmark.rpca import (  # noqa: E402
    RPCA_BENCHMARK_ROOT,
    RPCA_SEEDS,
    run_rpca_benchmark_suite,
    write_rpca_benchmark_docs,
)

RSCRIPT = REPO / "trackcell" / "tl" / "rscripts" / "export_rpca_benchmark_r.R"
ATA = Path("/Volumes/data/GOPATH/bin/ata")


def _resolve_rscript() -> list[str]:
    st_r = Path("/pmo/anaconda3/envs/st/bin/Rscript")
    if st_r.exists():
        return [str(st_r)]
    if shutil.which("Rscript"):
        return ["Rscript"]
    for candidate in ("/pmo/miniforge3/bin/Rscript",):
        if Path(candidate).exists():
            return [candidate]
    return ["conda", "run", "-n", "st", "Rscript"]


def export_r_commands(benchmark_root: Path, seeds: tuple[int, ...]) -> Path:
    cmd_file = benchmark_root / "steps" / "rpca_benchmark" / "export_r_commands.txt"
    cmd_file.parent.mkdir(parents=True, exist_ok=True)
    rscript = _resolve_rscript()
    lines: list[str] = []
    for seed in seeds:
        prefix = f"RPCA_BENCHMARK_ROOT={benchmark_root} RPCA_SEED={seed} "
        if rscript[0] == "conda":
            lines.append(prefix + " ".join(rscript + [str(RSCRIPT)]))
        else:
            lines.append(prefix + " ".join(rscript + [str(RSCRIPT)]))
    cmd_file.write_text("\n".join(lines) + "\n")
    return cmd_file


def main() -> None:
    parser = argparse.ArgumentParser(description="trackcell RPCA integration benchmark")
    parser.add_argument("--benchmark-root", type=Path, default=RPCA_BENCHMARK_ROOT)
    parser.add_argument("--export-r", action="store_true", help="Run R exports via ata")
    parser.add_argument("--run", action="store_true", help="Run Python comparisons")
    parser.add_argument(
        "--skip-e2e",
        action="store_true",
        help="Skip end-to-end native integrate_sct_rpca (slow)",
    )
    parser.add_argument("--threads", type=int, default=1, help="ata concurrency for R export")
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=list(RPCA_SEEDS),
        help="Integration seeds (R set.seed / Python random_state)",
    )
    args = parser.parse_args()

    seeds = tuple(args.seeds)
    root = args.benchmark_root

    if args.export_r:
        cmd_file = export_r_commands(root, seeds)
        if not ATA.exists():
            print(f"ata not found at {ATA}; run commands manually from {cmd_file}")
        else:
            print(f"Exporting R RPCA reference: {len(seeds)} job(s) …")
            subprocess.run(
                [str(ATA), "-i", str(cmd_file), "-t", str(min(args.threads, len(seeds)))],
                check=True,
                env={**os.environ, "RPCA_BENCHMARK_ROOT": str(root)},
            )

    if args.run:
        print("Running Python RPCA benchmark suite …")
        payload = run_rpca_benchmark_suite(
            benchmark_root=root, seeds=seeds, skip_e2e=args.skip_e2e
        )
        json_path, rst_path = write_rpca_benchmark_docs(payload)
        print(f"Wrote {json_path}")
        print(f"Wrote {rst_path}")
        print(f"Status: {payload['status']}")
        if payload.get("missing_r_exports"):
            print("Missing R exports:", payload["missing_r_exports"])
        print("Summary:", payload.get("summary"))

    if not args.export_r and not args.run:
        parser.print_help()


if __name__ == "__main__":
    main()
