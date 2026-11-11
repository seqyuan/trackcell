#!/usr/bin/env python
"""Measure what importing trackcell actually costs, and what it drags in.

Exists because "the import is slow" is not a testable claim. This runs each
import target in a **fresh subprocess** (``python -X importtime``), so a warm
``sys.modules`` from the parent cannot flatter the numbers, and reports both the
wall time and which heavyweight dependencies were loaded as a side effect.

``matplotlib.pyplot`` is the interesting one: importing it *locks the backend*, so
``matplotlib.use("Agg")`` must run before ``import trackcell`` to take effect —
a silent trap on headless CI that only shows up as a hang or an X11 error.

Usage
-----
    python scripts/measure_import.py                 # report
    python scripts/measure_import.py --check         # assert the budget
    python scripts/measure_import.py --json          # machine-readable

Budget (``--check``):
    import trackcell        < 1.0 s, and must not load matplotlib.pyplot or scanpy
    import trackcell.tl     < 1.0 s, and must not load matplotlib.pyplot
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Modules whose presence in a *core* import is a regression. ``matplotlib.pyplot``
#: locks the backend; ``scanpy`` pulls the whole single-cell stack.
HEAVY = {
    "matplotlib.pyplot": "pyplot",
    "scanpy": "scanpy",
    "geopandas": "geopandas",
    "zarr": "zarr",
    "imageio": "imageio",
}

TARGETS = ("trackcell", "trackcell.tl", "trackcell.io", "trackcell.pl")

#: Wall-time ceiling per target. ``trackcell.pl`` is excluded: it is the plotting
#: module, so it will always pay for pyplot.
BUDGET = {
    "trackcell": 1.0,
    "trackcell.tl": 1.0,
    "trackcell.io": 1.6,
}

_PROBE = r"""
import json, sys, time
t0 = time.perf_counter()
__import__({module!r})
elapsed = time.perf_counter() - t0
loaded = [m for m in {heavy!r} if m in sys.modules]
print("@RESULT@" + json.dumps({{"elapsed": elapsed, "loaded": loaded}}))
"""


def measure(module: str) -> dict:
    """Import ``module`` in a fresh interpreter and report cost and side effects."""
    code = _PROBE.format(module=module, heavy=sorted(HEAVY))
    proc = subprocess.run(
        [sys.executable, "-X", "importtime", "-c", code],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    result = None
    for line in proc.stdout.splitlines():
        if line.startswith("@RESULT@"):
            result = json.loads(line[len("@RESULT@") :])
    if result is None:
        raise RuntimeError(f"probe for {module} produced no result:\n{proc.stderr}")

    # Slowest self-times from ``-X importtime``, for the narrative report.
    slowest = []
    for line in proc.stderr.splitlines():
        parts = line.split("|")
        if len(parts) != 3:
            continue
        try:
            micros = int(parts[1])
        except ValueError:
            continue
        slowest.append((micros, parts[2].strip()))
    slowest.sort(reverse=True)
    result["slowest"] = slowest[:8]
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="enforce the budget")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args()

    results = {target: measure(target) for target in TARGETS}

    if args.json:
        print(json.dumps({k: v for k, v in results.items()}, indent=2))
    else:
        for target, info in results.items():
            loaded = sorted(info["loaded"])
            print(f"{target:16s} {info['elapsed']:5.2f}s  loads: {loaded or '-'}")
            for micros, name in info["slowest"][:5]:
                print(f"    {micros / 1000:7.1f} ms  {name}")

    if not args.check:
        return 0

    failures = []
    core = results["trackcell"]
    for module_name, label in HEAVY.items():
        if module_name in core["loaded"] and module_name in (
            "matplotlib.pyplot",
            "scanpy",
        ):
            failures.append(f"import trackcell pulls {label} ({module_name})")
    if "matplotlib.pyplot" in results["trackcell.tl"]["loaded"]:
        failures.append("import trackcell.tl pulls matplotlib.pyplot")

    for target, ceiling in BUDGET.items():
        elapsed = results[target]["elapsed"]
        if elapsed > ceiling:
            failures.append(f"{target} took {elapsed:.2f}s, budget {ceiling:.2f}s")

    if failures:
        print("\nFAILED:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print("\nimport budget met")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
