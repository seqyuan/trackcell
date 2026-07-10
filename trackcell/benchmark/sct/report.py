"""Generate Sphinx artifacts for SCT benchmark."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _fmt(x: float | None, nd: int = 4) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and (x != x):  # NaN
        return "—"
    return f"{x:.{nd}f}"


def render_sct_benchmark_table(payload: dict[str, Any]) -> str:
    """RST list-table for ref-seed stepwise runs."""
    ref = payload.get("ref_seed", 1448145)
    rows = [s for s in payload.get("stepwise", []) if s.get("seed") == ref]
    lines = [
        ".. list-table:: SCT stepwise parity (ref seed, native Python vs R)",
        "   :header-rows: 1",
        "   :widths: 12 8 10 10 10 10",
        "",
        "   * - Sample",
        "     - Native HVG J",
        "     - R fit→Py HVG J",
        "     - θ corr (step1)",
        "     - RV Spearman",
        "     - Residual corr",
    ]
    for s in rows:
        native = s["native_full_vst"]
        golden = s["r_step1_golden_path"]
        rfit = s["r_fit_to_py_residuals"]
        rv = native.get("gene_attr", {}).get("residual_variance", {})
        res = native.get("residuals", {})
        lines.extend(
            [
                f"   * - {s['sample']}",
                f"     - {_fmt(native.get('hvg_jaccard'))}",
                f"     - {_fmt(rfit.get('hvg_jaccard'))}",
                f"     - {_fmt(native.get('model_pars_step1', {}).get('theta', {}).get('corr'), 3)}",
                f"     - {_fmt(rv.get('spearman'), 3)}",
                f"     - {_fmt(res.get('corr'), 3)}",
            ]
        )

    rand = payload.get("randomness", {})
    r_inter = rand.get("r_inter_seed_hvg_jaccard", {})
    py_vs = rand.get("py_vs_r_ref_hvg_jaccard", {})
    lines.extend(
        [
            "",
            ".. list-table:: HVG randomness (6 R seeds × 3 samples)",
            "   :header-rows: 1",
            "   :widths: 30 20",
            "",
            "   * - Metric",
            "     - Value",
            f"   * - R inter-seed HVG Jaccard (median)",
            f"     - {_fmt(r_inter.get('median'))}",
            f"   * - R inter-seed HVG Jaccard (min – max)",
            f"     - {_fmt(r_inter.get('min'))} – {_fmt(r_inter.get('max'))}",
            f"   * - Python vs R (ref seed) HVG Jaccard (median)",
            f"     - {_fmt(py_vs.get('median'))}",
            f"   * - Py median within R seed min–max",
            f"     - {rand.get('py_median_within_r_seed_range')}",
            f"   * - Py median within R inter-seed IQR",
            f"     - {rand.get('py_median_within_r_inter_seed_iqr')}",
        ]
    )
    return "\n".join(lines) + "\n"


def write_sct_benchmark_docs(
    payload: dict[str, Any],
    *,
    output_dir: Path | None = None,
    docs_generated: Path | None = None,
) -> tuple[Path, Path]:
    """Write JSON + RST under ``docs/source/_generated/``."""
    repo = Path(__file__).resolve().parents[3]
    out = output_dir or (repo / "benchmark_out" / "sct")
    gen = docs_generated or (repo / "docs" / "source" / "_generated")
    out.mkdir(parents=True, exist_ok=True)
    gen.mkdir(parents=True, exist_ok=True)

    json_path = out / "sct_benchmark_results.json"
    json_path.write_text(json.dumps(payload, indent=2))
    (gen / "sct_benchmark_results.json").write_text(json_path.read_text())

    rst = render_sct_benchmark_table(payload)
    rst_path = out / "sct_benchmark_table.rst"
    rst_path.write_text(rst)
    (gen / "sct_benchmark_table.rst").write_text(rst)
    return json_path, rst_path
