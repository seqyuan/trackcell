"""Generate Sphinx artifacts for RPCA integration benchmark."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _fmt(x: float | None, nd: int = 4) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and (x != x):
        return "—"
    return f"{x:.{nd}f}"


def render_rpca_benchmark_table(payload: dict[str, Any]) -> str:
    ref = payload.get("ref_seed", 1448145)
    sct = next((s for s in payload.get("sct_stepwise", []) if s.get("seed") == ref), None)
    step = next((s for s in payload.get("stepwise", []) if s.get("seed") == ref), None)
    e2e = next((s for s in payload.get("e2e", []) if s.get("seed") == ref), None)

    lines = [
        ".. list-table:: Integration SCT stepwise (ref seed, vs R SCTransform)",
        "   :header-rows: 1",
        "   :widths: 14 12 10 24",
        "",
        "   * - Sample / merged",
        "     - Metric",
        "     - Value",
        "     - Notes",
    ]
    if sct:
        for row in sct.get("per_sample", []):
            lines.extend(
                [
                    f"   * - {row['sample']}",
                    "     - Native HVG J",
                    f"     - {_fmt(row['native_hvg_jaccard'])}",
                    "     - native vst vs R SCTransform",
                    f"   * - {row['sample']}",
                    "     - R fit→Py HVG J",
                    f"     - {_fmt(row['r_fit_to_py_hvg_jaccard'])}",
                    f"     - residual corr {_fmt(row.get('r_fit_to_py_residual_corr'), 3)}",
                ]
            )
        merged = sct.get("merged", {})
        prep = merged.get("prep_native_sct_r_features", {})
        lines.extend(
            [
                "   * - merged",
                "     - Features J (native)",
                f"     - {_fmt(merged.get('features_jaccard_native_vs_r'))}",
                "     - native SCT SelectIntegrationFeatures",
                "   * - merged",
                "     - Prep corr (R features)",
                f"     - {_fmt(prep.get('gene_corr_median'))}",
                "     - native SCT + forced R feature list",
            ]
        )

    lines.extend(
        [
            "",
            ".. list-table:: RPCA stepwise parity (ref seed, R input → Python step)",
            "   :header-rows: 1",
            "   :widths: 14 10 10 10 10",
            "",
            "   * - Step",
            "     - Metric",
            "     - Value",
            "     - Notes",
            "     -",
        ]
    )
    if step:
        st = step["steps"]
        rows = [
            ("03 features", "Jaccard", _fmt(st["03_features"]["jaccard"]), "R SCT → Py SelectIntegrationFeatures", ""),
            ("04 prep", "Gene corr (median)", _fmt(st["04_prep"].get("gene_corr_median")), "R SCT + R features → Py Prep", ""),
            (
                "05 anchors",
                "Local pair Jaccard",
                _fmt(st["05_anchors"]["local_pair_jaccard"]),
                f"count ratio {_fmt(st['05_anchors']['unique_pair_ratio_py_over_r'], 3)}",
                "",
            ),
            (
                "06 integrate",
                "Gene corr (median)",
                _fmt(st["06_integrated_r_anchors"].get("gene_corr_median")),
                "R prep + R anchors → Py IntegrateData",
                "",
            ),
            (
                "06 integrate",
                "Gene corr (Py anchors)",
                _fmt(st["06_integrated_py_anchors"].get("gene_corr_median")),
                "R prep + Py anchors → Py IntegrateData",
                "",
            ),
        ]
        for r in rows:
            lines.extend([f"   * - {r[0]}", f"     - {r[1]}", f"     - {r[2]}", f"     - {r[3]}", f"     - {r[4]}"])

    if e2e:
        summary = payload.get("summary", {})
        step_ref = summary.get("stepwise_ref", {})
        e2e_ref = summary.get("e2e_ref", {})
        lines.extend(
            [
                "",
                ".. list-table:: Python vs Seurat quick comparison (ref seed)",
                "   :header-rows: 1",
                "   :widths: 28 12 12",
                "",
                "   * - Metric",
                "     - Stepwise",
                "     - E2E",
                f"   * - Features Jaccard",
                f"     - {_fmt(step_ref.get('features_jaccard'))}",
                f"     - {_fmt(e2e_ref.get('features_jaccard'))}",
                f"   * - Prep gene corr",
                f"     - {_fmt(step_ref.get('prep_corr_median'))}",
                f"     - {_fmt(e2e_ref.get('prep_corr_median'))}",
                f"   * - Anchor local-pair Jaccard",
                f"     - {_fmt(step_ref.get('anchor_local_jaccard'))}",
                f"     - {_fmt(e2e_ref.get('anchor_local_jaccard'))}",
                f"   * - Integrate gene corr",
                f"     - {_fmt(step_ref.get('integrate_corr_median_r_anchors'))}",
                f"     - {_fmt(e2e_ref.get('integrate_corr_median'))}",
            ]
        )

    lines.extend(
        [
            "",
            ".. list-table:: RPCA end-to-end native Python vs R (ref seed)",
            "   :header-rows: 1",
            "   :widths: 22 12",
            "",
            "   * - Metric",
            "     - Value",
        ]
    )
    if e2e:
        lines.extend(
            [
                f"   * - SCT HVG Jaccard (median per batch)",
                f"     - {_fmt(payload.get('summary', {}).get('e2e_ref', {}).get('sct_hvg_jaccard_median'))}",
                f"   * - Integration features Jaccard",
                f"     - {_fmt(e2e['features_jaccard'])}",
                f"   * - Shared integration features",
                f"     - {e2e.get('n_shared_features', '—')}",
                f"   * - Prep gene corr (aligned columns)",
                f"     - {_fmt(e2e['prep'].get('gene_corr_median'))}",
                f"   * - Prep gene corr (shared features only)",
                f"     - {_fmt(e2e.get('prep_shared_features_only', {}).get('gene_corr_median'))}",
                f"   * - Anchor local pair Jaccard",
                f"     - {_fmt(e2e['anchors']['local_pair_jaccard'])}",
                f"   * - Integrated gene corr (aligned)",
                f"     - {_fmt(e2e['integrated'].get('gene_corr_median'))}",
                f"   * - Integrated gene corr (shared features)",
                f"     - {_fmt(e2e.get('integrated_shared_features_only', {}).get('gene_corr_median'))}",
            ]
        )
    return "\n".join(lines) + "\n"


def write_rpca_benchmark_docs(
    payload: dict[str, Any],
    *,
    output_dir: Path | None = None,
    docs_generated: Path | None = None,
) -> tuple[Path, Path]:
    repo = Path(__file__).resolve().parents[3]
    out = output_dir or (repo / "benchmark_out" / "rpca")
    gen = docs_generated or (repo / "docs" / "source" / "_generated")
    out.mkdir(parents=True, exist_ok=True)
    gen.mkdir(parents=True, exist_ok=True)

    json_path = out / "rpca_benchmark_results.json"
    json_path.write_text(json.dumps(payload, indent=2))
    (gen / "rpca_benchmark_results.json").write_text(json_path.read_text())

    rst = render_rpca_benchmark_table(payload)
    rst_path = out / "rpca_benchmark_table.rst"
    rst_path.write_text(rst)
    (gen / "rpca_benchmark_table.rst").write_text(rst)
    return json_path, rst_path
