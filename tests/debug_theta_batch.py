"""Debug batch theta parity: pyglmGamPoi vs R on Tier-B GSM8779707.

Run::

    PYTHONPATH=. python tests/debug_theta_batch.py

Findings (2026-07)
------------------
1. **Cell order:** ``umi[:, cells.isin(cells_step1)]`` misaligns counts vs
   ``cell_attr.loc[cells_step1]`` / offset — use ``cells.get_indexer(cells_step1)``.
2. **Gene names:** Seurat ``.N`` duplicate suffix order ≠ scanpy ``-N``; assign R
   rownames by position via ``ref/{sample}_genes_min_cells.txt`` for benchmarks.
3. After (1)+(2), pyglmGamPoi batch θ matches R ``fit_glmGamPoi_offset`` (raw inf
   224/224; after ``mark_suspicious_theta`` ~277).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from tests.sct_benchmark_data import BENCHMARK_ROOT, DEFAULT_SAMPLE, load_seurat_parity_with_cell_attr
from trackcell.tl._r_sctransform import fit_glmGamPoi_offset_r, load_r_vst_export
from trackcell.tl._sctransform import _build_design_matrix, vst
from trackcell.tl._sctransform_v2 import mark_suspicious_theta
from pyglmGamPoi.trackcell_compat import fit_offset_model as py_fit

ROOT = BENCHMARK_ROOT


def main() -> None:
    sample = DEFAULT_SAMPLE
    export = load_r_vst_export(ROOT / "steps/r_sct_stepwise" / sample)
    r_hvg = set((ROOT / "steps/r_sct_stepwise" / sample / "hvg_top3000.txt").read_text().strip().split("\n"))
    umi, genes, cells, cell_attr = load_seurat_parity_with_cell_attr(sample)
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)
    genes_step1 = pd.Index(export["genes_step1"]).intersection(genes)
    cells_idx = cells.get_indexer(cells_step1)
    umi_step1 = umi[genes.get_indexer(genes_step1)][:, cells_idx]
    reg = _build_design_matrix(cell_attr.loc[cells_step1], "y ~ log_umi")
    r_export = export["model_pars_step1"].reindex(genes_step1)

    py_raw = py_fit(umi_step1, reg, genes_step1, allow_inf_theta=True)
    py_ms = mark_suspicious_theta(py_raw, umi, genes).reindex(genes_step1)
    r_raw = fit_glmGamPoi_offset_r(umi_step1, reg, genes_step1, allow_inf_theta=True)
    r_ms = mark_suspicious_theta(r_raw, umi, genes).reindex(genes_step1)

    fin = np.isfinite(py_ms["theta"]) & np.isfinite(r_export["theta"])
    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=1448145,
    )
    py_hvg = set(out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str))
    report = {
        "sample": sample,
        "py_raw_inf": int((~np.isfinite(py_raw["theta"])).sum()),
        "r_raw_inf": int((~np.isfinite(r_raw["theta"])).sum()),
        "py_mark_inf": int((~np.isfinite(py_ms["theta"])).sum()),
        "r_export_inf": int((~np.isfinite(r_export["theta"])).sum()),
        "py_vs_r_raw_theta_corr": float(py_raw.loc[fin, "theta"].corr(r_ms.loc[fin, "theta"])) if fin.any() else None,
        "py_vs_export_log10_theta_corr": (
            float(
                np.log10(py_ms.loc[fin, "theta"].clip(1e-9)).corr(
                    np.log10(r_export.loc[fin, "theta"].clip(1e-9))
                )
            )
            if fin.any()
            else None
        ),
        "HVG_jaccard": len(py_hvg & r_hvg) / len(py_hvg | r_hvg),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
