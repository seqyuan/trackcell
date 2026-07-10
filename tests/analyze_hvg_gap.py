"""HVG Jaccard ~0.55 decomposition (GSM8779707 Tier-B, gene_symbols loader).

Run::

    PYTHONPATH=. python tests/analyze_hvg_gap.py

Summary
-------
Native Python HVG Jaccard ~0.99 vs R ``vst()`` export after fixing:

1. **Cell column order:** ``umi[:, cells.isin(cells_step1)]`` misaligns counts vs
   ``cell_attr.loc[cells_step1]`` offset — use ``cells.get_indexer(cells_step1)``.
2. **Gene rownames:** Seurat ``.N`` duplicate order ≠ scanpy ``-N``; benchmark loader
   assigns R rownames by position from ``ref/{sample}_genes_min_cells.txt``.

pyglmGamPoi batch θ matches R ``fit_glmGamPoi_offset`` once inputs align (raw Inf 224/224).

Injection controls
------------------
- ``r_step1`` (R step1 + Py ksmooth/residuals): HVG Jaccard **~0.997**
- ``r_fit`` (R model_pars_fit + Py residuals): HVG Jaccard **~0.999**

Native pipeline (Tier-B fixed cells/genes)
------------------------------------------
| Stage | Metric | Py vs R |
|-------|--------|---------|
| genes_vst membership | Jaccard | **1.0** |
| step1 Inf θ | count | **277 vs 277** |
| step1 finite θ | log10 Pearson | **~0.97** |
| model_pars_fit θ | Pearson | **~1.0** |
| residual_variance | Spearman | **~1.0** |
| HVG top-3000 | Jaccard | **~0.99** |

Root cause (fixed)
------------------
1. **Cell order bug** in ``vst()`` step1 slice (``cells.isin`` vs ``get_indexer``).
2. **Gene name order** for duplicate symbols (benchmark loader sidecar rownames).

Mitigation
----------
- Native path: fixed in ``trackcell/tl/_sctransform.py`` + ``load_seurat_parity_umi``.
- Benchmarks: export ``ref/{sample}_genes_min_cells.txt`` from Seurat rownames after filter.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    from .sct_benchmark_data import DEFAULT_SAMPLE, load_seurat_parity_with_cell_attr
except ImportError:
    from sct_benchmark_data import DEFAULT_SAMPLE, load_seurat_parity_with_cell_attr
from trackcell.tl._r_sctransform import align_r_gene_list, load_r_vst_export
from trackcell.tl._sctransform import vst

ROOT = Path("/Volumes/process/tmp/tcl_test")


def main() -> None:
    sample = DEFAULT_SAMPLE
    export = load_r_vst_export(ROOT / "steps/r_sct_stepwise" / sample)
    r_hvg = set((ROOT / "steps/r_sct_stepwise" / sample / "hvg_top3000.txt").read_text().strip().split("\n"))
    r_ga = pd.read_csv(ROOT / "steps/r_sct_stepwise" / sample / "gene_attr.csv", index_col=0)
    genes_vst_r = set((ROOT / "steps/r_sct_stepwise" / sample / "genes_vst.txt").read_text().strip().split("\n"))

    umi, genes, cells, cell_attr = load_seurat_parity_with_cell_attr(sample)
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))

    def hvg_j(out: dict) -> float:
        py = set(out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str))
        return len(py & r_hvg) / len(py | r_hvg)

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
    shared = genes_step1.intersection(export["model_pars_step1"].index)
    py_t = out["model_pars"].loc[shared, "theta"]
    rr_t = export["model_pars_step1"].loc[shared, "theta"]
    fin = np.isfinite(py_t) & np.isfinite(rr_t)
    fit = out["model_pars_fit"]
    rf = export["model_pars_fit"]
    sh = shared.intersection(fit.index).intersection(rf.index)
    ffin = np.isfinite(fit.loc[sh, "theta"]) & np.isfinite(rf.loc[sh, "theta"])
    py_vst = set(out["gene_attr"].index.astype(str))
    rv_sh = out["gene_attr"].index.intersection(r_ga.index)
    rv_fin = np.isfinite(out["gene_attr"].loc[rv_sh, "residual_variance"]) & np.isfinite(
        r_ga.loc[rv_sh, "residual_variance"]
    )

    out_r1 = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        model_pars_step1=export["model_pars_step1"],
        genes_log_gmean=export["genes_log_gmean"],
        genes_log_gmean_step1=export["genes_log_gmean_step1"],
        return_corrected_umi=False,
        seed=1448145,
    )
    out_fit = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        model_pars_fit_override=export["model_pars_fit"],
        return_corrected_umi=False,
        seed=1448145,
    )

    report = {
        "sample": sample,
        "loader": "seurat_parity (gene_symbols + CreateSeuratObject QC)",
        "native_py": {
            "HVG_jaccard": hvg_j(out),
            "genes_vst_jaccard": len(py_vst & genes_vst_r) / len(py_vst | genes_vst_r),
            "step1_inf_py": int((~np.isfinite(py_t)).sum()),
            "step1_inf_r": int((~np.isfinite(rr_t)).sum()),
            "step1_theta_log10_corr": float(np.log10(py_t[fin].clip(1e-9)).corr(np.log10(rr_t[fin].clip(1e-9)))),
            "fit_theta_corr": float(fit.loc[sh, "theta"][ffin].corr(rf.loc[sh, "theta"][ffin])),
            "RV_spearman": float(
                spearmanr(
                    out["gene_attr"].loc[rv_sh, "residual_variance"][rv_fin],
                    r_ga.loc[rv_sh, "residual_variance"][rv_fin],
                ).correlation
            ),
        },
        "r_step1_golden": {"HVG_jaccard": hvg_j(out_r1)},
        "r_fit_residuals": {"HVG_jaccard": hvg_j(out_fit)},
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
