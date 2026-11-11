#!/usr/bin/env python
"""Generate the tiny fixtures that let an agent verify trackcell code paths.

Why this exists
---------------
``tests/`` shipped assertions but no data, so an agent writing trackcell code had
no way to close the loop: it could read a docstring, emit a snippet, and never
learn whether the snippet runs. These fixtures are small enough to commit (~170 KB
total) and shaped like real platform output, so every reader, subsampler,
distance function and polygon plotter can be exercised in seconds.

Everything is deterministically generated from a fixed seed, so re-running this
script reproduces identical *content*. It does not reproduce identical *bytes*,
and should not be expected to: pyarrow, HDF5 and libpng all change their output
encoding between versions. Pinning any of them to keep a hash stable would make
the fixture an assertion about that library instead of about trackcell.

Practical consequence: regenerating with a different pyarrow version will show
binary diffs in ``tests/data/**/*.parquet`` even though nothing meaningful
changed. **Only re-run this script when the builders themselves change** — there
is no reason to regenerate otherwise.

Why not reuse ``tests/conftest.py``'s writers?
----------------------------------------------
It has ``write_10x_h5`` / ``write_scalefactors`` / ``write_cell_segmentations``,
and reusing them was the first instinct. Three of the four turned out to be too
rigid for what a *usable* fixture needs: ``write_scalefactors`` hardcodes its
values and omits ``microns_per_pixel`` (so ``resolve_microns_per_unit`` cannot be
exercised), ``write_cell_segmentations`` only emits axis-aligned squares on a
``10*i`` diagonal (no spatial grid, no cell-type blocks), and ``write_image``
emits a flat grey field (overlay plots then verify nothing).

Extending those helpers would mean editing shared test infrastructure, and the
result would couple a *product* artifact — the fixtures agents read — to *test
scaffolding* internals. A conftest refactor would then silently change what the
fixtures contain. Keeping the ~90 lines of writers local is the cheaper option,
and ``--check`` is what actually prevents rot: it re-validates the committed
files against the current readers, which is the contract that matters.

Layout
------
``tests/data/hd_cellseg_tiny``  SpaceRanger 4.0 ``segmented_outputs`` (pixels)
``tests/data/hd_bin_tiny``      SpaceRanger 4.0 ``binned_outputs/square_016um``
``tests/data/xenium_tiny``      Xenium ``outs`` with two tissue sections (µm)

Usage
-----
    python scripts/make_fixtures.py            # (re)build
    python scripts/make_fixtures.py --check    # build, then read back and report
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

SEED = 20260101
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"

# ── Genes ────────────────────────────────────────────────────────────────────
# Three coherent marker panels plus housekeeping fillers. The panels are what
# give the fixtures real structure: a cluster computed on them must recover the
# spatial blocks, so a broken pipeline becomes visible instead of silently
# returning noise.
MARKERS = {
    "Tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "MKI67", "TOP2A", "PCNA"],
    "Stroma": ["COL1A1", "COL1A2", "COL3A1", "DCN", "ACTA2", "TAGLN", "VIM"],
    "Immune": [
        "PTPRC",
        "CD3D",
        "CD3E",
        "CD2",
        "CD79A",
        "MS4A1",
        "CD37",
        "LYZ",
        "CD68",
        "CD14",
    ],
}
HOUSEKEEPING = [
    "ACTB",
    "GAPDH",
    "B2M",
    "TMSB4X",
    "RPLP0",
    "MALAT1",
    "EEF1A1",
    "RPL13A",
    "RPS18",
    "UBC",
    "HSP90AA1",
    "HSPA1A",
    "FTL",
    "FTH1",
    "S100A9",
    "CFL1",
]
GENES = [g for panel in MARKERS.values() for g in panel] + HOUSEKEEPING


def _gene_index(names):
    return [GENES.index(n) for n in names]


# ── 10x HDF5 writer ──────────────────────────────────────────────────────────
def _write_10x_h5(path: Path, matrix_cells_by_genes, barcodes, gene_id_prefix="ENSGT"):
    """Write a CellRanger-compatible ``.h5``.

    ``scanpy.read_10x_h5`` wants genes x cells in CSC layout, which is what the
    real pipeline stores, so the input is transposed back here.
    """
    genes_by_cells = matrix_cells_by_genes.T.tocsc()
    genes_by_cells.sort_indices()

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("data", data=genes_by_cells.data.astype(np.int32))
        m.create_dataset("indices", data=genes_by_cells.indices.astype(np.int32))
        m.create_dataset("indptr", data=genes_by_cells.indptr.astype(np.int32))
        m.create_dataset("shape", data=np.array(genes_by_cells.shape, dtype=np.int32))

        m.create_dataset(
            "barcodes",
            data=np.array(
                [b.encode() for b in barcodes], dtype=f"S{max(map(len, barcodes))}"
            ),
        )

        feats = m.create_group("features")
        ids = [f"{gene_id_prefix}{i:06d}" for i in range(len(GENES))]
        for key, values in (
            ("id", ids),
            ("name", GENES),
            ("genome", ["tiny"] * len(GENES)),
            ("feature_type", ["Gene Expression"] * len(GENES)),
        ):
            width = max(map(len, values))
            feats.create_dataset(
                key, data=np.array([v.encode() for v in values], dtype=f"S{width}")
            )
        # ``_all_tag_keys`` is present in real files; written for fidelity.
        feats.create_dataset("_all_tag_keys", data=np.array([b"genome"], dtype="S6"))


# ── Tissue images ────────────────────────────────────────────────────────────
def _synthetic_tissue_image(height: int, width: int, rng) -> np.ndarray:
    """Pink H&E-ish background with a darker patch, so overlays are legible.

    Shading is generated on a coarse grid and upsampled. Per-pixel noise would
    look marginally more realistic but PNG cannot compress it, and a fixture is
    only useful if it stays small enough to commit.
    """
    yy, xx = np.mgrid[0:height, 0:width]
    base = np.zeros((height, width, 3), dtype=np.float32)
    base[..., 0] = 246
    base[..., 1] = 226
    base[..., 2] = 232

    cy, cx = height / 2.0, width / 2.0
    radius = min(height, width) * 0.42
    dist = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    tissue = dist < radius
    base[tissue, 0] = 214
    base[tissue, 1] = 152
    base[tissue, 2] = 178

    # Coarse, smooth shading: a few hundred cells upsampled with repeat.
    gh = int(np.ceil(height / 16))
    gw = int(np.ceil(width / 16))
    coarse = rng.normal(0, 6.0, size=(gh, gw, 1))
    shading = np.repeat(np.repeat(coarse, 16, axis=0), 16, axis=1)[:height, :width]

    # Soft edge so the patch does not look like a hard disc.
    edge = np.clip((radius - dist) / 12.0, 0.0, 1.0)[..., None]
    base = base * edge + np.array([246, 226, 232], dtype=np.float32) * (1 - edge)

    return np.clip(base + shading, 0, 255).astype(np.uint8)


def _write_png(path: Path, array: np.ndarray):
    import imageio.v3 as iio

    path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(path, array)


# ── Expression ───────────────────────────────────────────────────────────────
def _counts_for(labels, rng, n_cells, n_genes=len(GENES)):
    """Poisson counts with panel-specific rates.

    A labelled cell expresses its own panel ~6x higher than the others, which is
    enough for PCA/neighbours to separate the blocks without any batch effect.
    """
    rate = np.full((n_cells, n_genes), 1.2)
    for label, panel in MARKERS.items():
        idx = _gene_index(panel)
        rows = np.array([lab == label for lab in labels])
        rate[np.ix_(rows, idx)] = 7.0
    # Housekeeping is flatter and higher, like real data.
    rate[:, _gene_index(HOUSEKEEPING)] = 12.0
    return rng.poisson(rate).astype(np.float32)


# ── Fixture 1: Visium HD cell segmentation ───────────────────────────────────
def build_hd_cellseg(out: Path):
    """SpaceRanger 4.0 ``segmented_outputs``: 36 cells on a 6x6 grid, 4 blocks.

    Coordinates are **full-resolution image pixels** (100 px pitch, 0.3 µm/px),
    which is what makes ``resolve_microns_per_unit`` and ``cells_in_radius``
    meaningful on this fixture.
    """
    rng = np.random.default_rng(SEED)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    n_side, pitch, origin = 6, 100.0, 1000.0
    microns_per_pixel = 0.3
    hires_scalef = 0.5

    blocks = [
        ("Tumor", (-1, -1)),
        ("Stroma", (-1, 1)),
        ("Immune", (1, -1)),
        ("Tumor", (1, 1)),
    ]
    rows = []
    cell_id = 1
    for r in range(n_side):
        for c in range(n_side):
            block_name = blocks[
                (0 if r < n_side // 2 else 2) + (0 if c < n_side // 2 else 1)
            ][0]
            cls_idx = sorted(MARKERS).index(block_name)
            rows.append((cell_id, r, c, block_name, cls_idx))
            cell_id += 1

    labels = [r[3] for r in rows]
    n_cells = len(rows)
    counts = _counts_for(labels, rng, n_cells)

    # Irregular octagon per cell so geometry is not trivially rectangular.
    geojson_features = []
    barcodes = []
    for cid, r, c, label, cls_idx in rows:
        cx = origin + c * pitch
        cy = origin + r * pitch
        radius = pitch * 0.36
        angles = np.linspace(0, 2 * np.pi, 8, endpoint=False) + rng.normal(0, 0.05, 8)
        radii = radius * rng.uniform(0.85, 1.15, 8)
        ring = np.column_stack(
            [cx + radii * np.cos(angles), cy + radii * np.sin(angles)]
        )
        geojson_features.append(
            {
                "type": "Feature",
                "properties": {
                    "cell_id": cid,
                    # SpaceRanger stores this as a stringified Python dict.
                    "classification": str({"name": label, "cls": cls_idx}),
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [list(map(float, p)) for p in ring]
                        + [list(map(float, ring[0]))]
                    ],
                },
            }
        )
        barcodes.append(f"cellid_{str(cid).zfill(9)}-1")

    with open(out / "graphclust_annotated_cell_segmentations.geojson", "w") as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f)

    from scipy import sparse as sp

    _write_10x_h5(
        out / "filtered_feature_cell_matrix.h5", sp.csr_matrix(counts), barcodes
    )

    scalefactors = {
        "microns_per_pixel": microns_per_pixel,
        "tissue_hires_scalef": hires_scalef,
        "tissue_lowres_scalef": 0.125,
        "fiducial_diameter_fullres": 700.0,
        "regist_target_img_scalef": hires_scalef,
    }
    (out / "spatial").mkdir(parents=True, exist_ok=True)
    (out / "spatial" / "scalefactors_json.json").write_text(
        json.dumps(scalefactors, indent=2)
    )

    fullres_span = int(n_side * pitch + 400)
    _write_png(
        out / "spatial" / "tissue_hires_image.png",
        _synthetic_tissue_image(
            int(fullres_span * hires_scalef), int(fullres_span * hires_scalef), rng
        ),
    )
    _write_png(
        out / "spatial" / "tissue_lowres_image.png",
        _synthetic_tissue_image(
            int(fullres_span * 0.125), int(fullres_span * 0.125), rng
        ),
    )
    return {"cells": n_cells, "genes": len(GENES), "blocks": sorted(set(labels))}


# ── Fixture 2: Visium HD 16 µm bins ──────────────────────────────────────────
def build_hd_bin(out: Path):
    """SpaceRanger ``binned_outputs/square_016um``: 100 bins, 32 px pitch."""
    rng = np.random.default_rng(SEED + 1)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    n_side, pitch, origin = 10, 32.0, 1000.0
    barcodes, coords, labels = [], [], []
    for r in range(n_side):
        for c in range(n_side):
            barcodes.append("".join("ACGT"[b] for b in rng.integers(0, 4, 16)) + "-1")
            coords.append((origin + c * pitch, origin + r * pitch))
            labels.append(
                "Stroma" if (r < n_side // 2) ^ (c < n_side // 2) else "Tumor"
            )

    counts = _counts_for(labels, rng, len(barcodes))

    from scipy import sparse as sp

    _write_10x_h5(
        out / "filtered_feature_bc_matrix.h5", sp.csr_matrix(counts), barcodes
    )

    positions = pd.DataFrame(
        {
            "barcode": barcodes,
            "in_tissue": np.ones(len(barcodes), dtype=np.int64),
            "array_row": [r for r in range(n_side) for _ in range(n_side)],
            "array_col": [c for _ in range(n_side) for c in range(n_side)],
            "pxl_row_in_fullres": [y for _, y in coords],
            "pxl_col_in_fullres": [x for x, _ in coords],
        }
    )
    (out / "spatial").mkdir(parents=True, exist_ok=True)
    positions.to_parquet(out / "spatial" / "tissue_positions.parquet", index=False)

    scalefactors = {
        "microns_per_pixel": 0.3,
        "tissue_hires_scalef": 0.5,
        "tissue_lowres_scalef": 0.125,
        "fiducial_diameter_fullres": 700.0,
    }
    (out / "spatial" / "scalefactors_json.json").write_text(
        json.dumps(scalefactors, indent=2)
    )

    span = int(n_side * pitch + 200)
    _write_png(
        out / "spatial" / "tissue_hires_image.png",
        _synthetic_tissue_image(int(span * 0.5), int(span * 0.5), rng),
    )
    _write_png(
        out / "spatial" / "tissue_lowres_image.png",
        _synthetic_tissue_image(int(span * 0.125), int(span * 0.125), rng),
    )
    return {"bins": len(barcodes), "genes": len(GENES), "binsize": 16}


# ── Fixture 3: Xenium with two tissue sections ───────────────────────────────
def build_xenium(out: Path):
    """Xenium ``outs``: 40 cells in two well-separated µm-space sections.

    The 2 000 µm gap is what makes ``spatial_slice_cluster(eps=80)`` return
    exactly two slices, so slice-splitting code can be asserted without guesswork.
    """
    rng = np.random.default_rng(SEED + 2)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    n_per = 20
    section_origins = np.array([[1000.0, 1000.0], [1000.0, 4000.0]])  # 3000 µm apart
    px = 0.2125  # µm per pixel, a real Xenium value

    cell_ids, xy, labels, boundaries, nuclei = [], [], [], [], []
    for sec, (sx, sy) in enumerate(section_origins):
        side = int(np.ceil(np.sqrt(n_per)))
        for k in range(n_per):
            r, c = divmod(k, side)
            cx = sx + c * 60.0 + rng.normal(0, 5)
            cy = sy + r * 60.0 + rng.normal(0, 5)
            cid = str(len(cell_ids) + 1)
            # Section 0 is Tumor/Stroma mix, section 1 is Immune — different
            # composition so per-slice summaries differ measurably.
            label = ("Tumor" if k % 3 else "Stroma") if sec == 0 else "Immune"
            cell_ids.append(cid)
            xy.append((cx, cy))
            labels.append(label)

            for scale, bucket in ((28.0, boundaries), (12.0, nuclei)):
                angles = np.linspace(0, 2 * np.pi, 10, endpoint=False)
                radii = scale * rng.uniform(0.85, 1.15, 10)
                for a, rad in zip(angles, radii):
                    bucket.append(
                        (cid, float(cx + rad * np.cos(a)), float(cy + rad * np.sin(a)))
                    )

    xy = np.asarray(xy)
    counts = _counts_for(labels, rng, len(cell_ids))

    cells = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "x_centroid": xy[:, 0],
            "y_centroid": xy[:, 1],
            "cell_area": rng.uniform(2400, 3200, len(cell_ids)),
            "nucleus_area": rng.uniform(450, 520, len(cell_ids)),
            "transcript_counts": counts.sum(axis=1).astype(int),
            "total_counts": counts.sum(axis=1).astype(int),
            "nucleus_count": rng.integers(30, 60, len(cell_ids)),
            "control_probe_counts": rng.integers(0, 3, len(cell_ids)),
            "genomic_control_counts": rng.integers(0, 2, len(cell_ids)),
            "control_codeword_counts": rng.integers(0, 2, len(cell_ids)),
            "unassigned_codeword_counts": rng.integers(0, 2, len(cell_ids)),
            "segmentation_method": ["multimodal"] * len(cell_ids),
        }
    )
    cells.to_parquet(out / "cells.parquet", index=False)

    for name, bucket in (
        ("cell_boundaries.parquet", boundaries),
        ("nucleus_boundaries.parquet", nuclei),
    ):
        pd.DataFrame(bucket, columns=["cell_id", "vertex_x", "vertex_y"]).to_parquet(
            out / name, index=False
        )

    from scipy import sparse as sp

    _write_10x_h5(
        out / "cell_feature_matrix.h5",
        sp.csr_matrix(counts),
        cell_ids,
        gene_id_prefix="ENSX",
    )

    (out / "experiment.xenium").write_text(
        json.dumps(
            {
                "run_name": "trackcell_tiny_xenium",
                "region": "tiny",
                "panel_tissue_type": "tissue",
                "pixel_size": px,
                "num_cells": len(cell_ids),
            },
            indent=2,
        )
    )

    (out / "gene_panel.json").write_text(
        json.dumps(
            {
                "targets": [
                    {
                        "type": {"descriptor": "gene"},
                        "name": g,
                        "gene_id": f"ENSX{i:06d}",
                    }
                    for i, g in enumerate(GENES)
                ],
            },
            indent=2,
        )
    )
    return {
        "cells": len(cell_ids),
        "genes": len(GENES),
        "sections": int(len(section_origins)),
    }


# ── Registry ─────────────────────────────────────────────────────────────────
BUILDERS = {
    "hd_cellseg_tiny": build_hd_cellseg,
    "hd_bin_tiny": build_hd_bin,
    "xenium_tiny": build_xenium,
}


def _size(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def check() -> list[str]:
    """Read every fixture back through trackcell and assert the contract.

    Returns the list of failure messages (empty on success) so ``pytest`` can
    call this directly and treat the fixtures themselves as tested — otherwise
    they rot silently the first time a reader's expectations change.
    """
    import trackcell as tcl

    failures: list[str] = []

    def expect(cond, message):
        if not cond:
            failures.append(message)
        print(f"  {'PASS' if cond else 'FAIL'}  {message}")

    print("\n[hd_cellseg_tiny]")
    adata = tcl.io.read_hd_cellseg(DATA / "hd_cellseg_tiny")
    expect(adata.n_obs == 36, f"36 cells (got {adata.n_obs})")
    expect("spatial" in adata.obsm, "obsm['spatial'] present")
    expect(
        "geometries" in adata.uns["spatial"][list(adata.uns["spatial"])[0]],
        "uns['spatial'][sample]['geometries'] present",
    )
    expect("classification" in adata.obs, "obs['classification'] present")
    expect("geometry" in adata.obs, "obs['geometry'] WKT present")
    mpp, family = tcl.pl.resolve_microns_per_unit(adata)
    expect(
        abs(mpp - 0.3) < 1e-9 and family == "px",
        f"resolve_microns_per_unit -> (0.3, 'px') got ({mpp}, {family!r})",
    )

    sub = adata[:10].copy()
    tcl.io.sync_geometries_after_subset(sub)
    expect(
        len(sub.uns["spatial"][list(sub.uns["spatial"])[0]]["geometries"]) == 10,
        "sync_geometries_after_subset keeps 10 geometries",
    )

    center = tuple(adata.obsm["spatial"][0])
    mask = tcl.tl.cells_in_radius(adata, center, radius=30, radius_unit="um")
    expect(mask.sum() >= 1, f"cells_in_radius(30um) selects {mask.sum()} cell(s)")

    tcl.tl.hd_labeldist(
        adata, groupby="classification", label="Tumor", coordinate_system="fullres"
    )
    expect("Tumor_dist" in adata.obs, "hd_labeldist writes obs['Tumor_dist'] (µm)")
    expect(adata.obs["Tumor_dist"].notna().all(), "hd_labeldist distances are finite")

    print("\n[hd_bin_tiny]")
    bins = tcl.io.read_hd_bin(DATA / "hd_bin_tiny", binsize=16)
    expect(bins.n_obs == 100, f"100 bins (got {bins.n_obs})")
    expect(bins.obsm["spatial"].shape == (100, 2), "obsm['spatial'] is (100, 2)")

    print("\n[xenium_tiny]")
    xad = tcl.io.read_xenium_cellseg(DATA / "xenium_tiny", verbose=False)
    expect(xad.n_obs == 40, f"40 cells (got {xad.n_obs})")
    expect(
        xad.uns["spatial"][list(xad.uns["spatial"])[0]]["metadata"]["pixel_size"]
        == 0.2125,
        "pixel_size 0.2125 carried in metadata",
    )
    tcl.tl.spatial_slice_cluster(xad, eps=80.0, min_samples=3, min_cells=5)
    n_slices = xad.obs["slice_id"].astype(str).nunique()
    expect(
        n_slices == 2,
        f"spatial_slice_cluster finds 2 sections (got {n_slices}: "
        f"{sorted(xad.obs['slice_id'].astype(str).unique())})",
    )

    print("\n" + "=" * 60)
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("  -", f)
    else:
        print("all fixture checks passed")
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check",
        action="store_true",
        help="read every fixture back through trackcell and verify",
    )
    args = ap.parse_args()

    total = 0
    for name, builder in BUILDERS.items():
        out = DATA / name
        info = builder(out)
        size = _size(out)
        total += size
        print(f"{name:20s} {info}  {size/1024:.1f} KB")
    print(f"{'total':20s} {total/1024:.1f} KB")

    if args.check:
        if check():
            raise SystemExit(1)


if __name__ == "__main__":
    main()
