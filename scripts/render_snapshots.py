#!/usr/bin/env python
"""Fingerprint what the plotting functions actually draw.

Built to make the ``spatial_cell`` split reviewable. Coverage cannot tell a
behaviour-preserving refactor from one that quietly drops a branch, and a
plotting function that still renders *something* will pass any test that only
checks "no exception". This renders a fixed set of scenarios and records the
artists produced — collection count, face colours, path counts, axis limits,
label text — so "before" and "after" can be compared mechanically.

Hashes are computed for the **artists**, not for PNG bytes: PNG output depends on
the matplotlib version, the freetype build and the Agg rasteriser, so a byte
comparison would fail on any environment change and be useless in a repository.

Usage
-----
    python scripts/render_snapshots.py > /tmp/before.json
    # ... refactor ...
    python scripts/render_snapshots.py > /tmp/after.json
    python scripts/render_snapshots.py --diff /tmp/before.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import trackcell as tcl  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "data" / "hd_cellseg_tiny"


def _round(values, digits: int = 4):
    return [round(float(v), digits) for v in np.ravel(values)]


def _fingerprint_axes(ax) -> dict:
    """Everything a reader would notice about the rendered panel."""
    collection_faces = []
    collection_face_counts = []
    paths = 0
    collection_style = []
    for coll in ax.collections:
        try:
            face = np.atleast_2d(coll.get_facecolor())
            unique_faces = np.unique(np.round(face, 4), axis=0)
            collection_faces.append(_round(unique_faces))
            collection_face_counts.append(int(len(unique_faces)))
        except Exception:  # noqa: BLE001 - a collection may have no colours
            pass
        try:
            paths += len(coll.get_paths())
        except Exception:  # noqa: BLE001
            pass
        # ``get_facecolor`` does **not** include the collection-level alpha — it is
        # applied at draw time. Without this, changing alpha 0.8 -> 0.79 left every
        # fingerprint identical, which is exactly the silent regression a plotting
        # refactor can introduce.
        try:
            edge = np.atleast_2d(coll.get_edgecolor())
            collection_style.append(
                {
                    "alpha": (
                        None
                        if coll.get_alpha() is None
                        else round(float(coll.get_alpha()), 6)
                    ),
                    "edge": _round(np.unique(np.round(edge, 4), axis=0)),
                    "lw": _round(
                        np.unique(np.round(np.atleast_1d(coll.get_linewidths()), 4))
                    ),
                    "zorder": float(coll.get_zorder()),
                }
            )
        except Exception:  # noqa: BLE001
            collection_style.append({"style": "unavailable"})

    legend = ax.get_legend()
    legend_labels = [t.get_text() for t in legend.get_texts()] if legend else None

    image_shapes = [
        list(im.get_array().shape[:2]) for im in ax.images if im.get_array() is not None
    ]

    # No "spines_visible": ``frameon=False`` goes through ``ax.axis("off")``,
    # which sets ``ax.axison`` and leaves each spine's own ``get_visible()``
    # true, so a per-spine check reports True for an invisible frame. Use
    # ``axis_on`` below.
    return {
        "collections": len(ax.collections),
        "unique_facecolors": sum(collection_face_counts),
        "paths": paths,
        "faces_hash": hashlib.sha256(
            json.dumps(collection_faces, sort_keys=True).encode()
        ).hexdigest()[:16],
        "style_hash": hashlib.sha256(
            json.dumps(collection_style, sort_keys=True).encode()
        ).hexdigest()[:16],
        "images": image_shapes,
        "image_alpha": [
            None if im.get_alpha() is None else round(float(im.get_alpha()), 6)
            for im in ax.images
        ],
        "xlim": _round(ax.get_xlim(), 3),
        "ylim": _round(ax.get_ylim(), 3),
        "xlabel": ax.get_xlabel(),
        "ylabel": ax.get_ylabel(),
        "legend": legend_labels,
        "xticks_shown": bool(any(t.get_visible() for t in ax.get_xticklabels())),
        "aspect": ax.get_aspect(),
        "axis_on": ax.axison,
    }


def _fixture():
    adata = tcl.io.read_hd_cellseg(FIXTURE)
    adata.obs["num"] = np.linspace(0, 1, adata.n_obs)
    return adata


def scenarios() -> dict:
    """Name -> callable rendering one panel, returning the Axes."""
    out = {}

    out["categorical"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", show=False
    )
    out["continuous"] = lambda a: tcl.pl.spatial_cell(a, color="num", show=False)
    out["no_color"] = lambda a: tcl.pl.spatial_cell(a, show=False)
    out["with_image"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", img_key="hires", show=False
    )
    out["no_legend"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", legend=False, show=False
    )
    out["no_frame"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", frameon=False, show=False
    )
    out["ticks_shown"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", show_ticks=True, show=False
    )
    out["two_panels"] = lambda a: tcl.pl.spatial_cell(
        a, color=["classification", "num"], show=False
    )
    out["edges_by_key"] = lambda a: tcl.pl.spatial_cell(
        a, color="num", edge_color="classification", show=False
    )
    out["scalebar"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", scalebar=True, show=False
    )
    out["img_key_none"] = lambda a: tcl.pl.spatial_cell(
        a, color="classification", img_key=None, show=False
    )
    out["squarebin_bin"] = lambda a: tcl.pl.spatial_bin(a, color="num", show=False)

    def _raw_hex(a):
        tcl.tl.multigene_blend(a, genes=["EPCAM", "COL1A1", "PTPRC"], key_added="blend")
        return tcl.pl.spatial_cell(a, color="blend", show=False)

    out["raw_hex"] = _raw_hex

    def _squarebin_markers(a):
        return tcl.pl.spatial_squarebin(a, color="num", shape="square", show=False)

    out["squarebin_square_marker"] = _squarebin_markers
    return out


def subsampled():
    """The fixture cut down, so a subset path is exercised too."""
    adata = _fixture()
    mask = adata.obs["classification"] == "Tumor"
    sub = adata[mask.to_numpy()].copy()
    tcl.io.sync_geometries_after_subset(sub)
    return sub


def collect() -> dict:
    results = {}
    adata = _fixture()
    for name, render in scenarios().items():
        plt.close("all")
        try:
            ax = render(adata)
        except Exception as exc:  # noqa: BLE001 - a failure is a result
            results[name] = {"error": f"{type(exc).__name__}: {exc}"}
            continue
        axes = np.ravel(np.atleast_1d(ax)).tolist()
        results[name] = [
            _fingerprint_axes(a) for a in axes if hasattr(a, "collections")
        ]

    plt.close("all")
    sub = subsampled()
    try:
        ax = tcl.pl.spatial_cell(sub, color="classification", show=False)
        results["subset_after_sync"] = [_fingerprint_axes(ax)]
    except Exception as exc:  # noqa: BLE001
        results["subset_after_sync"] = {"error": f"{type(exc).__name__}: {exc}"}
    plt.close("all")
    return results


def diff(reference_path: Path) -> int:
    reference = json.loads(reference_path.read_text())
    current = collect()

    failures = []
    for name in sorted(set(reference) | set(current)):
        if name not in reference:
            failures.append(f"{name}: new scenario (not in reference)")
            continue
        if name not in current:
            failures.append(f"{name}: scenario disappeared")
            continue
        if reference[name] != current[name]:
            failures.append(
                f"{name}: rendered output changed\n    before={reference[name]}\n    after ={current[name]}"
            )

    if failures:
        print("RENDER DIFFERS:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print(f"all {len(current)} scenarios render identically")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diff", type=Path, help="compare against a saved fingerprint")
    args = ap.parse_args()
    if args.diff:
        return diff(args.diff)
    print(json.dumps(collect(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
