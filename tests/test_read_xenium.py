"""Regression tests for the Xenium reader's polygon/boundary handling."""

import logging

import h5py
import numpy as np
import pandas as pd
import pytest

from trackcell.io.read_data import (
    _classification_name,
    _infer_sample_from_path,
    read_xenium_cellseg,
)


def _write_xenium_fixture(root):
    """Minimal 2-cell Xenium bundle (no cell_boundaries.parquet)."""
    genes = np.array([b"g1", b"g2", b"g3"])
    cells = np.array([b"c1-1", b"c2-1"])
    with h5py.File(root / "cell_feature_matrix.h5", "w") as handle:
        matrix = handle.create_group("matrix")
        # 10x HDF5 stores CSC (genes x cells): indptr length = n_cells + 1.
        matrix.create_dataset("shape", data=np.array([3, 2]))
        matrix.create_dataset("data", data=np.array([1, 2, 3], dtype=np.int32))
        matrix.create_dataset("indices", data=np.array([0, 1, 0], dtype=np.int32))
        matrix.create_dataset("indptr", data=np.array([0, 2, 3], dtype=np.int32))
        matrix.create_dataset("barcodes", data=cells)
        features = matrix.create_group("features")
        features.create_dataset("name", data=genes)
        features.create_dataset("id", data=np.array([b"E1", b"E2", b"E3"]))
        features.create_dataset("feature_type", data=np.array([b"Gene Expression"] * 3))

    pd.DataFrame(
        {
            "cell_id": ["c1-1", "c2-1"],
            "x_centroid": [1.0, 2.0],
            "y_centroid": [1.0, 2.0],
            "cell_area": [10.0, 10.0],
            "nucleus_area": [5.0, 5.0],
            "transcript_counts": [1, 1],
            "total_counts": [1, 1],
            "nucleus_count": [1, 1],
        }
    ).to_parquet(root / "cells.parquet")

    # Only nucleus boundaries exist: build_polygons=True must not crash.
    pd.DataFrame(
        {
            "cell_id": ["c1-1"] * 3 + ["c2-1"] * 3,
            "vertex_x": [0.0, 1.0, 0.5, 1.5, 2.5, 2.0],
            "vertex_y": [0.0, 0.0, 1.0, 1.5, 1.5, 2.5],
        }
    ).to_parquet(root / "nucleus_boundaries.parquet")


def test_read_xenium_without_cell_boundaries_keeps_nucleus(tmp_path):
    _write_xenium_fixture(tmp_path)

    adata = read_xenium_cellseg(
        tmp_path, sample="sample", verbose=False, include_nucleus=True
    )
    assert adata.shape == (2, 3)
    assert "nucleus_boundaries" in adata.uns
    assert "geometries" not in adata.uns["spatial"]["sample"]


def test_read_xenium_build_polygons_false_with_nucleus(tmp_path):
    _write_xenium_fixture(tmp_path)

    adata = read_xenium_cellseg(
        tmp_path, sample="sample", verbose=False, build_polygons=False
    )
    assert adata.shape == (2, 3)
    assert "cell_boundaries" not in adata.uns
    assert "nucleus_boundaries" in adata.uns


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("/data/Cse1/outs", "Cse1"),
        ("/data/Cse1/outs/binned_outputs", "Cse1"),
        ("/data/Cse1/outs/binned_outputs/square_016um", "Cse1"),
        ("/data/Cse1/outs/segmented_outputs", "Cse1"),
    ],
)
def test_infer_sample_from_spaceranger_layout(path, expected):
    assert _infer_sample_from_path(path) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ({"name": "B cell"}, "B cell"),
        ([{"name": "T cell"}, {"name": "other"}], "T cell"),
        (None, None),
    ],
)
def test_classification_name_handles_dicts_and_lists(value, expected):
    assert _classification_name(value) == expected


def test_read_xenium_rejects_no_matching_cell_ids(tmp_path):
    """Zero overlap between matrix barcodes and cells.parquet must raise.

    The reader used to emit a ``print`` that ``verbose=False`` swallowed, then
    return an AnnData with zero observations as if the load had succeeded.
    """
    _write_xenium_fixture(tmp_path)
    pd.DataFrame(
        {
            "cell_id": ["otherA-1", "otherB-1"],
            "x_centroid": [1.0, 2.0],
            "y_centroid": [1.0, 2.0],
        }
    ).to_parquet(tmp_path / "cells.parquet")

    with pytest.raises(ValueError, match="No cell IDs match"):
        read_xenium_cellseg(tmp_path, sample="sample", verbose=False)


def test_read_xenium_rejects_partial_overlap_below_threshold(tmp_path):
    _write_xenium_fixture(tmp_path)
    # 3 matrix cells, only 1 present in metadata -> 33% match.
    pd.DataFrame(
        {
            "cell_id": ["c1-1"],
            "x_centroid": [1.0],
            "y_centroid": [1.0],
        }
    ).to_parquet(tmp_path / "cells.parquet")

    with pytest.raises(ValueError, match="Refusing to continue below 90%"):
        read_xenium_cellseg(tmp_path, sample="sample", verbose=False)


# ---------------------------------------------------------------------------
# Progress output goes to the logger, not stdout
# ---------------------------------------------------------------------------


def test_progress_does_not_touch_stdout(tmp_path, capsys):
    """Library code must not print.

    ``read_xenium_cellseg`` used ``log = print if verbose else noop``, which
    broke the rule the 0.4.0 changelog set for library code: stdout gets
    polluted in notebooks and cannot be controlled by log level.
    """
    _write_xenium_fixture(tmp_path)

    read_xenium_cellseg(tmp_path, verbose=True)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_progress_reaches_the_logger(tmp_path, caplog):
    """``verbose=True`` must still be observable — through logging."""
    _write_xenium_fixture(tmp_path)

    with caplog.at_level(logging.INFO, logger="trackcell"):
        read_xenium_cellseg(tmp_path, verbose=True)

    messages = [r.getMessage() for r in caplog.records if r.levelno == logging.INFO]
    assert any("[Xenium]" in m for m in messages), messages
    assert any("Expression matrix" in m for m in messages), messages


def test_quiet_by_default_and_with_verbose_false(tmp_path, caplog):
    """The package logger carries a NullHandler, so nothing leaks either way."""
    _write_xenium_fixture(tmp_path)

    with caplog.at_level(logging.INFO, logger="trackcell"):
        read_xenium_cellseg(tmp_path, verbose=False)
        verbose_false = len(caplog.records)

        caplog.clear()
        read_xenium_cellseg(tmp_path, verbose=True)
        verbose_true = len(caplog.records)

    assert verbose_false == 0
    assert verbose_true > 0
