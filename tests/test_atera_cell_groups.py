"""Atera's own per-cell clustering (``cell_groups.csv``) must attach cleanly."""

import matplotlib  # noqa: F401  (keeps the backend import order consistent)

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.io.read_atera import _attach_cell_groups


def _adata(barcodes):
    obs = pd.DataFrame(
        {"barcode": list(barcodes)}, index=[f"i{i}" for i in range(len(barcodes))]
    )
    return AnnData(X=np.ones((len(barcodes), 1), dtype=np.float32), obs=obs)


def _write_groups(root, rows):
    pd.DataFrame(rows, columns=["cell_id", "group", "color"]).to_csv(
        root / "cell_groups.csv", index=False
    )


def test_attaches_groups_and_colors(tmp_path):
    adata = _adata(["c0", "c1", "c2"])
    _write_groups(
        tmp_path,
        [
            ("c0", "Tumour", "#ff0000"),
            ("c1", "Stroma", "#00ff00"),
            ("c2", "Tumour", "#ff0000"),
        ],
    )
    assert _attach_cell_groups(adata, tmp_path) is True

    assert isinstance(adata.obs["cell_group"].dtype, pd.CategoricalDtype)
    assert adata.obs["cell_group"].tolist() == ["Tumour", "Stroma", "Tumour"]
    # scanpy / squidpy convention: a list aligned with cat.categories (sorted: Stroma, Tumour)
    assert adata.uns["cell_group_colors"] == ["#00ff00", "#ff0000"]
    assert adata.uns["atera_cell_group_palette"] == {
        "Tumour": "#ff0000",
        "Stroma": "#00ff00",
    }
    meta = adata.uns["atera_cell_groups"]
    assert meta["n_groups"] == 2
    assert meta["n_unassigned"] == 0


def test_unassigned_cells_become_nan_and_are_counted(tmp_path):
    adata = _adata(["c0", "c1"])
    _write_groups(tmp_path, [("c0", "Tumour", "#ff0000")])
    assert _attach_cell_groups(adata, tmp_path) is True

    assert adata.obs["cell_group"].isna().tolist() == [False, True]
    assert adata.uns["atera_cell_groups"]["n_unassigned"] == 1


def test_duplicate_cell_ids_collapse(tmp_path):
    adata = _adata(["c0"])
    _write_groups(
        tmp_path,
        [("c0", "Tumour", "#ff0000"), ("c0", "Tumour", "#ff0000")],
    )
    assert _attach_cell_groups(adata, tmp_path) is True
    assert adata.obs["cell_group"].tolist() == ["Tumour"]


def test_missing_file_is_silent_and_returns_false(tmp_path):
    adata = _adata(["c0"])
    assert _attach_cell_groups(adata, tmp_path) is False
    assert "cell_group" not in adata.obs.columns


def test_foreign_ids_warn_and_do_not_attach(tmp_path):
    adata = _adata(["zzz"])
    _write_groups(tmp_path, [("c0", "Tumour", "#ff0000")])
    with pytest.warns(RuntimeWarning, match="do not match"):
        assert _attach_cell_groups(adata, tmp_path) is False
    assert "cell_group" not in adata.obs.columns


def test_wrong_columns_warn(tmp_path):
    pd.DataFrame({"id": ["c0"], "label": ["Tumour"]}).to_csv(
        tmp_path / "cell_groups.csv", index=False
    )
    adata = _adata(["c0"])
    with pytest.warns(RuntimeWarning, match="lacks"):
        assert _attach_cell_groups(adata, tmp_path) is False


def test_falls_back_to_cell_id_column(tmp_path):
    """The official format stores ``cell_id`` rather than ``barcode``."""
    obs = pd.DataFrame({"cell_id": [11, 12]}, index=["a", "b"])
    adata = AnnData(X=np.ones((2, 1), dtype=np.float32), obs=obs)
    _write_groups(tmp_path, [("11", "A", "#111111"), ("12", "B", "#222222")])
    assert _attach_cell_groups(adata, tmp_path) is True
    assert adata.obs["cell_group"].tolist() == ["A", "B"]


def test_custom_key(tmp_path):
    adata = _adata(["c0"])
    _write_groups(tmp_path, [("c0", "Tumour", "#ff0000")])
    assert _attach_cell_groups(adata, tmp_path, key="my_cluster") is True
    assert adata.obs["my_cluster"].tolist() == ["Tumour"]
    assert "cell_group" not in adata.obs.columns


def test_colors_are_stored_in_scanpy_convention(tmp_path):
    """``sq.pl.spatial_scatter(color=key)`` reads this key positionally.

    Storing ``{group: hex}`` here made squidpy raise
    ``Expected palette to be of length 19, found 1``; storing a list aligned with
    ``cat.categories`` is what both scanpy and squidpy expect.
    """
    adata = _adata(["c0", "c1", "c2"])
    _write_groups(
        tmp_path,
        [("c0", "B", "#222222"), ("c1", "A", "#111111"), ("c2", "B", "#222222")],
    )
    assert _attach_cell_groups(adata, tmp_path) is True

    colors = adata.uns["cell_group_colors"]
    assert isinstance(colors, list)
    assert len(colors) == len(adata.obs["cell_group"].cat.categories)
    # aligned by name with cat.categories (sorted here: A, B), not by file order
    assert list(adata.obs["cell_group"].cat.categories) == ["A", "B"]
    assert colors == ["#111111", "#222222"]
    assert adata.uns["atera_cell_group_palette"] == {"A": "#111111", "B": "#222222"}


def test_group_without_color_gets_gray_and_warns(tmp_path):
    adata = _adata(["c0", "c1"])
    _write_groups(tmp_path, [("c0", "A", "#111111"), ("c1", "B", "")])
    with pytest.warns(RuntimeWarning, match="has no color for 1 of 2 groups"):
        assert _attach_cell_groups(adata, tmp_path) is True

    assert adata.uns["cell_group_colors"] == ["#111111", "#d3d3d3"]
    assert adata.uns["atera_cell_group_palette"] == {"A": "#111111"}
