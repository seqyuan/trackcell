"""Tests for read_data helper functions."""

import pandas as pd
import pytest

from trackcell.io.read_data import (
    _read_tissue_positions,
    convert_classification_to_color_dict,
)


def test_convert_classification_to_color_dict_accepts_dicts_and_strings():
    df = pd.DataFrame(
        {
            "classification": [
                {"name": "A", "color": [255, 0, 0]},
                "{'name': 'B', 'color': [0, 128, 255]}",
                None,
            ]
        }
    )

    colors = convert_classification_to_color_dict(df)

    assert colors == {"A": "#ff0000", "B": "#0080ff"}


def _write_tissue_positions(tmp_path, sep: str, name: str = "tissue_positions.csv"):
    path = tmp_path / name
    path.write_text(
        sep.join(["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        + "\n"
        + sep.join(["AAAC", "1", "100", "200"])
        + "\n"
        + sep.join(["AAAG", "1", "101", "201"])
        + "\n"
    )
    return path


@pytest.mark.parametrize("sep", [",", "\t"])
def test_read_tissue_positions_handles_csv_and_tsv(tmp_path, sep):
    """A tab-separated table must parse correctly, not become one column.

    ``read_csv(sep=",")`` returns a single column for a TSV instead of raising,
    which made the old ``sep="\\t"`` fallback unreachable and pushed the failure
    into the barcode check ("barcodes are missing from tissue positions").
    """
    path = _write_tissue_positions(tmp_path, sep)
    frame = _read_tissue_positions(path)

    assert frame.shape == (2, 4)
    assert list(frame.columns) == [
        "barcode",
        "in_tissue",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    assert frame["barcode"].tolist() == ["AAAC", "AAAG"]
    assert frame["pxl_col_in_fullres"].tolist() == [200, 201]


def test_read_tissue_positions_raises_on_unparseable_file(tmp_path):
    """A file with no usable columns must say so, not silently pass garbage.

    ``sep=None`` (Sniffer) is the dangerous one: on a header-only barcode column
    it invents a delimiter and splits "barcode" into ``ba`` / ``code``, so any
    check based on column *count* would accept it.
    """
    path = tmp_path / "tissue_positions.csv"
    path.write_text("barcode\nAAAC\nAAAG\n")

    with pytest.raises(ValueError, match="no separator produced a barcode"):
        _read_tissue_positions(path)
