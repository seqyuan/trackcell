"""Tests for read_data helper functions."""

import pandas as pd

from trackcell.io.read_data import convert_classification_to_color_dict


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
