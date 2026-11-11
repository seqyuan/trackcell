"""Smoke-test every rendering path ``spatial_cell`` exposes.

Written while splitting ``spatial_cell`` (backlog C-1). Coverage can tell you a
branch executed; it cannot tell you the branch produced a picture. A plotting
refactor that drops a branch, double-inverts an axis or loses the legend still
"passes" any test that only checks for the absence of an exception — and the
15-scenario fingerprint in ``scripts/render_snapshots.py`` is the tool that
catches it, but its hashes are environment-sensitive and do not belong in CI.

So: this module asserts the *structural* facts that must hold in any environment,
one per scenario, and the fingerprint script remains the tool for before/after
comparison of an actual refactor.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "data" / "hd_cellseg_tiny"

pytestmark = pytest.mark.skipif(
    not FIXTURE.exists(),
    reason="committed fixture missing; run python scripts/make_fixtures.py",
)


def _load_snapshots():
    path = ROOT / "scripts" / "render_snapshots.py"
    spec = importlib.util.spec_from_file_location("_trackcell_render_snapshots", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def snapshots():
    return _load_snapshots()


@pytest.fixture(scope="module")
def fingerprints(snapshots):
    """One render of every scenario, shared across the assertions below."""
    return snapshots.collect()


def test_every_scenario_renders(fingerprints):
    """No scenario may raise — that is the whole point of the sweep."""
    failures = {
        name: info["error"]
        for name, info in fingerprints.items()
        if isinstance(info, dict) and "error" in info
    }
    assert not failures, failures
    assert len(fingerprints) >= 15


def test_categorical_panel_draws_polygons_and_a_legend(fingerprints):
    panels = fingerprints["categorical"]
    assert len(panels) == 1
    panel = panels[0]
    assert panel["collections"] == 1
    assert panel["legend"] == ["Immune", "Stroma", "Tumor"]
    # ``get_aspect()`` returns the numeric form of "equal".
    assert panel["aspect"] == 1.0


def test_continuous_panel_has_no_legend(fingerprints):
    """A gene is continuous; a legend listing its values would be nonsense."""
    panel = fingerprints["continuous"][0]
    assert panel["legend"] is None
    assert panel["collections"] == 1


def test_background_only_panel_has_no_collections(fingerprints):
    """``color=None`` draws the image and nothing else."""
    panel = fingerprints["no_color"][0]
    assert panel["collections"] == 0
    assert panel["xticks_shown"] is True


def test_raw_hex_column_uses_one_collection_with_many_colours(fingerprints):
    """Raw per-cell colours must not be collapsed into a categorical palette.

    One collection carries all facecolors. Grouping by unique colour made a
    typical blend (one colour per cell) create one collection per cell: 600
    cells took ~49 s and 5,000 cells effectively never finished.
    """
    panel = fingerprints["raw_hex"][0]
    assert panel["collections"] == 1
    assert panel["unique_facecolors"] > 3


def test_two_colour_request_produces_two_panels(fingerprints):
    panels = fingerprints["two_panels"]
    assert len(panels) == 2
    assert panels[0]["legend"] == ["Immune", "Stroma", "Tumor"]
    assert panels[1]["legend"] is None


def test_edge_colouring_adds_collections(fingerprints):
    """Fill + edges by a second key must produce more than the fill alone."""
    assert fingerprints["edges_by_key"][0]["collections"] > 1


def test_frameon_false_turns_the_axes_off(fingerprints):
    """``frameon=False`` calls ``ax.axis("off")``, which clears ``ax.axison``.

    It does **not** hide the individual spines: their ``get_visible()`` stays
    true and the frame disappears only because the axes is off. Asserting on the
    spines would report the frame as present.
    """
    assert fingerprints["no_frame"][0]["axis_on"] is False
    assert fingerprints["categorical"][0]["axis_on"] is True


def test_subset_still_renders_after_geometry_sync(fingerprints):
    panels = fingerprints["subset_after_sync"]
    assert isinstance(panels, list)
    assert panels[0]["legend"] == ["Tumor"]


def test_squarebin_paths_render(fingerprints):
    for name in ("squarebin_bin", "squarebin_square_marker"):
        assert fingerprints[name][0]["collections"] == 1
