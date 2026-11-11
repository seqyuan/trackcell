"""Guard the committed fixtures in ``tests/data`` against reader drift.

``scripts/make_fixtures.py`` builds three tiny vendor bundles that exist so an
agent (or a human) can run trackcell code without owning real data. Committed
fixtures rot silently: change what a reader stores, and the files on disk keep
looking fine until someone tries to use them.

This test is the anti-rot mechanism. It reads each committed bundle through the
real trackcell readers and asserts the documented contract — no rebuilding, so
a failure means either the fixtures are stale or a reader changed incompatibly.
Regenerate with::

    python scripts/make_fixtures.py --check
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FIXTURE_ROOT = ROOT / "tests" / "data"

pytestmark = pytest.mark.skipif(
    not (FIXTURE_ROOT / "hd_cellseg_tiny").exists(),
    reason="committed fixtures missing; run python scripts/make_fixtures.py",
)


def _load_make_fixtures():
    """Import ``scripts/make_fixtures.py`` by path.

    ``scripts/`` is not a package, so there is no importable module name.
    """
    path = ROOT / "scripts" / "make_fixtures.py"
    spec = importlib.util.spec_from_file_location("_trackcell_make_fixtures", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixture_contract():
    """Every fixture still satisfies the contract trackcell documents."""
    make_fixtures = _load_make_fixtures()
    assert (
        make_fixtures.DATA == FIXTURE_ROOT
    ), "make_fixtures.DATA and this test disagree about the fixture directory"

    failures = make_fixtures.check()
    assert not failures, "fixture contract drifted:\n  - " + "\n  - ".join(failures)


def test_fixture_builders_are_deterministic(tmp_path, monkeypatch):
    """Rebuilding must reproduce byte-identical files.

    Without this, ``tests/data`` becomes a set of unreviewable binary blobs: a
    regeneration would produce a diff nobody can audit, and the committed copy
    would drift from what the script actually emits.
    """
    import hashlib

    make_fixtures = _load_make_fixtures()
    monkeypatch.setattr(make_fixtures, "DATA", tmp_path)

    def digest(root: Path) -> dict[str, str]:
        return {
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*"))
            if p.is_file()
        }

    first: dict[str, dict[str, str]] = {}
    for name, builder in make_fixtures.BUILDERS.items():
        builder(tmp_path / name)
        first[name] = digest(tmp_path / name)

    for name, builder in make_fixtures.BUILDERS.items():
        builder(tmp_path / name)
        assert (
            digest(tmp_path / name) == first[name]
        ), f"{name} is not reproducible: rebuilding changed file contents"

    # And the committed copy must match what the script produces.
    committed = {
        name: digest(FIXTURE_ROOT / name)
        for name in make_fixtures.BUILDERS
        if (FIXTURE_ROOT / name).exists()
    }
    stale = [name for name in committed if committed[name] != first[name]]
    assert not stale, (
        "tests/data is out of date with scripts/make_fixtures.py for: "
        f"{stale}. Re-run: python scripts/make_fixtures.py"
    )
