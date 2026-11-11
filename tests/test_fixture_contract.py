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


def _snapshot(root: Path) -> dict[str, tuple[str, object]]:
    """Map relative path to a *content* fingerprint.

    Bytes are only used where the serializer is genuinely stable. Parquet is
    not: pyarrow's writer changes its byte layout between versions, so hashing
    the file turns this test into an assertion about pyarrow rather than about
    the fixtures. The first CI run failed exactly that way — ``hd_bin_tiny`` and
    ``xenium_tiny`` are the two bundles containing ``.parquet``, and they were
    the two reported as stale while the parquet-free bundle passed.

    HDF5 and PNG get the same treatment for the same reason: their encoders are
    versioned too, so their bytes are not a property of this repository.
    """
    import hashlib

    import h5py
    import imageio.v3 as iio
    import numpy as np
    import pandas as pd

    snapshot: dict[str, tuple[str, object]] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = str(path.relative_to(root))

        if path.suffix == ".parquet":
            snapshot[rel] = ("frame", pd.read_parquet(path))
        elif path.suffix == ".png":
            snapshot[rel] = ("image", np.asarray(iio.imread(path)))
        elif path.suffix == ".h5":
            datasets: dict[str, np.ndarray] = {}

            def collect(name, obj):
                if isinstance(obj, h5py.Dataset):
                    datasets[name] = obj[()]

            with h5py.File(path, "r") as handle:
                handle.visititems(collect)
            snapshot[rel] = ("hdf5", datasets)
        else:
            snapshot[rel] = ("bytes", hashlib.sha256(path.read_bytes()).hexdigest())
    return snapshot


def _assert_same_fixture(
    name: str,
    compare: dict[str, tuple[str, object]],
    reference: dict[str, tuple[str, object]],
) -> None:
    """Compare two snapshots and report the first meaningful difference."""
    import numpy as np
    import pandas as pd

    assert set(compare) == set(reference), (
        f"{name}: file set changed; only in one side: "
        f"{sorted(set(compare) ^ set(reference))}"
    )
    for rel in sorted(compare):
        kind_c, value_c = compare[rel]
        kind_r, value_r = reference[rel]
        assert kind_c == kind_r, f"{name}/{rel}: fingerprint kind {kind_c} != {kind_r}"

        if kind_c == "bytes":
            assert value_c == value_r, f"{name}/{rel}: file bytes differ"
        elif kind_c == "frame":
            # ``check_dtype=False``: the contract is the data, not the storage
            # width pyarrow happened to choose.
            pd.testing.assert_frame_equal(value_c, value_r, check_dtype=False)
        elif kind_c == "image":
            assert np.array_equal(value_c, value_r), f"{name}/{rel}: pixels differ"
        else:
            assert set(value_c) == set(
                value_r
            ), f"{name}/{rel}: HDF5 dataset set differs"
            for key in value_c:
                assert np.array_equal(
                    value_c[key], value_r[key]
                ), f"{name}/{rel}: HDF5 dataset {key!r} differs"


def test_fixture_builders_are_deterministic(tmp_path, monkeypatch):
    """Rebuilding must reproduce identical *content*.

    Without this, ``tests/data`` becomes a set of unreviewable binary blobs: a
    regeneration would produce a diff nobody can audit, and the committed copy
    would drift from what the script actually emits.
    """
    make_fixtures = _load_make_fixtures()
    monkeypatch.setattr(make_fixtures, "DATA", tmp_path)

    first: dict[str, dict[str, tuple[str, object]]] = {}
    for name, builder in make_fixtures.BUILDERS.items():
        builder(tmp_path / name)
        first[name] = _snapshot(tmp_path / name)

    for name, builder in make_fixtures.BUILDERS.items():
        builder(tmp_path / name)
        _assert_same_fixture(
            f"{name} (rebuild)", _snapshot(tmp_path / name), first[name]
        )

    # The committed copy must match what the script produces today. This is the
    # assertion that catches a stale fixture after a reader change.
    for name in make_fixtures.BUILDERS:
        committed_dir = FIXTURE_ROOT / name
        if not committed_dir.exists():
            continue
        _assert_same_fixture(
            f"{name} (committed vs regenerated)",
            _snapshot(committed_dir),
            first[name],
        )
