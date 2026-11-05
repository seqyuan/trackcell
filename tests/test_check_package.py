"""Regression tests for ``scripts/check_package.py``.

The packaging check must compare the built wheel's ``Version`` metadata
against the version declared in ``pyproject.toml``; a stale wheel left in
``dist/`` used to slip through as a false positive.
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import check_package  # noqa: E402  (path manipulation must come first)

EXTRA_PACKAGES = {
    "clustering": ("leidenalg", "igraph"),
    "integration": ("annoy", "harmonypy", "bbknn"),
    "sct": ("pyglmgampoi", "rrng"),
    "ml": ("xgboost",),
    "interactive": ("ipympl",),
}


def _write_fake_wheel(tmp_path: Path, version: str = "0.0.1") -> Path:
    """Build a minimal wheel whose METADATA satisfies the extras/license checks."""
    lines = [
        "Metadata-Version: 2.1",
        "Name: trackcell",
        f"Version: {version}",
        "License: MIT",
    ]
    for extra, packages in EXTRA_PACKAGES.items():
        lines.append(f"Provides-Extra: {extra}")
        for package in packages:
            lines.append(f'Requires-Dist: {package}; extra == "{extra}"')

    wheel = tmp_path / f"trackcell-{version}-py3-none-any.whl"
    metadata_name = f"trackcell-{version}.dist-info/METADATA"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(metadata_name, "\n".join(lines) + "\n")
    return wheel


def test_stale_wheel_is_rejected(tmp_path: Path) -> None:
    wheel = _write_fake_wheel(tmp_path)
    errors: list[str] = []

    check_package.check_wheel(wheel, errors, expected_version="9.9.9")

    assert any("declares Version" in error for error in errors), errors


def test_matching_wheel_passes(tmp_path: Path) -> None:
    wheel = _write_fake_wheel(tmp_path)
    errors: list[str] = []

    check_package.check_wheel(wheel, errors, expected_version="0.0.1")

    assert errors == []
