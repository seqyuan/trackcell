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


def test_expect_version_mismatch_fails(tmp_path: Path, monkeypatch, capsys) -> None:
    """A git tag that disagrees with pyproject must fail the release job.

    This is the only check that can catch a forgotten version bump: the wheel
    and pyproject versions always agree because they come from one source.
    """
    wheel = _write_fake_wheel(tmp_path, version="0.4.0")
    monkeypatch.setattr(check_package, "check_versions", lambda errors: "0.4.0")

    rc = check_package.main(
        ["check_package.py", str(wheel), "--expect-version", "0.5.0"]
    )

    assert rc == 1
    assert "release version '0.5.0'" in capsys.readouterr().out


def test_expect_version_match_passes(tmp_path: Path, monkeypatch) -> None:
    wheel = _write_fake_wheel(tmp_path, version="0.4.0")
    monkeypatch.setattr(check_package, "check_versions", lambda errors: "0.4.0")

    rc = check_package.main(
        ["check_package.py", str(wheel), "--expect-version", "v0.4.0"]
    )

    assert rc == 0


def test_expect_version_requires_a_value(tmp_path: Path) -> None:
    rc = check_package.main(["check_package.py", "--expect-version"])
    assert rc == 2
