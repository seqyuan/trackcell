#!/usr/bin/env python3
"""Validate package metadata before CI/release.

Checks
------
1. ``pyproject.toml`` version matches ``trackcell.__version__``.
2. The built wheel declares every optional extra with its dependencies.
3. The wheel carries the MIT license metadata.
4. Optional: a caller-supplied version (e.g. the git tag that triggered the
   release) matches ``pyproject.toml``.

Check 4 is the one that actually catches a forgotten version bump. Checks 1 and
3 compare artifacts built from the same source, so they agree by construction
and never fire in ``publish-to-pypi``; the tag is the only value that can
disagree.

Usage::

    python scripts/check_package.py [path/to/package.whl] [--expect-version X.Y.Z]

With no wheel argument the most recent wheel in ``dist/`` is used.
"""

from __future__ import annotations

import glob
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# extra -> required substrings in Requires-Dist lines
EXPECTED_EXTRAS: dict[str, tuple[str, ...]] = {
    "clustering": ("leidenalg", "igraph"),
    "integration": ("annoy", "harmonypy", "bbknn"),
    "sct": ("pyglmgampoi", "rrng"),
    "ml": ("xgboost",),
    "interactive": ("ipympl",),
}


def _fail(message: str, errors: list[str]) -> None:
    errors.append(message)
    print(f"ERROR: {message}")


def _metadata_field(lines: list[str], field: str) -> str | None:
    prefix = f"{field}:"
    for line in lines:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return None


def check_versions(errors: list[str]) -> str | None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.MULTILINE)
    if not match:
        _fail("could not read version from pyproject.toml", errors)
        return None
    pyproject_version = match.group(1)

    init = (ROOT / "trackcell" / "__init__.py").read_text(encoding="utf-8")
    init_match = re.search(r'__version__\s*=\s*"([^"]+)"', init)
    if not init_match:
        _fail("could not read __version__ from trackcell/__init__.py", errors)
        return None
    init_version = init_match.group(1)

    if pyproject_version != init_version:
        _fail(
            f"version mismatch: pyproject.toml={pyproject_version} "
            f"!= trackcell/__init__.py={init_version}",
            errors,
        )
    else:
        print(f"version OK: {pyproject_version}")

    return pyproject_version


def check_wheel(
    whl: Path, errors: list[str], expected_version: str | None = None
) -> None:
    with zipfile.ZipFile(whl) as archive:
        metadata_name = next(
            (n for n in archive.namelist() if n.endswith("METADATA")), None
        )
        if metadata_name is None:
            _fail(f"{whl.name} has no METADATA entry", errors)
            return
        metadata = archive.read(metadata_name).decode("utf-8")

    lines = metadata.splitlines()

    wheel_version = _metadata_field(lines, "Version")
    if expected_version is not None:
        if wheel_version == expected_version:
            print(f"wheel version OK: {wheel_version}")
        else:
            _fail(
                f"{whl.name} declares Version {wheel_version!r}, but pyproject.toml "
                f"says {expected_version!r}. Rebuild with `poetry build` "
                f"(or remove stale wheels in dist/).",
                errors,
            )

    provides = {
        line.split(":", 1)[1].strip()
        for line in lines
        if line.startswith("Provides-Extra:")
    }
    missing = set(EXPECTED_EXTRAS) - provides
    if missing:
        _fail(f"{whl.name} is missing extras: {sorted(missing)}", errors)
    else:
        print(f"extras present: {sorted(EXPECTED_EXTRAS)}")

    for extra, packages in EXPECTED_EXTRAS.items():
        deps = [
            line.lower()
            for line in lines
            if line.startswith("Requires-Dist:") and f'extra == "{extra}"' in line
        ]
        for package in packages:
            if not any(package in dep for dep in deps):
                _fail(
                    f"extra '{extra}' does not require '{package}' "
                    f"(Requires-Dist: {deps or 'none'})",
                    errors,
                )

    if not any(line.lower().startswith("license: mit") for line in lines):
        _fail(f"{whl.name} does not declare the MIT license", errors)
    else:
        print("license OK: MIT")


def main(argv: list[str]) -> int:
    args = list(argv[1:])
    expect_version: str | None = None
    if "--expect-version" in args:
        idx = args.index("--expect-version")
        if idx + 1 >= len(args):
            print("ERROR: --expect-version requires a value", file=sys.stderr)
            return 2
        expect_version = args[idx + 1].lstrip("v")
        del args[idx : idx + 2]
    args = [a for a in args if a.strip()]

    errors: list[str] = []
    source_version = check_versions(errors)

    if expect_version is not None and source_version is not None:
        if expect_version == source_version:
            print(f"tag/release version OK: {expect_version}")
        else:
            _fail(
                f"release version {expect_version!r} (git tag) != source version "
                f"{source_version!r}. Bump ``pyproject.toml`` and "
                "``trackcell/__init__.py`` before tagging.",
                errors,
            )

    if len(args) > 0:
        wheels = [Path(args[0])]
    else:
        wheels = [Path(p) for p in sorted(glob.glob(str(ROOT / "dist" / "*.whl")))]
        if not wheels:
            _fail("no wheel found in dist/; run `poetry build` first", errors)
            return 1
        wheels = [max(wheels, key=lambda p: p.stat().st_mtime)]

    for wheel in wheels:
        check_wheel(wheel, errors, expected_version=source_version)

    if errors:
        print(f"\n{len(errors)} packaging check(s) failed.", file=sys.stderr)
        return 1
    print("\nAll packaging checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
