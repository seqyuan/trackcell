"""Guard the module paths that A-1/0.5.0 retired.

``trackcell.tl.sctransform`` and ``trackcell.io.read_atera`` used to be modules
whose names were shadowed by same-named public functions bound in the package
``__init__``. That made ``import trackcell.tl.sctransform as m`` hand back a
*function*, so ``inspect``, ``mock.patch`` and any ``import a.b as c`` reached a
completely different object than they asked for — with an error that pointed
nowhere near the cause.

0.5.0 renamed the files to ``tl/_sct_entry.py`` and ``io/_atera.py``, keeping the
public API identical. The retired paths must now fail loudly instead of silently
resolving to a function. That only stays true if nothing re-introduces them,
which is easy to do by accident: the rename itself broke a *function-local*
import in ``pl/plot.py`` that no top-level grep would find.

This test therefore parses every file's AST and inspects real import statements.
It deliberately does not grep text, so prose in a docstring
(``trackcell/tl/__init__.py`` warns about this exact trap) is not a violation.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCANNED = ("trackcell", "tests", "scripts")

#: This file names the retired paths on purpose — in the assertions that they
#: raise, and in the pattern that recognises them. Scanning itself would only
#: ever report the guard.
SELF = Path(__file__).resolve()

#: Module basenames retired in 0.5.0. Importing a module by these names now
#: raises ModuleNotFoundError, which is the intended behaviour.
RETIRED = ("sctransform", "read_atera", "spatial_cluster")


def _python_files():
    for directory in SCANNED:
        for path in sorted((ROOT / directory).rglob("*.py")):
            if "__pycache__" in path.parts or path.resolve() == SELF:
                continue
            yield path


_RETIRED_TARGET = re.compile(r"^trackcell\.(?:tl\.sctransform|io\.read_atera)(?:\.|$)")


def _retired_references(tree: ast.AST):
    """Yield ``(lineno, what)`` for every reference to a retired module path.

    Two shapes are checked:

    * ``import`` statements, where the dotted name is the module;
    * string constants that *start* with a retired path. ``mock.patch(
      "trackcell.io.read_atera._subset_cells")`` and
      ``importlib.import_module("trackcell.tl.sctransform")`` never appear as
      import statements, so an AST walk over ``Import`` nodes alone missed them.

    ``from a.b import c`` yields ``a.b`` — ``c`` is a name, not a module — which
    is what keeps ``from ._sct_entry import sctransform`` legal. Note that
    ``from trackcell.tl import sctransform`` is also legitimate: it imports the
    *function* bound in the package ``__init__``, which is the recommended form.
    ``fullmatch``-style anchoring keeps docstring prose from tripping the check.
    """
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _is_retired_module(alias.name):
                    yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.module is not None and _is_retired_module(node.module):
                yield node.lineno, node.module
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            if _RETIRED_TARGET.match(node.value):
                yield node.lineno, f"string target {node.value!r}"


def _is_retired_module(dotted: str) -> bool:
    """True when the last dotted component is a retired module basename.

    Leading dots (relative imports) are irrelevant to the basename check:
    ``..io.read_atera`` and ``trackcell.io.read_atera`` are the same file.
    """
    return dotted.lstrip(".").rsplit(".", 1)[-1] in RETIRED


@pytest.mark.parametrize(
    "path", list(_python_files()), ids=lambda p: str(p.relative_to(ROOT))
)
def test_no_import_of_retired_module_path(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    # Every yield is already a violation: import statements are filtered by
    # module basename inside ``_retired_references``, and string targets are
    # matched against the full anchored pattern there. Re-filtering here by
    # basename silently dropped the string cases (``…read_atera._subset_cells``
    # has basename ``_subset_cells``), which is how the gap survived the first
    # version of this guard.
    violations = list(_retired_references(tree))

    assert not violations, (
        f"{path.relative_to(ROOT)} imports a module path retired in 0.5.0: "
        + ", ".join(f"line {n}: {m}" for n, m in violations)
        + ". Those files are now tl/_sct_entry.py and io/_atera.py; the public "
        "functions tcl.tl.sctransform / tcl.io.read_atera are unchanged."
    )


def test_retired_module_paths_fail_loudly():
    """The old paths must raise, not silently resolve to the shadowing function."""
    import importlib

    for dotted in (
        "trackcell.tl.sctransform",
        "trackcell.io.read_atera",
        "trackcell.tl.spatial_cluster",
    ):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(dotted)


def test_public_api_still_binds_functions():
    """The rename is internal: both attributes stay callable functions."""
    import inspect

    import trackcell as tcl

    assert inspect.isfunction(tcl.tl.sctransform)
    assert inspect.isfunction(tcl.io.read_atera)
    # Introspection must reach real source, which the shadowed version could not
    # do when the object was reached through the module path.
    assert "def sctransform" in inspect.getsource(tcl.tl.sctransform)
    assert "def read_atera" in inspect.getsource(tcl.io.read_atera)
