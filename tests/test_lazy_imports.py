"""Guard the lazy import graph.

``trackcell/__init__.py`` and the ``io`` / ``tl`` package inits resolve their
public names through PEP 562 ``__getattr__`` instead of importing everything up
front. That bought ``import trackcell`` going from 2.43 s to ~0 s, and — more
importantly — stopped ``matplotlib.pyplot`` from being imported, which locks the
backend and breaks ``matplotlib.use("Agg")`` on headless machines.

Lazy loading introduces two failure modes that do not exist in eager code, and
neither shows up as a normal test failure:

1. A name in ``__all__`` with no entry in the mapping. Eager imports would fail at
   import time; lazy ones only fail when someone touches that name.
2. A submodule whose name equals a public name. Importing the submodule (from
   anywhere) sets it as an attribute on the package, which *shadows*
   ``__getattr__`` — the caller then gets a module where a function was promised.
   That is the same defect A-1 fixed for ``sctransform``, and it resurfaced here
   for ``spatial_cluster`` the moment the eager binding stopped masking it.

The subprocess checks assert the contract that actually matters: in a *fresh*
interpreter, importing the core package must not load the heavy dependencies.
"""

from __future__ import annotations

import json
import subprocess
import sys
import types
from pathlib import Path

import pytest

import trackcell as tcl

ROOT = Path(__file__).resolve().parent.parent

PROBE = r"""
import json, sys
import {module}
print("@RESULT@" + json.dumps({{
    "loaded": [m for m in {watch!r} if m in sys.modules],
}}))
"""


def _fresh_import(module: str, watch: tuple[str, ...]) -> list[str]:
    """Import ``module`` in a new interpreter and report which of ``watch`` loaded."""
    code = PROBE.format(module=module, watch=sorted(watch))
    proc = subprocess.run(
        [sys.executable, "-c", code], cwd=str(ROOT), capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
    for line in proc.stdout.splitlines():
        if line.startswith("@RESULT@"):
            return json.loads(line[len("@RESULT@") :])["loaded"]
    raise AssertionError(f"probe for {module} produced no result")


# ── the contract: core imports stay light ───────────────────────────────────


@pytest.mark.parametrize("module", ["trackcell", "trackcell.tl", "trackcell.io"])
def test_core_imports_do_not_load_scanpy_or_pyplot(module):
    """The whole point of the change, asserted in a fresh interpreter.

    ``matplotlib.pyplot`` locks the backend on import, so a caller that does
    ``matplotlib.use("Agg")`` *after* ``import trackcell`` silently keeps the
    default backend. That is not a speed issue; it is a correctness trap on any
    headless machine.
    """
    loaded = _fresh_import(module, ("matplotlib.pyplot", "scanpy"))

    assert loaded == [], f"{module} pulled {loaded}"


def test_public_names_still_resolve():
    """Lazy failures are silent until touched, so touch every name."""
    unresolved = []
    for namespace in ("io", "tl", "pl"):
        module = getattr(tcl, namespace)
        for name in module.__all__:
            try:
                getattr(module, name)
            except Exception as exc:  # noqa: BLE001 - reporting, not handling
                unresolved.append(f"{namespace}.{name}: {type(exc).__name__}: {exc}")
    assert not unresolved, "unresolved public names:\n  " + "\n  ".join(unresolved)


@pytest.mark.parametrize("namespace", ["io", "tl"])
def test_lazy_mapping_covers_all(namespace):
    """A missing mapping entry is an AttributeError only at first use."""
    module = getattr(tcl, namespace)

    assert set(module._LAZY_EXPORTS) == set(module.__all__)
    for name in module.__all__:
        # ``dir()`` must advertise them, or tab-completion and ``help`` lie.
        assert name in dir(module)
        assert not isinstance(getattr(module, name), types.ModuleType), name


def test_cached_after_first_access():
    """``__getattr__`` caches into globals, so the indirection is paid once."""
    module = sys.modules["trackcell.tl"]

    assert "sctransform" in vars(module), "name was not cached after access"
    assert module.__getattr__("sctransform") is module.sctransform


# ── aliasing must survive the indirection ───────────────────────────────────


def test_io_tl_alias_identity_holds():
    """Both namespaces must hand back the same function object.

    ``tl`` re-exports two QC helpers from ``io``; with two independent lazy
    resolvers it would be easy to end up with two distinct wrappers.
    """
    assert tcl.tl.atera_cell_qc is tcl.io.atera_cell_qc
    assert tcl.tl.filter_atera_genes is tcl.io.filter_atera_genes


# ── the shadowing hazard ────────────────────────────────────────────────────


def test_submodule_import_cannot_shadow_a_public_function():
    """Importing ``tl._cluster`` must not replace the ``spatial_cluster`` function.

    Before the rename, ``trackcell/tl/spatial_cluster.py`` shared its name with
    the public function ``tl.spatial_cluster``. With eager binding the function
    was bound last and won; with PEP 562 the import system sets the module
    attribute directly, ``__getattr__`` is never consulted, and callers got a
    module. Assert on the *function* while the module is definitely imported.
    """
    import trackcell.tl._cluster  # noqa: F401  (deliberately load the submodule)

    assert not isinstance(tcl.tl.spatial_cluster, types.ModuleType)
    assert callable(tcl.tl.spatial_cluster)
    assert tcl.tl.spatial_cluster.__module__.endswith("._cluster")
