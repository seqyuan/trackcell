# Contributing to TrackCell

Thanks for helping improve TrackCell. This project is a solo-maintained
library, so small, focused, well-tested pull requests are easiest to review and
merge.

## Development setup

```bash
git clone https://github.com/seqyuan/trackcell.git
cd trackcell
poetry install --with dev --extras "clustering integration"
```

Python 3.10–3.12 is supported. Heavy parity/benchmark tests look for data under
`TRACKCELL_BENCHMARK_ROOT` (default `~/trackcell_benchmark`) and are skipped
when it is absent.

## Running tests

```bash
# Fast suite (skips tests needing external R/Seurat benchmark data)
poetry run pytest -m "not slow" -q

# Full suite, including golden-vector parity tests
TRACKCELL_BENCHMARK_ROOT=/path/to/tcl_test poetry run pytest
```

New behavior needs a test. IO readers are tested with synthetic vendor bundles
built in `tests/conftest.py`; prefer extending those fixtures over shipping real
data.

## Code style

```bash
poetry run black --check trackcell tests scripts   # formatting (line length 88)
poetry run flake8 trackcell                        # config in .flake8
poetry run mypy trackcell                          # scoped type checks
```

Conventions worth knowing:

- Library code **never** uses `print()` and never configures logging. Emit
  progress with `logger = trackcell._logging.get_logger(__name__)` at
  `INFO`; report non-fatal data problems with `warnings.warn(..., RuntimeWarning, stacklevel=2)`.
- Do not silently swallow exceptions: log at `debug` (optional work) or
  `warning` (visible failure) instead of `except Exception: pass`.

`pre-commit` runs black/flake8/autoflake on staged files; install the hooks with
`pre-commit install` after `pip install pre-commit`.

## Pull requests

- Keep the change focused; one concern per PR.
- Add or update tests, and update `docs/source/changelog.rst`.
- Make sure `black --check`, `flake8`, and the fast test suite pass.

## Releasing

Release steps (version bump in both `trackcell/__init__.py` and
`pyproject.toml`, changelog, tag, push) are documented in `AGENTS.md`. Releases
are published to PyPI by the `auto-release` GitHub Actions workflow when a `v*`
tag is pushed.
