"""R-compatible ``sample()`` for SCT v2 step-1 subsampling parity."""

from __future__ import annotations

from typing import Any, Optional, Protocol, runtime_checkable

import numpy as np

try:
    from rrng import RRNG

    _RRNG_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    RRNG = None  # type: ignore[misc, assignment]
    _RRNG_AVAILABLE = False

from ._r_density import r_density_bw_factor  # noqa: F401 — re-export


@runtime_checkable
class VstRNG(Protocol):
    def choice(
        self,
        items: np.ndarray,
        size: int,
        *,
        replace: bool = False,
        p: Optional[np.ndarray] = None,
    ) -> np.ndarray: ...


class _RSampleRNG:
    """Stateful RNG matching ``set.seed()`` + ``sample()`` in R (Mersenne-Twister)."""

    def __init__(self, seed: int) -> None:
        if not _RRNG_AVAILABLE:
            raise ImportError(
                "rrng is required for R-compatible SCT subsampling. "
                "Install with: pip install 'trackcell[sct]' or pip install rrng"
            )
        self._rng = RRNG(int(seed))

    def choice(
        self,
        items: np.ndarray,
        size: int,
        *,
        replace: bool = False,
        p: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        items = np.asarray(items)
        n = int(len(items))
        if size > n and not replace:
            raise ValueError("Cannot take a larger sample than population when replace=False")
        if p is None:
            idx = self._rng.sample(n, int(size), replace=replace)
        else:
            prob = np.asarray(p, dtype=np.float64)
            if prob.shape[0] != n:
                raise ValueError("p length must match items")
            idx = self._rng.sample(n, int(size), replace=replace, prob=prob)
        return items[np.asarray(idx, dtype=int)]


class _NumpyRNG:
    """Fallback when rrng is not installed."""

    def __init__(self, seed: int) -> None:
        self._rng = np.random.default_rng(int(seed))

    def choice(
        self,
        items: np.ndarray,
        size: int,
        *,
        replace: bool = False,
        p: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        return self._rng.choice(np.asarray(items), size=int(size), replace=replace, p=p)


def r_sample_available() -> bool:
    """True when the optional ``rrng`` package is installed."""
    return _RRNG_AVAILABLE


def make_vst_rng(seed: int, *, prefer_r: bool = True) -> VstRNG:
    """
    RNG for ``vst()`` step-1 subsampling.

    Uses R-compatible ``sample()`` when ``rrng`` is installed and ``prefer_r`` is
    True; otherwise falls back to ``numpy.random.Generator``.
    """
    if prefer_r and _RRNG_AVAILABLE:
        return _RSampleRNG(seed)
    return _NumpyRNG(seed)


