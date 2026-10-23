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

from ._density import r_density_bw_factor  # noqa: F401 — re-export


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
        items = np.asarray(items)
        if p is None:
            return self._rng.choice(items, size=int(size), replace=replace)

        prob = np.asarray(p, dtype=np.float64)
        if prob.shape[0] != len(items):
            raise ValueError("p length must match items")
        if not np.isfinite(prob).all() or (prob < 0).any():
            raise ValueError("p must contain finite non-negative values")
        total = prob.sum()
        if total <= 0:
            raise ValueError("p must sum to a positive value")
        # R's sample() normalizes probabilities internally; numpy.Generator
        # requires a sum of 1 (within tolerance), so normalize explicitly.
        prob = prob / total
        return self._rng.choice(items, size=int(size), replace=replace, p=prob)


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


