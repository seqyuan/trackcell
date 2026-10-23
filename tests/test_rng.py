"""Tests for the numpy fallback RNG (used when the optional ``rrng`` is absent)."""

from __future__ import annotations

import numpy as np
import pytest

from trackcell.tl._rng import _NumpyRNG, make_vst_rng


def test_numpy_fallback_normalizes_unnormalized_probabilities():
    """R's ``sample(prob=)`` normalizes internally; numpy requires sum == 1."""
    rng = _NumpyRNG(0)
    items = np.arange(10)
    prob = np.arange(1, 11, dtype=np.float64)  # sums to 55
    out = rng.choice(items, size=500, replace=True, p=prob)
    assert out.shape == (500,)
    assert np.isin(out, items).all()


def test_numpy_fallback_without_probabilities():
    rng = _NumpyRNG(0)
    out = rng.choice(np.arange(5), size=5, replace=False)
    assert sorted(out.tolist()) == [0, 1, 2, 3, 4]


def test_numpy_fallback_rejects_invalid_probabilities():
    rng = _NumpyRNG(0)
    with pytest.raises(ValueError, match="length"):
        rng.choice(np.arange(5), size=1, replace=True, p=np.ones(4))
    with pytest.raises(ValueError, match="non-negative"):
        rng.choice(np.arange(5), size=1, replace=True, p=np.array([-1.0, 1, 1, 1, 1]))
    with pytest.raises(ValueError, match="positive"):
        rng.choice(np.arange(5), size=1, replace=True, p=np.zeros(5))


def test_make_vst_rng_prefers_numpy_when_requested():
    rng = make_vst_rng(0, prefer_r=False)
    assert isinstance(rng, _NumpyRNG)
