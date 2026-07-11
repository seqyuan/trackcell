"""R ``stats::density.default`` (FFT, Gaussian kernel) for SCT parity."""

from __future__ import annotations

import numpy as np


def bw_nrd(x: np.ndarray) -> float:
    """R ``stats::bw.nrd``."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 2:
        return float("nan")
    sd = float(np.std(x, ddof=1))
    q75, q25 = np.percentile(x, [75, 25])
    h = float(q75 - q25) / 1.34
    return float(1.06 * min(sd, h) * np.power(n, -0.2))


def r_density_bw_factor() -> float:
    """Legacy ratio ``bw.nrd / bw.nrd0`` (1.06/0.9); kept for callers/tests."""
    return 1.06 / 0.9


def _bin_dist(
    x: np.ndarray,
    weights: np.ndarray,
    lo: float,
    hi: float,
    n: int,
) -> np.ndarray:
    """R ``stats`` ``C_BinDist``: linear binning onto ``2 * n`` points."""
    y = np.zeros(2 * n, dtype=np.float64)
    if n <= 1:
        return y
    xdelta = (hi - lo) / (n - 1)
    ixmin = 0
    ixmax = n - 2
    for xi, wi in zip(x, weights):
        if not np.isfinite(xi):
            continue
        xpos = (xi - lo) / xdelta
        ix = int(np.floor(xpos))
        fx = xpos - ix
        if ixmin <= ix <= ixmax:
            y[ix] += (1.0 - fx) * wi
            y[ix + 1] += fx * wi
        elif ix == -1:
            y[0] += fx * wi
        elif ix == ixmax + 1:
            y[ix] += (1.0 - fx) * wi
    return y


def r_density_fft(
    x: np.ndarray,
    *,
    bw: str = "nrd",
    adjust: float = 1.0,
    n: int = 512,
    cut: float = 3.0,
    ext: float = 4.0,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Port of R ``density.default`` (Gaussian kernel, FFT method).

    Returns ``(x_grid, density, bw)`` matching ``stats::density()``.
    """
    values = np.asarray(x, dtype=np.float64)
    finite = np.isfinite(values)
    values = values[finite]
    nx = values.size
    if nx == 0:
        return np.array([]), np.array([]), float("nan")

    if bw == "nrd":
        bandwidth = bw_nrd(values) * adjust
    elif isinstance(bw, (int, float)):
        bandwidth = float(bw) * adjust
    else:
        raise ValueError(f"unsupported bw rule: {bw!r}")

    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("non-finite or non-positive bandwidth")

    n_grid = max(int(n), 512)
    n_user = int(n)

    weights = np.full(nx, 1.0 / nx, dtype=np.float64)
    tot_mass = 1.0

    from_x = float(np.min(values)) - cut * bandwidth
    to_x = float(np.max(values)) + cut * bandwidth
    lo = from_x - ext * bandwidth
    up = to_x + ext * bandwidth

    y = _bin_dist(values, weights, lo, up, n_grid) * tot_mass

    # Installed R stats (≤4.3) uses legacy kernel coords: ``2 * (up - lo)``.
    kmax = 2.0 * (up - lo)
    kords = np.linspace(0.0, kmax, 2 * n_grid, dtype=np.float64)
    kords[n_grid + 1 : 2 * n_grid] = -kords[n_grid - 1 : 0 : -1]
    kords = np.exp(-0.5 * (kords / bandwidth) ** 2) / (bandwidth * np.sqrt(2.0 * np.pi))

    fy = np.fft.fft(y)
    fk = np.fft.fft(kords)
    conv = np.fft.ifft(fy * np.conj(fk))
    # NumPy ``ifft`` includes 1/N scaling; R divides ``fft(..., inverse=TRUE)`` by N once.
    dens_internal = np.maximum(0.0, np.real(conv[:n_grid]))

    xords = np.linspace(lo, up, n_grid, dtype=np.float64)
    x_out = np.linspace(from_x, to_x, n_user, dtype=np.float64)
    y_out = np.interp(x_out, xords, dens_internal)
    return x_out, y_out, bandwidth


def r_density_interp(
    x: np.ndarray,
    xout: np.ndarray,
    *,
    bw: str = "nrd",
    adjust: float = 1.0,
    n: int = 512,
) -> np.ndarray:
    """R ``approx(density(x)$x, density(x)$y, xout=xout)$y``."""
    grid_x, grid_y, _ = r_density_fft(x, bw=bw, adjust=adjust, n=n)
    xout = np.asarray(xout, dtype=np.float64)
    if grid_x.size == 0:
        return np.full(xout.shape, np.nan, dtype=np.float64)
    return np.interp(
        xout,
        grid_x,
        grid_y,
        left=float(grid_y[0]),
        right=float(grid_y[-1]),
    )
