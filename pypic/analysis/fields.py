"""Field analysis helpers: energy densities, magnetic flux function, front tracking."""
from __future__ import annotations

import numpy as np


def energy_density(*components) -> np.ndarray:
    """0.5 * sum(F_i^2) (normalised units: epsilon0 = 1)."""
    return 0.5 * sum(np.asarray(c, dtype=float) ** 2 for c in components)


def flux_function(Bx: np.ndarray, By: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Magnetic flux function A_z of an (x, y) plane with B = curl(A_z z): Bx = dA/dy, By = -dA/dx.
    Integrated along x on the first row and along y elsewhere (exact for the Yee staggering up to
    the half-cell offsets; good enough for X/O-point detection and reconnected flux)."""
    Bx = np.asarray(Bx, dtype=float)
    By = np.asarray(By, dtype=float)
    A = np.zeros_like(Bx)
    A[1:, 0] = -np.cumsum(By[:-1, 0]) * dx
    A[:, 1:] = A[:, :1] + np.cumsum(Bx[:, :-1], axis=1) * dy
    return A


def reconnected_flux(A: np.ndarray, axis: int = 1) -> float:
    """Flux difference between the O-line maximum and the X-line saddle along the current sheet
    (rows along x, sheet at the y index of the extremum): max(A) - min over x of max over y."""
    a = np.asarray(A)
    return float(a.max() - np.min(a.max(axis=axis)))


def front_position(x: np.ndarray, profile: np.ndarray, threshold: float, from_right: bool = True) -> float:
    """Position where a profile first crosses ``threshold`` (searching from the right by default),
    linearly interpolated: shock fronts, expansion fronts."""
    x = np.asarray(x, dtype=float)
    p = np.asarray(profile, dtype=float)
    idx = np.nonzero(p >= threshold)[0]
    if idx.size == 0:
        return float("nan")
    i = int(idx[-1] if from_right else idx[0])
    j = i + 1 if from_right else i - 1
    if 0 <= j < len(p) and p[j] != p[i]:
        return float(x[i] + (threshold - p[i]) * (x[j] - x[i]) / (p[j] - p[i]))
    return float(x[i])


def fit_velocity(t: np.ndarray, x: np.ndarray, t_min=None, t_max=None):
    """Least-squares front velocity over [t_min, t_max]. Returns (v, x0)."""
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    sel = np.isfinite(x)
    if t_min is not None:
        sel &= t >= t_min
    if t_max is not None:
        sel &= t <= t_max
    v, x0 = np.polyfit(t[sel], x[sel], 1)
    return float(v), float(x0)
