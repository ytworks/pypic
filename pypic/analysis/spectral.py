"""Spectral analysis: omega-k diagrams, k spectra, growth / damping rates of Fourier modes."""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------- omega - k
def omega_k(f_tx: np.ndarray, dt: float, dx: float, window: bool = True):
    """Power spectrum |F(omega, k)|^2 of a space-time signal f(t, x) (e.g. a lineout series).

    Convention: a wave exp(i(k x - omega t)) appears at (k, +omega); returned for omega >= 0
    with k in increasing order. Returns (k, omega, P[omega, k])."""
    f = np.asarray(f_tx, dtype=float)
    f = f - f.mean(axis=0, keepdims=True)
    nt, nx = f.shape
    if window:
        f = f * np.hanning(nt)[:, None]
    F = np.fft.fft2(f)  # sum f exp(-i(omega t + k x))  -> peak of exp(i(kx - wt)) at (-omega, k)
    F = np.fft.fftshift(F)
    w = np.fft.fftshift(np.fft.fftfreq(nt, d=dt)) * 2 * np.pi
    k = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2 * np.pi
    P = np.abs(F) ** 2
    # flip omega so that exp(i(kx - wt)) sits at positive omega
    P = P[::-1]
    w = -w[::-1]
    sel = w >= 0
    return k, w[sel], P[sel]


def ridge(k, w, P, k_range=None, w_range=None, both_signs=True):
    """Frequency of the strongest response at each k (parabolic peak interpolation).
    both_signs: fold +k and -k (standing / symmetric spectra). Returns (k, omega_peak, peak_power)."""
    k = np.asarray(k)
    if both_signs:
        kk = np.unique(np.abs(k))
        Pk = np.stack([P[:, np.isclose(np.abs(k), v)].sum(axis=1) for v in kk], axis=1)
        k, P = kk, Pk
    wsel = np.ones_like(w, dtype=bool) if w_range is None else (w >= w_range[0]) & (w <= w_range[1])
    ksel = np.ones_like(k, dtype=bool) if k_range is None else (k >= k_range[0]) & (k <= k_range[1])
    ws = w[wsel]
    out_k, out_w, out_p = [], [], []
    for j in np.nonzero(ksel)[0]:
        col = P[wsel, j]
        i = int(np.argmax(col))
        wp = ws[i]
        if 0 < i < len(col) - 1:
            a, b, c = np.log(col[i - 1] + 1e-300), np.log(col[i] + 1e-300), np.log(col[i + 1] + 1e-300)
            den = a - 2 * b + c
            if den < 0:
                wp = ws[i] + 0.5 * (a - c) / den * (ws[1] - ws[0])
        out_k.append(k[j])
        out_w.append(wp)
        out_p.append(col[i])
    return np.array(out_k), np.array(out_w), np.array(out_p)


def branch_ridge(k, w, P, k_th, w_th, rel_window=0.15, both_signs=True, bounds=None):
    """Measured frequency of one theoretical branch: at every k_th (with w_th(k) finite) the
    power maximum within w_th (1 +- rel_window), clipped to ``bounds`` = (w_lo, w_hi) if given
    (e.g. a Bernstein band between cyclotron harmonics), parabolically interpolated. A maximum on
    the window edge means the branch is not resolved there (NaN). Returns (k, w_meas, w_theory)."""
    k = np.asarray(k)
    if both_signs:
        kk = np.unique(np.abs(k))
        P = np.stack([P[:, np.isclose(np.abs(k), v)].sum(axis=1) for v in kk], axis=1)
        k = kk
    out_k, out_w, out_t = [], [], []
    for kt, wt in zip(np.atleast_1d(k_th), np.atleast_1d(w_th)):
        if not np.isfinite(wt) or wt <= 0:
            continue
        j = int(np.argmin(np.abs(k - kt)))
        lo, hi = wt * (1 - rel_window), wt * (1 + rel_window)
        if bounds is not None:
            lo, hi = max(lo, bounds[0]), min(hi, bounds[1])
        sel = np.nonzero((w >= lo) & (w <= hi))[0]
        if sel.size < 3:
            continue
        col = P[sel, j]
        i = int(np.argmax(col))
        wm = np.nan
        if 0 < i < len(col) - 1:
            a, b, c = (np.log(col[i - 1] + 1e-300), np.log(col[i] + 1e-300), np.log(col[i + 1] + 1e-300))
            den = a - 2 * b + c
            wm = w[sel[i]] + (0.5 * (a - c) / den * (w[1] - w[0]) if den < 0 else 0.0)
        out_k.append(k[j])
        out_w.append(wm)
        out_t.append(wt)
    return np.array(out_k), np.array(out_w), np.array(out_t)


def kspectrum_1d(f_x: np.ndarray, dx: float):
    """|F(k)|^2 of a profile (k >= 0)."""
    F = np.fft.rfft(np.asarray(f_x) - np.mean(f_x)) / len(f_x)
    return 2 * np.pi * np.fft.rfftfreq(len(f_x), d=dx), np.abs(F) ** 2


def mode_series(lineout: np.ndarray, m: int) -> np.ndarray:
    """Complex Fourier coefficient of mode m (k = 2 pi m / L) of every profile in a lineout series."""
    a = np.asarray(lineout)
    return np.fft.rfft(a, axis=1)[:, m] / a.shape[1]


# ------------------------------------------------------------ growth fits
def linear_window(t, amp, above_noise=2.0, sat_drop=1.0):
    """Index range of the linear phase: amplitude from `above_noise` e-folds above the initial
    noise up to `sat_drop` e-folds below the first saturation peak."""
    la = np.log(np.maximum(np.abs(np.asarray(amp, dtype=np.complex128)), 1e-300))
    n0 = max(3, len(la) // 20)
    noise = np.median(la[:n0])
    isat = int(np.argmax(la))
    lo, hi = noise + above_noise, la[isat] - sat_drop
    if hi <= lo:
        lo, hi = noise + 1.0, la[isat] - 0.5
    above = np.where(la[: isat + 1] >= lo)[0]
    i0 = int(above[0]) if len(above) else 0
    # end at the first time the amplitude exceeds `hi` (not the last time it is below: after a
    # first saturation the amplitude can dip and grow again in a secondary, nonlinear phase)
    over = np.where(la[i0: isat + 1] > hi)[0]
    i1 = i0 + int(over[0]) - 1 if len(over) and over[0] > 0 else isat
    return i0, i1


def fit_two_wave(ts, zs, gamma0, omega0):
    """Fit z(t) = exp(gamma t) (a exp(-i w t) + b exp(i w t)) to complex samples by variable
    projection (a, b linear; gamma, w nonlinear), relative residuals: weights 1/|z| in float64,
    floored at 1e-6 of the largest |z| so that exact zeros (a standing wave sampled at a node,
    complex64 data) keep a finite weight. Returns (gamma, w > 0, a, b), a and b for the time
    origin t = 0."""
    from scipy.optimize import minimize

    ts = np.asarray(ts, dtype=np.float64)
    zs = np.asarray(zs).astype(np.complex128)
    tt = ts - ts[0]
    amp = np.abs(zs)
    if not amp.max() > 0:
        raise ValueError("fit_two_wave: all samples are zero")
    wgt = 1.0 / np.maximum(amp, 1e-6 * amp.max())

    def solve(pv):
        gg, ww = pv
        M = np.exp(gg * tt)[:, None] * np.stack([np.exp(-1j * ww * tt), np.exp(1j * ww * tt)], 1)
        coef, *_ = np.linalg.lstsq(M * wgt[:, None], zs * wgt, rcond=None)
        return M, coef

    def resid(pv):
        M, coef = solve(pv)
        return float(np.sum(np.abs((M @ coef - zs) * wgt) ** 2))

    w0 = abs(omega0)
    best = min((minimize(resid, [gamma0, w], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-12})
                for w in (w0, 0.8 * w0, 1.25 * w0)), key=lambda r: r.fun)
    g, w = float(best.x[0]), float(best.x[1])
    _, coef = solve(best.x)
    # back from tt = t - t0 to t: exp(g tt) exp(-+i w tt) = exp(g t) exp(-+i w t) exp(-g t0) exp(+-i w t0)
    a = coef[0] * np.exp(-g * ts[0] + 1j * w * ts[0])
    b = coef[1] * np.exp(-g * ts[0] - 1j * w * ts[0])
    if w < 0:  # the same model with w > 0 has the coefficients swapped
        w, a, b = -w, b, a
    return g, w, a, b


def fit_growth(t, z, min_pts=8, omega_r=0.0):
    """Growth rate of one Fourier mode over its linear phase.

    omega_r == 0 (purely growing): straight-line fit of log|z|.
    omega_r != 0 (oscillating): the two counter-propagating growing waves make |z| beat with
    period pi/omega_r, so the complex coefficients are fitted with :func:`fit_two_wave`,
    started from the log fit and the theoretical omega_r. Returns (gamma, (t, fitted envelope))."""
    t = np.asarray(t, dtype=float)
    z = np.asarray(z)
    i0, i1 = linear_window(t, z)
    sel = np.arange(i0, i1 + 1)
    if len(sel) < min_pts:
        return float("nan"), None
    ts, zs = t[sel], z[sel]
    g = float(np.polyfit(ts, np.log(np.abs(zs)), 1)[0])
    if np.iscomplexobj(z) and abs(omega_r) > 1e-3:
        g = fit_two_wave(ts, zs, g, omega_r)[0]
    a = np.abs(zs)
    c0 = np.exp(np.mean(np.log(a)) - g * np.mean(ts))
    return g, (ts, c0 * np.exp(g * ts))


def fit_damped(t, z, omega0, gamma0=-0.05, t_min=None, t_max=None, noise_factor=5.0):
    """Frequency and damping rate of a decaying (standing) mode, e.g. a Landau-damped Langmuir
    wave. z: complex mode coefficients (or real field values). The fit window runs from t_min to
    the time the envelope falls to ``noise_factor`` times the noise floor (median of the last 20 %)
    unless t_max is given. Returns dict(gamma, omega, t_window)."""
    t = np.asarray(t, dtype=float)
    z = np.asarray(z)
    if not np.iscomplexobj(z):
        z = z.astype(complex)
    amp = np.abs(z)
    lo = 0 if t_min is None else int(np.searchsorted(t, t_min))
    if t_max is None:
        floor = np.median(amp[int(0.8 * len(amp)):])
        # running maximum over half a period removes the beating zeros of a standing wave
        nwin = max(1, int(np.pi / max(omega0, 1e-6) / max(t[1] - t[0], 1e-12)))
        env = np.array([amp[i: i + nwin].max() for i in range(len(amp))])
        above = np.nonzero(env[lo:] > noise_factor * floor)[0]
        hi = lo + (int(above[-1]) if len(above) else len(amp) - 1 - lo)
    else:
        hi = int(np.searchsorted(t, t_max))
    sel = slice(lo, max(hi, lo + 8))
    g, w, a, b = fit_two_wave(t[sel], z[sel], gamma0, omega0)
    return {"gamma": g, "omega": w, "t_window": (float(t[sel][0]), float(t[sel][-1])), "a": a, "b": b}
