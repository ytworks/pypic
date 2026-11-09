"""Particle analysis for acceleration studies: energy spectra, power-law fits, thermal fits,
maximum energy, and work done by each field component along tracked orbits."""
from __future__ import annotations

import numpy as np
from scipy.special import kve


# ---------------------------------------------------------------- spectra
def spectrum(edges: np.ndarray, counts: np.ndarray):
    """Differential spectrum dN/dE from histogram counts on (possibly logarithmic) bins.
    Returns (bin centres (geometric for log bins), dN/dE)."""
    e = np.asarray(edges, dtype=float)
    c = np.asarray(counts, dtype=float)
    width = np.diff(e)
    log = e[0] > 0 and np.allclose(e[2:] / e[1:-1], e[1] / e[0], rtol=1e-6) if len(e) > 2 else False
    centres = np.sqrt(e[1:] * e[:-1]) if log else 0.5 * (e[1:] + e[:-1])
    return centres, c / width


def local_slope(E, f):
    """d ln f / d ln E (central differences); a power law f ~ E^-p has slope -p."""
    lE, lf = np.log(E), np.log(np.maximum(f, 1e-300))
    return np.gradient(lf, lE)


def fit_power_law(E, f, lo=None, hi=None, min_decades=0.5, max_scatter=0.15):
    """Power-law index p of f ~ E^-p.

    With lo/hi given, a least-squares fit of ln f vs ln E on [lo, hi]. Otherwise the fit range is
    the widest interval (at least ``min_decades`` decades) where the local slope stays within
    ``max_scatter`` of its median, searched above the spectral peak. Returns dict(p, p_err, lo, hi)."""
    E = np.asarray(E, dtype=float)
    f = np.asarray(f, dtype=float)
    ok = (f > 0) & np.isfinite(f) & (E > 0)
    E, f = E[ok], f[ok]
    if lo is None or hi is None:
        s = local_slope(E, f)
        ipk = int(np.argmax(E * f))  # peak of E dN/dE
        best = None
        for i in range(ipk + 1, len(E)):
            for j in range(len(E) - 1, i + 2, -1):
                if np.log10(E[j] / E[i]) < min_decades:
                    break
                seg = s[i: j + 1]
                med = np.median(seg)
                if np.all(np.abs(seg - med) <= max_scatter * max(abs(med), 1.0)):
                    if best is None or np.log10(E[j] / E[i]) > np.log10(best[1] / best[0]):
                        best = (E[i], E[j])
                    break
        if best is None:
            return {"p": float("nan"), "p_err": float("nan"), "lo": float("nan"), "hi": float("nan")}
        lo, hi = best
    sel = (E >= lo) & (E <= hi)
    x, y = np.log(E[sel]), np.log(f[sel])
    A = np.stack([x, np.ones_like(x)], 1)
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    n = len(x)
    sigma2 = float(res[0]) / max(n - 2, 1) if len(res) else 0.0
    cov = sigma2 * np.linalg.inv(A.T @ A)
    return {"p": float(-coef[0]), "p_err": float(np.sqrt(cov[0, 0])), "lo": float(lo), "hi": float(hi),
            "norm": float(np.exp(coef[1]))}


def maxwell_juttner(gamma, theta):
    """Normalised Maxwell-Juttner distribution dN/dgamma (3-D) at temperature theta = T / m c^2."""
    g = np.asarray(gamma, dtype=float)
    b = np.sqrt(np.maximum(1.0 - 1.0 / g ** 2, 0.0))
    # exp(-g/theta) / (theta K2(1/theta)) with the scaled Bessel function kve(2, x) = K2(x) e^x
    return g ** 2 * b * np.exp(-(g - 1.0) / theta) / (theta * kve(2, 1.0 / theta))


def fit_temperature(gamma_centres, dNdg, n_total=None):
    """Temperature of the thermal part of a spectrum: matches the Maxwell-Juttner shape below the
    peak region (least squares in log space around the peak). Returns (theta, nonthermal_fraction)."""
    from scipy.optimize import minimize_scalar

    g = np.asarray(gamma_centres, dtype=float)
    f = np.asarray(dNdg, dtype=float)
    dg = np.gradient(g)
    N = n_total if n_total is not None else float(np.sum(f * dg))
    ipk = int(np.argmax(f))
    sel = (np.arange(len(g)) <= ipk + max(3, ipk // 2)) & (f > 0)

    def cost(lt):
        th = np.exp(lt)
        m = N * maxwell_juttner(g[sel], th)
        return float(np.sum((np.log(f[sel]) - np.log(np.maximum(m, 1e-300))) ** 2))

    r = minimize_scalar(cost, bounds=(np.log(1e-5), np.log(1e3)), method="bounded")
    th = float(np.exp(r.x))
    thermal = N * maxwell_juttner(g, th)
    excess = np.clip(f - thermal, 0, None)
    return th, float(np.sum(excess * dg) / max(N, 1e-300))


def energy_cutoff(E, f, fraction=1e-3):
    """Energy above which the integral N(>E) falls below ``fraction`` of the total."""
    E = np.asarray(E, dtype=float)
    f = np.asarray(f, dtype=float)
    dE = np.gradient(E)
    n = f * dE
    above = np.cumsum(n[::-1])[::-1]
    tot = above[0] if above.size else 0.0
    idx = np.nonzero(above >= fraction * tot)[0]
    return float(E[idx[-1]]) if idx.size else float("nan")


# --------------------------------------------------------------- tracks
def work_by_component(tracks, charge: float, mass: float):
    """Energy gained by each tracked particle from each E component, W_i(t) = int q E_i v_i dt
    (fields at the particle, v from the recorded momenta; explicit runs record u at t - dt/2,
    this is the time-centred estimate). Returns dict t, W_x, W_y, W_z (n_t, n), gamma (n_t, n)."""
    t = tracks.t
    u = np.stack([tracks["ux"], tracks["uy"], tracks["uz"]])
    g = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    v = u / g
    E = np.stack([tracks["Ex"], tracks["Ey"], tracks["Ez"]])
    p = charge * E * v / mass  # d gamma / dt contributions
    dt = np.diff(t)[:, None]
    W = np.zeros_like(p)
    W[:, 1:] = np.nancumsum(0.5 * (p[:, 1:] + p[:, :-1]) * dt[None], axis=1)
    return {"t": t, "W_x": W[0], "W_y": W[1], "W_z": W[2], "gamma": g}


def parallel_perpendicular_work(tracks, charge: float, mass: float):
    """Work split into E parallel and perpendicular to the local B (reconnection / shock studies)."""
    t = tracks.t
    u = np.stack([tracks["ux"], tracks["uy"], tracks["uz"]])
    g = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    v = u / g
    E = np.stack([tracks["Ex"], tracks["Ey"], tracks["Ez"]])
    B = np.stack([tracks["Bx"], tracks["By"], tracks["Bz"]])
    b = B / np.maximum(np.sqrt(np.sum(B ** 2, axis=0)), 1e-300)
    Epar = np.sum(E * b, axis=0) * b
    Eperp = E - Epar
    dt = np.diff(t)[:, None]
    out = {"t": t, "gamma": g}
    for k, Ek in (("W_par", Epar), ("W_perp", Eperp)):
        p = charge * np.sum(Ek * v, axis=0) / mass
        W = np.zeros_like(p)
        W[1:] = np.nancumsum(0.5 * (p[1:] + p[:-1]) * dt, axis=0)
        out[k] = W
    return out
