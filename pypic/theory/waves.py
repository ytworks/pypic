"""Linear waves of uniform plasmas (normalised units: c = 1, omega_pe = 1 for the reference
density, lengths c/omega_pe, thermal speeds v_t = sqrt(T/m) in c).

Unmagnetised:  Bohm-Gross, kinetic Langmuir (Landau damping) and ion-acoustic roots of the
               Maxwellian dielectric function, EM waves, Yee numerical dispersion.
Magnetised:    cold-plasma branches (Stix R, L, O, X; whistler; upper hybrid) and electron
               Bernstein waves.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import ive, wofz


# ------------------------------------------------------------------ Z function
def Z(zeta):
    """Plasma dispersion function Z(zeta) = i sqrt(pi) w(zeta) (Fried & Conte)."""
    return 1j * math.sqrt(math.pi) * wofz(zeta)


def dZ(zeta):
    """Z'(zeta) = -2 (1 + zeta Z(zeta))."""
    return -2.0 * (1.0 + zeta * Z(zeta))


@dataclass
class Maxwellian:
    """A (non-drifting) Maxwellian component: plasma frequency^2, thermal speed sqrt(T/m)."""

    wp2: float
    vt: float


def epsilon(w, k, comps):
    """Longitudinal dielectric function of Maxwellian components (Landau contour)."""
    eps = 1.0 + 0j
    for c in comps:
        zeta = w / (math.sqrt(2.0) * k * c.vt)
        eps -= c.wp2 / (k * k * c.vt * c.vt) * 0.5 * dZ(zeta)
    return eps


def kinetic_root(k, comps, w0: complex, tol=1e-12, maxit=200):
    """Complex root omega of epsilon(omega, k) = 0 near w0 (secant iteration)."""
    w1 = w0 * (1 + 1e-3) + 1e-4j
    f0, f1 = epsilon(w0, k, comps), epsilon(w1, k, comps)
    for _ in range(maxit):
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0)
        if abs(w2 - w1) < tol * max(1.0, abs(w2)):
            return complex(w2)
        w0, f0, w1, f1 = w1, f1, w2, epsilon(w2, k, comps)
    raise RuntimeError(f"kinetic root did not converge at k={k}")


def bohm_gross(k, vt, wp=1.0):
    return np.sqrt(wp * wp + 3.0 * (np.asarray(k) * vt) ** 2)


def langmuir(k, vt, wp=1.0):
    """Kinetic Langmuir wave of Maxwellian electrons (ions immobile): complex omega (Im < 0 =
    Landau damping). k and vt in normalised units; k lambda_De = k vt / wp."""
    comps = [Maxwellian(wp * wp, vt)]
    out = []
    for kk in np.atleast_1d(k):
        w0 = complex(bohm_gross(kk, vt, wp), -0.01)
        out.append(kinetic_root(kk, comps, w0))
    return np.array(out)


def ion_acoustic(k, Te, Ti, mass_ratio, Z_i=1, wp=1.0):
    """Kinetic ion-acoustic root (electrons + ions, Maxwellian): complex omega."""
    ni = 1.0 / Z_i
    comps = [Maxwellian(wp * wp, math.sqrt(Te)), Maxwellian(wp * wp * ni * Z_i * Z_i / mass_ratio,
                                                             math.sqrt(Ti / mass_ratio))]
    out = []
    for kk in np.atleast_1d(k):
        cs = math.sqrt((Z_i * Te + 3.0 * Ti) / mass_ratio)
        w0 = complex(kk * cs / math.sqrt(1.0 + kk * kk * Te), -0.05 * kk * cs)
        out.append(kinetic_root(kk, comps, w0))
    return np.array(out)


def ion_acoustic_fluid(k, Te, Ti, mass_ratio, Z_i=1):
    cs = math.sqrt((Z_i * Te + 3.0 * Ti) / mass_ratio)
    k = np.asarray(k, dtype=float)
    return k * cs / np.sqrt(1.0 + k * k * Te)


def em_wave(k, wp=1.0):
    """Electromagnetic wave in a cold unmagnetised plasma: omega^2 = omega_p^2 + k^2."""
    return np.sqrt(wp * wp + np.asarray(k) ** 2)


def yee_vacuum(k, dx, dt):
    """Numerical dispersion of the Yee scheme for a wave along a grid axis."""
    return 2.0 / dt * np.arcsin(dt / dx * np.sin(np.asarray(k) * dx / 2.0))


def em_wave_yee(k, dx, dt, wp=1.0):
    """EM wave in a cold plasma with Yee differencing in space and leapfrog in time:
    (2/dt sin(omega dt/2))^2 = omega_p^2 + (2/dx sin(k dx/2))^2."""
    kt = 2.0 / dx * np.sin(np.asarray(k) * dx / 2.0)
    return 2.0 / dt * np.arcsin(np.minimum(dt / 2.0 * np.sqrt(wp * wp + kt * kt), 1.0))


# ---------------------------------------------------------- cold magnetised
def stix(w, wpe2, wce, ions=()):
    """Stix S, D, P for electrons (omega_ce = |e| B / m_e > 0) and ions [(wp2, wc), ...]."""
    w = np.asarray(w, dtype=float)
    R = 1.0 - wpe2 / (w * (w - wce))
    L = 1.0 - wpe2 / (w * (w + wce))
    P = 1.0 - wpe2 / w ** 2
    for wp2, wc in ions:
        R -= wp2 / (w * (w + wc))
        L -= wp2 / (w * (w - wc))
        P -= wp2 / w ** 2
    return R, L, P


def cold_parallel(k, wpe=1.0, wce=1.0):
    """k parallel to B (electrons only): R branches (incl. whistler) and L branch.
    Returns dict name -> omega(k) (NaN where a branch does not exist)."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    out = {"R_upper": [], "whistler": [], "L": []}
    for kk in k:
        # R: (k^2 - w^2)(w - wce) + wpe^2 w = 0  ->  -w^3 + wce w^2 + (k^2 + wpe^2) w - k^2 wce = 0
        rr = np.roots([-1.0, wce, kk * kk + wpe * wpe, -kk * kk * wce])
        rr = np.sort(rr[(np.abs(rr.imag) < 1e-9) & (rr.real > 0)].real)
        out["whistler"].append(rr[0] if len(rr) >= 2 and rr[0] < wce else np.nan)
        out["R_upper"].append(rr[-1] if len(rr) else np.nan)
        # L: (k^2 - w^2)(w + wce) + wpe^2 w = 0  ->  -w^3 - wce w^2 + (k^2 + wpe^2) w + k^2 wce = 0
        rl = np.roots([-1.0, -wce, kk * kk + wpe * wpe, kk * kk * wce])
        rl = rl[(np.abs(rl.imag) < 1e-9) & (rl.real > 0)].real
        out["L"].append(rl.max() if len(rl) else np.nan)
    return {kname: np.array(v) for kname, v in out.items()}


def cold_perpendicular(k, wpe=1.0, wce=1.0):
    """k perpendicular to B (electrons only): O mode (E parallel B) and the two X-mode branches.
    X: (w^2 - k^2)(w^2 - w_UH^2) - wpe^2 (w^2 - wpe^2) = 0, a quadratic in w^2."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    wuh2 = wpe * wpe + wce * wce
    a = 1.0
    b = -(k * k + wuh2 + wpe * wpe)
    c = k * k * wuh2 + wpe ** 4
    disc = np.sqrt(b * b - 4 * a * c)
    return {"O": np.sqrt(wpe * wpe + k * k), "X_upper": np.sqrt((-b + disc) / 2), "X_lower": np.sqrt((-b - disc) / 2),
            "upper_hybrid": np.full_like(k, math.sqrt(wuh2))}


def bernstein(k, vt, wpe=1.0, wce=1.0, n_bands=4, n_terms=40):
    """Electron Bernstein waves (k perpendicular to B, Maxwellian electrons, electrostatic):
    1 = (2 wpe^2 / lambda) exp(-lambda) sum_n n^2 I_n(lambda) / (omega^2 - n^2 wce^2),
    lambda = (k vt / wce)^2. Returns an array (n_bands, len(k)): the root in each band
    (n wce, (n+1) wce), NaN if none."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    out = np.full((n_bands, len(k)), np.nan)
    ns = np.arange(1, n_terms + 1)
    for j, kk in enumerate(k):
        lam = (kk * vt / wce) ** 2
        if lam <= 0:
            continue
        coef = 2.0 * wpe * wpe / lam * ns ** 2 * ive(ns, lam)  # ive = I_n e^{-lambda}

        def eps(w, coef=coef):
            return 1.0 - np.sum(coef / (w * w - (ns * wce) ** 2))

        for b in range(n_bands):
            lo, hi = (b + 1) * wce * (1 + 1e-9), (b + 2) * wce * (1 - 1e-9)
            try:
                if eps(lo) * eps(hi) < 0:
                    out[b, j] = brentq(eps, lo, hi, xtol=1e-12)
            except ValueError:
                pass
    return out
