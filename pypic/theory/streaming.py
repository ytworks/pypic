"""Streaming instabilities: multi-beam electrostatic two-stream / Buneman (k parallel to the
drift) and transverse filamentation (Weibel) of counter-streaming beams.

Electrostatic two-stream / Buneman (k || x, drift along x):

Relativistic warm-fluid dispersion relation (k || x, drift along x):

    1 = sum_s  w_ps^2 / gamma_s^3 / [ (w - k v_s)^2 - 3 k^2 v_ts^2 ]

with w_ps^2 = n_s q_s^2 / m_s (normalised: w_pe = 1) and v_ts^2 = theta_s / gamma_s (fluid
thermal correction; set warm=False for cold beams). Multiplying out gives a polynomial in w
whose complex roots are found with numpy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import polynomial as Pn
from scipy.optimize import minimize_scalar


@dataclass
class Component:
    wp2: float  # omega_p^2 / gamma^3
    v: float
    vt2: float


def components(species, warm: bool = True) -> list[Component]:
    comps: dict = {}
    for s in species:
        for b in s.beams:
            u = b.drift[0]
            g = math.sqrt(1.0 + sum(x * x for x in b.drift))
            v = u / g
            wp2 = b.density * s.charge ** 2 / s.mass / g ** 3
            vt2 = 3.0 * b.theta / g if warm else 0.0
            key = (round(v, 14), round(vt2, 14))
            if key in comps:
                comps[key].wp2 += wp2
            else:
                comps[key] = Component(wp2, v, vt2)
    return list(comps.values())


def roots(comps: list[Component], k: float) -> np.ndarray:
    # factor_s(w) = (w - k v)^2 - k^2 vt2 = w^2 - 2 k v w + k^2 (v^2 - vt2)
    facs = [np.array([k * k * (c.v * c.v - c.vt2), -2 * k * c.v, 1.0]) for c in comps]
    full = np.array([1.0])
    for f in facs:
        full = Pn.polymul(full, f)
    poly = full.copy()
    for i, c in enumerate(comps):
        prod = np.array([1.0])
        for j, f in enumerate(facs):
            if j != i:
                prod = Pn.polymul(prod, f)
        poly = Pn.polysub(poly, c.wp2 * prod)
    return Pn.polyroots(poly)


def growth_rate(comps, k: float) -> tuple[float, float]:
    r = roots(comps, k)
    i = int(np.argmax(r.imag))
    return float(r[i].imag), float(r[i].real)


@dataclass
class LinearResult:
    k_max: float
    gamma_max: float
    omega_r: float
    k_grid: np.ndarray
    gamma_grid: np.ndarray

    @property
    def wavelength(self) -> float:
        return 2.0 * math.pi / self.k_max


def fastest_mode(species, warm: bool = True, kmin: float = 1e-3, kmax: float | None = None,
                 nk: int = 400) -> LinearResult:
    comps = components(species, warm)
    vmax = max(abs(c.v) for c in comps)
    if kmax is None:
        kmax = 10.0 / max(vmax, 1e-3)
    ks = np.geomspace(kmin, kmax, nk)
    gs = np.array([growth_rate(comps, k)[0] for k in ks])
    i = int(np.argmax(gs))
    if gs[i] <= 1e-12:
        raise ValueError("no electrostatic two-stream instability for these parameters")
    lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, nk - 1)]
    res = minimize_scalar(lambda k: -growth_rate(comps, k)[0], bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-10 * hi})
    k = float(res.x)
    g, wr = growth_rate(comps, k)
    return LinearResult(k, g, wr, ks, gs)


def buneman_cold(mass_ratio: float, Z: float = 1.0) -> float:
    """Cold Buneman growth rate (sqrt(3)/2) (Z m_e / 2 m_i)^(1/3) omega_pe."""
    return math.sqrt(3.0) / 2.0 * (Z / (2.0 * mass_ratio)) ** (1.0 / 3.0)


def filamentation_cold(k, u0: float, wp: float = 1.0) -> np.ndarray:
    """Growth rate of the purely transverse (k perpendicular to the drift) filamentation /
    current-filamentation (Weibel) mode of two cold, symmetric, counter-streaming beams of
    proper velocity +-u0 and total plasma frequency wp (each beam wp^2/2; electrons of a pair or
    electron-electron system with the same response). From the cold relativistic fluid theory
    (e.g. Bret, Gremillet & Dieckmann, Phys. Plasmas 17, 120501 (2010)):

        Gamma^2 = [sqrt((k^2 + a)^2 + 4 k^2 beta^2 a) - (k^2 + a)] / 2,   a = wp^2 / gamma0,

    which tends to beta wp / sqrt(gamma0) for k >> wp/c (c = 1)."""
    k = np.asarray(k, dtype=float)
    g0 = math.sqrt(1.0 + u0 * u0)
    beta = u0 / g0
    a = wp * wp / g0
    s = k * k + a
    return np.sqrt(np.maximum(0.5 * (np.sqrt(s * s + 4.0 * k * k * beta * beta * a) - s), 0.0))


def filamentation_estimate(species) -> float:
    """Rough growth rate of the transverse (filamentation/Weibel) mode of counter-streaming
    cold beams, gamma_F ~ beta * omega_p,eff / sqrt(gamma) (for comparison in 3-D runs)."""
    num = 0.0
    beta = 0.0
    for s in species:
        for b in s.beams:
            g = math.sqrt(1.0 + sum(x * x for x in b.drift))
            num += b.density * s.charge ** 2 / s.mass / g
            beta = max(beta, abs(b.drift[0]) / g)
    return beta * math.sqrt(num)
