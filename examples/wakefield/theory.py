"""Fields of the drivers of a plasma wake.

Units: c = eps0 = 1, lengths in c/omega_pe, densities in n0, charges in e.

Gaussian bunch. A rigid bunch of the charge Q with the density

    rho = Q / ((2 pi)^(3/2) sigma_x sigma_r^2) exp(-xi^2 / (2 sigma_x^2) - r^2 / (2 sigma_r^2)),   xi = x - beta t,

has in vacuum the potential phi with (k_perp^2 + k_x^2 / gamma^2) phi = rho, the fields
E_r = -d phi / dr, B_theta = beta E_r. In its rest frame it is a Gaussian at rest with the length
gamma sigma_x. With 1 / |x| = (2 / sqrt(pi)) int_0^inf exp(-s^2 x^2) ds the convolution with a
Gaussian is closed, and in the mid-plane xi = 0 (derived here, no paper was consulted)

    E_r(r) = gamma Q / (pi^(3/2)) r int_0^inf s^2 exp(-s^2 r^2 / (1 + 2 s^2 sigma_r^2))
                                        / ((1 + 2 s^2 sigma_r^2)^2 sqrt(1 + 2 s^2 gamma^2 sigma_x^2)) ds.

Limits: sigma -> 0 and gamma = 1 give Q / (4 pi r^2); gamma sigma_x >> r, sigma_r gives the line
charge, lambda / (2 pi r) (1 - exp(-r^2 / (2 sigma_r^2))) with lambda = Q / (sqrt(2 pi) sigma_x).

Start of the bunch (PIC-D21). The run loads the bunch on top of a partner of the opposite charge at
rest, with E = B = 0, and the bunch moves from t = 0 on. Mode by mode the transverse potential
obeys a'' + w^2 a = J_T exp(-i W t), a(0) = a'(0) = 0, with w = |k| and W = k_x beta; the
longitudinal field follows the charge at once. In the mid-plane of the bunch, relative to the
steady field of the same mode,

    B_theta:  1 - [(1 + W/w) / 2] exp(-i (w - W) t) - [(1 - W/w) / 2] exp(i (w + W) t)
    E_r:      the transverse part beta^2 k_x^2 / k^2 of the steady field with W/w -> w/W in the
              bracket, the longitudinal part unchanged, and the field of the partner, which
              is -(1 / k^2) exp(i W t) against 1 / (k_perp^2 + k_x^2 / gamma^2).

The part of the start that travels with the bunch (small k_perp) loses its phase at the rate
k_perp^2 / (2 k_x) + k_x / (2 gamma^2): the field is there after a time of the order
2 gamma^2 sigma_x, not gamma r.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.special import j1


def bunch_charge(density: float, sigma_x: float, sigma_r: float) -> float:
    """Charge (in e n0 (c/omega_pe)^3, positive number) of a Gaussian bunch with the peak density."""
    return density * (2.0 * math.pi) ** 1.5 * sigma_x * sigma_r ** 2


def bunch_field(r, sigma_x: float, sigma_r: float, gamma: float, charge: float = 1.0):
    """E_r in the mid-plane of the steady Gaussian bunch at the radii r (B_theta = beta E_r)."""
    def integrand(s, rr):
        a = 1.0 + 2.0 * (s * sigma_r) ** 2
        return s * s * math.exp(-(s * rr) ** 2 / a) / (a * a * math.sqrt(1.0 + 2.0 * (s * gamma * sigma_x) ** 2))

    out = [gamma * charge / math.pi ** 1.5 * rr * quad(integrand, 0.0, math.inf, args=(rr,), epsabs=0.0,
                                                    epsrel=1e-10, limit=400)[0] for rr in np.atleast_1d(r)]
    return np.asarray(out)


def bunch_startup(r, t: float, sigma_x: float, sigma_r: float, gamma: float, nk: int = 1500, nx: int = 1201):
    """(E_r, B_theta) in the mid-plane at the time t after the start on the partner, each relative
    to its steady value, at the radii r. Trapezoid rule over k_perp <= 8 / sigma_r and
    |k_x| <= 7 / sigma_x."""
    beta = math.sqrt(1.0 - 1.0 / gamma ** 2)
    kp = np.linspace(1e-4 / sigma_r, 8.0 / sigma_r, nk)[:, None]
    kx = np.linspace(-7.0 / sigma_x, 7.0 / sigma_x, nx)[None, :]
    k2 = kp ** 2 + kx ** 2
    w, W = np.sqrt(k2), kx * beta
    steady = 1.0 / (kp ** 2 + (kx / gamma) ** 2)
    near, far = np.exp(-1j * (w - W) * t), np.exp(1j * (w + W) * t)
    b = steady * (1.0 - 0.5 * (1.0 + W / w) * near - 0.5 * (1.0 - W / w) * far)
    # beta^2 kx^2 / k^2 times the bracket with w / W, written without the division by W
    et = steady * beta / k2 * (0.5 * (W * kx + w * kx) * near + 0.5 * (W * kx - w * kx) * far)
    e = steady - et - np.exp(1j * W * t) / k2
    rho = np.exp(-0.5 * (kp * sigma_r) ** 2 - 0.5 * (kx * sigma_x) ** 2)

    def mid(f, rr):
        return np.trapezoid(np.trapezoid(kp * kp * j1(kp * rr) * rho * f, kx[0], axis=1), kp[:, 0]).real

    out = [(mid(e, rr) / mid(steady, rr), mid(b, rr) / mid(steady, rr)) for rr in np.atleast_1d(r)]
    return np.asarray(out).T
