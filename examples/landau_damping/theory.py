"""Langmuir oscillations: Landau damping (kinetic theory) and the exact cold-plasma solutions.

Landau damping: complex root of the Maxwellian dielectric function, pypic.theory.waves.langmuir
(k lambda_De = 0.5 -> omega = 1.4157 - 0.1534 i omega_pe).

Cold plasma (Dawson 1959): an electron sheet displaced by xi from its equilibrium x0 over a uniform
ion background feels E = n0 e xi / eps0 exactly as long as no two sheets have crossed, so every
sheet is an independent oscillator d(gamma v)/dt = -omega_pe^2 xi. Started at equilibrium with
velocity v(x0) = v0 sin(k x0) (no initial field):
    non-relativistic: x(t) = x0 + (v0 / omega_pe) sin(k x0) sin(omega_pe t)  (exact while no crossing)
    sheets cross (wave breaking) when 1 + (k v0 / omega_pe) cos(k x0) sin(omega_pe t) = 0 first happens:
        omega_pe t_cross = arcsin(omega_pe / (k v0)) for k v0 >= omega_pe, never for k v0 < omega_pe
    relativistic: gamma + omega_pe^2 xi^2 / 2 = gamma_0 (energy invariant), amplitude
        a = sqrt(2 (gamma_0 - 1)) / omega_pe, period T = 4 int_0^a dxi / v(xi)
        = (2 pi / omega_pe) (1 + 3 beta_max^2 / 16 + ...)
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad

from pypic.theory.waves import langmuir


def landau(k_lambda_D: float, vt: float) -> complex:
    """Complex Langmuir frequency at k lambda_De for thermal speed vt (omega_pe = 1)."""
    return complex(langmuir(np.array([k_lambda_D / vt]), vt)[0])


def dawson_orbit(t, x0, v0, k):
    """Non-relativistic cold sheet started at equilibrium x0 with velocity v0 sin(k x0)."""
    return x0 + v0 * np.sin(k * x0) * np.sin(np.asarray(t, dtype=float))


def crossing_time(k_v0: float) -> float:
    """First sheet crossing (wave breaking) of the cold standing oscillation started with
    velocity v0 sin(k x0) (k v0 in units of omega_pe), in 1/omega_pe."""
    return math.asin(1.0 / k_v0) if k_v0 >= 1.0 else math.inf


def amplitude_from_u0(u0: float) -> float:
    """Oscillation amplitude of a relativistic sheet started at equilibrium with proper velocity u0."""
    return math.sqrt(2.0 * (math.sqrt(1.0 + u0 * u0) - 1.0))


def relativistic_period(a: float) -> float:
    """Exact period of a cold relativistic sheet oscillation of amplitude a (c/omega_pe units)."""
    g_m = 1.0 + 0.5 * a * a

    def inv_v(x):
        g = g_m - 0.5 * x * x
        return g / math.sqrt(g * g - 1.0)

    # substitution x = a sin(s) removes the turning-point singularity
    val, _ = quad(lambda s: inv_v(a * math.sin(s)) * a * math.cos(s), 0.0, math.pi / 2, limit=200)
    return 4.0 * val


def beta_max(a: float) -> float:
    g = 1.0 + 0.5 * a * a
    return math.sqrt(1 - 1 / g ** 2)
