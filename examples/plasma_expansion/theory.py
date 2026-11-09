"""Isothermal plasma expansion into vacuum (P. Mora, Phys. Rev. Lett. 90, 185002 (2003)).

Normalised units (c = 1, omega_pe = 1 for the initial electron density n_e0 = Z n_i0 = 1,
lengths c/omega_pe): T_e in m_e c^2, ion mass m_i in m_e, charge number Z. The plasma initially
fills x < 0 (edge at x = 0); hot isothermal electrons, cold ions.

    omega_pi = sqrt(Z / m_i),  c_s = sqrt(Z T_e / m_i),  lambda_D0 = sqrt(T_e),
    tau = omega_pi t / sqrt(2 e)                      (e = Euler's number)
    front velocity   v_f = 2 c_s ln(tau + sqrt(tau^2 + 1))
    front position   x_f = 2 sqrt(2 e) lambda_D0 [tau ln(tau + sqrt(tau^2+1)) - sqrt(tau^2+1) + 1]
    maximum energy   E_max = 2 Z T_e [ln(tau + sqrt(tau^2 + 1))]^2
    front field      E_f = 2 E_0 / sqrt(2 e + omega_pi^2 t^2),  E_0 = sqrt(n_e0 T_e) = sqrt(T_e)
    self-similar     n_i = n_i0 exp(-(x / (c_s t) + 1)),  v_i = c_s + x / t   (-c_s t < x < x_f)
    ion spectrum     dN/dE = n_i0 c_s t / sqrt(2 E Z T_e) exp(-sqrt(2 E / (Z T_e)))   (per area, E < E_max)
"""
from __future__ import annotations

import math

import numpy as np

EULER = math.e


def scales(Te: float, mass_ratio: float, Z: float = 1.0) -> dict:
    return {"omega_pi": math.sqrt(Z / mass_ratio), "cs": math.sqrt(Z * Te / mass_ratio), "lambda_D0": math.sqrt(Te),
            "E0": math.sqrt(Te)}


def tau(t, Te, mass_ratio, Z=1.0):
    return scales(Te, mass_ratio, Z)["omega_pi"] * np.asarray(t, dtype=float) / math.sqrt(2 * EULER)


def front_velocity(t, Te, mass_ratio, Z=1.0):
    s = scales(Te, mass_ratio, Z)
    tt = tau(t, Te, mass_ratio, Z)
    return 2 * s["cs"] * np.log(tt + np.sqrt(tt * tt + 1))


def front_position(t, Te, mass_ratio, Z=1.0):
    s = scales(Te, mass_ratio, Z)
    tt = tau(t, Te, mass_ratio, Z)
    q = np.sqrt(tt * tt + 1)
    return 2 * math.sqrt(2 * EULER) * s["lambda_D0"] * (tt * np.log(tt + q) - q + 1)


def max_energy(t, Te, mass_ratio, Z=1.0):
    tt = tau(t, Te, mass_ratio, Z)
    return 2 * Z * Te * np.log(tt + np.sqrt(tt * tt + 1)) ** 2


def front_field(t, Te, mass_ratio, Z=1.0):
    s = scales(Te, mass_ratio, Z)
    return 2 * s["E0"] / np.sqrt(2 * EULER + (s["omega_pi"] * np.asarray(t, dtype=float)) ** 2)


def ion_spectrum(E, t, Te, mass_ratio, Z=1.0, n_i0=None):
    """dN/dE per unit area (normalised: n0 (c/omega_pe) per m_e c^2); zero above E_max."""
    s = scales(Te, mass_ratio, Z)
    ni0 = 1.0 / Z if n_i0 is None else n_i0
    E = np.asarray(E, dtype=float)
    f = ni0 * s["cs"] * t / np.sqrt(2 * E * Z * Te) * np.exp(-np.sqrt(2 * E / (Z * Te)))
    return np.where(E <= max_energy(t, Te, mass_ratio, Z), f, 0.0)


def ion_velocity_profile(x, t, Te, mass_ratio, Z=1.0):
    """Self-similar ion velocity v = c_s + x / t (x measured from the initial edge)."""
    s = scales(Te, mass_ratio, Z)
    return s["cs"] + np.asarray(x, dtype=float) / t
