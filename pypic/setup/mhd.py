"""Ideal-MHD scales of a magnetised plasma for the planners of the examples.

Units: c = eps0 = 1, masses in m_e, densities in the reference density (omega_pe = 1). The Alfven
speed uses the ion mass density, v_A = B / sqrt(n m_i); beta = 2 n (T_e + T_i) / B^2.

Perpendicular shock. The Rankine-Hugoniot conditions of an MHD shock with the field perpendicular
to the flow give the density compression X = n2 / n1 as the positive root of (e.g. Priest 1982;
written here with the Alfven Mach number)

    (2 - gamma) X^2 + [gamma (beta + 1) + (gamma - 1) M_A^2] X - (gamma + 1) M_A^2 = 0.

X = 1 at the fast magnetosonic Mach number 1 and X -> (gamma + 1) / (gamma - 1) for a strong shock.
"""
from __future__ import annotations

import math

GAMMA = 5.0 / 3.0


def alfven_speed(b: float, mass_ratio: float, n: float = 1.0) -> float:
    return b / math.sqrt(n * mass_ratio)


def plasma_beta(vte: float, ti_te: float, b: float, n: float = 1.0) -> float:
    """beta of electrons and ions together; vte = sqrt(T_e / m_e c^2)."""
    return 2.0 * n * vte * vte * (1.0 + ti_te) / (b * b)


def fast_speed(v_a: float, beta: float, gamma: float = GAMMA) -> float:
    """Fast magnetosonic speed across the field: sqrt(v_A^2 + c_s^2), c_s^2 = gamma beta v_A^2 / 2."""
    return v_a * math.sqrt(1.0 + 0.5 * gamma * beta)


def fast_mach(mach_alfven: float, beta: float, gamma: float = GAMMA) -> float:
    """Fast magnetosonic Mach number of a flow perpendicular to the field."""
    return mach_alfven / math.sqrt(1.0 + 0.5 * gamma * beta)


def compression(mach_alfven: float, beta: float, gamma: float = GAMMA) -> float:
    """Density compression n2 / n1 of a perpendicular MHD shock (1 if the flow is sub-magnetosonic)."""
    m2 = mach_alfven * mach_alfven
    a = 2.0 - gamma
    b = gamma * (beta + 1.0) + (gamma - 1.0) * m2
    x = (-b + math.sqrt(b * b + 4.0 * a * (gamma + 1.0) * m2)) / (2.0 * a)
    return max(x, 1.0)
