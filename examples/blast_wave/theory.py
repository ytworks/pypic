"""Scales of a debris cloud that expands into a magnetised plasma (Le et al. 2021; Keenan et al. 2021).

Units: c = eps0 = 1, lengths in c/omega_pe (d_e), masses in m_e, the ambient density n = 1.

A cloud of N_d ions with the radial speed V_d sweeps up the ambient plasma and pushes the field out
of a cavity. Two lengths follow from its mass and its energy:

    equal-mass radius      (4 pi / 3) R_m^3 n m_i = N_d m_i
    stopping radius        (4 pi / 3) R_B^3 B^2 / 2 = N_d m_i V_d^2 / 2,   R_B = R_m M_A^(2/3)

with M_A = V_d / v_A. The debris ions gyrate with rho_d = V_d / Omega_ci = M_A d_i; the ambient
plasma is coupled to the cloud, and a shock can leave it, when rho_d < R_m.

The cloud of this example is a shell, n(r) = n_pk exp(-((r - r_s) / w)^2). Its particle number is
4 pi n_pk sqrt(pi) w (r_s^2 + w^2 / 2) for w << r_s (the part of the Gaussian at r < 0 is dropped).
"""
from __future__ import annotations

import math


def stopping_radius(r_m: float, mach_alfven: float) -> float:
    """R_B: the debris energy equals the energy of the ambient field inside R_B."""
    return r_m * mach_alfven ** (2.0 / 3.0)


def shell_peak_density(r_m: float, r_s: float, w: float) -> float:
    """Peak density of a Gaussian shell (radius r_s, width w) that holds the ambient mass of R_m."""
    return r_m ** 3 / (3.0 * math.sqrt(math.pi) * w * (r_s * r_s + 0.5 * w * w))


def fast_speed(v_a: float, beta: float, gamma: float = 5.0 / 3.0) -> float:
    """Fast magnetosonic speed across the field: sqrt(v_A^2 + c_s^2), c_s^2 = gamma beta v_A^2 / 2."""
    return v_a * math.sqrt(1.0 + 0.5 * gamma * beta)
