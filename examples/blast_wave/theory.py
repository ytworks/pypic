"""Scales of a debris cloud that expands into a magnetised plasma (Le et al. 2021; Keenan et al. 2022).

Units: c = eps0 = 1, lengths in c/omega_pe (d_e), masses in m_e, the ambient density n = 1.

A cloud of N_d ions with the radial speed V_d sweeps up the ambient plasma and pushes the field out
of a cavity. Two lengths follow from its mass and its energy:

    equal-mass radius      (4 pi / 3) R_m^3 n m_i = N_d m_i
    stopping radius        (4 pi / 3) R_B^3 B^2 / 2 = N_d m_i V_d^2 / 2,   R_B = R_m M_A^(2/3)

with M_A = V_d / v_A. The debris ions gyrate with rho_d = V_d / Omega_ci = M_A d_i; the ambient
plasma is coupled to the cloud, and a shock can leave it, when rho_d < R_m.

Extent of the debris (Le et al. 2021, arXiv:2109.00583). Eq. 3 measures it at a time t on the cut
through the centre perpendicular to the field, from the density n_d of the debris ions,

    int_0^D98 n_d dz = 0.98 int_0^inf n_d dz,

and D_max is the largest D98 over time. Eq. 4 fits it for M_A > 1 and rho_d < R_m:

    D_max / R_m ~ (1 + c0 M_A) / (1 + c1 / M_A^2)^(1/D),      c0 = 0.03, c1 = 2,

D the number of dimensions. The denominator is the energy balance (the debris energy against the
field and the swept-up plasma inside the radius), the numerator a correction for the ions that
decouple. c0 and c1 were fitted to 2-D hybrid runs (D = 2) with m_d / m_b = 3, Z_d / Z_b = 1 and a
Gaussian cloud; for 3-D the paper gives no fit and expects coefficients of similar magnitude.

The cloud of this example is a shell, n(r) = n_pk exp(-((r - r_s) / w)^2). Its particle number is
4 pi n_pk sqrt(pi) w (r_s^2 + w^2 / 2) for w << r_s (the part of the Gaussian at r < 0 is dropped).
"""
from __future__ import annotations

import math

LE_C0, LE_C1 = 0.03, 2.0  # Le et al. 2021, Eq. 4 (2-D fit)


def stopping_radius(r_m: float, mach_alfven: float) -> float:
    """R_B: the debris energy equals the energy of the ambient field inside R_B."""
    return r_m * mach_alfven ** (2.0 / 3.0)


def shell_peak_density(r_m: float, r_s: float, w: float) -> float:
    """Peak density of a Gaussian shell (radius r_s, width w) that holds the ambient mass of R_m."""
    return r_m ** 3 / (3.0 * math.sqrt(math.pi) * w * (r_s * r_s + 0.5 * w * w))


def debris_extent(r_m: float, mach_alfven: float, dim: int = 3, c0: float = LE_C0, c1: float = LE_C1) -> float:
    """Largest extent D_max of the debris ions across the field: Eq. 4 of Le et al. 2021 with D = dim.
    The coefficients are the paper's 2-D fit; with dim = 3 this is an extrapolation of that fit."""
    return r_m * (1.0 + c0 * mach_alfven) / (1.0 + c1 / mach_alfven ** 2) ** (1.0 / dim)


def shell_fraction_inside(r_s: float, w: float, radius: float) -> float:
    """Fraction of the particles of the Gaussian shell (radius r_s, width w) closer to its centre than
    ``radius``: the integral of exp(-u^2) (r_s + w u)^2 over u = (r - r_s) / w, in closed form."""
    def part(u):  # antiderivative; it tends to (r_s^2 + w^2 / 2) sqrt(pi) / 2 for u -> inf
        e = math.exp(-u * u)
        return (r_s * r_s + 0.5 * w * w) * 0.5 * math.sqrt(math.pi) * math.erf(u) - (r_s * w + 0.5 * w * w * u) * e

    lo, total = part(-r_s / w), (r_s * r_s + 0.5 * w * w) * 0.5 * math.sqrt(math.pi)
    return (part((radius - r_s) / w) - lo) / (total - lo)
