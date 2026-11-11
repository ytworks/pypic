"""Scales of a magnetised planet in the solar wind: magnetopause distance and bow shock.

Units: c = eps0 = 1, lengths in c/omega_pe (d_e), masses in m_e, the wind density n = 1.

Wind. The Alfven speed of the interplanetary field (IMF) with the ion mass density is
v_A = B / sqrt(n m_i); the wind moves at v = M_A v_A and carries the motional field E = -v x B.

Magnetopause. The wind stops where its ram pressure equals the magnetic pressure of the compressed
dipole field (Chapman-Ferraro; Spreiter et al. 1966):

    (f B_d)^2 / 2 = kappa n m_i v^2,      B_d = M / R^3  (equatorial dipole field at the distance R)

with f = 2 (the field in front of a plane current sheet is twice the dipole field) and kappa = 1
(specular reflection of a cold stream gives 2, a hypersonic gas with gamma = 5/3 gives 0.88).
The planner turns this around: it takes the distance R_mp and sets the moment M = B_d R_mp^3.
R depends on f^2 / kappa only as its 1/6 power.

Bow shock. On the subsolar line the IMF is perpendicular to the flow, and the Rankine-Hugoniot
conditions of a perpendicular MHD shock give the density compression X = n2 / n1 as the positive
root of (e.g. Priest 1982; written here with the Alfven Mach number)

    (2 - gamma) X^2 + [gamma (beta + 1) + (gamma - 1) M_A^2] X - (gamma + 1) M_A^2 = 0,

beta = 2 n (T_e + T_i) / B^2. X = 1 at the fast magnetosonic Mach number 1 and X -> (gamma + 1) /
(gamma - 1) = 4 for a strong shock. The shock stands off the obstacle by the thickness of the
sheath, which gas-dynamic models tie to the compression (Spreiter et al. 1966):

    (R_bs - R_mp) / R_mp = 1.1 / X.

Kinetic condition (known, not analytic): a shock needs a magnetosphere larger than the ion scales,
Dp = R_mp / d_i >~ 1, and is Earth-like from Dp ~ 20 (Omidi et al. 2004); the PIC model of Mercury of
Lavorenti et al. (2022) runs at a reduced Dp ~ 8.
"""
from __future__ import annotations

import math

GAMMA = 5.0 / 3.0


def alfven_speed(b: float, mass_ratio: float, n: float = 1.0) -> float:
    return b / math.sqrt(n * mass_ratio)


def plasma_beta(vte: float, ti_te: float, b: float, n: float = 1.0) -> float:
    """beta of electrons and ions together; vte = sqrt(T_e / m_e c^2)."""
    return 2.0 * n * vte * vte * (1.0 + ti_te) / (b * b)


def fast_mach(mach_alfven: float, beta: float, gamma: float = GAMMA) -> float:
    """Fast magnetosonic Mach number of a flow perpendicular to the field."""
    return mach_alfven / math.sqrt(1.0 + 0.5 * gamma * beta)


def standoff_field(v: float, mass_ratio: float, n: float = 1.0, f: float = 2.0, kappa: float = 1.0) -> float:
    """Undisturbed equatorial dipole field B_d at the magnetopause: (f B_d)^2 / 2 = kappa n m_i v^2."""
    return math.sqrt(2.0 * kappa * n * mass_ratio) * v / f


def dipole_moment(r_mp: float, v: float, mass_ratio: float, **kw) -> float:
    """Moment that puts the magnetopause of the wind (v, n, m_i) at the distance r_mp."""
    return standoff_field(v, mass_ratio, **kw) * r_mp ** 3


def standoff_distance(moment: float, v: float, mass_ratio: float, **kw) -> float:
    return (moment / standoff_field(v, mass_ratio, **kw)) ** (1.0 / 3.0)


def compression(mach_alfven: float, beta: float, gamma: float = GAMMA) -> float:
    """Density compression n2 / n1 of a perpendicular MHD shock (1 if the flow is sub-magnetosonic)."""
    m2 = mach_alfven * mach_alfven
    a = 2.0 - gamma
    b = gamma * (beta + 1.0) + (gamma - 1.0) * m2
    x = (-b + math.sqrt(b * b + 4.0 * a * (gamma + 1.0) * m2)) / (2.0 * a)
    return max(x, 1.0)


def downstream(mach_alfven: float, beta: float, gamma: float = GAMMA) -> dict:
    """State behind the perpendicular shock for an upstream n = B = 1 (so v1 = M_A, p1 = beta / 2):
    v2 = v1 / X and B2 = X from mass and flux conservation, p2 from the momentum flux."""
    x = compression(mach_alfven, beta, gamma)
    v1, p1 = mach_alfven, 0.5 * beta
    return {"compression": x, "v": v1 / x, "B": x, "p": p1 + v1 * v1 * (1.0 - 1.0 / x) + 0.5 * (1.0 - x * x)}


def shock_standoff(r_mp: float, x: float) -> float:
    """Distance of the bow shock from the centre of the obstacle on the subsolar line (Spreiter et al. 1966)."""
    return r_mp * (1.0 + 1.1 / x)


# Region of the undisturbed wind where the run is compared with the wind (gate "upstream"), in R_mp:
# range of x from the upstream face, and distance from the four side faces.
UPSTREAM_X, UPSTREAM_SIDE = (0.2, 0.8), 0.4


def upstream_blocks(n_cells, block, dx, r_mp) -> list:
    """Per axis, the indices of the blocks of the coarse output (``block`` cells of size ``dx`` each)
    that lie wholly inside the upstream region. An empty list: the blocks are too large to measure it."""
    lim = [(UPSTREAM_X[0] * r_mp, UPSTREAM_X[1] * r_mp)]
    lim += [(UPSTREAM_SIDE * r_mp, n_cells[a] * dx[a] - UPSTREAM_SIDE * r_mp) for a in (1, 2)]
    eps, out = 1e-9 * r_mp, []
    for a, (lo, hi) in enumerate(lim):
        w = block[a] * dx[a]
        out.append([k for k in range(-(-n_cells[a] // block[a])) if k * w >= lo - eps and (k + 1) * w <= hi + eps])
    return out
