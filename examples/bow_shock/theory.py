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
conditions of a perpendicular MHD shock give the density compression X = n2 / n1
(``pypic.setup.mhd.compression``): X = 1 at the fast magnetosonic Mach number 1 and X -> 4 for a
strong shock. The shock stands off the obstacle by the thickness of the sheath, which gas-dynamic
models tie to the compression (Spreiter et al. 1966):

    (R_bs - R_mp) / R_mp = 1.1 / X.

Oblique IMF. With the angle theta between the IMF and the shock normal (the flow direction on the
subsolar line), m = M_A^2, c = cos^2 theta, s = sin^2 theta and c_s^2 = gamma beta / 2 (speeds in
units of the Alfven speed of the whole field), the compression of the fast shock is the root of

    (m - c X)^2 [X c_s^2 + m (X (gamma - 1) - (gamma + 1)) / 2]
        + m s X [(gamma + X (2 - gamma)) m - c X ((gamma + 1) - X (gamma - 1))] / 2 = 0

with 1 < X <= (gamma + 1) / (gamma - 1) and m > c X (the flow behind a fast shock stays faster than
the Alfven speed of the normal field). It follows from the conservation of mass, of both components
of the momentum, of the tangential electric field and of the energy across a plane shock (the form
given by Priest 1982, "Solar Magnetohydrodynamics", sec. 5.4, written with the normal Mach number).
At theta = 90 deg it is the quadratic of the perpendicular shock; at theta = 0 it is the gas-dynamic
shock, or the switch-on shock X = m where the gas-dynamic one would leave the flow slower than the
Alfven speed. The tangential field grows by X (m - c) / (m - c X).

Kinetic condition (known, not analytic): a shock needs a magnetosphere larger than the ion scales,
Dp = R_mp / d_i >~ 1, and is Earth-like from Dp ~ 20 (Omidi et al. 2004); the PIC model of Mercury of
Lavorenti et al. (2022) runs at a reduced Dp ~ 8.
"""
from __future__ import annotations

import math

import numpy as np
from numpy.polynomial import Polynomial

from pypic.setup.mhd import GAMMA, compression


def standoff_field(v: float, mass_ratio: float, n: float = 1.0, f: float = 2.0, kappa: float = 1.0) -> float:
    """Undisturbed equatorial dipole field B_d at the magnetopause: (f B_d)^2 / 2 = kappa n m_i v^2."""
    return math.sqrt(2.0 * kappa * n * mass_ratio) * v / f


def dipole_moment(r_mp: float, v: float, mass_ratio: float, **kw) -> float:
    """Moment that puts the magnetopause of the wind (v, n, m_i) at the distance r_mp."""
    return standoff_field(v, mass_ratio, **kw) * r_mp ** 3


def standoff_distance(moment: float, v: float, mass_ratio: float, **kw) -> float:
    return (moment / standoff_field(v, mass_ratio, **kw)) ** (1.0 / 3.0)


def downstream(mach_alfven: float, beta: float, gamma: float = GAMMA) -> dict:
    """State behind the perpendicular shock for an upstream n = B = 1 (so v1 = M_A, p1 = beta / 2):
    v2 = v1 / X and B2 = X from mass and flux conservation, p2 from the momentum flux."""
    x = compression(mach_alfven, beta, gamma)
    v1, p1 = mach_alfven, 0.5 * beta
    return {"compression": x, "v": v1 / x, "B": x, "p": p1 + v1 * v1 * (1.0 - 1.0 / x) + 0.5 * (1.0 - x * x)}


def _parallel_shock(m: float, beta: float, gamma: float) -> tuple:
    """Compression and tangential field (over the upstream field) behind a shock along the field,
    m = M_A^2: the gas-dynamic shock while the flow behind it stays faster than the Alfven speed
    (X < m), the switch-on shock X = m otherwise; (1, 0) without a shock."""
    x_gas = (gamma + 1.0) * m / ((gamma - 1.0) * m + gamma * beta)
    if x_gas <= 1.0 or m <= 1.0:
        return 1.0, 0.0
    if x_gas < m:
        return x_gas, 0.0
    return m, math.sqrt((m - 1.0) * ((gamma + 1.0) - gamma * beta - (gamma - 1.0) * m))


def oblique_compression(mach_alfven: float, beta: float, cos2: float, gamma: float = GAMMA) -> float:
    """Density compression of the fast MHD shock whose normal makes the angle theta with the upstream
    field, cos2 = cos^2 theta; mach_alfven is the normal flow speed over the Alfven speed of the
    whole field. 1 when there is no fast shock."""
    if cos2 == 0.0:
        return compression(mach_alfven, beta, gamma)
    m, c, s, cs2 = mach_alfven ** 2, cos2, 1.0 - cos2, 0.5 * gamma * beta
    if s == 0.0:  # the cubic has a double root at X = m there
        return _parallel_shock(m, beta, gamma)[0]
    x = Polynomial([0.0, 1.0])
    poly = ((m - c * x) ** 2 * (x * cs2 + 0.5 * m * (x * (gamma - 1.0) - (gamma + 1.0)))
            + 0.5 * m * s * x * ((gamma + x * (2.0 - gamma)) * m - c * x * ((gamma + 1.0) - x * (gamma - 1.0))))
    roots = poly.roots()
    real = roots[np.abs(roots.imag) <= 1e-6 * np.maximum(1.0, np.abs(roots.real))].real
    x_max = (gamma + 1.0) / (gamma - 1.0)
    fast = real[(real > 1.0 + 1e-12) & (real <= x_max * (1.0 + 1e-12)) & (m - c * real > 0.0)]
    return float(fast.max()) if fast.size else 1.0


def oblique_downstream(mach_alfven: float, beta: float, cos2: float, gamma: float = GAMMA) -> dict:
    """State behind the oblique fast shock for an upstream n = |B| = 1, the normal along x and the
    tangential field along z (v1 = M_A along x, p1 = beta / 2): v_x and B_z from the mass flux and
    the tangential electric field, v_z from the tangential momentum, p from the normal momentum."""
    x = oblique_compression(mach_alfven, beta, cos2, gamma)
    v1, bx, bz1 = mach_alfven, math.sqrt(cos2), math.sqrt(1.0 - cos2)
    m = v1 * v1
    bz = _parallel_shock(m, beta, gamma)[1] if bz1 == 0.0 else bz1 * x * (m - cos2) / (m - cos2 * x)
    return {"compression": x, "v_x": v1 / x, "v_z": bx * (bz - bz1) / v1, "B_x": bx, "B_z": bz,
            "p": 0.5 * beta + m * (1.0 - 1.0 / x) + 0.5 * (bz1 * bz1 - bz * bz)}


def shock_standoff(r_mp: float, x: float) -> float:
    """Distance of the bow shock from the centre of the obstacle on the subsolar line (Spreiter et al. 1966)."""
    return r_mp * (1.0 + 1.1 / x)


# Region of the undisturbed wind where the run is compared with the wind (gate "upstream"), in R_mp:
# range of x from the upstream face, and distance from the four side faces.
UPSTREAM_X, UPSTREAM_SIDE = (0.2, 0.8), 0.4


def upstream_blocks(n_cells, block, dx, r_mp, z_range=None) -> list:
    """Per axis, the indices of the blocks of the coarse output (``block`` cells of size ``dx`` each)
    that lie wholly inside the upstream region. An empty list: the blocks are too large to measure it.
    ``z_range = (low, high)`` narrows the region along z (an oblique IMF: the side without foreshock)."""
    lim = [(UPSTREAM_X[0] * r_mp, UPSTREAM_X[1] * r_mp)]
    lim += [(UPSTREAM_SIDE * r_mp, n_cells[a] * dx[a] - UPSTREAM_SIDE * r_mp) for a in (1, 2)]
    if z_range is not None:
        lim[2] = (max(lim[2][0], z_range[0]), min(lim[2][1], z_range[1]))
    eps, out = 1e-9 * r_mp, []
    for a, (lo, hi) in enumerate(lim):
        w = block[a] * dx[a]
        out.append([k for k in range(-(-n_cells[a] // block[a])) if k * w >= lo - eps and (k + 1) * w <= hi + eps])
    return out
