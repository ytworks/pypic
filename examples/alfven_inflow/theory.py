"""Theory of an Alfvénic wave carried by a flow: the cold two-fluid parallel R and L waves, the
fluctuation table of an inflow face (PIC-D26), the loaded wave (PIC-D24) and the jump conditions of
the perpendicular shock that a reflecting wall drives into the flow.

Units: c = eps0 = 1, lengths in c/omega_pe, masses in m_e, the density of each species n = 1.

Parallel waves. The plasma (species s: charge q_s, mass m_s, density n_s, at rest) sits in B0 = B y.
A wave E e^{i (k y - w t)} moves every species with the cold response

    -i w v = (q / m) (E + v x B0),

and the plasma current sigma(w) E = sum_s q_s n_s v_s closes Maxwell's equations,

    sigma(w) E = i (w - k^2 / w) E        (E transverse: in the plane (z, x)).

On a circular vector of that plane the two roots below the cyclotron resonances are the R wave (it
rotates with the electrons; the whistler branch, w > k v_A) and the L wave (with the ions; the ion
cyclotron branch, w < k v_A). Both are the textbook relations of Stix,

    k^2 / w^2 = 1 - sum_s w_ps^2 / (w (w + s Omega_s)),    s = +1 (R), -1 (L),  Omega_s = q_s B / m_s,

which tests/cpu check against this solver. For k d_i -> 0 both tend to w = k v_A / sqrt(1 + v_A^2)
(v_A = B / sqrt(sum n m)); to first order in k d_i (d_i = sqrt(m_i) for n = 1) the Hall term splits
them: w / (k v_A) = 1 +/- k d_i / 2 (R: +, L: -).

Box frame. The plasma flows at v_d = -V0 x; the wave has k along y, normal to the flow, so its phase
k y - w t is the same in both frames (to order V0^2). The fields of the box frame are B and
E = E_rest - v_d x B (the motional field of the wave, -V0 dB_z along y, comes in here).

Perpendicular shock. A wall at rest stops the flow; the downstream is at rest in the wall frame and
the shock moves out at V_sh. With the upstream speed U1 = V0 + V_sh in the shock frame, mass
conservation gives r = n2 / n1 = U1 / V_sh, so V_sh = V0 / (r - 1), and r is the compression of a
perpendicular MHD shock at the Alfvén Mach number M_A = U1 / v_A (pypic.setup.mhd.compression,
gamma = 5/3). r -> 4 for M_A -> infinity.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from pypic.setup import mhd

BRANCHES = ("R", "L")


def two_fluid(mass_ratio: float) -> list:
    """Electrons and ions of density 1: (name, charge, mass, density)."""
    return [("e", -1.0, 1.0, 1.0), ("i", 1.0, float(mass_ratio), 1.0)]


def alfven_speed(b: float, species) -> float:
    """B / sqrt(sum n m) (non-relativistic)."""
    return b / math.sqrt(sum(n * m for _, _, m, n in species))


def ion_skin_depth(species) -> float:
    """c / omega_pi of the heaviest species."""
    _, q, m, n = max(species, key=lambda s: s[2])
    return math.sqrt(m / (q * q * n))


def cold_velocity(w: float, q: float, m: float, b: float, E) -> np.ndarray:
    """Cold response v of a species at rest in B0 = b y to the field E e^{i(k y - w t)}."""
    om = q * b / m
    A = np.array([[-1j * w, 0, om], [0, -1j * w, 0], [-om, 0, -1j * w]], dtype=complex)
    return np.linalg.solve(A, (q / m) * np.asarray(E, dtype=complex))


def conductivity(w: float, b: float, species) -> np.ndarray:
    s = np.zeros((3, 3), dtype=complex)
    for _, q, m, n in species:
        for c in range(3):
            e = np.zeros(3)
            e[c] = 1.0
            s[:, c] += q * n * cold_velocity(w, q, m, b, e)
    return s


def _circular(sign: float) -> np.ndarray:
    """A circular unit vector of the plane (z, x): components (z, x) = (1, sign i) / sqrt 2."""
    return np.array([1.0, sign * 1j]) / math.sqrt(2.0)


def _residual(w: float, k: float, b: float, species, sign: float) -> float:
    sub = conductivity(w, b, species)[np.ix_([2, 0], [2, 0])]  # (z, x)
    e = _circular(sign)
    return float(np.imag(e.conj() @ sub @ e - 1j * (w - k * k / w)))


def stix(w: float, k: float, b: float, species, branch: str) -> float:
    """k^2 / w^2 - (1 - sum w_ps^2 / (w (w + s Omega_s))) for the R (s = +1) or L (s = -1) wave."""
    s = 1.0 if branch == "R" else -1.0
    n2 = 1.0 - sum(q * q * n / m / (w * (w + s * q * b / m)) for _, q, m, n in species)
    return k * k / (w * w) - n2


def parallel_roots(k: float, b: float, species) -> dict:
    """Frequencies of the R and L waves at k (below the cyclotron resonances): {"R": (w, sign),
    "L": (w, sign)}, sign the circular vector (z, x) = (1, sign i) of each."""
    w_hi = min(abs(q) * b / m for _, q, m, _ in species if q < 0)  # |Omega_ce|
    roots = []
    for sign in (1.0, -1.0):
        ws = np.geomspace(1e-4 * k * alfven_speed(b, species), w_hi * (1 - 1e-9), 4000)
        f = np.array([_residual(w, k, b, species, sign) for w in ws])
        for a in np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]:
            w = brentq(_residual, ws[a], ws[a + 1], args=(k, b, species, sign), xtol=1e-15, rtol=1e-14)
            if abs(_residual(w, k, b, species, sign)) < 1e-6 * k * k / w:  # a root, not a pole (there: ~1e12)
                roots.append((w, sign))
    if len(roots) != 2 or roots[0][1] == roots[1][1]:
        raise ValueError(f"expected one R and one L root at k = {k:g}, found {roots}")
    roots.sort()
    return {"L": roots[0], "R": roots[1]}


def parallel_wave(k: float, b: float, species, branch: str = "R", amplitude: float = 0.1) -> dict:
    """The cold R or L wave at k along B0 = b y with |delta B| = amplitude * b, in the rest frame of the
    plasma: omega, the complex amplitudes E and B (x, y, z) and the velocity u of every species."""
    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}")
    w, sign = parallel_roots(k, b, species)[branch]
    E = np.array([sign * 1j, 0.0, 1.0])  # (x, y, z)
    Bh = np.cross([0.0, k, 0.0], E) / w
    E = E * amplitude * b / np.linalg.norm(Bh)
    Bh = np.cross([0.0, k, 0.0], E) / w
    u = {name: cold_velocity(w, q, m, b, E) for name, q, m, _ in species}
    return {"omega": float(w), "k": float(k), "E": E, "B": Bh, "u": u, "branch": branch}


def _cos_sin(z) -> dict:
    """Re[z e^{i th}] = Re z cos th - Im z sin th."""
    return {"cos": [float(v) for v in np.real(z)], "sin": [float(v) for v in -np.imag(z)]}


def fluctuation(wave: dict, v0: float, face: str = "x+", phase: float = 0.0) -> dict:
    """``injection.fluctuation`` of PIC-D26 for a plasma that flows at -v0 along x and carries ``wave``
    (from parallel_wave): box-frame E = E_rest - v_d x B."""
    vd = np.array([-v0, 0.0, 0.0])
    Eb = wave["E"] - np.cross(vd, wave["B"])
    mode = {"k": [0.0, wave["k"], 0.0], "omega": wave["omega"], "phase": float(phase), "B": _cos_sin(wave["B"]),
            "E": _cos_sin(Eb), "u": {name: _cos_sin(v) for name, v in wave["u"].items()}}
    return {"face": face, "modes": [mode]}


def incident(table: dict, comp: str, x, y, z, t: float, E0=(0.0, 0.0, 0.0), B0=(0.0, 0.0, 0.0)):
    """The incident field component ``comp`` (Ex..Bz) of the table at the points (x, y, z) (NumPy,
    broadcastable) and time t: the uniform E0 / B0 plus the modes, c cos th + s sin th with
    th = k . r - w t + phase (the engine's incident field, pypic.simulation)."""
    kind, c = comp[0], "xyz".index(comp[1])
    out = np.full(np.broadcast(x, y, z).shape, float((E0 if kind == "E" else B0)[c]))
    for m in table["modes"]:
        k = m["k"]
        th = k[0] * x + k[1] * y + k[2] * z - m["omega"] * t + m.get("phase", 0.0)
        cs = m[kind]
        out = out + cs.get("cos", [0, 0, 0])[c] * np.cos(th) + cs.get("sin", [0, 0, 0])[c] * np.sin(th)
    return out


def _cosine(C: float, S: float):
    """(scale, phase) with C cos(k y) + S sin(k y) = scale cos(k y + phase)."""
    return math.hypot(C, S), math.atan2(-S, C)


def _ycos(ly: float, phase: float) -> dict:
    return {"shape": "cos", "mode": 1, "length": ly, "phase": phase}


def envelope(lx: float, wall: tuple | None, face: float | None, n: int = 32, gap: float = 1.0) -> dict:
    """Envelope h(x) of the loaded wave as a table (sin^2 ramps of n points each): 0 below wall[0],
    rising to 1 at wall[1]; with ``face``, falling to 0 over ``face`` c/omega_pe and exactly 0 over the
    last ``gap`` before x = lx (no displacement may reach the face, PIC-D24)."""
    xs, vs = [0.0], [0.0 if wall else 1.0]

    def ramp(x0, x1, up):
        for s in np.linspace(0.0, 1.0, n):
            xs.append(x0 + s * (x1 - x0))
            v = math.sin(0.5 * math.pi * s) ** 2
            vs.append(v if up else 1.0 - v)

    if wall:
        ramp(wall[0], wall[1], True)
    if face:
        ramp(lx - gap - face, lx - gap, False)
    xs.append(lx)
    vs.append(vs[-1])
    # strictly increasing positions: of equal neighbours (a ramp that starts at 0) the later value stays
    px, pv = [xs[0]], [vs[0]]
    for x, v in zip(xs[1:], vs[1:]):
        if x > px[-1] + 1e-12:
            px.append(float(x))
            pv.append(float(v))
        else:
            pv[-1] = float(v)
    return {"shape": "table", "x": px, "values": pv}


def preload(table: dict, ly: float, dy: float, h: dict | None = None) -> tuple:
    """Initial fields and beam perturbations that load the (single) wave of ``table`` at t = 0.

    Without ``h``: B_x, B_z and E_x as profiles along y (uniform in x), E_z from a displacement of
    the electrons along z, every velocity component as a shaped perturbation along y.
    With the envelope ``h`` (a 1-D shape along x, pypic.profiles), every part is multiplied by h(x)
    (PIC-D24): B from the vector potential, A_z = h a_x(y) for B_x = dA_z/dy and A_x = -h a_z(y) for
    B_z = -dA_x/dy (div B = 0 bitwise; the ramps of h add B_y = -h' a_x); E_x and E_z from
    displacements of the electrons along x and z (xi = -E / (q n), q = -1, n = 1). h must vanish on
    a reflecting wall (the wall zeroes the normal B_x) and on an open face (an A tangential to it and
    a displacement along its normal are refused there).

    In both cases the box-frame E_y of the wave (-V0 dB_z, whose charge density v0 dv / c^2 the
    particles do not carry) is left out; the plasma sets it up."""
    m = table["modes"][0]
    k = m["k"][1]
    gph = float(m.get("phase", 0.0))  # the mode's own phase, th = k y - w t + phase (review of PR #50)
    kd = 2.0 * math.sin(0.5 * k * dy) / dy  # the discrete derivative of the potential along y
    fields, pert = [], {name: [] for name in m["u"]}
    hx = {"x": h} if h else {}
    for c in "xz":
        a = "xyz".index(c)
        R, ph = _cosine(m["B"]["cos"][a], m["B"]["sin"][a])
        ph += gph
        if R == 0.0:
            continue
        if h:
            # B_c = R cos(k y + ph) h(x) from A ∝ sin(k y + ph) = cos(k y + ph - pi/2)
            comp, sgn = ("Az", 1.0) if c == "x" else ("Ax", -1.0)
            fields.append({"kind": "potential", "component": comp, "scale": sgn * R / kd,
                           "shape": {"y": _ycos(ly, ph - 0.5 * math.pi), **hx}})
        else:
            fields.append({"kind": "profile", "component": "B" + c, "axis": "y", "scale": R, **_ycos(ly, ph)})
    for c in "xz":
        a = "xyz".index(c)
        R, ph = _cosine(m["E"]["cos"][a], m["E"]["sin"][a])
        ph += gph
        if R == 0.0:
            continue
        if c == "x" and not h:
            fields.append({"kind": "profile", "component": "Ex", "axis": "y", "scale": R, **_ycos(ly, ph)})
        else:  # xi_c = -E_c / (q n) = E_c for the electrons
            pert["e"].append({"shape": {"y": _ycos(ly, ph), **hx}, "displacement": R, "displacement_component": c})
    for name, u in m["u"].items():
        for a, c in enumerate("xyz"):
            R, ph = _cosine(u["cos"][a], u["sin"][a])
            ph += gph
            if R > 0.0:
                pert[name].append({"shape": {"y": _ycos(ly, ph), **hx}, "velocity": R, "velocity_component": c})
    return fields, pert


def piston_shock(v0: float, v_a: float, beta: float, mass: float = 1.0, gamma: float = mhd.GAMMA) -> dict:
    """The perpendicular shock that a wall drives into a flow at v0 (wall = downstream frame):
    compression r, shock speed V_sh = v0 / (r - 1) and Alfvén Mach number M_A = (v0 + V_sh) / v_A
    of the shock frame, with r = mhd.compression(M_A, beta); and the downstream temperature
    T2 = p2 / n2 (the sum T_e + T_i of an electron-ion plasma) from the momentum balance

        p1 + rho1 U1^2 + B1^2 / 2 = p2 + rho1 U1 U2 + (r B1)^2 / 2,    U2 = U1 / r,

    with B1^2 = rho1 v_A^2, p1 = beta B1^2 / 2 and ``mass`` = rho1 / n1 (m_i + m_e per pair)."""
    if v0 <= 0.0:
        raise ValueError("the flow speed must be positive")

    def f(r):
        return mhd.compression(v0 * r / (r - 1.0) / v_a, beta, gamma) - r

    r_max = (gamma + 1.0) / (gamma - 1.0)
    lo, hi = 1.0 + 1e-3, r_max * (1 - 1e-12)  # M_A from 1e3 v0 / v_A down to 3 v0 / (4 v_A)
    if f(lo) * f(hi) > 0:
        raise ValueError(f"no shock for v0 = {v0:g}, v_A = {v_a:g}, beta = {beta:g}")
    r = brentq(f, lo, hi, xtol=1e-14)
    v_sh = v0 / (r - 1.0)
    u1 = v0 + v_sh
    p2 = 0.5 * beta * v_a ** 2 + u1 * u1 * (1.0 - 1.0 / r) + 0.5 * v_a ** 2 * (1.0 - r * r)  # per rho1
    return {"compression": r, "v_shock": v_sh, "mach_alfven": u1 / v_a,
            "mach_fast": mhd.fast_mach(u1 / v_a, beta, gamma), "temperature": mass * p2 / r}
