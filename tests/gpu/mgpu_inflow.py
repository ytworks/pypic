"""An inflow face that carries a wave in (PIC-D26), and the helper for its multi-GPU test
(``python -m pypic.launch -n N mgpu_inflow.py [piston]``).

Electrons and ions (m_i = 16) flow at -v0 along x through B0 = B y (E0 = v0 B z). The face x+ is a
paired inflow face that carries one circularly polarized R wave along B0 (k = k y), the exact
mode of the cold two-fluid plasma: its frequency is the root of the parallel dispersion relation
with the displacement current, its particle velocities the cold response. Each species enters
through a flux face with its own flux (a paired face ties the electrons' birth rate to the ions'
flux, which differs from theirs in a wave by the wave current); the Marder zone takes the charge that
the entering particles leave on the face nodes. The box starts without
the wave unless ``preload``: then the box holds the same wave at t = 0 (B and E_x as profiles along y,
E_z from a displacement of the electrons along z, the velocities as shaped perturbations; the charge
of E_y, of order v0 dv / c^2, is left out), and the face has to keep it.
- ``open`` (default): x- is an open face where the plasma leaves.
- ``piston``: x- is a reflecting wall; the flow piles up into a shock that moves into the plasma
  that carries the wave.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import CurrentFilter, Injection

NX, NY, NZ, DX, DT = 256, 128, 2, 0.25, 0.1
LX, LY = NX * DX, NY * DX
MI, B, V0, AMP = 16.0, 0.5, 0.05, 0.2  # AMP: |delta B| / B
THETA = 1e-3  # T_e / m c^2: lambda_De = 0.13 = dx / 2
K = 2 * math.pi / LY
SPECIES = (("e", -1.0, 1.0), ("i", 1.0, MI))  # name, charge, mass; density 1 each


def _velocity(w, q, m, E):
    """Cold response v of a species at rest in B0 = B y to the field E e^{i(k y - w t)}."""
    om = q * B / m
    A = np.array([[-1j * w, 0, om], [0, -1j * w, 0], [-om, 0, -1j * w]], dtype=complex)
    return np.linalg.solve(A, (q / m) * np.asarray(E, dtype=complex))


def _sigma(w):
    s = np.zeros((3, 3), dtype=complex)
    for _, q, m in SPECIES:
        for c in range(3):
            e = np.zeros(3)
            e[c] = 1.0
            s[:, c] += q * _velocity(w, q, m, e)
    return s


def r_mode():
    """(omega, E-hat) of the R wave (rotating with the electrons) at k: sigma E = i (w - k^2 / w) E
    on the plane (z, x) normal to B0 = B y."""
    def branch(w, sign):
        sub = _sigma(w)[np.ix_([2, 0], [2, 0])]  # (z, x)
        e = np.array([1.0, sign * 1j]) / math.sqrt(2)  # a circular vector of the plane
        lam = (e.conj() @ sub @ e)
        return float(np.imag(lam - 1j * (w - K * K / w)))

    va = B / math.sqrt(MI)
    out = None
    for sign in (1.0, -1.0):
        ws = np.linspace(0.3 * K * va, 5 * K * va, 400)
        f = [branch(w, sign) for w in ws]
        for a in range(len(ws) - 1):
            if f[a] * f[a + 1] < 0:
                w = brentq(branch, ws[a], ws[a + 1], args=(sign,), xtol=1e-14)
                if out is None or w > out[0]:
                    out = (w, sign)
    w, sign = out
    Ez, Ex = 1.0, sign * 1j
    E = np.array([Ex, 0.0, Ez])
    Bh = np.cross([0.0, K, 0.0], E) / w
    E *= AMP * B / np.linalg.norm(Bh)
    return w, E


def fluctuation():
    w, Er = r_mode()
    k = np.array([0.0, K, 0.0])
    Bh = np.cross(k, Er) / w
    vd = np.array([-V0, 0.0, 0.0])
    Eb = Er - np.cross(vd, Bh)  # box frame: E = E_rest - v_d x B
    u = {name: _velocity(w, q, m, Er) for name, q, m in SPECIES}

    def cs(z):  # Re[z e^{i th}] = Re z cos th - Im z sin th
        return {"cos": [float(v) for v in np.real(z)], "sin": [float(v) for v in -np.imag(z)]}

    mode = {"k": [0.0, K, 0.0], "omega": float(w), "phase": 0.0, "B": cs(Bh), "E": cs(Eb),
            "u": {name: cs(v) for name, v in u.items()}}
    return {"face": "x+", "modes": [mode]}


def _cosine(C, S):
    """(scale, phase) with C cos(k y) + S sin(k y) = scale cos(k y + phase)."""
    return math.hypot(C, S), math.atan2(-S, C)


def _shape(C, S):
    R, ph = _cosine(C, S)
    return {"shape": "cos", "mode": 1, "length": LY, "phase": ph}, R


def preload(fl):
    """Initial fields and beam perturbations of the wave of ``fl`` at t = 0."""
    m = fl["modes"][0]
    fields, pert = [], {name: [] for name, _, _ in SPECIES}
    for kind, comps in (("B", "xz"), ("E", "x")):
        for c in comps:
            a = "xyz".index(c)
            R, ph = _cosine(m[kind]["cos"][a], m[kind]["sin"][a])
            if R > 0:
                fields.append({"kind": "profile", "component": kind + c, "axis": "y", "shape": "cos", "mode": 1,
                               "length": LY, "phase": ph, "scale": R})
    shp, R = _shape(m["E"]["cos"][2], m["E"]["sin"][2])  # xi_z = -E_z / (q n) with q = -1, n = 1
    pert["e"].append({"shape": {"y": shp}, "displacement": R, "displacement_component": "z"})
    for name, _, _ in SPECIES:
        for a, c in enumerate("xyz"):
            shp, R = _shape(m["u"][name]["cos"][a], m["u"][name]["sin"][a])
            if R > 0:
                pert[name].append({"shape": {"y": shp}, "velocity": R, "velocity_component": c})
    return fields, pert


def config(dims=None, steps=2000, piston=False, faces=None, preloaded=False):
    u0 = V0 / math.sqrt(1 - V0 * V0)
    fl = fluctuation()
    fields, pert = preload(fl) if preloaded else ([], {"e": [], "i": []})
    bc = (("reflecting" if piston else "open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    return SimConfig(
        n_cells=(NX, NY, NZ), dx=(DX,) * 3, dt=DT, n_steps=steps, boundaries=bc, ppc=32, npml=12,
        current_filter=CurrentFilter(4, True),
        B0=(0.0, B, 0.0), E0=(0.0, 0.0, V0 * B), dims=dims, diag_interval=50, gauss_interval=50,
        injection=Injection(faces=faces or {"x+": "flux"}, fluctuation=fl), initial_fields=fields, marder_zone=0.125,
        species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (-u0, 0, 0), THETA, perturbation=pert["e"])]),
                 SpeciesConfig("i", 1, MI, [Beam(1.0, (-u0, 0, 0), THETA / MI, perturbation=pert["i"])])])


if __name__ == "__main__":
    from pypic.simulation import Simulation

    n = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    piston = len(sys.argv) > 1 and sys.argv[1] == "piston"
    sim = Simulation(config(dims=(2, 1, 1) if n == 2 else None, steps=400, piston=piston), verbose=False)
    sim.run()
    if sim.rank == 0:
        print(json.dumps({k: sim.history[k] for k in ("kinetic", "electric", "magnetic", "gauss", "divB",
                                                      "n_injected")}))
