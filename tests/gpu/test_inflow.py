"""An inflow face that carries a prescribed fluctuation (PIC-D26)."""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys

import cupy as cp
import numpy as np
import pytest

from pypic import SimConfig
from pypic.output.names import stagger_of
from pypic.parallel.comm import SelfComm
from pypic.simulation import Simulation

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mgpu_inflow as M  # noqa: E402


def _mode(sim, comp, rows=None):
    """Complex amplitude of e^{i k y} per x row of a field and of its incident value."""
    g = sim.grid
    st = stagger_of(comp)
    F = cp.asnumpy(getattr(g, comp)[g.physical_slice()])[:, :, 0]
    x = (np.arange(F.shape[0]) + st[0]) * M.DX
    y = (np.arange(F.shape[1]) + st[1]) * M.DX
    inc = cp.asnumpy(sim._incident(comp, (cp.asarray(x[:, None, None]), cp.asarray(y[None, :, None]),
                                          cp.full((1, 1, 1), st[2] * M.DX)), sim.time))[:, :, 0]
    e = np.exp(-1j * M.K * y)[None, :]
    return x, 2 * (F * e).mean(axis=1), 2 * (inc * e).mean(axis=1)


@pytest.mark.parametrize("side", ["hi", "lo"])
def test_tfsf_face_injects_a_vacuum_wave_without_scattering(side):
    """The total-field / scattered-field correction alone, in vacuum: a plane wave that runs into the
    box through an open face leaves a scattered field below 1 % of its amplitude outside the face
    (measured 0.26 %, 2026-11-13)."""
    nx, ny, dx, dt = 128, 8, 0.25, 0.1
    cfg = SimConfig(n_cells=(nx, ny, 2), dx=(dx,) * 3, dt=dt, npml=12, species=[], gauss_interval=0,
                    diag_interval=0, boundaries=(("open", "open"), ("periodic", "periodic"), ("periodic", "periodic")))
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    k = 2 * math.pi / 8.0
    sgn = -1.0 if side == "hi" else 1.0
    sim.fluct = {"face": 1 if side == "hi" else 0,
                 "modes": [{"k": [sgn * k, 0, 0], "omega": k, "phase": 0.0, "E": ([0, 1, 0], [0, 0, 0]),
                            "B": ([0, 0, sgn], [0, 0, 0])}]}
    sim.fl_face = sim.fluct["face"]
    g = sim.grid
    for _ in range(300):
        t0 = sim.time
        g.push_B(0.5 * dt)
        sim._tfsf_B(0.5 * dt, t0)
        g.push_E(dt)
        sim._tfsf_E(dt, t0 + 0.5 * dt)
        g.push_B(0.5 * dt)
        sim._tfsf_B(0.5 * dt, t0 + dt)
        sim.time += dt
    Ey = cp.asnumpy(g.Ey)[:, g.g, g.g]
    pos = g.global_positions(0, 0.0)
    outside = pos > nx + 1 if side == "hi" else pos < -1
    assert np.max(np.abs(Ey[outside])) < 1e-2
    sim.close()


def test_inflow_face_keeps_the_wave_in_a_flow():
    """A cold R wave loaded in a flow and carried in through a flux face keeps its pattern: after
    t = 300 (15 c/omega_pe of new plasma) the k-mode of B_x and B_z differs from the incident one by
    at most 12 % on any x row and by 4 % on average away from the outflow face (measured 8 % and
    2.6 %); the density stays within 10 %; Gauss's law holds; the energy budget, with the Poynting
    flux of the total field through the face, closes within 1 % of the total energy (measured 0.31 %;
    4.2 % without the incident field in that flux)."""
    sim = Simulation(M.config(steps=3000, preloaded=True), comm=SelfComm(), verbose=False)
    sim.run()
    assert np.nanmax(sim.history["gauss"]) < 1e-10
    budget = np.asarray(sim.history["budget"])
    assert np.max(np.abs(budget - budget[0])) < 1e-2 * sim.history["total"][0]
    for comp in ("Bx", "Bz"):
        x, a, b = _mode(sim, comp)
        sel = x > 10.0  # the outflow face x- carries no incident field
        err = np.abs(a - b)[sel] / np.abs(b[sel]).mean()
        assert err.max() < 0.12, comp
        assert err[x[sel] < M.LX - 4].mean() < 0.04, comp
    for sp in sim.species:
        xs = cp.asnumpy(sp.rows()[0])
        h, _ = np.histogram(xs, np.arange(10.0, M.LX + 1e-9, 2.0))
        n = h / (2.0 / M.DX * M.NY * M.NZ * 32)
        assert np.all(np.abs(n - 1.0) < 0.1), sp.name
    sim.close()


def test_inflow_wave_meets_a_piston_shock():
    """With a reflecting wall opposite the inflow face (E0 tangential to it, PIC-D26) the flow piles up
    into a shock and the wave reaches it: downstream density above 2, the k-mode of B_z ahead of the
    shock between 0.5 and 2 of the incident one (the shock reflects part of the wave and the foot
    compresses it, so the upstream pattern is not the incident one; measured 0.78-1.50 on
    2026-11-14), Gauss's law."""
    sim = Simulation(M.config(steps=10000, piston=True), comm=SelfComm(), verbose=False)
    sim.run()
    assert np.nanmax(sim.history["gauss"]) < 1e-10
    xs = cp.asnumpy(sim.species[1].rows()[0])
    h, e = np.histogram(xs, np.arange(0.0, M.LX + 1e-9, 2.0))
    n = h / (2.0 / M.DX * M.NY * M.NZ * 32)
    assert n[:12].max() > 2.0
    shock = e[np.nonzero(n > 1.5)[0].max() + 1]
    x, a, b = _mode(sim, "Bz")
    ahead = (x > shock + 4) & (x < M.LX - 4)
    r = np.abs(a[ahead]) / np.abs(b[ahead]).mean()
    assert ahead.any() and r.min() > 0.5 and r.max() < 2.0
    sim.close()


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_inflow_is_decomposition_independent():
    """The fluctuation is evaluated from global coordinates and the injection keys do not depend on
    the layout: 1 GPU and 2 GPUs split along x give the same run."""
    script = os.path.join(HERE, "mgpu_inflow.py")
    outs = {}
    for n in (1, 2):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script],
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    assert outs[1]["n_injected"] == outs[2]["n_injected"]
    for k in ("kinetic", "electric", "magnetic"):
        a, b = np.asarray(outs[1][k]), np.asarray(outs[2][k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10, k
    assert np.nanmax(np.asarray(outs[2]["gauss"], dtype=float)) < 1e-10


def _entering_modes(precision):
    cfg = M.config(steps=50, preloaded=True)
    cfg.precision = precision
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    out = {}
    for sp in sim.species:
        P = cp.asnumpy(sp.rows())
        sel = P[0] > M.LX - 0.3  # entered in the last steps
        for a in (0, 2):
            out[sp.name, a] = 2 * np.mean(P[3 + a][sel] * np.exp(-1j * M.K * P[1][sel]))
    t = sim.time
    sim.close()
    return out, t


def test_entering_particles_take_the_mode_at_their_physical_place():
    """The compact storage works in cell units; the mode is evaluated at the physical place of the
    entering particles with both storages (review of PR #50): the k-mode of their velocity is the
    same with fp64 and compact within 2 % (measured: equal to the printed digits; before the fix the
    compact particles took the mode at 1/dy times their place), and the ions' z velocity is the
    incident one within 20 % (measured 8 %; the electrons entering through the face carry the
    flux-weighted bias of a thermal speed close to the flow)."""
    f64, t = _entering_modes("fp64")
    cmp_, _ = _entering_modes("compact")
    for key, v in f64.items():
        assert abs(cmp_[key] / v - 1) < 0.02, key
    m = M.fluctuation()["modes"][0]
    C, S = m["u"]["i"]["cos"][2], m["u"]["i"]["sin"][2]
    ref = (C - 1j * S) * np.exp(-1j * m["omega"] * t)
    assert abs(f64["i", 2] / ref - 1) < 0.2


def test_particles_next_to_the_face_feel_the_total_field():
    """A particle within half a cell of the face interpolates between the total field inside and the
    scattered field outside; the push gathers the total field on both sides (review of PR #50). With a
    table of zero amplitude the flow is a balanced E0 = -v x B0: the electrons in the last quarter cell
    keep u_z = 0 within four standard errors (before the fix: -0.0026, seven of them)."""
    from pypic.config import Injection

    fl = M.fluctuation()
    for m in fl["modes"]:
        m["B"] = m["E"] = {"cos": [0, 0, 0], "sin": [0, 0, 0]}
        m["u"] = {}
    cfg = M.config(steps=100)
    cfg.injection = Injection(faces={"x+": "flux"}, fluctuation=fl)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    P = cp.asnumpy(sim.species[0].rows())
    sel = P[0] >= M.LX - 0.25
    err = P[5][sel].std() / np.sqrt(sel.sum())
    assert abs(P[5][sel].mean()) < 4 * err
    sim.close()


def test_tracks_next_to_the_face_report_the_total_field(tmp_path):
    """The tracks gather the fields at the particles like the push: across the fluctuation face they
    see the total field (review of PR #50). A table of zero amplitude in a uniform B0 = B y: a test
    particle a quarter of a cell inside the face reports B_y = B0 (before the fix: 0.75 B0)."""
    from pypic import SpeciesConfig
    from pypic.analysis import open_run
    from pypic.config import Injection, OutputConfig

    fl = M.fluctuation()
    for m in fl["modes"]:
        m["B"] = m["E"] = {"cos": [0, 0, 0], "sin": [0, 0, 0]}
        m["u"] = {}
    cfg = M.config(steps=1)
    cfg.injection = Injection(faces={"x+": "flux"}, fluctuation=fl)
    u0 = M.V0 / math.sqrt(1 - M.V0 ** 2)
    cfg.species = list(cfg.species) + [SpeciesConfig("t", -1, 1, [], test=True,
                                                     particles=[[M.LX - 0.25 * M.DX, 3.0, 0.2, -u0, 0, 0]])]
    cfg.output = OutputConfig(dir=str(tmp_path), diagnostics=[
        dict(name="tr", kind="tracks", species=["t"], stride=1, every=1, with_fields=True, dtype="float64")])
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    sim.close()
    tr = open_run(str(tmp_path)).tracks("tr", "t")
    assert abs(float(np.asarray(tr["By"])[0, 0]) / M.B - 1.0) < 1e-9


def test_a_fluctuation_without_an_injecting_face_is_refused():
    """A fluctuation table that no face carries in is refused, not dropped (review of PR #50)."""
    from pypic.config import Injection

    cfg = M.config(steps=1)
    cfg.injection = Injection(faces={}, fluctuation=M.fluctuation())
    with pytest.raises(ValueError, match="injecting faces"):
        Simulation(cfg, comm=SelfComm(), verbose=False)


def test_marder_zone_leaves_no_incident_field_on_the_outer_conductor():
    """The Marder zone acts on the total field beyond the face; taking the incident field away again
    must leave the tangential E on the conductor behind the absorbing layer at zero, as a restore
    onto another layout rebuilds it (review of PR #50; before: -E0 there)."""
    sim = Simulation(M.config(steps=2, preloaded=True), comm=SelfComm(), verbose=False)
    sim.run()
    g = sim.grid
    pos = g.global_positions(0, 0.0)
    wall = np.nonzero(np.abs(pos - (M.NX + sim.cfg.npml)) < 1e-9)[0]
    assert wall.size == 1
    assert float(cp.max(cp.abs(g.Ez[int(wall[0])]))) == 0.0
    assert float(cp.max(cp.abs(g.Ey[int(wall[0])]))) == 0.0
    sim.close()


def test_a_fluctuation_alone_staggers_the_momenta(monkeypatch):
    """Without E0, B0, initial fields or displacements, a fluctuation table still acts on the particles
    next to its face from t = 0: the loaded momenta are taken back half a step (review of PR #50)."""
    calls = []
    orig = Simulation._initial_half_kick

    def spy(self):
        calls.append(True)
        return orig(self)

    monkeypatch.setattr(Simulation, "_initial_half_kick", spy)
    cfg = M.config(steps=1)
    cfg.B0, cfg.E0 = (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.close()
    assert calls
