"""GPU tests of the engine (single GPU; the multi-GPU test uses the launcher).
Skipped as a whole without CuPy or a CUDA device (tests/gpu/conftest.py)."""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys

import cupy as cp
import numpy as np
import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.parallel.comm import SelfComm
from pypic.simulation import Simulation

HERE = os.path.dirname(__file__)


def plasma(solver="explicit", bc=None, n=(16, 8, 8), dt=None, drift=0.2, theta=1e-3, steps=0, **kw):
    dt = dt if dt is not None else (0.08 if solver == "explicit" else 0.3)
    bc = bc or (("periodic", "periodic"),) * 3
    electrons = [Beam(0.5, (drift, 0, 0), theta), Beam(0.5, (-drift, 0, 0), theta)]
    cfg = SimConfig(n_cells=n, dx=(0.2, 0.2, 0.2), dt=dt, n_steps=steps, solver=solver, boundaries=bc,
                    species=[SpeciesConfig("e", -1, 1, electrons),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), theta / 25)])],
                    ppc=4, diag_interval=5, gauss_interval=5, npml=6, **kw)
    return Simulation(cfg, comm=SelfComm(), verbose=False)


def inject(sp, pos, u, w=1e-20):
    pos, u = np.atleast_2d(pos), np.atleast_2d(u)
    n = len(pos)
    sp.ensure_capacity(n)
    sp.P[:3, :n] = cp.asarray(pos.T)
    sp.P[3:6, :n] = cp.asarray(u.T)
    sp.P[6, :n] = w
    sp.n = n


# ----------------------------------------------------------------------- pusher
def test_boris_gyration_conserves_gamma_and_period():
    B = 0.5
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.5,) * 3, dt=0.05, species=[SpeciesConfig("e", -1, 1, [])],
                    B0=(0, 0, B), diag_interval=0, gauss_interval=0, sort_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sp = sim.species[0]
    u0 = 2.0  # relativistic: gamma = sqrt(5)
    inject(sp, [2.0, 2.0, 2.0], [u0, 0.0, 0.0])
    g0 = math.sqrt(1 + u0 * u0)
    wc = B / g0  # relativistic gyro-frequency
    period = 2 * math.pi / wc
    nst = int(round(period / cfg.dt))
    ang = []
    for _ in range(nst):
        sim.step()
        ux, uy = float(sp.P[3, 0]), float(sp.P[4, 0])
        assert abs(math.hypot(ux, uy) - u0) < 1e-13  # |u| exactly conserved by the rotation
        ang.append(math.atan2(uy, ux))
    # Boris phase error: the rotation angle per step is 2 atan(wc dt/2)
    expected = nst * 2 * math.atan(0.5 * wc * cfg.dt)
    total = np.unwrap(np.array(ang))[-1] - 0.0
    assert abs(abs(total) - expected) < 2 * 2 * math.atan(0.5 * wc * cfg.dt)


def test_exb_drift():
    E, B = 0.1, 0.5
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.5,) * 3, dt=0.05, species=[SpeciesConfig("e", -1, 1, [])],
                    B0=(0, 0, B), diag_interval=0, gauss_interval=0, sort_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.grid.Ey[...] = E  # uniform Ey
    sp = sim.species[0]
    inject(sp, [2.0, 2.0, 2.0], [0.0, 0.0, 0.0])
    x0 = 2.0
    nst = 4000
    xs = []
    for _ in range(nst):
        sim.step()
        xs.append(float(sp.P[0, 0]))
    xs = np.unwrap(np.array(xs) / sim.grid.L[0] * 2 * np.pi) * sim.grid.L[0] / (2 * np.pi)
    vd = (xs[-1] - x0) / (nst * cfg.dt)
    assert abs(vd - E / B) < 0.01 * E / B


# ------------------------------------------------------------ conservation laws
@pytest.mark.parametrize("bc", ["periodic", "reflecting", "open"])
def test_explicit_gauss_and_divB(bc):
    sim = plasma("explicit", bc=((bc, bc), ("periodic", "periodic"), (bc if bc != "open" else "periodic",) * 2))
    sim.run(40)
    g = np.nanmax(sim.history["gauss"])
    b = np.nanmax(sim.history["divB"])
    assert g < 1e-12, g
    assert b < 1e-12, b


@pytest.mark.parametrize("bc", ["periodic", "reflecting", "open"])
def test_implicit_exact_energy_and_gauss(bc):
    sim = plasma("implicit", bc=((bc, bc), ("periodic", "periodic"), ("periodic", "periodic")), drift=0.3,
                 theta=0.01)
    sim.run(30)
    h = sim.history
    e = np.asarray(h["budget"])
    assert np.max(np.abs(e - e[0])) / e[0] < 1e-11
    assert np.nanmax(h["gauss"]) < 1e-11
    assert np.nanmax(h["divB"]) < 1e-12


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
def test_mean_current_removal_keeps_conservation(solver):
    """Net electron current (Buneman set-up) with remove_mean_current: the box-averaged E stays
    zero (no k = 0 plasma oscillation), Gauss's law is exact and the implicit energy stays exact."""
    cfg = SimConfig(n_cells=(16, 8, 8), dx=(0.2,) * 3, dt=0.08 if solver == "explicit" else 0.3, solver=solver,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.3, 0, 0), 0.01)]),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                    ppc=4, diag_interval=5, gauss_interval=5, remove_mean_current=True)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(30)
    o = sim.grid.owned_slice()
    assert abs(float(cp.mean(sim.grid.Ex[o]))) < 1e-13  # without the removal <E_x> ~ -<J_x> t ~ 0.3
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    if solver == "implicit":
        e = np.asarray(sim.history["budget"])
        assert np.max(np.abs(e - e[0])) / e[0] < 1e-11


def test_reflecting_keeps_particles():
    bc = (("reflecting", "reflecting"),) * 3
    sim = plasma("explicit", bc=bc, drift=0.5, theta=0.01)
    n0 = sim.total_particles()
    sim.run(60)
    assert sim.total_particles() == n0
    L = sim.grid.L
    for sp in sim.species:
        x = cp.asnumpy(sp.P[:3, : sp.n])
        for a in range(3):
            assert x[a].min() >= 0.0 and x[a].max() <= L[a]


def test_open_boundary_removes_particles():
    bc = (("open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    sim = plasma("explicit", bc=bc, drift=0.5, theta=0.01)
    n0 = sim.total_particles()
    sim.run(60)
    assert sim.total_particles() < n0
    assert sim.history["lost_kinetic"][-1] > 0


# --------------------------------------------------------------------- fields
def test_vacuum_wave_dispersion():
    """A y-polarised plane wave along x follows the Yee dispersion relation."""
    N, dx, dt = 64, 0.1, 0.05
    cfg = SimConfig(n_cells=(N, 4, 4), dx=(dx,) * 3, dt=dt, species=[], diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    g = sim.grid
    k = 2 * math.pi * 4 / (N * dx)
    w = 2 / dt * math.asin(dt / dx * math.sin(k * dx / 2))  # Yee dispersion
    xi = g.global_positions(0, 0.0) * dx  # Ey at integer x
    xh = g.global_positions(0, 0.5) * dx  # Bz at half x
    g.Ey[...] = cp.asarray(np.cos(k * xi))[:, None, None]
    g.Bz[...] = cp.asarray(np.cos(k * xh))[:, None, None]  # B^0 for a +x wave (approximately in phase)
    nst = 400
    for _ in range(nst):
        sim.step()
    s = g.physical_slice()
    ey = cp.asnumpy(g.Ey[s][:, 0, 0])
    x = xi[s[0]]
    phase = math.atan2(-np.sum(ey * np.sin(k * x)), np.sum(ey * np.cos(k * x)))
    expected = (-w * nst * dt + math.pi) % (2 * math.pi) - math.pi
    assert abs((phase - expected + math.pi) % (2 * math.pi) - math.pi) < 0.05


def test_cpml_absorbs_pulse():
    N, dx = 96, 0.1
    cfg = SimConfig(n_cells=(N, 4, 4), dx=(dx,) * 3, dt=0.05, species=[], npml=12,
                    boundaries=(("open", "open"), ("periodic", "periodic"), ("periodic", "periodic")),
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    g = sim.grid
    xi = g.global_positions(0, 0.0) * dx
    xh = g.global_positions(0, 0.5) * dx
    x0, wid = N * dx / 2, 0.5
    g.Ey[...] = cp.asarray(np.exp(-((xi - x0) / wid) ** 2))[:, None, None]
    g.Bz[...] = cp.asarray(np.exp(-((xh - x0) / wid) ** 2))[:, None, None]
    g.apply_bc_E()
    e0 = sum(g.field_energy())
    sim.run(400)  # pulse travels 20 c/omega_pe >> half box
    e1 = sum(g.field_energy())
    assert e1 / e0 < 1e-4, e1 / e0


def test_two_stream_growth_matches_theory(tmp_path):
    """Small 3-D electron two-stream run: measured growth of the k_max mode (read back from the
    modes output at the cadence the planner derived) vs linear theory."""
    from pypic.analysis import open_run
    from pypic.analysis.spectral import fit_growth
    from pypic.setup.planner import TwoStreamParams, plan_two_stream

    p = TwoStreamParams(mixture="e-i", mass_ratio=1836.0, u_drift=0.2, T_e=1e-4, n_wavelengths=2,
                        cells_per_wavelength=32, transverse_ratio=0.125, ppc=32, growth_times=14,
                        n_records=300)
    cfg, _, info = plan_two_stream(p)
    cfg.output.diagnostics = [d for d in cfg.output.diagnostics if d.name == "modes"]
    cfg.output.dir = str(tmp_path)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    sim.close()
    md = open_run(str(tmp_path)).series("modes")
    gfit, _ = fit_growth(md.t, md["Ex"][:, 1])
    assert abs(gfit - info["gamma_max"]) / info["gamma_max"] < 0.15, (gfit, info["gamma_max"])


# ------------------------------------------------------------------ multi-GPU
@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
@pytest.mark.parametrize("solver", ["explicit", "implicit"])
def test_multi_gpu_matches_single(solver, tmp_path):
    script = os.path.join(HERE, "mgpu_case.py")
    outs = {}
    for n in (1, 2):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, solver],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    for k in ("kinetic", "electric", "magnetic"):
        a, b = np.asarray(outs[1][k]), np.asarray(outs[2][k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10


# --------------------------------------------------------- review regressions
def test_cg_handles_exact_convergence_between_checks():
    """A curl-free uniform field is solved in one CG step (A = I on that mode); the following
    iteration must not divide 0/0 before the next convergence check."""
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.5,) * 3, dt=0.5, solver="implicit", species=[],
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    b = cp.zeros((3,) + sim.grid.shape)
    for c, v in enumerate((1.0, -2.0, 0.5)):
        b[c][sim.own] = v
    x = cp.zeros_like(b)
    sim._cg(x, b)
    assert bool(cp.all(cp.isfinite(x)))
    assert float(cp.max(cp.abs(x[(slice(None),) + sim.own] - b[(slice(None),) + sim.own]))) < 1e-12


@pytest.mark.parametrize("theta", [1e-3, 0.05, 0.3, 3.0])
def test_maxwell_juttner_sampling(theta):
    from scipy.special import kve

    cfg = SimConfig(n_cells=(32, 32, 32), dx=(0.5,) * 3, dt=0.1, ppc=8,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), theta)])],
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sp = sim.species[0]
    u = cp.asnumpy(sp.P[3:6, : sp.n])
    eps = np.sqrt(1 + (u ** 2).sum(0)) - 1
    expected = kve(1, 1 / theta) / kve(2, 1 / theta) + 3 * theta - 1
    assert abs(eps.mean() - expected) < 5 * eps.std() / math.sqrt(eps.size)


def test_unconverged_implicit_step_is_not_committed():
    from pypic.simulation import ImplicitStepError

    sim = plasma("implicit", drift=0.3, theta=0.01)
    sim.run(2)
    sim.cfg.implicit.picard_maxit = 1
    P0 = [sp.P[:, : sp.n].copy() for sp in sim.species]
    E0 = sim.grid.Ex.copy()
    with pytest.raises(ImplicitStepError):
        sim.step()
    for sp, p in zip(sim.species, P0):
        assert bool(cp.all(sp.P[:, : sp.n] == p))
    assert bool(cp.all(sim.grid.Ex == E0))


def test_ghost_wider_than_local_grid_is_rejected():
    cfg = SimConfig(n_cells=(4, 4, 4), dx=(1.0,) * 3, dt=3.5, solver="implicit", species=[],
                    diag_interval=0, gauss_interval=0)
    with pytest.raises(ValueError, match="ghost width"):
        Simulation(cfg, comm=SelfComm(), verbose=False)




def test_implicit_reflection_with_displacement_beyond_half_box_conserves_charge():
    """Orbit folded at the lower wall that ends past the box centre (|d| > L/2): the deposit must
    fold about the wall the unfolded path crossed, or Gauss's law breaks (review, PR #1)."""
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(1.0,) * 3, dt=5.1, solver="implicit",
                    boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    species=[SpeciesConfig("e", -1, 1, []), SpeciesConfig("p", 1, 1, [])],
                    diag_interval=0, gauss_interval=0, sort_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    e, p = sim.species
    inject(e, [0.1, 2.3, 2.3], [-4.9, 0.04, 0.03], w=1e-3)  # d_x ~ -5 cells: reflects, ends at x ~ 4.9
    inject(p, [0.1, 2.3, 2.3], [0.0, 0.0, 0.0], w=1e-3)  # neutral start: rho = 0, E = 0
    assert sim.gauss_error()[0] < 1e-15
    sim.step()
    x = float(e.P[0, 0])
    assert 4.0 < x < 8.0 and float(e.P[3, 0]) > 0  # reflected past the box centre
    assert sim.gauss_error()[0] < 1e-12
