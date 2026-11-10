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
def test_esirkepov_current_moment_equals_the_displacement():
    """Without fields, sum J dV dt = sum q w (x1 - x0) per component: the deposit of moves along
    every axis and across cell faces carries the exact current moment (issue #6)."""
    sim = plasma("explicit", B0=(0.0, 0.0, 0.0))
    sp = sim.species[0]
    rng = np.random.default_rng(5)
    n = 2000
    pos = rng.uniform(0.5, 1.1, (n, 3)) * np.asarray(sim.grid.L) / 1.6
    u = rng.uniform(-3.0, 3.0, (n, 3))
    inject(sp, pos, u, w=1.0)
    sim.species[1].n = 0
    sim.step()
    g = sim.grid
    o = g.owned_slice()
    J = [float(cp.sum(F[o])) * g.cell_volume * sim.cfg.dt for F in (g.Jx, g.Jy, g.Jz)]
    gam = np.sqrt(1 + np.sum(u * u, axis=1))
    ref = sp.q * np.sum(u / gam[:, None], axis=0) * sim.cfg.dt
    np.testing.assert_allclose(J, ref, rtol=1e-12, atol=1e-12 * np.abs(ref).max())


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


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
def test_perturbation_is_gauss_exact_with_the_requested_density(solver):
    """Density mode n = n0 (1 + alpha cos kx) set up by a charge-conserving virtual move (PIC-D9):
    Gauss's law exact at t = 0 and after stepping, the density mode and the field energy are the
    requested ones, the implicit energy stays exact."""
    alpha, m = 0.3, 2
    cfg = SimConfig(n_cells=(32, 8, 8), dx=(0.1,) * 3, dt=0.05 if solver == "explicit" else 0.2, solver=solver,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.001, perturbation=[
                        {"axis": "x", "mode": m, "density": alpha}])]),
                             SpeciesConfig("i", 1, 100, [Beam(1.0, (0, 0, 0), 1e-5)])],
                    ppc=16, diag_interval=5, gauss_interval=5)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    assert sim.gauss_error()[0] < 1e-13
    sp = sim.species[0]
    x, w = cp.asnumpy(sp.P[0, : sp.n]), cp.asnumpy(sp.P[6, : sp.n])
    L = sim.grid.L[0]
    k = 2 * math.pi * m / L
    nk = abs(np.sum(w * np.exp(-1j * k * x)) / np.sum(w))
    assert abs(nk - alpha / 2) < 4 / math.sqrt(sp.n)  # random positions: noise ~ 1/sqrt(N)
    V = L * 0.8 * 0.8
    assert abs(sim.grid.field_energy()[0] / (alpha ** 2 / (4 * k * k) * V) - 1) < 0.02
    sim.run(20)
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    if solver == "implicit":
        e = np.asarray(sim.history["budget"])
        assert np.max(np.abs(e - e[0])) / e[0] < 1e-11


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
@pytest.mark.parametrize("solver", ["explicit", "implicit"])
@pytest.mark.parametrize("variant", ["plain", "filter"])
def test_perturbation_across_ranks_matches_single_gpu(solver, variant):
    """A displacement of 21 cells (more than a rank's 16) is moved in rounds with migration: the
    run is the same on 1 and 2 GPUs and Gauss's law holds; `filter` adds reflecting walls on the
    decomposed axis and the current filter."""
    script = os.path.join(HERE, "mgpu_perturb.py")
    outs = {}
    for n in (1, 2):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, solver, variant],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    for k in ("kinetic", "electric", "magnetic"):
        a, b = np.asarray(outs[1][k]), np.asarray(outs[2][k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10
    assert np.nanmax(np.asarray(outs[2]["gauss"], dtype=float)) < 1e-12


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
def test_test_species_do_not_change_the_plasma(solver):
    """A test species (loaded beam + explicit particles) feels the fields but deposits nothing, is
    outside the budget and does not shift the plasma's random numbers (PIC-D9). Its beam has a
    density perturbation: the virtual move that sets it up deposits nothing either (review of
    PR #2: it used the real charge and changed E)."""
    plasma = [SpeciesConfig("e", -1, 1, [Beam(1.0, (0.1, 0, 0), 0.01)]),
              SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 1e-3)])]
    test = SpeciesConfig("t", -1, 1, [Beam(0.3, (0.5, 0, 0), 0.0, perturbation=[{"axis": "x", "density": 0.5}])],
                         test=True, particles=[[0.5, 0.5, 0.5, 2.0, 0, 0], [1.0, 0.2, 0.3, 0, 1.0, 0]])
    hist = []
    for sp in (plasma, [test] + plasma):  # declared first: must not consume the plasma's RNG streams
        cfg = SimConfig(n_cells=(16, 8, 8), dx=(0.2,) * 3, dt=0.08 if solver == "explicit" else 0.3, solver=solver,
                        species=sp, ppc=4, diag_interval=5, gauss_interval=5, B0=(0, 0, 0.3))
        sim = Simulation(cfg, comm=SelfComm(), verbose=False)
        sim.run(30)
        hist.append(sim.history)
    for k in ("kinetic", "electric", "budget"):
        a, b = np.asarray(hist[0][k]), np.asarray(hist[1][k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-12, k
    assert np.nanmax(hist[1]["gauss"]) < 1e-12


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
def test_plane_wave_is_an_exact_eigenmode(solver):
    """initial_fields plane_wave is the discrete eigenmode of the solver (Yee with the stored-B
    half-step average, or Crank-Nicolson): after 400 steps it equals the analytic discrete wave."""
    from pypic.initial_fields import omega_of

    N, dx = 64, 0.1
    dt = 0.05 if solver == "explicit" else 0.2
    cfg = SimConfig(n_cells=(N, 4, 4), dx=(dx,) * 3, dt=dt, solver=solver, species=[], diag_interval=0,
                    gauss_interval=0, initial_fields=[{"kind": "plane_wave", "axis": "x", "mode": 3,
                                                       "polarization": "y", "amplitude": 0.7, "phase": 0.4}])
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    k = 2 * math.pi * 3 / (N * dx)
    w = omega_of(k, dx, dt, solver == "implicit")
    sim.run(400)
    g = sim.grid
    s = g.physical_slice()
    x = g.global_positions(0, 0.0)[s[0]] * dx
    assert np.max(np.abs(cp.asnumpy(g.Ey[s][:, 0, 0]) - 0.7 * np.cos(k * x - w * 400 * dt + 0.4))) < 1e-12


def test_density_profile_is_exact_and_neutral():
    """A step profile shared by electrons and ions: only particles inside are kept, rho = 0."""
    prof = {"axis": "x", "shape": "step", "x0": 0.0, "x1": 3.0}
    cfg = SimConfig(n_cells=(32, 4, 4), dx=(0.25,) * 3, dt=0.1, ppc=8,
                    boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.01, profile=prof)]),
                             SpeciesConfig("i", 1, 100, [Beam(1.0, (0, 0, 0), 0.0, profile=prof)])],
                    diag_interval=5, gauss_interval=5)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    e = sim.species[0]
    assert e.n == 12 * 4 * 4 * 8 and float(cp.max(e.P[0, : e.n])) < 3.0
    assert sim.gauss_error()[0] < 1e-13
    sim.run(20)
    assert np.nanmax(sim.history["gauss"]) < 1e-12


def test_current_filter_transfer_function():
    """(1,2,1)/4 x N on cos(kx): amplitude cos^2N(k dx/2), times 1 + (N/2)(1 - cos k dx) with the compensator."""
    cfg = SimConfig(n_cells=(32, 4, 4), dx=(0.1,) * 3, dt=0.05, species=[], diag_interval=0, gauss_interval=0)
    g = Simulation(cfg, comm=SelfComm(), verbose=False).grid
    x = g.global_positions(0, 0.0) * 0.1
    for N, comp in ((4, False), (4, True)):
        for m in (1, 4, 8):
            k = 2 * math.pi * m / g.L[0]
            F = cp.asarray(np.broadcast_to(np.cos(k * x)[:, None, None], g.shape).copy())
            g.filter([F], [(0, 0, 0)], N, comp)
            T = math.cos(k * 0.05) ** (2 * N) * ((1 + N / 2 * (1 - math.cos(k * 0.1))) if comp else 1)
            assert abs(float(cp.max(F[g.physical_slice()])) - T) < 1e-12


def _filter_reference(g, arrays, staggers, passes, compensate):
    """The array-expression filter that the fused kernels replaced (issue #9)."""
    w = [(0.25, 0.5, 0.25)] * passes + ([(-passes / 4.0, 1.0 + passes / 2.0, -passes / 4.0)] if compensate else [])
    for a in range(3):
        ax = g.d.axes[a]
        lo_wall = g.d.bc[a][0] != "periodic" and ax.lo is None
        hi_wall = g.d.bc[a][1] != "periodic" and ax.hi is None
        f = g.g + ax.n_int
        for wl, wc, wr in w:
            g.fill_ghosts(arrays, axes=(a,))
            for F, st in zip(arrays, staggers):
                node = st[a] == 0
                if lo_wall:
                    g._slab(F, a, g.g - 1, g.g)[...] = g._slab(F, a, g.g, g.g + 1) if node else 0.0
                if hi_wall:
                    if node:
                        g._slab(F, a, f + 1, f + 2)[...] = g._slab(F, a, f, f + 1)
                    else:
                        g._slab(F, a, f, f + 1)[...] = 0.0
                hi = f + 1 if (hi_wall and node) else f
                new = (wl * g._slab(F, a, g.g - 1, hi - 1) + wc * g._slab(F, a, g.g, hi)
                       + wr * g._slab(F, a, g.g + 1, hi + 1))
                g._slab(F, a, g.g, hi)[...] = new


@pytest.mark.parametrize("xbc", ["periodic", "reflecting"])
def test_fill_ghosts_gives_the_periodic_images(xbc):
    """Every ghost cell of a periodic axis, edges and corners included, holds exactly the value of its
    periodic image; cells of a wall axis are left alone. All axes in one call, and one axis alone."""
    cfg = SimConfig(n_cells=(7, 5, 3), dx=(0.1,) * 3, dt=0.05, species=[], diag_interval=0, gauss_interval=0,
                    boundaries=((xbc, xbc), ("periodic", "periodic"), ("periodic", "periodic")), npml=2)
    g = Simulation(cfg, comm=SelfComm(), verbose=False).grid
    rs = cp.random.RandomState(11)
    for axes in ((0, 1, 2), (1,)):
        arrays = [rs.standard_normal(g.shape) for _ in range(4)]  # two launches: 3 arrays, then 1
        before = [cp.asnumpy(F) for F in arrays]
        g.fill_ghosts(arrays, axes=axes)
        image = []
        for a in range(3):
            i, n = np.arange(g.shape[a]), g.d.axes[a].n_int
            wraps = a in axes and g.d.bc[a][0] == "periodic"
            image.append(g.g + (i - g.g) % n if wraps else i)
        for F, F0 in zip(arrays, before):
            assert np.array_equal(cp.asnumpy(F), F0[np.ix_(*image)])


@pytest.mark.parametrize("bc", ["periodic", "reflecting"])
def test_current_filter_kernels_match_the_array_expressions(bc):
    """The fused filter kernels give bitwise the result of the former array expressions, for J
    (three staggers) and node data, at periodic and conducting faces."""
    cfg = SimConfig(n_cells=(12, 10, 8), dx=(0.1,) * 3, dt=0.05, species=[], boundaries=((bc, bc),) * 3,
                    diag_interval=0, gauss_interval=0, npml=2)
    g = Simulation(cfg, comm=SelfComm(), verbose=False).grid
    rng = np.random.default_rng(7)
    for staggers in (Simulation.J_STAGGER, ((0, 0, 0),)):
        for comp in (True, False):  # an even and an odd number of passes per axis
            A = [cp.asarray(rng.standard_normal(g.shape)) for _ in staggers]
            B = [a.copy() for a in A]
            g.filter(A, staggers, 3, comp)
            _filter_reference(g, B, staggers, 3, comp)
            for a, b in zip(A, B):
                assert bool(cp.all(a == b))


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
@pytest.mark.parametrize("bc", ["periodic", "reflecting"])
def test_current_filter_keeps_gauss_and_implicit_energy(solver, bc):
    """With a 4-pass compensated current filter (PIC-D10) and particles bouncing at the walls,
    Gauss's law (against the identically filtered rho) and the implicit energy stay exact."""
    from pypic.config import CurrentFilter

    cfg = SimConfig(n_cells=(16, 8, 8), dx=(0.2,) * 3, dt=0.08 if solver == "explicit" else 0.3, solver=solver,
                    boundaries=((bc, bc), ("periodic",) * 2, (bc, bc)),
                    species=[SpeciesConfig("e", -1, 1, [Beam(0.5, (0.5, 0.1, 0), 0.05),
                                                        Beam(0.5, (-0.5, 0, 0.1), 0.05)]),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.01)])],
                    ppc=4, diag_interval=5, gauss_interval=5, current_filter=CurrentFilter(4, True))
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(40)
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    if solver == "implicit":
        e = np.asarray(sim.history["budget"])
        assert np.max(np.abs(e - e[0])) / e[0] < 1e-11


# ------------------------------------------------------------------- storage
BIG_CAP = 360_000_000  # 6 * cap > 2**31: the weight row starts beyond int32 offsets (issue #11)


@pytest.mark.skipif(cp.cuda.Device().mem_info[0] < 32 * 2**30, reason="needs ~24 GB of free GPU memory")
def test_particle_rows_beyond_int32_offsets():
    """A capacity whose row offsets (row * cap) exceed 2**31 gives the same run as a small one."""
    ref = plasma("explicit", sort_interval=5)
    sim = plasma("explicit", sort_interval=5)
    sp = sim.species[0]
    ids = cp.asnumpy(sp.P[7, : sp.n]).copy()
    sp.ensure_capacity(BIG_CAP)
    assert 6 * sp.cap > 2**31
    ref.run(10)
    sim.run(10)
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    np.testing.assert_allclose(sim.history["total"], ref.history["total"], rtol=1e-12)
    assert np.array_equal(np.sort(cp.asnumpy(sp.P[7, : sp.n])), np.sort(ids))
    del sim, sp
    cp.get_default_memory_pool().free_all_blocks()


def test_capacity_beyond_int32_particle_index_is_rejected():
    from pypic.particles import Species

    with pytest.raises(ValueError, match="2\\*\\*31"):
        Species("e", -1, 1, 2**31)


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


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
def test_locally_computed_ghost_layer_of_B_equals_the_neighbours_values():
    """The explicit update computes B on one ghost layer instead of exchanging it. With open x on 4
    ranks the edge ranks run the kernel with the absorbing layer and the middle ranks the one
    without: both must give the same bits, or the ghosts drift from the neighbours' values."""
    r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", "4", os.path.join(HERE, "mgpu_regress.py"),
                        "b_ghost"], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    for rank in range(4):
        assert f"RESULT rank={rank} ranks_with_a_difference=0 b_is_nonzero=True" in r.stdout, r.stdout


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
def test_migration_across_faces_and_corners_keeps_every_particle():
    """Two species cross faces and corners of a 2x2 rank grid every step (issue #7): the run equals
    the 1-GPU run, every particle id survives, and an overflowing flag list (full-compaction path
    on some ranks) gives the same run. With one species at rest, steps that move particles of the
    other species only work too (PR #15 review: CuPy 14.0's bincount raised on an empty input)."""
    script = os.path.join(HERE, "mgpu_migrate.py")
    outs = {}
    for n, dims, variant in ((1, "1,1,1", "plain"), (4, "2,2,1", "plain"), (4, "2,2,1", "overflow"),
                             (1, "1,1,1", "idle"), (4, "2,2,1", "idle")):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, dims, variant],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[(n, variant)] = json.loads(r.stdout.strip().splitlines()[-1])
    for key in ((4, "plain"), (4, "overflow"), (4, "idle")):
        ref = outs[(1, "idle" if key[1] == "idle" else "plain")]
        o = outs[key]
        assert o["ids"] == ref["ids"], key
        for k in ("kinetic", "electric"):
            a, b = np.asarray(ref[k]), np.asarray(o[k])
            assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10, (key, k)
        assert np.nanmax(np.asarray(o["gauss"], dtype=float)) < 1e-12


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
def test_mixed_gather_is_decomposition_independent():
    """The float32 gather (PIC-D11) takes its weights from global coordinates: 1 GPU and a 2x2
    rank grid give the same run, Gauss's law holds, and it stays close to the FP64 run."""
    script = os.path.join(HERE, "mgpu_migrate.py")
    outs = {}
    for n, dims, variant in ((1, "1,1,1", "mixed"), (4, "2,2,1", "mixed"), (1, "1,1,1", "plain")):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, dims, variant],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[(n, variant)] = json.loads(r.stdout.strip().splitlines()[-1])
    a, b = (np.asarray(outs[(n, "mixed")]["kinetic"]) for n in (1, 4))
    assert np.max(np.abs(a - b) / np.abs(a)) < 1e-10
    assert np.nanmax(np.asarray(outs[(4, "mixed")]["gauss"], dtype=float)) < 1e-12
    f = np.asarray(outs[(1, "plain")]["kinetic"])
    assert np.max(np.abs(a - f) / np.abs(f)) < 1e-5


def test_mixed_precision_is_explicit_only():
    with pytest.raises(ValueError, match="explicit"):
        plasma("implicit", precision="mixed")


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


def test_cg_inner_products_match_the_owned_unknowns():
    """The inner products that the CG takes with its matrix product (sums per thread block, then
    of the blocks) are those over the owned, unconstrained E unknowns: (r, v), (A v, v), (r, r)
    agree with the sums on the host (fixed seed; mixed signs, and the norm of the convergence test)."""
    bc = (("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))
    cfg = SimConfig(n_cells=(48, 32, 24), dx=(0.5,) * 3, dt=0.5, solver="implicit", species=[], boundaries=bc,
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    rng = np.random.default_rng(3)
    v, r = (cp.asarray(rng.standard_normal((3,) + sim.grid.shape)) for _ in range(2))
    Av, tmp = cp.empty_like(v), cp.empty_like(v)
    got = cp.asnumpy(sim._matvec(v, Av, tmp, r))
    W, v, r, Av = (cp.asnumpy(a) for a in (sim.Wmask, v, r, Av))
    for value, (a, b) in zip(got, ((r, v), (Av, v), (r, r))):
        assert abs(value - np.sum(W * a * b)) <= 1e-13 * np.sum(np.abs(W * a * b))


def test_cg_matrix_product_in_one_launch_matches_the_kernel_sequence():
    """On a periodic grid the CG takes its matrix product and inner products in one launch: they are those
    of the general sequence (curl E, curl B, diagonal term) to rounding, on every unknown (unequal cell
    sizes, fixed seed)."""
    cfg = SimConfig(n_cells=(12, 10, 8), dx=(0.5, 0.4, 0.3), dt=0.7, solver="implicit", species=[],
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    assert sim.cg_fused
    rs = cp.random.RandomState(11)
    v, r = (rs.standard_normal((3,) + sim.grid.shape) for _ in range(2))
    result = []
    for fused in (True, False):
        sim.cg_fused = fused
        Av, tmp = cp.zeros_like(v), cp.zeros_like(v)
        sums = cp.asnumpy(sim._matvec(v.copy(), Av, tmp, r))
        result.append((cp.asnumpy(Av[(slice(None),) + sim.own]), sums))
    (Av1, s1), (Av3, s3) = result
    assert np.abs(Av3).max() > 1 and np.abs(Av1 - Av3).max() <= 1e-14 * np.abs(Av3).max()
    assert np.all(np.abs(s1 - s3) <= 1e-13 * np.abs(Av3).sum())


@pytest.mark.parametrize("x_faces", ["reflecting", "open"])
def test_cg_solves_the_field_equation_to_its_tolerance(x_faces):
    """The CG keeps a recurrence for the residual and for the matrix product of its search
    direction; the residual of its result, taken with the matrix itself, is at the level of
    cg_tol (conducting walls, and the absorbing layer of open faces; fixed seed)."""
    bc = ((x_faces, x_faces), ("periodic", "periodic"), ("periodic", "periodic"))
    cfg = SimConfig(n_cells=(48, 16, 12), dx=(0.5,) * 3, dt=1.0, solver="implicit", species=[], boundaries=bc,
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    b = cp.asarray(np.random.default_rng(5).standard_normal((3,) + sim.grid.shape)) * sim.Wmask
    x = cp.zeros_like(b)
    its = sim._cg(x, b)
    Ax, tmp = cp.empty_like(x), cp.empty_like(x)
    sim._matvec(x, Ax, tmp, b)
    res = float(cp.linalg.norm((b - Ax) * sim.Wmask)) / float(cp.linalg.norm(b))
    assert 0 < its < cfg.implicit.cg_maxit
    assert res <= 10 * cfg.implicit.cg_tol


@pytest.mark.parametrize("x_faces", ["reflecting", "open"])
def test_half_step_fields_and_picard_change_match_the_array_expressions(x_faces):
    """The kernels of the half-step fields give, on every point (ghosts included), the bits of
    Eh = (E0 + E1) / 2 with a ghost fill of its own and of Bh = D (B0 - dt/2 curl Eh). The change over a
    Picard sweep is the largest |Enew - E1| and |Enew| over the unknowns, and a NaN unknown is kept."""
    bc = ((x_faces, x_faces), ("periodic", "periodic"), ("periodic", "periodic"))
    cfg = SimConfig(n_cells=(24, 6, 5), dx=(0.5, 0.4, 0.3), dt=0.2, solver="implicit", species=[], boundaries=bc,
                    diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    g = sim.grid
    rs = cp.random.RandomState(3)
    E0, B0 = (g.Ex, g.Ey, g.Ez), (g.Bx, g.By, g.Bz)
    for F in (*E0, *B0, sim.E1):
        F[...] = rs.standard_normal(F.shape)
    for E in (E0, sim.E1):  # as inside a step: boundary values set, ghosts filled
        g.apply_bc_E(E)
        g.fill_ghosts(E)
    g.apply_bc_B(B0)
    g.fill_ghosts(B0)
    sim._half_fields(E0, B0)
    Eh = cp.stack([0.5 * (e0 + e1) for e0, e1 in zip(E0, sim.E1)])
    g.apply_bc_E(Eh)
    g.fill_ghosts(Eh)
    Bh = cp.zeros_like(Eh)
    g.curl_E(Eh, Bh, -0.5 * cfg.dt)
    for c in range(3):
        Bh[c] += B0[c]
    sim._applyD(Bh)
    g.apply_bc_B(Bh)
    g.fill_ghosts(Bh)
    assert bool(cp.all(sim.Eh == Eh)) and bool(cp.all(sim.Bh == Bh))

    Enew = sim.E1 + 1e-3 * rs.standard_normal(sim.E1.shape)
    expected = (float(cp.max(cp.abs(Enew - sim.E1) * sim.Wmask)), float(cp.max(cp.abs(Enew) * sim.Wmask)))
    assert sim._picard_change(Enew) == expected and expected[0] > 0
    Enew[1][sim.own][-1, -1, -1] = np.nan
    assert np.isnan(sim._picard_change(Enew)[0])


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


def test_implicit_momentum_is_the_fixed_point_of_the_midpoint_rotation():
    """In uniform E and B the implicit push returns u1 = 2 ub - u0 with ub = u0 + eps E + ub x t,
    t = 2 eps B / (g0 + g1) and g1 = sqrt(1 + u1^2). Compared with that fixed point iterated in extended
    precision, for half rotation angles |t| from 1e-9 to 0.14 and gamma from 1.1 to 30: the sweeps
    over g1 must not stop before the direction of the momentum is converged to rounding. One evaluation
    of the particle map (no iteration over the path, which would add sweeps): the fields are uniform."""
    from pypic.backend import launch

    cfg = SimConfig(n_cells=(8, 8, 8), dx=(1.0,) * 3, dt=0.1, solver="implicit", diag_interval=0,
                    gauss_interval=0, species=[SpeciesConfig("e", -1, 1, [])])
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sp, g = sim.species[0], sim.grid
    eps = sp.q * cfg.dt / (2.0 * sp.m)
    rs = np.random.RandomState(5)
    n = 2048

    def vectors(length):
        v = rs.standard_normal((n, 3))
        return v / np.linalg.norm(v, axis=1)[:, None] * np.asarray(length).reshape(-1, 1)

    for bmag in (1e-6, 1e-3, 0.3, 3.0):
        u0 = vectors(10 ** rs.uniform(-0.3, 1.5, n))
        E, B = vectors(2.0)[0], vectors(bmag)[0]
        inject(sp, rs.uniform(0.0, 8.0, (n, 3)), u0)
        F = [cp.full(g.shape, v) for v in (*E, *B)]
        sim.err.fill(0)
        launch(sim.k_ipush, n, (sp.P, sp.Pn, np.int32(n), np.int64(sp.cap), *F, g.Jx, g.Jy, g.Jz, np.float64(eps),
                                np.float64(sp.q / sp.m), np.float64(0.0), np.int32(1),
                                np.float64(cfg.implicit.particle_tol), np.int32(0), sim.err))
        assert int(sim.err[0]) == 0
        x = np.longdouble
        a0, av, t0 = u0.astype(x), u0.astype(x) + x(eps) * E.astype(x), 2 * x(eps) * B.astype(x)
        g0 = np.sqrt(1 + (a0 * a0).sum(1))
        g1 = g0
        for _ in range(100):
            t = t0 / (g0 + g1)[:, None]
            ub = (av + np.cross(av, t) + (av * t).sum(1)[:, None] * t) / (1 + (t * t).sum(1))[:, None]
            u1 = 2 * ub - a0
            g1 = np.sqrt(1 + (u1 * u1).sum(1))
        # directions: the energy closure rescales |u1| to the work along the path, here not iterated
        got, ref = cp.asnumpy(sp.Pn[3:6, :n]).T, u1.astype(np.float64)
        unit = [v / np.linalg.norm(v, axis=1)[:, None] for v in (got, ref)]
        assert np.abs(unit[0] - unit[1]).max() < 2e-15  # 4.4e-16 measured; 9e-15 with a bound 100 times looser


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
