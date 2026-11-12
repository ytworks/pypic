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
from pypic.body import relax_edges
from pypic.config import Body, CurrentFilter, ExternalFields, Injection
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
    blk = np.zeros((8, len(pos)))
    blk[:3], blk[3:6], blk[6] = pos.T, u.T, w
    sp.n = 0
    sp.append(blk)


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
        ux, uy = float(sp.mom(0)[0]), float(sp.mom(1)[0])
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
        xs.append(float(sp.pos(0)[0]))
    xs = np.unwrap(np.array(xs) / sim.grid.L[0] * 2 * np.pi) * sim.grid.L[0] / (2 * np.pi)
    vd = (xs[-1] - x0) / (nst * cfg.dt)
    assert abs(vd - E / B) < 0.01 * E / B


@pytest.mark.parametrize("precision", ["fp64", "mixed"])
def test_external_uniform_fields_give_the_orbits_of_grid_fields(precision):
    """Test particles in external uniform E and B (PIC-D12) move bit for bit as in the same
    fields held by the grid (B0, E0): the initial half kick and the gather use the total field."""
    E, B = (0.0, 0.03, 0.01), (0.02, 0.0, 0.5)
    rng = np.random.default_rng(3)
    pts = [[*rng.uniform(0.5, 3.5, 3), *rng.uniform(-2.0, 2.0, 3)] for _ in range(50)]
    orbits = []
    for fields in (dict(B0=B, E0=E), dict(external_fields=ExternalFields(E=E, B=B))):
        cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.5,) * 3, dt=0.05, precision=precision, diag_interval=0,
                        gauss_interval=0, species=[SpeciesConfig("t", -1, 1, [], test=True, particles=pts)], **fields)
        sim = Simulation(cfg, comm=SelfComm(), verbose=False)
        for _ in range(200):
            sim.step()
        orbits.append(cp.asnumpy(sim.species[0].rows()[:6, :50]))
        sim.close()
    assert np.array_equal(*orbits)
    assert np.max(np.abs(orbits[0][3:] - np.asarray(pts).T[3:])) > 0.1  # the fields did act


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
    x, w = cp.asnumpy(sp.pos(0)), cp.asnumpy(sp.weight())
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
    assert e.n == 12 * 4 * 4 * 8 and float(cp.max(e.pos(0))) < 3.0
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


@pytest.mark.parametrize("bc", ["periodic", "reflecting", "mixed"])
def test_current_filter_kernels_match_the_array_expressions(bc):
    """The filter kernels give the result of the array expressions pass by pass, for J (three
    staggers) and node data: bitwise at conducting faces, and to rounding on the owned cells of
    periodic axes, whose passes are one composite stencil (here wider than the 3 cells of z). The
    mixed case has a conducting axis between two periodic ones."""
    kinds = ("periodic", "reflecting", "periodic") if bc == "mixed" else (bc,) * 3
    cfg = SimConfig(n_cells=(12, 10, 3), dx=(0.1,) * 3, dt=0.05, species=[], diag_interval=0, gauss_interval=0,
                    boundaries=tuple((k, k) for k in kinds), npml=2)
    g = Simulation(cfg, comm=SelfComm(), verbose=False).grid
    rng = np.random.default_rng(7)
    own = tuple(slice(g.g, g.g + ax.n_int) for ax in g.d.axes)
    for staggers in (Simulation.J_STAGGER, ((0, 0, 0),)):
        for comp in (True, False):  # an even and an odd number of passes per axis
            A = [cp.asarray(rng.standard_normal(g.shape)) for _ in staggers]
            B = [a.copy() for a in A]
            g.filter(A, staggers, 3, comp)
            _filter_reference(g, B, staggers, 3, comp)
            for a, b in zip(A, B):
                if bc == "reflecting":
                    assert bool(cp.all(a == b))
                else:
                    assert float(cp.abs(a[own] - b[own]).max()) < 1e-14 * float(cp.abs(b[own]).max())


@pytest.mark.parametrize("solver", ["explicit", "implicit"])
@pytest.mark.parametrize("bc", ["periodic", "reflecting"])
def test_current_filter_keeps_gauss_and_implicit_energy(solver, bc):
    """With a 4-pass compensated current filter (PIC-D10) and particles bouncing at the walls,
    Gauss's law (against the identically filtered rho) and the implicit energy stay exact."""
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


@pytest.mark.parametrize("marder", [0.0, 0.125])
def test_filter_and_marder_zone_at_open_faces_keep_gauss(marder):
    """A drifting plasma leaves through open faces with the current filter on (PIC-D15), without
    and with the Marder zone (PIC-D16): the charge left on the face nodes spreads through the
    filter and the zone works on it, and Gauss's law stays exact at every node behind them."""
    cfg = SimConfig(n_cells=(16, 16, 12), dx=(0.5,) * 3, dt=0.2, npml=6,
                    boundaries=(("open", "open"), ("open", "open"), ("periodic", "periodic")),
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.1, 0.05, 0), 0.02)]),
                             SpeciesConfig("i", 1, 16, [Beam(1.0, (0.1, 0.05, 0), 0.001)])],
                    ppc=4, diag_interval=5, gauss_interval=5, current_filter=CurrentFilter(1, True),
                    marder_zone=marder)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(40)
    assert sim.history["lost_charge"][-1] != 0.0
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    if marder:  # the zone deposits the charge of the particles near the open faces only: the same there
        g = sim.grid
        i = np.arange(g.shape[0]) - (g.g + cfg.npml)  # node index from the face x = 0 (y alike)
        near = cp.asarray((i <= sim.zone_depth) | (i >= 16 - sim.zone_depth))
        zone = near[:, None, None] | near[None, :, None]
        assert float((cp.abs(sim._charge_density(zone=True) - sim._charge_density()) * zone).max()) < 1e-12


def test_marder_zone_takes_the_charge_left_at_an_open_face():
    """A static field with div E - rho != 0 around one node of an open face and no charge
    (PIC-D16): one step of the zone diffuses it with the Marder number, the conducting wall behind
    the layer takes it all in the end, and div E never changes outside the zone."""
    nu, npml = 0.125, 3
    cfg = SimConfig(n_cells=(8, 6, 6), dx=(0.5, 0.4, 0.25), dt=0.1, npml=npml, species=[], diag_interval=0,
                    gauss_interval=0, boundaries=(("open", "open"), ("periodic",) * 2, ("periodic",) * 2),
                    marder_zone=nu)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    g = sim.grid
    face = g.g + npml  # local index of the node plane x = 0
    node = (face, g.g + 2, g.g + 3)
    # E = -grad phi with phi at one node: curl-free, div E = 0 away from the node and its neighbours
    for a, E in enumerate((g.Ex, g.Ey, g.Ez)):
        lo = list(node)
        lo[a] -= 1
        E[node] = 1.0 / g.dx[a]
        E[tuple(lo)] = -1.0 / g.dx[a]
    g.fill_ghosts((g.Ex, g.Ey, g.Ez))
    inside = (slice(face + sim.zone_depth + 1, face + 8 - sim.zone_depth), slice(None), slice(None))
    own = tuple(slice(g.g, g.g + ax.n_int) for ax in g.d.axes)
    d0 = g.div_E()
    peak = float(d0[node])
    assert not bool(cp.any(d0[inside]))
    g.marder_zone(cp.zeros(g.shape), sim.zone_depth, nu)
    h2 = 3.0 / sum(1.0 / d ** 2 for d in g.dx)
    d1 = g.div_E()
    expect = peak + nu * h2 * sum(
        (float(d0[tuple(n + (b == a) for b, n in enumerate(node))]) - 2.0 * peak
         + float(d0[tuple(n - (b == a) for b, n in enumerate(node))])) / g.dx[a] ** 2 for a in range(3))
    assert abs(float(d1[node]) - expect) < 1e-12 * abs(peak)
    for _ in range(1500):
        g.marder_zone(cp.zeros(g.shape), sim.zone_depth, nu)
    d = g.div_E()
    assert float(cp.abs(d[own]).max()) < 1e-6 * float(cp.abs(d0[own]).max())
    assert not bool(cp.any(d[inside]))


# ------------------------------------------------------------------- storage
@pytest.mark.parametrize("passes", [0, 1])
def test_absorbing_sphere_keeps_the_charge_of_the_particles_it_took(passes):
    """A drifting plasma meets an absorbing sphere (PIC-D14), without and with the current filter.
    The load leaves no particle inside and stays neutral; every particle that ends a step inside is
    removed and counted; the charge the sphere holds (div E - rho summed over the nodes left out of
    the Gauss check) is the absorbed charge; and Gauss's law is exact at every other node."""
    sphere = Body(center=(2.0, 2.1, 1.9), radius=0.5)
    kw = {"current_filter": CurrentFilter(passes, True)} if passes else {}
    cfg = SimConfig(n_cells=(20, 20, 20), dx=(0.2,) * 3, dt=0.08, body=sphere, ppc=4, diag_interval=5,
                    gauss_interval=5, species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.2, 0, 0), 0.01)]),
                                               SpeciesConfig("i", 1, 25, [Beam(1.0, (0.2, 0, 0), 0.0004)])], **kw)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)

    def inside():
        r2 = [sum((sp.pos(a) - sphere.center[a]) ** 2 for a in range(3)) for sp in sim.species]
        return sum(int((v < sphere.radius**2).sum()) for v in r2)

    n0 = sim.total_particles()
    assert inside() == 0 and 0 < 2 * 4 * 20**3 - n0 < 0.01 * 2 * 4 * 20**3  # the sphere holds 0.8 % of the box
    assert float(cp.abs(sim._charge_density()).max()) < 1e-13
    sim.run(41)  # the record of step 40 holds what this last step absorbed
    h, g = sim.history, sim.grid
    assert h["n_absorbed"][-1] == n0 - sim.total_particles() > 100 and inside() == 0
    excess = (g.div_E() - sim._charge_density())[g.gauss_slice(sim.gauss_retreat)]
    assert float(cp.abs(excess).max()) > 0.1  # the sphere is charged
    assert h["absorbed_charge"][-1] != 0.0
    assert abs(float(excess.sum()) * 0.2**3 - h["absorbed_charge"][-1]) < 1e-12
    assert np.nanmax(h["gauss"]) < 1e-12


def test_body_relax_removes_the_field_inside_the_sphere_and_keeps_its_charge():
    """``body.relax`` (PIC-D14) lets the field decay on the edges inside the sphere: the inside
    becomes an equipotential, while the sphere as a whole still holds exactly the absorbed charge
    and Gauss's law stays exact at every node outside."""
    inner = {}
    for rate in (0.0, 30.0):  # exp(-rate dt) = 0.09 per step
        sphere = Body(center=(2.0, 2.1, 1.9), radius=0.7, relax=rate)
        cfg = SimConfig(n_cells=(20, 20, 20), dx=(0.2,) * 3, dt=0.08, body=sphere, ppc=4, diag_interval=5,
                        gauss_interval=5, species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.2, 0, 0), 0.01)]),
                                                   SpeciesConfig("i", 1, 25, [Beam(1.0, (0.2, 0, 0), 0.0004)])])
        sim = Simulation(cfg, comm=SelfComm(), verbose=False)
        sim.run(41)
        h, g = sim.history, sim.grid
        edges = relax_edges(cp, [g.global_positions(a, 0.0) * 0.2 for a in range(3)], sphere, g.dx)
        inner[rate] = sum(float((e[m] ** 2).sum()) for e, m in zip((g.Ex, g.Ey, g.Ez), edges))
        excess = (g.div_E() - sim._charge_density())[g.gauss_slice(sim.gauss_retreat)]
        assert h["n_absorbed"][-1] > 100
        assert abs(float(excess.sum()) * 0.2**3 - h["absorbed_charge"][-1]) < 1e-12
        assert np.nanmax(h["gauss"]) < 1e-12
    assert 0.0 < inner[30.0] < 0.01 * inner[0.0]


def test_absorbed_particle_is_booked_with_its_energy_and_charge():
    """One electron and one test particle fly into the sphere (PIC-D14): each is counted once, the
    electron's kinetic energy and charge go to the absorbed totals, the test particle adds no
    charge and no energy, and a particle that passes outside the sphere stays."""
    cfg = SimConfig(n_cells=(16, 16, 16), dx=(0.5,) * 3, dt=0.2, diag_interval=1, gauss_interval=0,
                    body=Body(center=(4.0, 4.0, 4.0), radius=1.0),
                    species=[SpeciesConfig("e", -1, 1, []), SpeciesConfig("t", 1, 4, [], test=True)])
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    e, t = sim.species
    w, u = 1e-20, 0.5
    inject(e, [[2.6, 4.0, 4.1], [2.6, 5.3, 4.0]], [[u, 0, 0], [u, 0, 0]], w=w)  # the second one misses
    inject(t, [4.1, 2.7, 3.9], [0, u, 0], w=w)
    ke = w * (math.sqrt(1.0 + u * u) - 1.0)
    sim.run(6)  # 0.089 per step: both are inside after 5 steps
    h = sim.history
    assert (e.n, t.n) == (1, 0) and abs(float(e.pos(1)[0]) - 5.3) < 1e-12
    assert h["n_absorbed"][0] == 0 and h["n_absorbed"][-1] == 2
    assert abs(h["absorbed_kinetic"][-1] - ke) < 1e-12 * ke and abs(h["absorbed_charge"][-1] + w) < 1e-12 * w
    assert abs(h["kinetic"][-1] - ke) < 1e-12 * ke  # of the electron that is left


def test_sort_orders_by_cell_and_keeps_every_particle():
    """The sort moves whole particles: every row of the particle array and of the implicit work
    array follows the same order, and the cells come out ascending."""
    from pypic.backend import launch

    sim = plasma("implicit", sort_interval=0)
    sim.run(2)
    for sp in sim.species:
        n = sp.n
        before = [cp.asnumpy(A) for A in (sp.rows(), sp.Pn[:, :n])]
        sp.sort(sim.k_cell)
        after = [cp.asnumpy(A) for A in (sp.rows(), sp.Pn[:, :n])]
        assert not np.array_equal(before[0][7], after[0][7])  # the particles were out of order
        by_id0, by_id1 = np.argsort(before[0][7]), np.argsort(after[0][7])
        for b, a in zip(before, after):
            assert np.array_equal(b[:, by_id0], a[:, by_id1])
        key = cp.empty(n, dtype=cp.int64)
        launch(sim.k_cell, n, (*sp.store, np.int32(n), np.int64(sp.cap), key))
        assert bool((key[1:] >= key[:-1]).all())


def test_sort_of_a_growing_species_does_not_grow_the_memory_pool():
    """Regression: CuPy's pool keeps a freed block only for requests that fit, so the temporaries of
    the sorts of a species that gains particles (injection) piled up until the device was full."""
    sim = plasma("explicit", sort_interval=0)
    sp = sim.species[0]
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    n, held = sp.n, []
    for part in (0.5, 0.75, 1.0):
        sp.n = int(part * n)
        sp.sort(sim.k_cell)
        held.append(pool.total_bytes())
    assert held[1] == held[0] and held[2] == held[0], held


BIG_CAP = 360_000_000  # 6 * cap > 2**31: the weight row starts beyond int32 offsets (issue #11)


@pytest.mark.skipif(cp.cuda.Device().mem_info[0] < 32 * 2**30, reason="needs ~24 GB of free GPU memory")
def test_particle_rows_beyond_int32_offsets():
    """A capacity whose row offsets (row * cap) exceed 2**31 gives the same run as a small one."""
    ref = plasma("explicit", sort_interval=5)
    sim = plasma("explicit", sort_interval=5)
    sp = sim.species[0]
    ids = cp.asnumpy(sp.ids()).copy()
    sp.ensure_capacity(BIG_CAP)
    assert 6 * sp.cap > 2**31
    ref.run(10)
    sim.run(10)
    assert np.nanmax(sim.history["gauss"]) < 1e-12
    np.testing.assert_allclose(sim.history["total"], ref.history["total"], rtol=1e-12)
    assert np.array_equal(np.sort(cp.asnumpy(sp.ids())), np.sort(ids))
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
        x = cp.asnumpy(sp.rows()[:3])
        for a in range(3):
            assert x[a].min() >= 0.0 and x[a].max() <= L[a]


def test_open_boundary_removes_particles():
    bc = (("open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    sim = plasma("explicit", bc=bc, drift=0.5, theta=0.01)
    n0 = sim.total_particles()
    sim.run(60)
    assert sim.total_particles() < n0
    assert sim.history["lost_kinetic"][-1] > 0


# ------------------------------------------------------------------ injection
def inflow_moments(u_d, theta, n=3001):
    """A drifting Maxwell-Juttner beam crossing a plane normal to its drift, by quadrature over
    (u_n, u_perp): ``E[max(v_n, 0)]`` of the beam, and the mean and the variance of ``gamma - 1``
    of the particles that cross (weight ``v_n f``)."""
    gd, s = math.sqrt(1 + u_d * u_d), math.sqrt(theta)
    un, up = np.meshgrid(u_d + np.linspace(-14, 14, n) * s * gd, np.linspace(0, 14, n // 2) * s, indexing="ij")
    g = np.sqrt(1 + un ** 2 + up ** 2)
    e = (gd * g - u_d * un - 1.0) / theta
    f = np.exp(e.min() - e) * up
    fw = np.maximum(un / g, 0.0) * f
    mean = (fw * (g - 1)).sum() / fw.sum()
    return fw.sum() / f.sum(), mean, (fw * (g - 1) ** 2).sum() / fw.sum() - mean * mean


def test_injected_flux_is_the_inflow_of_the_drifting_maxwell_juttner_beams():
    """An injecting open face (PIC-D13) lets in ``n E[max(v_n, 0)]`` particles per area and time of
    each species, and Gauss's law holds while particles enter and leave."""
    bc = (("open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    drift, thetas, steps = 0.3, (0.01, 0.001), 50
    cfg = SimConfig(n_cells=(16, 12, 12), dx=(0.5,) * 3, dt=0.2, boundaries=bc, ppc=8, npml=6,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (drift, 0, 0), thetas[0])]),
                             SpeciesConfig("i", 1, 16, [Beam(1.0, (drift, 0, 0), thetas[1])])],
                    injection=Injection(faces={"x-": "flux"}), diag_interval=10, gauss_interval=10)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(steps)
    trials = 12 * 12 * cfg.ppc * steps  # one per slot of the cell layer outside the face and step
    for d, theta in zip(sim.sdiag, thetas):
        p = inflow_moments(drift, theta)[0] * cfg.dt / 0.5
        assert abs(float(d[8]) - trials * p) < 5 * math.sqrt(trials * p * (1 - p))
    assert np.nanmax(sim.history["gauss"]) < 1e-10


def test_cold_diagonal_inflow_through_two_faces_fills_the_edge_once():
    """A cold electron-ion beam along the x-y diagonal enters through x- and y-: the entering
    number is the volume swept into the box, the cells outside both faces are sampled once, the
    ids stay unique, and the two species enter at the same points (no field appears)."""
    v, steps = 0.5, 30  # 30 steps: nothing that entered through the edge has crossed the box yet
    u = v / math.sqrt(1 - 2 * v * v)
    bc = (("open", "open"), ("open", "open"), ("periodic", "periodic"))
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.5,) * 3, dt=0.2, boundaries=bc, ppc=8, npml=6,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (u, u, 0), 0.0)]),
                             SpeciesConfig("i", 1, 16, [Beam(1.0, (u, u, 0), 0.0)])],
                    injection=Injection(faces={"x-": "flux", "y-": "flux"}), diag_interval=10, gauss_interval=10)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(steps)
    a, (lx, ly, lz) = v * cfg.dt, sim.grid.L
    swept = steps * (lx * ly - (lx - a) * (ly - a)) * lz * cfg.ppc / sim.grid.cell_volume
    edge_trials = 8 * cfg.ppc * steps  # slots of the cells outside both faces
    p = (a / 0.5) ** 2
    for sp, d in zip(sim.species, sim.sdiag):
        assert abs(float(d[8]) - swept) < 5 * math.sqrt(swept)
        # straight orbits: floor(x / a) counts the steps since the particle crossed x = 0
        kx, ky = (np.floor(cp.asnumpy(sp.pos(c)) / a) for c in (0, 1))
        entered = np.minimum(kx, ky) < steps
        assert abs(np.count_nonzero(entered & (kx == ky)) - edge_trials * p) < 5 * math.sqrt(edge_trials * p * (1 - p))
    ids = cp.concatenate([sp.ids() for sp in sim.species])
    assert int(cp.unique(ids).size) == int(ids.size)
    assert max(sim.history["electric"]) < 1e-20


def test_paired_face_keeps_the_charge_that_crossed_it_at_zero():
    """At a paired face the neutraliser never leaves, a carrier enters with a neutraliser born at
    the same point, and a carrier that left is paid for by the next one entering without a partner:
    the charge that crossed the face is the unpaid rest in the bank, and Gauss's law holds. The
    partners have the energy of the neutraliser's inflow."""
    bc = (("open", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))
    drift, theta = 0.1, 0.01
    cfg = SimConfig(n_cells=(16, 12, 12), dx=(0.5,) * 3, dt=0.2, boundaries=bc, ppc=8, npml=6,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (drift, 0, 0), theta)]),
                             SpeciesConfig("i", 1, 4, [Beam(1.0, (drift, 0, 0), theta)])],
                    injection=Injection(faces={"x-": "paired"}, neutralizer="e"),
                    diag_interval=10, gauss_interval=10)
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run(50)
    h, (e, i) = sim.history, (cp.asnumpy(d) for d in sim.sdiag)
    bank = float(sim.bank.sum())
    crossed = e[7] + i[7] - e[5] - i[5]  # injected - lost charge, up to the last step
    assert e[5] == 0.0
    assert i[5] > 0.0 and bank < 0.25 * i[5]  # ions left through the face and most were paid for
    assert float(sim.bank.min()) >= 0.0  # and none was paid for twice
    assert abs(crossed + bank) < 1e-9
    assert abs(sum(float(sp.q * sp.weight().sum()) for sp in sim.species) - crossed) < 1e-9
    assert np.nanmax(h["gauss"]) < 1e-10
    _, mean, var = inflow_moments(drift, theta)
    w = sim.grid.cell_volume / cfg.ppc
    assert abs(e[6] / (w * e[8]) - mean) < 5 * math.sqrt(var / e[8])


# ---------------------------------------------------------------- radial beams
def shell_plasma(electrons, ions, steps=0, seed=1):
    cfg = SimConfig(n_cells=(32, 32, 32), dx=(0.5,) * 3, dt=0.2, n_steps=steps, ppc=8, seed=seed,
                    species=[SpeciesConfig("e", -1, 1, electrons), SpeciesConfig("i", 1, 16, ions)],
                    diag_interval=5, gauss_interval=5)
    return Simulation(cfg, comm=SelfComm(), verbose=False)


def test_radial_shell_has_the_mass_and_the_momentum_of_its_shape():
    """A Gaussian shell with a radial drift (PIC-D17), fixed seed. Per species: the mass is the
    integral of the shape, within 5 sigma of a bound of the sampling error (ppc uniform positions
    per cell: variance <= w0 n0 int f^2 dV); the mean proper velocity along the radius is that of a
    drifting Maxwellian, u_r (1 + 5 theta / 2), within 5 sigma (variance theta sum w^2 / (sum w)^2);
    the particles are those of a uniform load of the box where the shape is not cut, so the block
    of cells holds the whole shell. Gauss's law holds at the load and while the shell expands."""
    origin, c, wd, n0, ur, theta = (8.1, 7.9, 8.0), 3.0, 0.6, 2.0, 0.1, 1e-2

    def shell(th):
        return [Beam(n0, theta=th, origin=origin, drift_radial=ur,
                     profile={"axis": "r", "shape": "gauss", "center": c, "width": wd})]

    sim = shell_plasma(shell(theta), shell(theta / 16), steps=40)
    ref = shell_plasma([Beam(n0, theta=theta)], [Beam(n0, theta=theta / 16)])
    r = (np.arange(400000) + 0.5) * (16.0 / 400000)  # midpoint rule; the shape ends at r = 5.9
    f = np.exp(-((r - c) / wd) ** 2)
    mass = n0 * np.sum(4 * np.pi * r * r * f) * (r[1] - r[0])
    w0 = n0 * sim.grid.cell_volume / sim.cfg.ppc
    sigma = math.sqrt(w0 * n0 * np.sum(4 * np.pi * r * r * f * f) * (r[1] - r[0]))
    for sp, full, th in zip(sim.species, ref.species, (theta, theta / 16)):
        P, Q = cp.asnumpy(sp.rows()), cp.asnumpy(full.rows())
        d = P[:3] - np.asarray(origin)[:, None]
        w = P[6]
        assert abs(w.sum() - mass) < 5 * sigma
        u_r = (P[3:6] * d).sum(0) / np.sqrt((d * d).sum(0))
        mean = (w * u_r).sum() / w.sum()
        assert abs(mean - ur * (1 + 2.5 * th)) < 5 * math.sqrt(th * (w * w).sum()) / w.sum()
        rq = np.sqrt(((Q[:3] - np.asarray(origin)[:, None]) ** 2).sum(0))
        assert np.array_equal(np.sort(P[7]), np.sort(Q[7][np.exp(-((rq - c) / wd) ** 2) >= 1e-10]))
    assert sim.gauss_error()[0] < 1e-10
    sim.run()
    assert np.nanmax(sim.history["gauss"]) < 1e-10


def test_a_block_that_is_not_neutral_is_rejected():
    """Beams that load the same block of cells share positions and must cancel there (PIC-D17): a
    shell of electrons alone, with the ions uniform, leaves a charge that the initial E = 0 does not
    hold."""
    shell = [Beam(1.0, origin=(8.0, 8.0, 8.0), profile={"axis": "r", "shape": "gauss", "center": 3.0, "width": 0.6})]
    with pytest.raises(ValueError, match="not charge neutral"):
        shell_plasma(shell, [Beam(1.0)])


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


@pytest.mark.parametrize("name, dtype", [("mig_rank32", np.int32), ("mig_rank64", np.int64)])
def test_migration_lists_are_ordered_like_a_sort(name, dtype):
    """The lists of a migration are ordered by a counting kernel up to RANK_SORT_MAX keys and by a
    sort above: both give the sorted keys, with equal keys and at lengths off the warp size too."""
    from pypic.particles import RANK_SORT_MAX, ascending

    kernel = plasma().pmod.get_function(name)
    rng = np.random.default_rng(3)
    for n in (2, 33, 1000, RANK_SORT_MAX, RANK_SORT_MAX + 1):
        for bound in (n // 2 + 1, np.iinfo(dtype).max):  # with and without equal keys
            a = cp.asarray(rng.integers(0, bound, n).astype(dtype))
            ref = cp.sort(a)
            assert bool((ascending(a, kernel) == ref).all()), (n, bound)


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


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
def test_external_dipole_is_decomposition_independent():
    """A dipole and uniform external fields (PIC-D12) with the float32 gather: 1 GPU and a 2x2 rank
    grid give the same run (the dipole is computed from global indices), and Gauss's law holds."""
    script = os.path.join(HERE, "mgpu_migrate.py")
    outs = {}
    for n, dims in ((1, "1,1,1"), (4, "2,2,1")):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, dims, "dipole"],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    a, b = (np.asarray(outs[n]["kinetic"]) for n in (1, 4))
    assert np.max(np.abs(a - b) / np.abs(a)) < 1e-10
    assert np.nanmax(np.asarray(outs[4]["gauss"], dtype=float)) < 1e-12


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
@pytest.mark.parametrize("variant", ["inject", "inject_filter", "body", "mixed"])
def test_open_flow_is_decomposition_independent(variant):
    """A magnetised oblique flow through six open faces, one paired and four flux (PIC-D13): what
    enters depends on (seed, step, global cell, slot) only, so 1 GPU and a 2x2 rank grid give the
    same run, the same particles and the same charge in the bank, and Gauss's law holds. With the
    current filter and the Marder zone too (``inject_filter``, PIC-D15 and PIC-D16), with an absorbing
    sphere on the edge shared by the ranks (``body``, PIC-D14), and with all of them, a dipole in the
    sphere and the float32 gather (``mixed``)."""
    script = os.path.join(HERE, "mgpu_open_flow.py")
    outs = {}
    for n, dims in ((1, "1,1,1"), (4, "2,2,1")):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, dims, variant],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    ref, o = outs[1], outs[4]
    assert o["ids"] == ref["ids"] and o["n_injected"] == ref["n_injected"] and o["n_absorbed"] == ref["n_absorbed"]
    assert ref["n_injected"][-1] > 0 and ref["bank"] > 0.0
    assert (ref["n_absorbed"][-1] > 0) == (variant in ("body", "mixed"))
    for k in ("kinetic", "electric", "lost_charge", "injected_charge", "absorbed_kinetic", "absorbed_charge", "bank"):
        a, b = np.asarray(ref[k]), np.asarray(o[k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10, k
    assert np.nanmax(np.asarray(o["gauss"], dtype=float)) < 1e-12


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
def test_radial_shell_is_decomposition_independent():
    """A shell loaded into a block of cells (PIC-D17) around an absorbing sphere, in an ambient plasma
    that enters through six open faces: the block lies across a rank face of the 2x2 rank grid and
    two of the four ranks load no shell particle. 1 GPU and the 2x2 grid give the same run and the
    same particles, shell particles reach the ranks that loaded none, and Gauss's law holds."""
    script = os.path.join(HERE, "mgpu_radial.py")
    outs = {}
    for n, dims in ((1, "1,1,1"), (4, "2,2,1")):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), script, dims],
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr[-2000:]
        outs[n] = json.loads(r.stdout.strip().splitlines()[-1])
    ref, o = outs[1], outs[4]
    for k in ("ids", "n_injected", "n_absorbed", "n_particles", "far", "shell"):
        assert o[k] == ref[k], k
    assert ref["shell"] > 0 and ref["far"] > 0 and ref["n_absorbed"][-1] > 0
    for k in ("kinetic", "electric", "magnetic", "lost_charge", "injected_charge", "absorbed_kinetic",
              "absorbed_charge"):
        a, b = np.asarray(ref[k]), np.asarray(o[k])
        assert np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)) < 1e-10, k
    assert np.nanmax(np.asarray(o["gauss"], dtype=float)) < 1e-10


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
    u = cp.asnumpy(sp.rows()[3:6])
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
    P0 = [sp.rows().copy() for sp in sim.species]
    E0 = sim.grid.Ex.copy()
    with pytest.raises(ImplicitStepError):
        sim.step()
    for sp, p in zip(sim.species, P0):
        assert bool(cp.all(sp.rows() == p))
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
    x = float(e.pos(0)[0])
    assert 4.0 < x < 8.0 and float(e.mom(0)[0]) > 0  # reflected past the box centre
    assert sim.gauss_error()[0] < 1e-12
