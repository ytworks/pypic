"""Helper for the moving-window tests (PIC-D22): run with ``python -m pypic.launch -n N``.

    mgpu_window.py <dims> <case> [<file>]    dims like 2,2,1; case thermal | tracers | laser | static |
                                             pulse | pulse_static

thermal, tracers: a warm plasma fills a box of 64 x 24 x 24 cells up to the vacuum margin before the x+ wall; the
window moves from t = 0 on and makes 400 shifts (more than six box lengths: every particle left
at the end was loaded by a shift).

thermal: electrons and ions. Prints the histories and, per slab of eight cells of the window, the
         electron density, the mean of u_x and the mean of u^2.
tracers: one species of test particles and no field at all: the particles move on straight lines,
         so their bits do not depend on the layout of the ranks. Prints the number of particles and
         a sum over the particles of a hash of (id, position, momentum), exact in any order. At the
         end every rank whose arrays hold the nodes around the middle of the window along x (two
         ranks, when x is split) deposits the charge of 32 single particles there, one at a time:
         prints the largest difference between two ranks at a node, the sum of the charge and the
         number of ranks. The cells of this case are 0.0999 long: with 0.1 the inverse is the round
         number 10 and the index coordinate of a particle has no rounding to differ by.

laser:   a focused pulse in a cold plasma, open faces across the window (absorbing layers that
         the shift moves through every rank), 200 shifts. Writes E on the cells of the window
         farther than 32 cells from x- to <file> (.npy), and output (a field dump, a slice, a
         lineout, particles and a histogram, with regions in the window) to the directory
         ``laser`` next to it.
static:  the same pulse and plasma in a box at rest that ends where the window ends after its
         200 shifts. Writes E on the same laboratory cells, and the same output with the regions
         at the laboratory places that the window has at the end, to the directory ``static``.
pulse:   a plane wave in vacuum, periodic across the window. The antenna emits the same pulse
         towards -x: it runs into the lower face of the window, which comes towards it. Writes
         field dumps of E and B every eight steps to the directory ``pulse`` next to <file>.
pulse_static: the same in a box at rest that holds every place of the window (directory
         ``pulse_static``).
"""
import json
import os
import sys

import cupy as cp
import numpy as np

from pypic import Beam, SimConfig, SpeciesConfig, window
from pypic.backend import launch
from pypic.config import Diagnostic, Laser, OutputConfig, Window
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
case = sys.argv[2]


def laser_run():
    n, dx, shifts, skip, margin = (128, 32, 32), (0.05, 0.1, 0.1), 200, 32, 8
    moving = case == "laser"
    end = (n[0] - margin + shifts) * dx[0]  # the plasma of the window ends here after the last shift
    prof = None if moving else {"shape": "step", "x0": -1.0, "x1": end}
    pulse = Laser(plane_x=60.5 * dx[0], a0=0.1, omega=10.0, waist=0.8, focus_x=8.0,
                  envelope={"shape": "gaussian", "duration": 0.3, "t_peak": 1.5})
    lab = 0 if moving else shifts  # cells from the x of the diagnostics to the one of the window at the end

    def cells(c0, c1):
        return ((c0 + lab) * dx[0], (c1 + lab) * dx[0])

    kw = {"every": 256, "dtype": "float64"}  # outputs at the steps 0, 256 (85 shifts) and 512 (200 shifts)
    inner = (cells(40, 100), None, None)
    diags = [Diagnostic("dump", "fields", fields=("Ey",), **kw,  # static: every place the window dump has had
                        region=(cells(skip, n[0]) if moving else (skip * dx[0], (n[0] + shifts) * dx[0]), None, None)),
             Diagnostic("plane", "slice", fields=("Ey",), axis="z", **kw),
             Diagnostic("line", "lineout", fields=("n:e",), axis="x", average=False, **kw),
             Diagnostic("part", "particles", species=("e",), region=inner, **kw),
             Diagnostic("hist", "histogram", species=("e",), quantities=("x",), bins=(24,),
                        ranges=(cells(skip, n[0]),), region=inner, **kw)]
    cfg = SimConfig(n_cells=n if moving else (n[0] + shifts,) + n[1:], dx=dx, dt=0.0225, n_steps=512, dims=dims,
                    npml=8, ppc=4, lasers=[pulse], diag_interval=8, gauss_interval=8,
                    output=OutputConfig(dir=os.path.join(os.path.dirname(sys.argv[3]), case), diagnostics=diags),
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, profile=prof)]),
                             SpeciesConfig("i", 1, 100, [Beam(1.0, profile=prof)])],
                    boundaries=(("open", "reflecting"), ("open", "open"), ("open", "open")),
                    window=Window(start=1.5, margin=margin * dx[0]) if moving else None)
    sim = Simulation(cfg, verbose=False)
    sim.run()
    # E on the cells i0 .. i0 + nx of the box: every rank puts its part into an array of zeros
    i0, nx = (skip if moving else skip + shifts), n[0] - skip
    g, ax = sim.grid, sim.decomp.axes
    full = cp.zeros((3, nx) + n[1:])
    x0, x1 = max(i0, ax[0].phys_lo), min(i0 + nx, ax[0].phys_lo + ax[0].n_phys)
    if x1 > x0:
        tr = tuple(slice(a.phys_lo, a.phys_lo + a.n_phys) for a in ax[1:])
        for c, E in enumerate((g.Ex, g.Ey, g.Ez)):
            full[(c, slice(x0 - i0, x1 - i0)) + tr] = E[g.physical_slice()][x0 - ax[0].phys_lo:x1 - ax[0].phys_lo]
    full = sim.comm.allreduce(full.ravel()).reshape(full.shape)
    if sim.rank == 0:
        np.save(sys.argv[3], cp.asnumpy(full))
        h = sim.history
        print(json.dumps({"shift": sim.window_shift, "n": h["n_particles"][-1], "gauss": h["gauss"],
                          "divB": h["divB"]}))


def pulse_run():
    moving = case == "pulse"
    n, dx = (128 if moving else 328, 8, 8), (0.05, 0.1, 0.1)
    pulse = Laser(plane_x=60.5 * dx[0], a0=0.1, omega=10.0,
                  envelope={"shape": "gaussian", "duration": 0.3, "t_peak": 1.5})
    dump = Diagnostic("dump", "fields", fields=("E", "B"), every=8, dtype="float64")
    cfg = SimConfig(n_cells=n, dx=dx, dt=0.0225, n_steps=512, dims=dims, npml=8, lasers=[pulse], species=[],
                    boundaries=(("open", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    output=OutputConfig(dir=os.path.join(os.path.dirname(sys.argv[3]), case), diagnostics=[dump]),
                    window=Window(start=1.5, margin=8 * dx[0]) if moving else None, diag_interval=8,
                    gauss_interval=8)
    sim = Simulation(cfg, verbose=False)
    sim.run()
    if sim.rank == 0:
        print(json.dumps({"shift": sim.window_shift}))


if case in ("laser", "static", "pulse", "pulse_static"):
    pulse_run() if case.startswith("pulse") else laser_run()
    sys.exit(0)
n, dx, theta, slab = (64, 24, 24), 0.1 if case == "thermal" else 0.0999, 4e-3, 8
species = ([SpeciesConfig("e", -1, 1, [Beam(1.0, theta=theta)]), SpeciesConfig("i", 1, 100, [Beam(1.0, theta=theta)])]
           if case == "thermal" else [SpeciesConfig("t", -1, 1, [Beam(1.0, theta=theta)], test=True)])
cfg = SimConfig(n_cells=n, dx=(dx,) * 3, dt=0.045, n_steps=889, dims=dims, npml=8, ppc=16, species=species,
                boundaries=(("open", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                window=Window(start=0.0, margin=1.0), diag_interval=7, gauss_interval=7)
sim = Simulation(cfg, verbose=False)
sim.run()
e = sim.species[0]
P = e.P[:, : e.n]
out = {"shift": sim.window_shift, "n": sim.total_particles(), "front": window.front_cell(cfg)}
if case == "thermal":
    cell = cp.floor(P[0] / dx).astype(cp.int64) - sim.window_shift
    rows = cp.stack([cp.bincount(cell // slab, weights=w, minlength=n[0] // slab)
                     for w in (P[6], P[6] * P[3], P[6] * (P[3] ** 2 + P[4] ** 2 + P[5] ** 2))])
    rows = cp.asnumpy(sim.comm.allreduce(rows.ravel())).reshape(3, -1)
    out.update({k: sim.history[k] for k in ("kinetic", "electric", "magnetic", "window_out", "budget", "gauss",
                                            "divB", "n_particles")})
    w = rows[0].clip(1e-300)
    out.update(density=(rows[0] / (slab * dx ** 3 * n[1] * n[2])).tolist(), ux=(rows[1] / w).tolist(),
               u2=(rows[2] / w).tolist(), theta=theta, slab=slab)
else:
    h = cp.zeros(e.n, dtype=cp.uint64)
    for r in (7, 0, 1, 2, 3, 4, 5):  # splitmix64 over the bits of the rows
        h = (h ^ P[r].view(cp.uint64)) + cp.uint64(0x9E3779B97F4A7C15)
        h = (h ^ (h >> cp.uint64(30))) * cp.uint64(0xBF58476D1CE4E5B9)
        h = (h ^ (h >> cp.uint64(27))) * cp.uint64(0x94D049BB133111EB)
        h = h ^ (h >> cp.uint64(31))
    # sums of the four 16-bit pieces: exact in float64 whatever the order
    parts = cp.stack([((h >> cp.uint64(16 * k)) & cp.uint64(0xFFFF)).astype(cp.float64).sum() for k in range(4)])
    out["hash"] = cp.asnumpy(sim.comm.allreduce(parts)).tolist()
    out["id_max"] = sim.comm.allreduce_scalar(float(P[7].max()) if e.n else 0.0, "max")
    # one particle at a time (no order of a sum), in the two cells at the middle plane of the
    # window along x, on the nodes of these cells
    g, K, mid, cy, cz = sim.grid, 32, n[0] // 2, 5, 6  # the cell (cy, cz): far from the rank faces along y and z
    lo = [c + g.g + a.pml_lo - a.phys_lo for c, a in zip((mid - 1, cy, cz), sim.decomp.axes)]
    sl = tuple(slice(i, i + m) for i, m in zip(lo, (3, 2, 2)))
    held = all(i >= 0 and i + m <= s for i, m, s in zip(lo, (3, 2, 2), g.shape))
    vals = cp.zeros((K, 3, 2, 2))
    one = cp.zeros((8, 1))
    one[1], one[2], one[6] = (cy + 0.3) * dx, (cz + 0.7) * dx, 1.0
    for k in range(K if held else 0):
        one[0] = (sim.window_shift + mid - 1.0 + 2.0 * ((0.137 + 0.6180339887 * k) % 1.0)) * dx
        rho = cp.zeros(g.shape)
        launch(sim.k_rho, 1, (one, np.int32(1), np.int64(1), np.float64(1.0), np.float64(0.0), rho))
        vals[k] = rho[sl]
    top = sim.comm.allreduce(vals if held else cp.full(vals.shape, -cp.inf), "max")
    low = sim.comm.allreduce(vals if held else cp.full(vals.shape, cp.inf), "min")
    out["probe_spread"], out["probe_charge"] = float(cp.max(top - low)), float(top.sum())
    out["probe_ranks"] = int(sim.comm.allreduce_scalar(float(held)))
if sim.rank == 0:
    print(json.dumps(out))
