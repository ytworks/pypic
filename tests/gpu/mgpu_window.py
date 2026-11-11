"""Helper for the moving-window tests (PIC-D22): run with ``python -m pypic.launch -n N``.

    mgpu_window.py <dims> <case>    dims like 2,2,1; case thermal | tracers

A warm plasma fills a box of 64 x 24 x 24 cells up to the vacuum margin before the x+ wall; the
window moves from t = 0 on and makes 400 shifts (more than six box lengths: every particle left
at the end was loaded by a shift).

thermal: electrons and ions. Prints the histories and, per slab of eight cells of the window, the
         electron density, the mean of u_x and the mean of u^2.
tracers: one species of test particles and no field at all: the particles move on straight lines,
         so their bits do not depend on the layout of the ranks. Prints the number of particles and
         a sum over the particles of a hash of (id, position, momentum), exact in any order.
"""
import json
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig, window
from pypic.config import Window
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
case = sys.argv[2]
n, dx, theta, slab = (64, 24, 24), 0.1, 4e-3, 8
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
if sim.rank == 0:
    print(json.dumps(out))
