"""Helper for the bunch test (PIC-D21): run with ``python -m pypic.launch -n N``.

    mgpu_bunch.py <dims>    dims like 2,2,1

A Gaussian electron bunch with gamma = 10 on a partner of ions at rest, in a vacuum box of
64 x 24 x 24 cells with six open faces. The bunch is loaded into a block of cells in the low half
of x, its axis on the middle of the box across x: with the rank grid 2,2,1 the two ranks of the
high half of x load no particle, and the whole bunch crosses the rank face along x during the run.
Prints the histories, the global id sums and the number of bunch particles in the high half of x.
"""
import json
import math
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
gamma = 10.0


def bunch(u):
    shape = {"shape": "gauss", "width": 0.8}  # cut 7.7 cells from its centre
    return [Beam(1.0, drift=(u, 0.0, 0.0), profile=[{**shape, "axis": a, "center": 6.0} for a in "xyz"])]


cfg = SimConfig(n_cells=(64, 24, 24), dx=(0.5,) * 3, dt=0.25, n_steps=60, dims=dims, npml=6,
                boundaries=(("open", "open"),) * 3,
                species=[SpeciesConfig("b", -1, 1, bunch(math.sqrt(gamma ** 2 - 1.0))),
                         SpeciesConfig("p", 1, 1836, bunch(0.0))],
                ppc=4, diag_interval=5, gauss_interval=5)
sim = Simulation(cfg, verbose=False)
loaded = float(sim.species[0].n)
sim.run()
b = sim.species[0]
ids = b.ids()
# sums that are exact in float64 in any order (the order of the particles is not fixed): the ids
# and the squares of their two low 16-bit halves
low, high = (((ids.astype(cp.int64) >> s) & 0xFFFF).astype(cp.float64) for s in (0, 16))
far = float(cp.count_nonzero(b.pos(0) >= 16.0))
s = sim.comm.allreduce(cp.stack([ids.sum(), (low * low).sum(), (high * high).sum(), cp.asarray(float(ids.size)),
                                 cp.asarray(far), cp.asarray(loaded)]))
if sim.rank == 0:
    h = sim.history
    keys = ("kinetic", "electric", "magnetic", "n_particles", "gauss", "divB")
    print(json.dumps({**{k: h[k] for k in keys}, "ids": [float(c) for c in cp.asnumpy(s[:4])],
                      "far": float(s[4]), "bunch": float(s[5])}))
