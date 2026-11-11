"""Helper for the radial-beam test (PIC-D17): run with ``python -m pypic.launch -n N``.

    mgpu_radial.py <dims>    dims like 2,2,1

The set-up of a blast wave in a small box of 16 x 24 x 12 cells with six open faces: an ambient
plasma at rest that enters through every face (``flux``), an absorbing sphere, and around it a
shell of two more species that drifts outwards. The shell is loaded into a block of cells in the
low half of y: with the rank grid 2,2,1 the block lies across the rank face along x, the sphere
on that face, and the two ranks of the high half load no shell particle and receive them later.
Prints the histories, the global id sums and the number of shell particles in the high half of y.
"""
import json
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Body, Injection
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
origin = (4.0, 3.0, 3.0)
shell = {"axis": "r", "shape": "gauss", "center": 1.2, "width": 0.3}


def debris(theta):
    return [Beam(2.0, theta=theta, profile=shell, origin=origin, drift_radial=0.3) for _ in range(2)]


cfg = SimConfig(n_cells=(16, 24, 12), dx=(0.5,) * 3, dt=0.2, n_steps=40, dims=dims, npml=6,
                boundaries=(("open", "open"),) * 3, body=Body(center=origin, radius=0.6),
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, theta=0.01)]),
                         SpeciesConfig("i", 1, 16, [Beam(1.0, theta=0.01 / 16)]),
                         SpeciesConfig("de", -1, 1, debris(0.01)),
                         SpeciesConfig("di", 1, 16, debris(0.01 / 16))],
                ppc=4, diag_interval=5, gauss_interval=5,
                injection=Injection(faces={f: "flux" for f in ("x-", "x+", "y-", "y+", "z-", "z+")},
                                    species=("e", "i")))
sim = Simulation(cfg, verbose=False)
loaded = [float(sp.n) for sp in sim.species[2:]]
sim.run()
ids = cp.concatenate([sp.P[7, : sp.n] for sp in sim.species])
# sums that are exact in float64 in any order (the order of the particles is not fixed): the ids
# and the squares of their two low 16-bit halves
low, high = (((ids.astype(cp.int64) >> b) & 0xFFFF).astype(cp.float64) for b in (0, 16))
far = sum(int(cp.count_nonzero(sp.P[1, : sp.n] >= 6.0)) for sp in sim.species[2:])
s = sim.comm.allreduce(cp.stack([ids.sum(), (low * low).sum(), (high * high).sum(), cp.asarray(float(ids.size)),
                                 cp.asarray(float(far)), cp.asarray(sum(loaded))]))
if sim.rank == 0:
    h = sim.history
    keys = ("kinetic", "electric", "magnetic", "lost_charge", "injected_charge", "n_injected", "absorbed_kinetic",
            "absorbed_charge", "n_absorbed", "n_particles", "gauss")
    print(json.dumps({**{k: h[k] for k in keys}, "ids": [float(c) for c in cp.asnumpy(s[:4])],
                      "far": float(s[4]), "shell": float(s[5])}))
