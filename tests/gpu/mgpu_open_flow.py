"""Helper for the open-flow test: run with ``python -m pypic.launch -n N``.

    mgpu_open_flow.py <dims> <variant>    dims like 2,2,1; variant inject

A magnetised, oblique, thermal electron-ion flow through a box with six open faces (PIC-D13):
x- is the paired inflow face, the y and z faces inject on their own (``flux``), x+ only removes.
Particles enter, leave and are returned on every rank and at the edges between injecting faces.
Prints the histories, the charge kept in the bank of the paired face and the global id sums.
"""
import json
import math
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import ExternalFields, Injection
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
variant = sys.argv[2]
assert variant == "inject"
u = (0.1, 0.1, 0.05)  # slow along x: some ions leave through the paired face and are paid for
gamma = math.sqrt(1.0 + sum(c * c for c in u))
v, bz = [c / gamma for c in u], 0.2
external = ExternalFields(E=(-v[1] * bz, v[0] * bz, 0.0), B=(0.0, 0.0, bz))  # E = -v x B: the flow is force-free
cfg = SimConfig(n_cells=(16, 16, 8), dx=(0.5,) * 3, dt=0.2, n_steps=40, dims=dims, npml=6,
                boundaries=(("open", "open"),) * 3,
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, u, 0.01)]),
                         SpeciesConfig("i", 1, 16, [Beam(1.0, u, 0.01)])],
                ppc=4, diag_interval=5, gauss_interval=5, external_fields=external,
                injection=Injection(faces={"x-": "paired", "y-": "flux", "y+": "flux", "z-": "flux", "z+": "flux"},
                                    neutralizer="e"))
sim = Simulation(cfg, verbose=False)
sim.run()
ids = cp.concatenate([sp.P[7, : sp.n] for sp in sim.species])
s = sim.comm.allreduce(cp.stack([ids.sum(), (ids * ids).sum(), cp.asarray(float(ids.size)), sim.bank.sum()]))
if sim.rank == 0:
    h = sim.history
    keys = ("kinetic", "electric", "lost_charge", "injected_charge", "n_injected", "gauss")
    print(json.dumps({**{k: h[k] for k in keys}, "ids": [float(c) for c in cp.asnumpy(s[:3])], "bank": float(s[3])}))
