"""Helper for the open-flow test: run with ``python -m pypic.launch -n N``.

    mgpu_open_flow.py <dims> <variant>    dims like 2,2,1; variant inject | inject_filter | body | mixed

A magnetised, oblique, thermal electron-ion flow through a box with six open faces (PIC-D13):
x- is the paired inflow face, the y and z faces inject on their own (``flux``), x+ only removes.
Particles enter, leave and are returned on every rank and at the edges between injecting faces.
``inject_filter`` adds the current filter and the Marder zone (PIC-D15, PIC-D16) in a box of
8 cells along x: with two ranks along x the zone ends at the rank face and the filter reaches across.
``body`` puts an absorbing sphere (PIC-D14) into the box of ``inject``, centred on the edge shared by
the four ranks, with ``body.relax`` on the edges inside it. ``mixed`` is the combination of a bow
shock: filter, Marder zone, sphere, a dipole inside it (PIC-D12) and the float32 gather (PIC-D11),
in a box of 20 cells per axis that holds the nodes left out of the Gauss check around the sphere.
Prints the histories, the charge kept in the bank of the paired face and the global id sums.
"""
import json
import math
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Body, CurrentFilter, Dipole, ExternalFields, Injection
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
variant = sys.argv[2]
assert variant in ("inject", "inject_filter", "body", "mixed")
extra = {"n_cells": (16, 16, 8)}
if variant == "inject_filter":
    extra = {"n_cells": (8, 16, 12), "current_filter": CurrentFilter(1, True), "marder_zone": 0.125}
dipoles = []
if variant == "body":
    extra["body"] = Body(center=(4.0, 4.0, 2.0), radius=0.6, relax=1.0)
if variant == "mixed":
    extra = {"n_cells": (20, 20, 20), "current_filter": CurrentFilter(1, True), "marder_zone": 0.125,
             "body": Body(center=(5.0, 5.0, 5.0), radius=0.5), "precision": "mixed"}
    dipoles = [Dipole(moment=(0.0, 0.01, 0.03), center=(5.05, 4.95, 5.0), core=0.4)]
u = (0.1, 0.1, 0.05)  # slow along x: some ions leave through the paired face and are paid for
gamma = math.sqrt(1.0 + sum(c * c for c in u))
v, bz = [c / gamma for c in u], 0.2
# E = -v x B: the uniform flow is force-free
external = ExternalFields(E=(-v[1] * bz, v[0] * bz, 0.0), B=(0.0, 0.0, bz), dipoles=dipoles)
cfg = SimConfig(dx=(0.5,) * 3, dt=0.2, n_steps=40, dims=dims, npml=6, **extra,
                boundaries=(("open", "open"),) * 3,
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, u, 0.01)]),
                         SpeciesConfig("i", 1, 16, [Beam(1.0, u, 0.01)])],
                ppc=4, diag_interval=5, gauss_interval=5, external_fields=external,
                injection=Injection(faces={"x-": "paired", "y-": "flux", "y+": "flux", "z-": "flux", "z+": "flux"},
                                    neutralizer="e"))
sim = Simulation(cfg, verbose=False)
sim.run()
ids = cp.concatenate([sp.ids() for sp in sim.species])
# sums that are exact in float64 in any order (the order of the particles is not fixed): the ids
# and the squares of their two low 16-bit halves
low, high = (((ids.astype(cp.int64) >> b) & 0xFFFF).astype(cp.float64) for b in (0, 16))
s = sim.comm.allreduce(cp.stack([ids.sum(), (low * low).sum(), (high * high).sum(), cp.asarray(float(ids.size)),
                                 sim.bank.sum()]))
if sim.rank == 0:
    h = sim.history
    keys = ("kinetic", "electric", "lost_charge", "injected_charge", "n_injected", "absorbed_kinetic",
            "absorbed_charge", "n_absorbed", "gauss")
    print(json.dumps({**{k: h[k] for k in keys}, "ids": [float(c) for c in cp.asnumpy(s[:4])], "bank": float(s[4])}))
