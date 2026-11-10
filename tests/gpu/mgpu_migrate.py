"""Helper for the migration test: run with ``python -m pypic.launch -n N``.

    mgpu_migrate.py <dims> <variant>    dims like 2,2,1; variant plain | overflow | idle | mixed | dipole

Hot, drifting electrons and positrons on two decomposed axes, so particles leave through faces and
corners every step. ``overflow`` shrinks the flag lists so that some ranks overflow and every rank
takes the full-compaction path. ``idle`` makes the positrons cold and at rest, so some steps
flag particles of one species only, with the empty-input semantics of CuPy 14.0's bincount
(``max`` of an empty array raises) in force. ``mixed`` uses the float32 field gather (PIC-D11).
``dipole`` adds external fields (PIC-D12: an off-centre dipole, uniform E and B) to ``mixed``.
Prints the energy history and the global count and id sums.
"""
import json
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Dipole, ExternalFields
from pypic.particles import Species
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
variant = sys.argv[2]
if variant == "overflow":
    def tiny_flag_list(self):
        if getattr(self, "_tiny", None) is None:
            self._tiny = cp.empty(64, dtype=cp.int32)
        return self._tiny

    Species.flag_list = tiny_flag_list
if variant == "idle":
    _bincount = cp.bincount

    def bincount_cupy14_0(x, weights=None, minlength=None):
        if x.size == 0:
            cp.max(x)  # CuPy 14.0 takes the maximum first, which raises on an empty array
        return _bincount(x, weights, minlength)

    cp.bincount = bincount_cupy14_0
positron_u = (0.0, 0.0, 0.0) if variant == "idle" else (-0.5, 0.6, 0.0)
positron_theta = 1e-8 if variant == "idle" else 0.05
external = ExternalFields()
if variant == "dipole":
    external = ExternalFields(E=(0.0, 0.05, 0.0), B=(0.0, 0.0, 0.2),
                              dipoles=[Dipole(moment=(0.0, 0.02, 0.05), center=(1.3, 2.1, 0.7), core=0.4)])
cfg = SimConfig(n_cells=(16, 16, 8), dx=(0.2,) * 3, dt=0.1, n_steps=30, dims=dims,
                boundaries=(("periodic", "periodic"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.6, 0.5, 0.1), 0.05)]),
                         SpeciesConfig("p", 1, 1, [Beam(1.0, positron_u, positron_theta)])],
                ppc=4, diag_interval=5, gauss_interval=5, external_fields=external,
                precision="mixed" if variant in ("mixed", "dipole") else "fp64")
sim = Simulation(cfg, verbose=False)
sim.run()
ids = cp.concatenate([sp.P[7, : sp.n] for sp in sim.species])
s = sim.comm.allreduce(cp.stack([ids.sum(), (ids * ids).sum(), cp.asarray(float(ids.size))]))
if sim.rank == 0:
    h = sim.history
    print(json.dumps({"kinetic": h["kinetic"], "electric": h["electric"], "gauss": h["gauss"],
                      "ids": [float(v) for v in cp.asnumpy(s)]}))
