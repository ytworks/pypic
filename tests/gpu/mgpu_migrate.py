"""Helper for the migration test: run with ``python -m pypic.launch -n N``.

    mgpu_migrate.py <dims> <variant>    dims like 2,2,1; variant plain | overflow | idle | mixed | dipole | compact

Hot, drifting electrons and positrons on two decomposed axes, so particles leave through faces and
corners every step. ``overflow`` shrinks the flag lists so that some ranks overflow and every rank
takes the full-compaction path. ``idle`` makes the positrons cold and at rest, so some steps
flag particles of one species only, with the empty-input semantics of CuPy 14.0's bincount
(``max`` of an empty array raises) in force. ``mixed`` uses the float32 field gather (PIC-D11).
``dipole`` adds external fields (PIC-D12: an off-centre dipole, uniform E and B) to ``mixed``, and the
absorbing sphere that a dipole needs around its core (PIC-D14). ``compact`` uses the compact
particle storage (PIC-D19) and also prints a sum over the loaded lattice points.
Prints the energy history and the global count and id sums.
"""
import json
import sys

import cupy as cp

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Body, Dipole, ExternalFields
from pypic.particles import TWO32, Species
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
external, body = ExternalFields(), None
if variant == "dipole":
    external = ExternalFields(E=(0.0, 0.05, 0.0), B=(0.0, 0.0, 0.2),
                              dipoles=[Dipole(moment=(0.0, 0.02, 0.05), center=(1.3, 2.1, 0.7), core=0.3)])
    body = Body(center=(1.3, 2.1, 0.8), radius=0.45)
cfg = SimConfig(n_cells=(16, 16, 8), dx=(0.2,) * 3, dt=0.1, n_steps=30, dims=dims,
                boundaries=(("periodic", "periodic"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.6, 0.5, 0.1), 0.05)]),
                         SpeciesConfig("p", 1, 1, [Beam(1.0, positron_u, positron_theta)])],
                ppc=4, diag_interval=5, gauss_interval=5, external_fields=external, body=body,
                precision={"mixed": "mixed", "dipole": "mixed", "compact": "compact"}.get(variant, "fp64"))
sim = Simulation(cfg, verbose=False)
lattice = None
if variant == "compact":
    # the loaded positions as a sum that is exact in any order: the index coordinates in 2^-32 of a
    # cell are integers; each is multiplied by odd numbers made of the id and the axis and added
    # modulo 2^64. The ranks add the four 16-bit pieces of their sums as float64.
    h = cp.zeros((), dtype=cp.int64)
    for sp in sim.species:
        key = sp.ids().astype(cp.int64) * cp.int64(2) + cp.int64(1)
        for a in range(3):
            h = h + ((sp.index_pos(a) * TWO32).astype(cp.int64) * key * cp.int64(2 * a + 3)).sum(dtype=cp.int64)
    pieces = cp.stack([((h >> cp.int64(16 * b)) & cp.int64(0xFFFF)).astype(cp.float64) for b in range(4)])
    pieces = cp.asnumpy(sim.comm.allreduce(pieces))
    lattice = sum(int(p) << (16 * b) for b, p in enumerate(pieces)) % 2**64
sim.run()
ids = cp.concatenate([sp.ids() for sp in sim.species])
s = sim.comm.allreduce(cp.stack([ids.sum(), (ids * ids).sum(), cp.asarray(float(ids.size))]))
if sim.rank == 0:
    h = sim.history
    print(json.dumps({"kinetic": h["kinetic"], "electric": h["electric"], "gauss": h["gauss"],
                      "ids": [float(v) for v in cp.asnumpy(s)], "lattice": lattice}))
