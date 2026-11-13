"""Helper for the restart tests (PIC-D25, OUT-D6): run with ``python -m pypic.launch -n N``.

    mgpu_restart.py <dims> <case> <mode> <dir> [seed]

case: two_stream (periodic counter-streaming electrons, every output diagnostic) | open_flow (six open
faces, x- paired: the bank) | compact (open flow with filter, Marder zone, sphere, dipole and the
compact storage) | window (a warm plasma in a moving window) | implicit (the implicit solver with an
absorbing layer).
mode: full (2N steps, a checkpoint every N) | part (N + 3 steps, a checkpoint at N: a run that stops
after its checkpoint) | resume (continue to 2N from the latest checkpoint in <dir>/checkpoint).
``seed`` changes the seed (a configuration that a restart refuses). Output goes to <dir>. Prints
the histories, the particle numbers and exact sums of their ids, and the charge in the bank; rank 0
writes E and B on the physical cells of the box to <dir>/fields.npy.
"""
import json
import math
import os
import sys

import cupy as cp
import numpy as np

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Body, Checkpoint, CurrentFilter, Dipole, ExternalFields, Injection, OutputConfig, Window
from pypic.grid import FIELD_NAMES
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
case, mode, out = sys.argv[2], sys.argv[3], sys.argv[4]
assert mode in ("full", "part", "resume")
N = {"implicit": 10, "compact": 10, "window": 30}.get(case, 20)
diags = []
if case == "two_stream":
    diags = [
        dict(name="fields", kind="fields", fields=["E", "B", "J", "rho", "n:e"], every=10),
        dict(name="slice", kind="slice", fields=["Ex", "n:e"], axis="x", every=5),
        dict(name="lineout", kind="lineout", fields=["Ex", "Ey"], axis="x", every=6, average=True),
        dict(name="modes", kind="modes", fields=["Ex"], modes=[[1, 0, 0], [2, 0, 0], [1, 1, 0]], every=6,
             average=True),
        dict(name="kspec", kind="kspectrum", fields=["Ex"], every=10),
        dict(name="phase", kind="histogram", quantities=["x", "ux"], bins=[16, 8], every=10),
        dict(name="sample", kind="particles", every=20),
        dict(name="tracks", kind="tracks", species=["e"], stride=13, every=2, with_fields=True),
    ]
common = dict(dims=dims, n_steps=2 * N if mode != "part" else N + 3, checkpoint=Checkpoint(every=N, keep=2),
              restart=os.path.join(out, "checkpoint") if mode == "resume" else None,
              output=OutputConfig(dir=out, diagnostics=diags, scalars=True))
beams_e = [Beam(0.5, (0.3, 0.05, 0), 0.01), Beam(0.5, (-0.3, 0, 0.05), 0.01)]
if case == "two_stream":
    cfg = SimConfig(n_cells=(32, 8, 8), dx=(0.2,) * 3, dt=0.08, ppc=4, diag_interval=2, gauss_interval=4,
                    species=[SpeciesConfig("e", -1, 1, beams_e), SpeciesConfig("i", 1, 25, [Beam(1.0, theta=0.001)])],
                    **common)
elif case == "implicit":
    cfg = SimConfig(n_cells=(32, 16, 16), dx=(0.2,) * 3, dt=0.3, solver="implicit", npml=6, ppc=4,
                    boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("open", "open")),
                    species=[SpeciesConfig("e", -1, 1, beams_e), SpeciesConfig("i", 1, 25, [Beam(1.0, theta=0.001)])],
                    diag_interval=2, gauss_interval=4, **common)
elif case == "window":
    th = 4e-3
    cfg = SimConfig(n_cells=(32, 8, 8), dx=(0.1,) * 3, dt=0.045, npml=8, ppc=8,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, theta=th)]),
                             SpeciesConfig("i", 1, 100, [Beam(1.0, theta=th)])],
                    boundaries=(("open", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    window=Window(start=0.0, margin=1.0), diag_interval=3, gauss_interval=3, **common)
else:  # open flow (tests/gpu/mgpu_open_flow.py: inject, compact)
    extra = {"n_cells": (16, 16, 8)}
    dipoles = []
    if case == "compact":
        extra = {"n_cells": (20, 20, 20), "current_filter": CurrentFilter(1, True), "marder_zone": 0.125,
                 "body": Body(center=(5.0, 5.0, 5.0), radius=0.5), "precision": "compact"}
        dipoles = [Dipole(moment=(0.0, 0.01, 0.03), center=(5.05, 4.95, 5.0), core=0.4)]
    u = (0.1, 0.1, 0.05)
    gamma = math.sqrt(1.0 + sum(c * c for c in u))
    v, bz = [c / gamma for c in u], 0.2
    external = ExternalFields(E=(-v[1] * bz, v[0] * bz, 0.0), B=(0.0, 0.0, bz), dipoles=dipoles)
    cfg = SimConfig(dx=(0.5,) * 3, dt=0.2, npml=6, **extra, boundaries=(("open", "open"),) * 3,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, u, 0.01)]),
                             SpeciesConfig("i", 1, 16, [Beam(1.0, u, 0.01)])],
                    ppc=4, diag_interval=2, gauss_interval=4, external_fields=external,
                    injection=Injection(faces={"x-": "paired", "y-": "flux", "y+": "flux", "z-": "flux",
                                               "z+": "flux"}, neutralizer="e"), **common)
if sys.argv[5:] == ["seed"]:
    cfg.seed += 1
sim = Simulation(cfg, verbose=False)
start = sim.step_count
sim.run()
sim.close()
res = {"start": start, "end": sim.step_count, "shift": sim.window_shift}
ids = []
for sp in sim.species:
    i = sp.ids()
    low, high = (((i.astype(cp.int64) >> b) & 0xFFFF).astype(cp.float64) for b in (0, 16))
    ids.append(cp.stack([i.sum(), (low * low).sum(), (high * high).sum(), cp.asarray(float(i.size))]))
s = cp.asnumpy(sim.comm.allreduce(cp.concatenate(ids + [sim.bank.sum()[None]])))
# E and B on the physical cells of the box: every rank puts its block into an array of zeros
g, ax = sim.grid, sim.decomp.axes
full = cp.zeros((6,) + tuple(cfg.n_cells))
blk = tuple(slice(a.phys_lo, a.phys_lo + a.n_phys) for a in ax)
for c, name in enumerate(FIELD_NAMES):
    full[(c,) + blk] = getattr(g, name)[g.physical_slice()]
full = sim.comm.allreduce(full.ravel()).reshape(full.shape)
if sim.rank == 0:
    np.save(os.path.join(out, "fields.npy"), cp.asnumpy(full))
    res.update({k: sim.history[k] for k in ("step", "kinetic", "electric", "magnetic", "total", "budget",
                                            "injected_charge", "lost_charge", "n_injected", "absorbed_charge",
                                            "n_absorbed", "boundary_flux", "dissipated", "window_out",
                                            "n_particles", "gauss")})
    res["ids"] = s[:-1].tolist()
    res["bank"] = float(s[-1])
    print(json.dumps(res))
