"""Multi-GPU regression cases (run with ``python -m pypic.launch -n 2``).

    edge          a particle exactly on the upper reflecting wall (x = L) of a decomposed axis must
                  stay with the edge rank (it was dropped by the migration, PR #2)
    implicit_err  an implicit-step error raised on one rank only must be raised on every rank, so a
                  caller that catches it can continue with collectives (it used to hang)
"""
import sys

import cupy as cp

from pypic import SimConfig, SpeciesConfig
from pypic.simulation import ImplicitStepError, Simulation

case = sys.argv[1]
bc = (("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))


def put(sp, pos, u, w=1e-3):
    sp.ensure_capacity(1)
    sp.P[:, 0] = 0
    sp.P[0:3, 0] = cp.asarray(pos)
    sp.P[3:6, 0] = cp.asarray(u)
    sp.P[6, 0] = w
    sp.n = 1


if case == "edge":
    cfg = SimConfig(n_cells=(16, 4, 4), dx=(0.25,) * 3, dt=0.1, species=[SpeciesConfig("e", -1, 1, [])],
                    boundaries=bc, diag_interval=0, gauss_interval=0, sort_interval=0)
    sim = Simulation(cfg, verbose=False)
    L = sim.grid.L[0]
    if sim.decomp.axes[0].hi is None:  # the rank owning the upper wall
        put(sim.species[0], [L, 0.5, 0.5], [0.0, 0.0, 0.0])
    sim.run(3)
    print(f"RESULT particles={sim.total_particles()}", flush=True)
elif case == "implicit_err":
    cfg = SimConfig(n_cells=(16, 4, 4), dx=(0.25,) * 3, dt=0.2, solver="implicit",
                    species=[SpeciesConfig("e", -1, 1, []), SpeciesConfig("p", 1, 1, [])], boundaries=bc,
                    diag_interval=0, gauss_interval=0, sort_interval=0)
    sim = Simulation(cfg, verbose=False)
    for sp in sim.species:
        put(sp, [1.0 + 2.0 * sim.rank, 0.5, 0.5], [0.1, 0.0, 0.0])
    # fault injection: rank 0 reports an energy-closure failure (err = 3), rank 1 does not
    real = sim.k_ipush

    def faulty(grid, block, args):
        real(grid, block, args)
        if sim.rank == 0:
            sim.err[0] = 3

    sim.k_ipush = faulty
    try:
        sim.step()
        raised = False
    except ImplicitStepError:
        raised = True
    n = sim.total_particles()  # a collective after the (caught) error
    print(f"RESULT rank={sim.rank} raised={raised} particles={n}", flush=True)
