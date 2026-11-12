"""Multi-GPU regression cases (run with ``python -m pypic.launch -n 2``; ``b_ghost`` with ``-n 4``).

    edge          a particle exactly on the upper reflecting wall (x = L) of a decomposed axis must
                  stay with the edge rank (it was dropped by the migration, PR #2)
    implicit_err  an implicit-step error raised on one rank only must be raised on every rank, so a
                  caller that catches it can continue with collectives (it used to hang)
    b_ghost       the ghost layer of B that the explicit update computes locally must equal the
                  neighbour's values bit for bit, also between a rank with an absorbing layer and
                  one without (their kernels differ; contracted products gave 1 ulp, issue #34)
"""
import sys

import cupy as cp
import numpy as np

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.simulation import ImplicitStepError, Simulation

case = sys.argv[1]
bc = (("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))


def put(sp, pos, u, w=1e-3):
    blk = np.zeros((8, 1))
    blk[0:3, 0], blk[3:6, 0], blk[6, 0] = pos, u, w
    sp.n = 0
    sp.append(blk)


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
elif case == "b_ghost":
    periodic = ("periodic", "periodic")
    electrons = [Beam(0.5, (0.3, 0.05, 0), 0.01), Beam(0.5, (-0.3, 0, 0.05), 0.01)]
    ions = [Beam(1.0, (0, 0, 0), 0.001)]
    cfg = SimConfig(n_cells=(32, 16, 16), dx=(0.2,) * 3, dt=0.08, n_steps=20, ppc=4, npml=6, dims=(4, 1, 1),
                    species=[SpeciesConfig("e", -1, 1, electrons), SpeciesConfig("i", 1, 25, ions)],
                    boundaries=(("open", "open"), periodic, periodic), diag_interval=0, gauss_interval=0)
    sim = Simulation(cfg, verbose=False)
    sim.run()
    g = sim.grid
    B = (g.Bx, g.By, g.Bz)
    local = [b.copy() for b in B]
    g.fill_ghosts(B)
    read = tuple(slice(g.g - 1, g.g + a.n_int + 1) for a in g.d.axes)  # the points of B a step reads
    differ = any(bool((a[read] != b[read]).any()) for a, b in zip(local, B))
    n = sim.comm.allreduce_scalar(float(differ))
    print(f"RESULT rank={sim.rank} ranks_with_a_difference={int(n)} b_is_nonzero={bool(cp.abs(B[2]).max() > 0)}",
          flush=True)
