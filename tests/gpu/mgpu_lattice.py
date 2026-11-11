"""Helper for the lattice test (PIC-D23): run with ``python -m pypic.launch -n N``.

    mgpu_lattice.py <dims>    dims like 2,2,1

A warm electron-ion plasma in a periodic box of 16 x 8 x 8 cells of three different sizes, loaded
on a 2 x 4 x 2 lattice. Prints the number of electron coordinates that differ from the place that
the id of the particle gives (global cell and slot), the largest deviation of the charge density
of the electrons alone from -1 at the nodes of the Gauss check, and the number of particles.
"""
import json
import sys

import cupy as cp
import numpy as np

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.backend import launch
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
n_cells, dx, lattice, ppc = (16, 8, 8), (0.5, 0.25, 0.4), (2, 4, 2), 16
cfg = SimConfig(n_cells=n_cells, dx=dx, dt=0.1, n_steps=1, dims=dims, ppc=ppc, lattice=lattice,
                species=[SpeciesConfig("e", -1, 1, [Beam(1.0, theta=1e-3)]),
                         SpeciesConfig("i", 1, 1836, [Beam(1.0, theta=1e-6)])])
sim = Simulation(cfg, verbose=False)
e = sim.species[0]
P = e.P[:, : e.n]
ids = P[7].astype(cp.int64)  # the first stream: linear global cell * ppc + slot
slot, lin = ids % ppc, ids // ppc
cell = (lin // (n_cells[1] * n_cells[2]), lin // n_cells[2] % n_cells[1], lin % n_cells[2])
sub = (slot % lattice[0], slot // lattice[0] % lattice[1], slot // (lattice[0] * lattice[1]))
wrong = sum(int(cp.count_nonzero(P[a] != (cell[a] + (sub[a] + 0.5) / lattice[a]) * dx[a])) for a in range(3))
g = sim.grid
rho = cp.zeros(g.shape)
launch(sim.k_rho, e.n, (e.P, np.int32(e.n), np.int64(e.cap), np.float64(e.q), np.float64(0.0), rho))
g.reduce_ghosts((rho,))
dev = float(cp.max(cp.abs(rho[g.gauss_slice(sim.gauss_retreat)] + 1.0)))
s = sim.comm.allreduce(cp.asarray([float(wrong), float(e.n)]))
dev = sim.comm.allreduce_scalar(dev, "max")
if sim.rank == 0:
    print(json.dumps({"wrong": float(s[0]), "n": float(s[1]), "density_deviation": dev}))
