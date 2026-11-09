"""Helper for test_multi_gpu_matches_single: run with ``python -m pypic.launch -n N``."""
import json
import sys

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.simulation import Simulation

solver = sys.argv[1]
cfg = SimConfig(n_cells=(32, 16, 16), dx=(0.2,) * 3, dt=0.08 if solver == "explicit" else 0.3, n_steps=20,
                solver=solver, boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("open", "open")),
                npml=6, species=[SpeciesConfig("e", -1, 1, [Beam(0.5, (0.3, 0.05, 0), 0.01),
                                                            Beam(0.5, (-0.3, 0, 0.05), 0.01)]),
                                 SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                ppc=4, diag_interval=5, gauss_interval=5, mode_diag=0)
sim = Simulation(cfg, verbose=False)
sim.run()
if sim.rank == 0:
    print(json.dumps({k: sim.history[k] for k in ("kinetic", "electric", "magnetic")}))
