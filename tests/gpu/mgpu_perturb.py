"""Helper for multi-GPU tests (``python -m pypic.launch -n N mgpu_perturb.py <solver> [filter]``):
a displacement of several ranks' widths (moved in rounds with migration) and, with ``filter``,
reflecting walls and a 4-pass compensated current filter must give the same run on 1 and 2 GPUs."""
import json
import os
import sys

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import CurrentFilter
from pypic.simulation import Simulation

solver = sys.argv[1]
filt = len(sys.argv) > 2 and sys.argv[2] == "filter"
pert = [{"axis": "x", "mode": 1, "displacement": 2.1}, {"axis": "x", "mode": 3, "density": 0.2, "phase": 0.4}]
bc_y = ("reflecting", "reflecting") if filt else ("periodic", "periodic")
cfg = SimConfig(n_cells=(32, 8, 8), dx=(0.1,) * 3, dt=0.05 if solver == "explicit" else 0.2, n_steps=10,
                solver=solver, species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.001, perturbation=pert)]),
                                        SpeciesConfig("i", 1, 100, [Beam(1.0, (0, 0, 0), 1e-5)])],
                ppc=8, diag_interval=2, gauss_interval=2, B0=(0.0, 0.0, 0.3),
                boundaries=(("periodic", "periodic"), bc_y, ("periodic", "periodic")),
                current_filter=CurrentFilter(4, True) if filt else CurrentFilter(),
                # with the filter, decompose the wall axis (y) to exercise walls and rank faces together
                dims=(1, 2, 1) if filt and int(os.environ.get("PYPIC_WORLD_SIZE", "1")) == 2 else None)
sim = Simulation(cfg, verbose=False)
sim.run()
if sim.rank == 0:
    print(json.dumps({k: sim.history[k] for k in ("kinetic", "electric", "magnetic", "gauss", "budget")}))
