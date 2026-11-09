"""Helper for the multi-GPU tests: run with ``python -m pypic.launch -n N``.

    mgpu_case.py <solver>                 history of a small mixed-boundary run (JSON on stdout)
    mgpu_case.py <solver> <output_dir>    the same run writing every output diagnostic
"""
import json
import sys

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import OutputConfig
from pypic.simulation import Simulation

solver = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else None
diags = []
if out:
    diags = [
        dict(name="fields", kind="fields", fields=["E", "B", "J", "rho", "n:e", "flux:e", "ekin:i"], every=10,
             coarsen=[3, 2, 2]),
        dict(name="slice", kind="slice", fields=["Ex", "Bz", "n:e"], axis="x", position=3.3, every=5),
        dict(name="lineout", kind="lineout", fields=["Ex", "Ey"], axis="x", every=5, average=True),
        dict(name="modes", kind="modes", fields=["Ex"], modes=[[1, 0, 0], [2, 0, 0], [1, 1, 0], [0, 1, 1]], every=5,
             average=True),
        dict(name="kspec", kind="kspectrum", fields=["Ex"], every=10),
        dict(name="phase", kind="histogram", quantities=["x", "ux"], bins=[32, 16], every=10),
        dict(name="spectrum", kind="histogram", quantities=["ekin"], bins=[40], log=[True], ranges=[[1e-6, 1]],
             every=10),
        dict(name="sample", kind="particles", every=20),  # stride 1: every particle, sorted by id
        dict(name="tracks", kind="tracks", species=["e"], stride=50, every=2, with_fields=True),
    ]
cfg = SimConfig(n_cells=(32, 16, 16), dx=(0.2,) * 3, dt=0.08 if solver == "explicit" else 0.3, n_steps=20,
                solver=solver, boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("open", "open")),
                npml=6, species=[SpeciesConfig("e", -1, 1, [Beam(0.5, (0.3, 0.05, 0), 0.01),
                                                            Beam(0.5, (-0.3, 0, 0.05), 0.01)]),
                                 SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                ppc=4, diag_interval=5, gauss_interval=5, output=OutputConfig(dir=out, diagnostics=diags))
sim = Simulation(cfg, verbose=False)
sim.run()
sim.close()
if sim.rank == 0:
    print(json.dumps({k: sim.history[k] for k in ("kinetic", "electric", "magnetic")}))
