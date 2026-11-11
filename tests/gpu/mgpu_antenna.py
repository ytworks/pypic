"""Helper for the laser antenna test (PIC-D20): run with ``python -m pypic.launch -n N``.

    mgpu_antenna.py <dims> <polarization>    dims like 2,2,1; polarization y | z

A focused pulse with a Gaussian envelope in a vacuum box with six open faces. The two node planes
of the antenna lie on either side of the face between the ranks along x, and the axis of the beam
on the faces between the ranks across it. Prints the histories.
"""
import json
import math
import sys

from pypic import SimConfig
from pypic.config import Laser
from pypic.simulation import Simulation

dims = tuple(int(v) for v in sys.argv[1].split(","))
pol = sys.argv[2]
lam = 2.0 * math.pi  # omega = 1
n, dx = (48, 32, 32), (lam / 12, lam / 6, lam / 6)
laser = Laser(plane_x=(n[0] // 2 - 0.5) * dx[0], a0=1.0, omega=1.0, polarization=pol, waist=lam,
              focus_x=32 * dx[0], phase=0.4, envelope={"shape": "gaussian", "duration": 2 * lam, "t_peak": 5 * lam})
cfg = SimConfig(n_cells=n, dx=dx, dt=lam / 16, n_steps=192, dims=dims, npml=6, species=[],
                boundaries=(("open", "open"),) * 3, lasers=[laser], diag_interval=8, gauss_interval=8)
sim = Simulation(cfg, verbose=False)
sim.run()
if sim.rank == 0:
    h = sim.history
    print(json.dumps({k: h[k] for k in ("electric", "magnetic", "antenna", "boundary_flux", "gauss", "divB")}))
