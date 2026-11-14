"""A frozen-in Alfvénic fluctuation in a flow between reflecting walls (PIC-D24), and the helper for
its multi-GPU test (``python -m pypic.launch -n N mgpu_frozen.py``).

Electrons and ions (m_i = 16) drift at -v0 along x in B0 = B y with E0 = v0 B z (PIC-D18). On top:
- B from A_z = -(b/k) cos(k y) h(x): delta B_x = b sin(k y) h(x), div B = 0 bitwise;
- the Alfvénic velocity delta v_x = -(v_A / B) delta B_x on both species (a wave running along +y);
- the current J_z = -d(delta B_x)/dy = -b k cos(k y) h(x) of Ampère's law, carried by the ions
  (their polarization drift in the wave): delta v_z = J_z / (e n);
- the motional field delta E_z = -delta v_x B (where h is flat) from a displacement of the electrons
  along z by xi = -delta E_z / (q n) (q = -1, n = 1).
h is zero within 8 c/omega_pe of the walls, so A_z vanishes on them and nothing is displaced there.
"""
import json
import math
import os

from pypic import Beam, SimConfig, SpeciesConfig

NX, NY, NZ, DX, DT = 384, 256, 2, 0.25, 0.1
LX, LY = NX * DX, NY * DX
MI, B, V0, BFRAC = 16.0, 0.5, 0.05, 0.3
VA = B / math.sqrt(MI)  # n = 1, c = 1
K = 2 * math.pi / LY
H = {"shape": "table", "x": [0.0, 8.0, 16.0, LX - 16.0, LX - 8.0, LX], "values": [0, 0, 1, 1, 0, 0]}
SIN_Y = {"shape": "cos", "mode": 1, "length": LY, "phase": -math.pi / 2}  # sin(k y)
COS_Y = {"shape": "cos", "mode": 1, "length": LY}


def config(dims=None, steps=240, displace=True, a_on_wall=False, xi_x=False):
    b = BFRAC * B
    u0 = V0 / math.sqrt(1 - V0 * V0)
    dvx = -VA / B * b  # amplitude of delta v_x (proper velocity; nonrelativistic here)
    pert = [{"shape": {"y": SIN_Y, "x": H}, "velocity": dvx, "velocity_component": "x"}]
    e_pert = list(pert)
    i_pert = pert + [{"shape": {"y": COS_Y, "x": H}, "velocity": -b * K, "velocity_component": "z"}]
    if displace:
        # delta E_z = -delta v_x B = -(dvx B) sin(k y) h; xi_z = -delta E_z / q with q = -1, n = 1
        e_pert.append({"shape": {"y": SIN_Y, "x": H}, "displacement": -dvx * B, "displacement_component": "z"})
    if xi_x:  # refused: a displacement along the wall normal that reaches the wall
        e_pert.append({"shape": {"y": SIN_Y}, "displacement": 0.01, "displacement_component": "x"})
    hx = {"shape": "uniform"} if a_on_wall else H
    bc = (("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))
    return SimConfig(
        n_cells=(NX, NY, NZ), dx=(DX,) * 3, dt=DT, n_steps=steps, boundaries=bc, ppc=8, lattice=(2, 2, 2),
        B0=(0.0, B, 0.0), E0=(0.0, 0.0, V0 * B), dims=dims, diag_interval=20, gauss_interval=20,
        initial_fields=[{"kind": "potential", "component": "Az", "shape": {"y": COS_Y, "x": hx}, "scale": -b / K}],
        species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (-u0, 0, 0), 1e-4, perturbation=e_pert)]),
                 SpeciesConfig("i", 1, MI, [Beam(1.0, (-u0, 0, 0), 1e-4 / MI, perturbation=i_pert)])])


if __name__ == "__main__":
    from pypic.simulation import Simulation

    n = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    sim = Simulation(config(dims=(2, 1, 1) if n == 2 else None, steps=60), verbose=False)
    sim.run()
    if sim.rank == 0:
        print(json.dumps({k: sim.history[k] for k in ("kinetic", "electric", "magnetic", "gauss", "divB")}))
