"""Writes the configs of the dipole runs (dipole-design.md). python dipole_configs.py

D1: a debris shell expands from a magnetised planet through its dipole. V_d = 0.25 c; the equatorial dipole field
gives M_A = 3 at r = 70 d_e (B_eq = 0.333 there); planet radius 32 d_e; shell at 48 d_e, width 5 d_e; R_m = 90 d_e;
B0 negligible (dipole radius 576 d_e = 3 half-widths). D1-null: the same without the debris.
"""
import math
import os

from campaign_a import outputs

HERE = os.path.dirname(os.path.abspath(__file__))
MI, V, R_REF, MA_REF = 16.0, 0.25, 70.0, 3.0
R_DIP, R_BODY, CORE, R_M = 576.0, 32.0, 16.0, 90.0
SHELL, WIDTH, HALF, T_RUN, BOX = 48.0, 5.0, 192.0, 450.0, 384.0


def config(name, debris):
    b_eq = V / MA_REF * math.sqrt(MI)
    moment = b_eq * R_REF ** 3
    b0 = moment / R_DIP ** 3
    ma0 = V / (b0 / math.sqrt(MI))
    r_b = R_M * ma0 ** (2.0 / 3.0)
    head = (f"# Dipole run {name} (topic spherical-shock-injection; design: ../dipole-design.md).\n"
            f"# Written by dipole_configs.py; run with ../run_dipole.py. B0 is negligible: mach_alfven, wpe_wce,\n"
            f"# half_width and run refer to it; the local field is the dipole (equatorial M_A = 3 at r = 70 d_e).\n")
    body = f"""mass_ratio: {MI}
mach_alfven: {ma0:.4f}
wpe_wce: {1.0 / b0:.4f}
vte: 0.1
ti_te: 1.0
rm: {R_M / math.sqrt(MI)}
shell_radius: {SHELL / R_M:.6f}
shell_width: {WIDTH / R_M:.6f}
debris_beams: 2
body: true
dipole_radius: {R_DIP / R_M:.6f}
body_radius: {R_BODY / R_DIP:.6f}
core: {CORE / R_DIP:.6f}
dipole: [0.0, 0.0, 1.0]
"""
    if not debris:
        body += "debris: false\n"
    body += f"""half_width: {HALF / r_b:.6f}
cells_per_de: 2.0
ppc: 2
run: {T_RUN * V / r_b:.6f}
filter_passes: 4
seed: 12345
overrides:
  precision: compact
"""
    return head + body + outputs(BOX, debris)


DISK = """    - name: disk_e             # every electron in an equatorial slab (|z - c| < 20, |x - c|, |y - c| < 130 d_e):
      # the electrons that collect in the equatorial disk of D1, followed by id from dump to dump (D1-disk)
      kind: particles
      species: [e]
      region: [[62.0, 322.0], [62.0, 322.0], [172.0, 212.0]]
      every: 36
"""


def config_disk():
    """D1-disk: D1 with the same physics and seed; the tracks and hot_e are replaced by the slab dump."""
    text = config("D1-disk", True)
    head, rest = text.split("    - name: hot_e", 1)
    return head.replace("Dipole run D1-disk", "Dipole run D1-disk (D1 rerun with a dump of the equatorial slab)") + DISK


if __name__ == "__main__":
    for name, debris in (("D1", True), ("D1-null", False)):
        with open(os.path.join(HERE, name + ".yaml"), "w") as f:
            f.write(config(name, debris))
        print(name)
    with open(os.path.join(HERE, "D1-disk.yaml"), "w") as f:
        f.write(config_disk())
    print("D1-disk")
