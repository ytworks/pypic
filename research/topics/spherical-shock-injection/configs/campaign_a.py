"""Writes the configs of campaign A (campaign-a.md): electron temperature and Alfven Mach number varied around P5.

python campaign_a.py      (writes A5-*.yaml and A7-*.yaml next to this file)

A5: the box, grid, field, debris and length of P5 (strong-P5.yaml), v_te 0.14 / 0.20 / 0.28 (P5: 0.10).
A7: M_A 7.5 (V_d of P5, B0 lower: wpe_wce 5.25), R_m 45 d_i = 6 rho_d, cells of 0.75 d_e (allowed by the longer
    Debye length of the hotter electrons), box 576 d_e, shell width 0.2 R_m (omega_pe dt at the shell peak),
    1600 steps; v_te 0.14 / 0.20 / 0.28.
T_i = 0.01 m_e c^2 (that of P5) in every run: ti_te = 0.01 / v_te^2. Seeds 12345 and 54321. Null runs: the same
plasma, field and grid without the debris (v_te 0.10 at A5 is L-null).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = {"": 12345, "-s2": 54321}
TEMPS = {"T14": 0.14, "T20": 0.20, "T28": 0.28}

GEOM = {
    # name: (mach_alfven, wpe_wce, rm, shell_width, cells_per_de, half_width [R_B], run [R_B / V_d], box [d_e])
    "A5": (5.0, 3.5, 35.0, 0.1, 2.0, 0.46902, 0.31408, 384.0),
    "A7": (7.5, 5.25, 45.0, 0.2, 1.3333333, 0.41758, 0.31070, 576.0),
}


def outputs(box, debris):
    c = box / 2.0
    sp = "[e, i, de, di]" if debris else "[e, i]"
    dens = '"n:e", "n:i", "n:de", "n:di", ' if debris else '"n:e", "n:i", '
    coarse = ('["n:e", "flux:e", "ekin:e", "n:i", "flux:i", "ekin:i", "n:de", "flux:de", "ekin:de",\n'
              '               "n:di", "flux:di", "ekin:di", E, B, rho]') if debris else \
        '["n:e", "flux:e", "ekin:e", "n:i", "flux:i", "ekin:i", E, B, rho]'
    nb = int(round(box / 2.0))
    hot_sp = "[e, de]" if debris else "[e]"
    stride = 947 if debris else 4735
    sl = ""
    for name, axis in (("slice_z", "z"), ("slice_y", "y"), ("slice_x", "x")):
        sl += (f"    - name: {name}\n      kind: slice\n      fields: [{dens}\"flux:i\", \"flux:e\", \"ekin:e\", "
               f"\"ekin:i\", E, B, J]\n      axis: {axis}\n      position: {c}\n      coarsen: [1, 1, 1]\n"
               f"      every: 36\n")
    return f"""output:
  # As strong-P5.yaml, for a box of {box:g} d_e (centre {c:g}); hot_e from gamma 1.5 every 36 steps (at the
  # times of the slices, as `analyze.py census` and `counts` need).
  diagnostics:
    - name: coarse
      kind: fields
      fields: {coarse}
      coarsen: [13, 13, 13]
      every: 144
{sl}    - name: spec_meridional
      kind: histogram
      species: {sp}
      quantities: [x, z, ekin]
      bins: [{nb}, {nb}, 96]
      ranges: [[0.0, {box}], [0.0, {box}], [0.0001, 100.0]]
      log: [false, false, true]
      region: [null, [{c - 8.0}, {c + 8.0}], null]
      every: 72
    - name: spec_equatorial
      kind: histogram
      species: {sp}
      quantities: [x, y, ekin]
      bins: [{nb}, {nb}, 96]
      ranges: [[0.0, {box}], [0.0, {box}], [0.0001, 100.0]]
      log: [false, false, true]
      region: [null, null, [{c - 8.0}, {c + 8.0}]]
      every: 72
    - name: spec_all
      kind: histogram
      species: {sp}
      quantities: [ekin]
      bins: [224]
      ranges: [[0.000001, 100.0]]
      log: [true]
      every: 36
    - name: spec_up
      kind: histogram
      species: [e, i]
      quantities: [ekin]
      bins: [224]
      ranges: [[0.000001, 10.0]]
      log: [true]
      region: [[16.0, 48.0], [16.0, 48.0], [16.0, {box - 16.0}]]
      every: 36
    - name: hot_e
      kind: particles
      species: {hot_sp}
      select: {{gamma_min: 1.5}}
      every: 36
    - name: tracks_e
      kind: tracks
      species: [e]
      stride: {stride}
      with_fields: true
      every: 2
"""


def config(geom, vte, seed, debris, name):
    ma, wpe_wce, rm, sw, cpd, hw, run, box = GEOM[geom]
    head = (f"# Campaign A run {name} (topic spherical-shock-injection; design: ../campaign-a.md).\n"
            f"# Written by campaign_a.py; run with ../run_piston.py.\n")
    body = f"""mass_ratio: 16.0
mach_alfven: {ma}
wpe_wce: {wpe_wce}
vte: {vte}
ti_te: {0.01 / vte ** 2:.6f}
rm: {rm}
shell_radius: 0.33
shell_width: {sw}
debris_beams: 2
body: false
"""
    if not debris:
        body += "debris: false\n"
    body += f"""half_width: {hw}
cells_per_de: {cpd}
ppc: 2
run: {run}
filter_passes: 4
seed: {seed}
overrides:
  precision: compact
"""
    return head + body + outputs(box, debris)


def main():
    names = []
    for geom in GEOM:
        for tname, vte in TEMPS.items():
            for suffix, seed in SEEDS.items():
                names.append((f"{geom}-{tname}{suffix}", geom, vte, seed, True))
            names.append((f"{geom}-{tname}-null", geom, vte, 12345, False))
    names.append(("A5-T10-s2", "A5", 0.1, 54321, True))
    for name, geom, vte, seed, debris in names:
        with open(os.path.join(HERE, name + ".yaml"), "w") as f:
            f.write(config(geom, vte, seed, debris, name))
    print("\n".join(n[0] for n in names))


if __name__ == "__main__":
    main()
