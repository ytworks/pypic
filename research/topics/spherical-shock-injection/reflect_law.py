"""Reflection law at the front of a dipole run, and returns to the other hemisphere (D-A) on the tracks.

python reflect_law.py <run>   (after analyze.py front, sites, campaign); writes analysis/reflect_law.json
"""
import sys

sys.path.insert(0, "/home/y-takagi/pypic/research/topics/spherical-shock-injection")
import json

import hdf5plugin  # noqa: F401
import numpy as np

from pypic.analysis import open_run

R = sys.argv[1]
run = open_run(R)
c = np.array(run.meta["planner"]["center"], float)
E = np.load(R + "/analysis/campaign_entries.npy")
col, tin, pol, tbn, vht, v, vpar, ret, gin, gout, ins = E[:, :11].T
comp = E[:, 13]
ref2 = (ret == 1) & (gout > 1.02 * gin)
# persistent: after the return the electron stays upstream for one gyro-period or gets d_i + 2 rho_e ahead, with
# the gyro-period and rho_e from the static dipole |B| at the front radius of that polar angle and time
fz = np.load(R + "/analysis/front.npz")
moment = float(run.config["external_fields"]["dipoles"][0]["moment"][2])
d_i = float(run.meta["planner"]["d_i"])
i_f = np.array([int(np.argmin(np.abs(fz["t"] - q))) for q in tin])
r_f = np.array([np.interp(np.radians(pp), fz["theta"], fz["r"][k]) for pp, k in zip(pol, i_f)])
b_loc = moment / r_f ** 3 * np.sqrt(1 + 3 * np.cos(np.radians(pol)) ** 2)
t_up, ahead = E[:, 14], E[:, 15]
persist = (t_up >= 2 * np.pi * gin / b_loc) | (ahead >= d_i + 2 * np.sqrt(gin ** 2 - 1) / b_loc)
ref = ref2 & persist
out = {"entries": int(len(E)), "reflections_2pc": int(ref2.sum()), "reflections_persistent": int(ref.sum()),
       "definition": "persistent: gamma_out > 1.02 gamma_in, and one local gyro-period upstream or d_i + 2 rho_e ahead"}
bins = {"V_HT [c]": (vht, [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 3.0]),
        "local theta_Bn [deg]": (tbn, [0, 15, 30, 45, 60, 75, 90]),
        "polar angle [deg]": (pol, [0, 15, 30, 45, 60, 75, 90]),
        "V_HT / v": (vht / v, [0, 1, 1.5, 2, 3, 5, 10, 100]),
        "compression B / B_static": (comp, [0, 1.5, 2, 3, 4, 6, 20])}
for name, (x, edges) in bins.items():
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x < b) & (ins == 0)
        if m.sum():
            g = (gout[m & ref] / gin[m & ref]) if (m & ref).any() else np.array([np.nan])
            rows.append({"bin": [a, b], "entries": int(m.sum()), "P_ref": float(ref[m].mean()),
                         "gain_ratio_p50_p90": [float(np.nanpercentile(g, 50)), float(np.nanpercentile(g, 90))]})
    out[name] = rows
# D-A: hemisphere of each entry from the tracks; an electron reflected in one hemisphere that later enters the layer
# in the other hemisphere (after travelling upstream)
tr = run.tracks("tracks_e", "e")
t = np.asarray(tr.t, float)
cols = np.unique(col[ref]).astype(int)
zz = np.asarray(tr.f["z"][:, cols], float) - c[2]
pos = {cc: j for j, cc in enumerate(cols)}
hemi_changes, n_el, gains_second = 0, 0, []
by_el = {}
for k in np.flatnonzero(np.isin(col, cols)):
    by_el.setdefault(int(col[k]), []).append(k)
for cc, ks in by_el.items():
    ks = sorted(ks, key=lambda q: tin[q])
    j = pos[cc]
    hem = [np.sign(zz[int(np.argmin(np.abs(t - tin[q]))), j]) for q in ks]
    first_ref = next((i for i, q in enumerate(ks) if ref[q]), None)
    if first_ref is None:
        continue
    n_el += 1
    for i in range(first_ref + 1, len(ks)):
        if hem[i] != hem[first_ref] and hem[i] != 0:
            hemi_changes += 1
            if ref[ks[i]]:
                gains_second.append(gout[ks[i]] / gin[ks[i]])
            break
out["D-A: reflected electrons (tracked)"] = n_el
out["D-A: of these, a later entry in the other hemisphere"] = hemi_changes
med2 = float(np.median(gains_second)) if gains_second else None
out["D-A: of these, reflected again there (gain ratio p50)"] = [len(gains_second), med2]
json.dump(out, open(R + "/analysis/reflect_law.json", "w"), indent=1)
print(json.dumps(out, indent=1))
