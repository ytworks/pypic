"""Analysis of the topic's runs; writes the run records of ``runs/``.

    python research/topics/shock-alfvenic-turbulence/analyze.py pilot <run dir> [--reference <laminar run dir>]
                                                                      [--block 1.0] [--record]

``pilot``: a periodic box of the pilot P1a (design/setup.md, Sect. 8). Measured on the stored plane z = 0,
frame by frame, after an average over blocks of ``--block`` d_i (it removes the particle noise of a cell):

- thermal energy per particle of each species, eth = (ekin - m |flux|^2 / (2 n)) / n in a block, averaged
  over the blocks with the weight n, and its change relative to the first frame (``dT_e``, ``dT_i``). The
  first frame is the reference, not 3 T / 2: the loaded Maxwell-Jüttner distribution holds more.
- density: std / mean over the blocks (``dn_e``, ``dn_i``);
- |B|^2: std / mean over the blocks (``b2_spread``);
- dB_x: the rms over y of the x average of B_x, over B0 (``dbx_rms``: the part with k_x = 0), the share of
  its power in the loaded mode numbers (``loaded_share``), and the share of the block-averaged B_x^2 that is
  not in k_x = 0 (``oblique_share``);
- field lines: peak-to-peak of xi(y) = int <B_x>_x dy / B0 (``xi_pp``, d_i) and Xi = xi_pp / L_sh;
- from the scalars: the change of the energy budget over the first total energy (``budget_drift``);
- from the histogram (x, ekin): the share of electrons above 10 and 30 T_e0 in the box frame
  (``share_10``, ``share_30``), T_e0 = theta_e m_e c^2.

The criteria C1 to C3 and the comparison with the linear table are in design/setup.md (Sect. 8): they were
fixed before the runs. ``--record`` writes runs/<date>-<commit>-<config>.json and its entry in runs/index.json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
TOPIC = os.path.basename(HERE)
sys.path.insert(0, HERE)

from plan import b2_spread  # noqa: E402

from pypic.analysis import open_run  # noqa: E402

# criteria of the pilot, fixed before the runs (design/setup.md, Sect. 8)
C1_DT_E = 0.05  # laminar: |dT_e / T_e| over 40 / Omega_ci
C2_DN = 0.05  # one wave: std / mean of the density over the blocks, at every frame
C2_B2 = 0.05  # one wave: std / mean of |B|^2 over the blocks, at every frame
C3_DBX = 0.9  # dbx_rms at the end over the loaded value
LINEAR_FACTOR = 2.0  # several modes: measured |B|^2 spread within this factor of the linear one
T_FULL = 40.0  # 1 / Omega_ci: the length the criteria refer to


def blocks(a: np.ndarray, n: int) -> np.ndarray:
    """Mean of a over blocks of n x n cells (the rest of a division is cut off)."""
    nx, ny = a.shape[0] // n, a.shape[1] // n
    return a[:nx * n, :ny * n].reshape(nx, n, ny, n).mean(axis=(1, 3), dtype=np.float64)


def share_above(counts: np.ndarray, edges: np.ndarray, e_min: float) -> tuple:
    """Share of the counts in the bins whose lower edge is at or above e_min, and the edge used."""
    j = int(np.searchsorted(edges[:-1], e_min))
    total = counts.sum()
    return (float(counts[j:].sum() / total) if total > 0 else 0.0), float(edges[min(j, edges.size - 1)])


def pilot(run_dir: str, block: float = 1.0) -> dict:
    r = open_run(run_dir)
    pl = r.meta["planner"]
    par = pl["params"]
    if par["wall"]:
        raise ValueError("pilot: a periodic box is expected (wall: false)")
    d_i, b0, w_ci = pl["d_i"], pl["b0"], pl["omega_ci"]
    nb = int(round(block * d_i / r.dx[0]))
    masses = {"e": 1.0, "i": par["mass_ratio"]}
    plane = r.series("plane")
    t = np.asarray(plane.t) * w_ci
    numbers = [int(round(m["k_di"] * par["ly"] / (2.0 * math.pi))) for m in pl["modes"]] if par["db"] > 0.0 else []
    out = {k: [] for k in ("eth_e", "eth_i", "dn_e", "dn_i", "b2_spread", "dbx_rms", "loaded_share",
                           "oblique_share", "xi_pp")}
    for f in range(t.size):
        for s, m in masses.items():
            n = blocks(plane.read(f"n:{s}", f), nb)
            flux2 = sum(blocks(plane.read(f"flux{c}:{s}", f), nb) ** 2 for c in "xyz")
            eth = (blocks(plane.read(f"ekin:{s}", f), nb) - 0.5 * m * flux2 / n) / n
            out[f"eth_{s}"].append(float((eth * n).sum() / n.sum()))
            out[f"dn_{s}"].append(float(n.std() / n.mean()))
        bx = plane.read("Bx", f).astype(np.float64)
        bxb = blocks(bx, nb)
        b2 = bxb ** 2 + blocks(plane.read("By", f), nb) ** 2 + blocks(plane.read("Bz", f), nb) ** 2
        out["b2_spread"].append(float(b2.std() / b2.mean()))
        line = bx.mean(axis=0)  # the part with k_x = 0, as a function of y
        power = np.abs(np.fft.rfft(line)) ** 2
        power[0] = 0.0
        out["dbx_rms"].append(float(math.sqrt((line ** 2).mean()) / b0))
        out["loaded_share"].append(float(power[numbers].sum() / power.sum()) if numbers else 0.0)
        total = float((bxb ** 2).mean())
        out["oblique_share"].append(1.0 - float((bxb.mean(axis=0) ** 2).mean()) / total if total > 0.0 else 0.0)
        xi = np.cumsum(line - line.mean()) * r.dx[1] / b0 / d_i
        out["xi_pp"].append(float(xi.max() - xi.min()))
    res = {"t": t.tolist(), **out}
    for s in masses:
        e = np.asarray(out[f"eth_{s}"])
        res[f"dT_{s}"] = (e / e[0] - 1.0).tolist()
    res["Xi"] = (np.asarray(out["xi_pp"]) / pl["l_sh"]).tolist()
    if par["db"] > 0.0:  # |B|^2 spread of the loaded modes with their linear phases, at the frame times
        kdi = np.array([m["k_di"] for m in pl["modes"]])
        amp = np.array([m["dbx"] for m in pl["modes"]])
        ph = np.array([m["phase"] for m in pl["modes"]])
        om = np.array([m["omega_over_omega_ci"] for m in pl["modes"]])
        sg = np.array([1.0 if m["branch"] == "R" else -1.0 for m in pl["modes"]])
        res["b2_spread_linear"] = [b2_spread(kdi, amp, ph - om * ti, sg) for ti in t]
    sc = r.scalars
    ts = np.asarray(sc.t) * w_ci
    budget, total = np.asarray(sc["budget"]), np.asarray(sc["total"])
    res["scalars"] = {"t_end": float(ts[-1]), "budget_drift": float((budget[-1] - budget[0]) / total[0]),
                      "kinetic_e_change": float(sc["kinetic_e"][-1] / sc["kinetic_e"][0] - 1.0),
                      "kinetic_i_change": float(sc["kinetic_i"][-1] / sc["kinetic_i"][0] - 1.0),
                      "magnetic_change": float(sc["magnetic"][-1] / sc["magnetic"][0] - 1.0),
                      "gauss_max": float(np.nanmax(np.abs(sc["gauss"]))),
                      "divB_max": float(np.nanmax(np.abs(sc["divB"])))}
    sp = r.series("spectrum")
    edges = np.asarray(sp.edges[1])
    counts = np.asarray(sp.read("e"), dtype=np.float64).sum(axis=1)  # (times, energy bins)
    res["spectrum"] = {"t": (np.asarray(sp.t) * w_ci).tolist()}
    for mult in (10, 30):
        shares = [share_above(c, edges, mult * pl["theta_e"]) for c in counts]
        res["spectrum"][f"share_{mult}"] = [s[0] for s in shares]
        res["spectrum"][f"edge_{mult}_over_T_e0"] = shares[0][1] / pl["theta_e"]
    res["block_cells"] = nb
    res["block_particles"] = nb * nb * r.config["ppc"]
    return {"run": r, "result": res}


def at(t, y, t0: float) -> float:
    return float(np.interp(t0, t, y))


def criteria(par: dict, res: dict, ref: dict | None) -> dict:
    """Verdicts of the pilot criteria that apply to this run (design/setup.md, Sect. 8)."""
    t, c = np.asarray(res["t"]), {}
    full = t[-1] >= T_FULL - 0.5
    dte = np.asarray(res["dT_e"])
    rate = float(np.polyfit(t, dte, 1)[0])
    c["dT_e_end"], c["dT_i_end"], c["dT_e_rate"] = float(dte[-1]), float(res["dT_i"][-1]), rate
    c["dT_e_first_half"] = at(t, dte, 0.5 * t[-1])
    if par["db"] == 0.0:
        c["C1"] = {"limit": C1_DT_E, "value": float(dte[-1]) if full else rate * T_FULL,
                   "extrapolated": not full, "passed": bool(abs(dte[-1]) <= C1_DT_E) if full else None}
    else:
        if par["state"] == "constant_b" and par["zeros"] == 0:
            dn, b2 = max(res["dn_e"]), max(res["b2_spread"])
            c["C2"] = {"limit_dn": C2_DN, "limit_b2": C2_B2, "dn_e_max": dn, "b2_spread_max": b2,
                       "passed": bool(dn <= C2_DN and b2 <= C2_B2) if full else None}
        c["C3"] = {"limit": C3_DBX * par["db"], "value": res["dbx_rms"][-1],
                   "passed": bool(res["dbx_rms"][-1] >= C3_DBX * par["db"]) if full else None}
        if len(res.get("b2_spread_linear", [])) and par["state"] == "constant_b" and par["zeros"] > 0:
            rows = []
            for t0 in (10.0, 20.0, 40.0):
                if t0 <= t[-1] + 0.5:
                    m, lin = at(t, res["b2_spread"], t0), at(t, res["b2_spread_linear"], t0)
                    rows.append({"t": t0, "measured": m, "linear": lin, "ratio": m / lin})
            c["linear"] = {"factor": LINEAR_FACTOR, "rows": rows,
                           "passed": bool(all(1.0 / LINEAR_FACTOR <= x["ratio"] <= LINEAR_FACTOR for x in rows))}
    if ref is not None:
        tr = np.asarray(ref["t"])
        tc = min(t[-1], tr[-1])
        c["relative_to_laminar"] = {
            "t": float(tc),
            "dT_e": at(t, res["dT_e"], tc) - at(tr, ref["dT_e"], tc),
            "dT_i": at(t, res["dT_i"], tc) - at(tr, ref["dT_i"], tc),
            "dn_e_floor": float(np.mean(ref["dn_e"])), "b2_spread_floor": float(np.mean(ref["b2_spread"])),
            "share_10_ratio": res["spectrum"]["share_10"][-1] / ref["spectrum"]["share_10"][-1]
            if ref["spectrum"]["share_10"][-1] > 0 and t[-1] == tr[-1] else None}
    return c


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()


def dir_bytes(path: str) -> int:
    return sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(path) for f in fs)


def write_record(run_dir: str, r, res: dict, crit: dict, reference: str | None) -> str:
    runs = os.path.join(HERE, "runs")
    os.makedirs(runs, exist_ok=True)
    info, pl, meta = r.info, r.meta["planner"], r.meta["run"]
    commit = info["git_commit"]
    dirty = git("status", "--porcelain", "--", ".", f":!research/topics/{TOPIC}/runs")
    config = os.path.splitext(os.path.basename(r.meta["config_file"]))[0]
    date = str(info["created"])[:10]
    name = f"{date}-{commit.removesuffix('-dirty')[:9]}-{config}.json"
    record = {
        "topic": TOPIC, "kind": "pilot", "config": config, "date": date, "commit": commit,
        "analysis_commit": git("rev-parse", "HEAD") + ("-dirty" if dirty else ""),
        "config_file": os.path.relpath(r.meta["config_file"], REPO), "params": pl["params"],
        "seed": pl["params"]["seed"], "gpus": meta["gpus"], "gpu_name": info.get("gpu_name"),
        "steps": meta["steps"], "wall_time_s": meta["wall_time_s"],
        "ns_per_particle_step": meta["ns_per_particle_step"],
        "output": {"path": os.path.abspath(run_dir), "bytes": dir_bytes(run_dir)}, "reference": reference,
        "planner": {k: pl[k] for k in ("v_a", "v0", "omega_ci", "theta_e", "theta_i", "modes", "xi_pp", "Xi", "l_sh",
                                       "b2_spread", "b_total", "debye_cells", "particles", "dims")},
        "criteria": crit, "analysis": res,
    }
    with open(os.path.join(runs, name), "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1, allow_nan=False)
        fh.write("\n")
    index = os.path.join(runs, "index.json")
    entries = []
    if os.path.exists(index):
        with open(index, encoding="utf-8") as fh:
            entries = json.load(fh)["entries"]
    entries = [e for e in entries if e["file"] != name]
    entries.append({"file": name, "kind": "pilot", "config": config, "date": date, "commit": commit,
                    "analysis_commit": record["analysis_commit"]})
    with open(index, "w", encoding="utf-8") as fh:
        json.dump({"entries": sorted(entries, key=lambda e: e["file"])}, fh, indent=1)
        fh.write("\n")
    return os.path.join(runs, name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["pilot"])
    ap.add_argument("run")
    ap.add_argument("--reference", help="laminar run of the same resolution (noise floors, heating)")
    ap.add_argument("--block", type=float, default=1.0, help="block size of the averages (d_i)")
    ap.add_argument("--record", action="store_true", help="write the record into runs/")
    a = ap.parse_args()
    got = pilot(a.run, a.block)
    ref = pilot(a.reference, a.block)["result"] if a.reference else None
    r, res = got["run"], got["result"]
    crit = criteria(r.meta["planner"]["params"], res, ref)
    print(json.dumps({"config": os.path.basename(r.meta["config_file"]), "commit": r.info["git_commit"],
                      "t_end": res["t"][-1], "criteria": crit, "scalars": res["scalars"],
                      "dn_e_end": res["dn_e"][-1], "b2_spread_end": res["b2_spread"][-1],
                      "dbx_rms": [res["dbx_rms"][0], res["dbx_rms"][-1]],
                      "loaded_share_end": res["loaded_share"][-1], "oblique_share_end": res["oblique_share"][-1],
                      "Xi": [res["Xi"][0], res["Xi"][-1]],
                      "share_10": [res["spectrum"]["share_10"][0], res["spectrum"]["share_10"][-1]]}, indent=1))
    if a.record:
        print("record:", write_record(a.run, r, res, crit, os.path.abspath(a.reference) if a.reference else None))


if __name__ == "__main__":
    main()
