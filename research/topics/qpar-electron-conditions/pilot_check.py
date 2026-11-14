"""Numbers of one pilot run (plan step 4): conservation, timing, temperatures, shock position, seeds.

    python research/topics/qpar-electron-conditions/pilot_check.py <run dir> [--json <file>]

CPU only (reads the output with ``pypic.analysis``). It prints numbers and judges nothing.

- Temperature of a species from the moments of the ``profile`` lineout, per cell column:
  T = (2/3) (ekin / n - m (gamma_b - 1)), v_b = flux / n. It is the non-relativistic estimate (the error is of
  relative order v^2); the first output gives its value for the loaded plasma, so read changes against that.
- Periodic box (N1) and flow-through control: the temperatures averaged over x against time.
- Flow-through control (``through``): B_y / B0y - 1, E_z / E0z - 1 and B_z / B0y, averaged over the box without
  ``margin`` d_i at each end, against time; the slope of the first two per 1 / Omega_ci with its standard error
  (straight-line fit that takes the outputs as independent) over t >= 1 / Omega_ci and over the second half.
- Wall: the shock position is the largest x where the ion density, averaged over ``smooth`` d_i, is above 2;
  its speed is the slope over the second half of the run. The upstream numbers are taken between the shock
  + ``margin`` d_i and the face - ``margin`` d_i.
- Tail: the share of a species' histogram above |u| = u0 + k sqrt(theta), k = 5 and 8, per output.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from pypic.analysis.run import open_run


def temperature(prof, name: str, mass: float):
    """(n, T) of species ``name``, arrays (n_t, n_x); T in m_e c^2, NaN where the density is zero."""
    n = prof[f"n:{name}"]
    with np.errstate(divide="ignore", invalid="ignore"):
        v2 = sum((prof[f"flux{c}:{name}"] / n) ** 2 for c in "xyz")
        gb = 1.0 / np.sqrt(1.0 - np.minimum(v2, 0.999999))
        t = (2.0 / 3.0) * (prof[f"ekin:{name}"] / n - mass * (gb - 1.0))
    return n, np.where(n > 0, t, np.nan)


def shock_position(x, n_i, cells: int, level: float = 2.0):
    """Largest x where the running mean of n_i over ``cells`` is above ``level`` (0 if nowhere)."""
    out = np.zeros(len(n_i))
    kern = np.ones(cells) / cells
    for j, row in enumerate(n_i):
        idx = np.nonzero(np.convolve(row, kern, mode="same") > level)[0]
        out[j] = x[idx[-1]] if len(idx) else 0.0
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--json")
    ap.add_argument("--smooth", type=float, default=0.5, help="d_i")
    ap.add_argument("--margin", type=float, default=3.0, help="d_i")
    a = ap.parse_args()
    run = open_run(a.run)
    info, par = run.meta["planner"], run.meta["params"]
    d_i, w_ci = info["d_i"], info["omega_ci"]
    theta = {"e": info["theta_e"], "i": info["theta_i"]} | {s["name"]: s["theta"] for s in par.get("seeds", [])}
    mass = {n: (par.get("mass_ratio", 100.0) if n == "i" else 1.0) for n in theta}
    out = {"run": a.run, "commit": run.info.get("git_commit"), "stats": run.meta.get("run")}

    sc = run.scalars
    out["scalars"] = {
        "gauss_max": float(np.nanmax(np.abs(sc["gauss"]))), "divB_max": float(np.nanmax(np.abs(sc["divB"]))),
        "n_particles": [float(sc["n_particles"][0]), float(sc["n_particles"][-1])],
        "n_injected": float(sc["n_injected"][-1]), "fill_max": float(np.nanmax(sc["fill"])),
        "kinetic_e": [float(sc["kinetic_e"][0]), float(sc["kinetic_e"][-1])],
        "kinetic_i": [float(sc["kinetic_i"][0]), float(sc["kinetic_i"][-1])]}

    prof = run.series("profile")
    x, t = prof.x("n:e"), prof.t * w_ci
    dx = x[1] - x[0]
    temps = {n: temperature(prof, n, mass[n]) for n in ("e", "i")}
    out["t_loaded_over_theta"] = {n: float(np.nanmean(temps[n][1][0]) / (theta[n] * mass[n])) for n in temps}
    pick = np.unique(np.linspace(0, len(t) - 1, 11).astype(int))
    if not info["wall"]:
        out["uniform"] = {"omega_ci_t": t[pick].tolist()} | {
            f"T_{n}_over_loaded": (np.nanmean(temps[n][1], axis=1)[pick] / np.nanmean(temps[n][1][0])).tolist()
            for n in temps}
        if info.get("through") and info["by"]:
            zone = (x > a.margin * d_i) & (x < x[-1] - a.margin * d_i)
            e0z = info["v0"] * info["by"]
            series = {"By": prof["By"][:, zone].mean(axis=1) / info["by"] - 1.0,
                      "Ez": prof["Ez"][:, zone].mean(axis=1) / e0z - 1.0,
                      "Bz": prof["Bz"][:, zone].mean(axis=1) / info["by"]}
            fld = {"zone_di": [float(x[zone][0] / d_i), float(x[zone][-1] / d_i)], "omega_ci_t": t[pick].tolist()}
            fld |= {f"{k}_dev": v[pick].tolist() for k, v in series.items()}
            for label, sel in (("from_1", t >= 1.0), ("second_half", t >= 0.5 * t[-1])):
                if sel.sum() > 3:
                    for k in ("By", "Ez"):
                        c, cov = np.polyfit(t[sel], series[k][sel], 1, cov=True)
                        res = series[k][sel] - np.polyval(c, t[sel])
                        fld[f"{k}_slope_{label}"] = {"per_inv_omega_ci": float(c[0]),
                                                     "stderr": float(np.sqrt(cov[0, 0])),
                                                     "rms_about_fit": float(res.std())}
            out["through"] = fld
    else:
        n_i = temps["i"][0]
        xs = shock_position(x, n_i, max(1, int(round(a.smooth * d_i / dx))))
        half = len(t) // 2
        speed = float(np.polyfit(prof.t[half:], xs[half:], 1)[0]) if len(t) - half > 2 else float("nan")
        down = (x > 0.2 * xs[-1]) & (x < 0.8 * xs[-1])
        up = (x > xs[-1] + a.margin * d_i) & (x < x[-1] - a.margin * d_i)
        sh = {"omega_ci_t": t[pick].tolist(), "x_shock_di": (xs[pick] / d_i).tolist(), "speed_c": speed,
              "compression": float(n_i[-1][down].mean()) if down.any() else float("nan"),
              "expected": {"r": info["shock"]["r_mhd"], "speed_c": info["v0"] / (info["shock"]["r_mhd"] - 1.0)}}
        if up.any():
            sh["upstream"] = {
                "x_di": [float(x[up][0] / d_i), float(x[up][-1] / d_i)],
                "n_e": float(temps["e"][0][-1][up].mean()),
                "v_x_over_v0": float((prof["fluxx:e"][-1][up] / temps["e"][0][-1][up]).mean() / -info["v0"]),
                "T_e_over_loaded": float(np.nanmean(temps["e"][1][-1][up]) / np.nanmean(temps["e"][1][0])),
                "T_i_over_loaded": float(np.nanmean(temps["i"][1][-1][up]) / np.nanmean(temps["i"][1][0])),
                "By_over_B0y": float(prof["By"][-1][up].mean() / info["by"]) if info["by"] else None,
                "Ez_over_E0z": float(prof["Ez"][-1][up].mean() / (info["v0"] * info["by"])) if info["by"] else None}
        out["shock"] = sh

    tails = {}
    for n in theta:
        if f"spectrum_{n}" not in run.names("histogram"):
            continue
        h = run.series(f"spectrum_{n}")
        xe, ue = h.edges
        cnt = h[n]  # (n_t, n_x, n_u)
        tot = cnt.sum(axis=(1, 2))
        rec = {"omega_ci_t": [float(h.t[0] * w_ci), float(h.t[-1] * w_ci)], "total": [float(tot[0]), float(tot[-1])]}
        for k in (5, 8):
            above = cnt[:, :, ue[:-1] >= info["u0"] + k * np.sqrt(theta[n])].sum(axis=(1, 2))
            rec[f"share_above_{k}"] = [float(above[0] / tot[0]), float(above[-1] / tot[-1])]
        rec["per_x_last"] = cnt[-1].sum(axis=1).tolist() if len(xe) - 1 <= 100 else None
        tails[n] = rec
    out["tails"] = tails

    text = json.dumps(out, indent=1)
    print(text)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
