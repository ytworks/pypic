"""Numbers of one pilot run (plan step 4): conservation, timing, temperatures, shock position, seeds.

    python research/topics/qpar-electron-conditions/pilot_check.py <run dir> [--json <file>]

CPU only (reads the output with ``pypic.analysis``). It prints numbers. It judges only the flow-through control,
by the criteria fixed in notebook.md (2026-11-15, third review) before that run was read (``control``, ``slab``,
``velocity`` below; each holds its numbers and ``pass_*`` flags).

- Temperature of a species from the moments of the ``profile`` lineout, per cell column:
  T = (2/3) (ekin / n - m (gamma_b - 1)), v_b = flux / n. It is the non-relativistic estimate (the error is of
  relative order v^2); the first output gives its value for the loaded plasma, so read changes against that.
- Periodic box (N1) and flow-through control: the temperatures averaged over x against time.
- Flow-through control (``through``): B_y / B0y - 1, E_z / E0z - 1 and B_z / B0y, averaged over the box without
  ``margin`` d_i at each end, against time; the slope of the first two per 1 / Omega_ci with its standard error
  (straight-line fit that takes the outputs as independent) over t >= 1 / Omega_ci and over the second half.
- ``control`` (flow-through): the slope of <B_y> / B0y in the zone ``--zone`` (30-37 d_i) over the outputs of
  ``--window`` (t = 4-40 / Omega_ci), with the standard error of a straight-line fit that takes the outputs as
  independent. Criteria: slope + 2 standard errors <= 3.3e-4 per 1 / Omega_ci; the slopes of the two halves of the
  window differ by at most 2 standard errors of their difference.
- ``slab`` (no shock; ``--slab``, 20-30 d_i; outputs from the start of the window): the electron-ion drift
  (flux_e - flux_i) / (n0 v0) per component, mean over slab and outputs, standard error from the scatter between
  outputs; criterion |mean| <= 0.002. dB_perp: the rms over the slab of B_y - <B_y> and of B_z, in units of B0,
  mean over the outputs (reported; the periodic N1 run gives the level to compare with).
- ``velocity`` (histogram ``velocity_e`` of the electrons in (u_x, u_y, u_z), bin centres): second moments about
  the mean, along B (T_par) and the mean of the two directions across it (T_perp), each relative to the first
  output, which holds the loaded plasma in the same bins; criterion: within 2 % at every output of the window.
  Shares with |u - u_drift| > k sqrt(theta_e), k = 2, 3, 4, pooled over the outputs of the window, against the
  first output; z = difference / sigma with sigma from the counts of both; criterion |z| <= 3.
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


def line(t, y):
    """(slope, standard error) of a straight-line fit."""
    c, cov = np.polyfit(t, y, 1, cov=True)
    return float(c[0]), float(np.sqrt(cov[0, 0]))


def control(prof, info, t, zone, window, limit: float = 3.3e-4) -> dict:
    x = prof.x("By") / info["d_i"]
    z = (x >= zone[0]) & (x <= zone[1])
    by = prof["By"][:, z].mean(axis=1) / info["by"]
    sel = (t >= window[0]) & (t <= window[1])
    mid = 0.5 * (window[0] + window[1])
    s, e = line(t[sel], by[sel])
    (s1, e1), (s2, e2) = (line(t[h], by[h]) for h in (sel & (t <= mid), sel & (t >= mid)))
    tol = 2.0 * float(np.hypot(e1, e2))
    return {"zone_di": list(zone), "window": list(window), "outputs": int(sel.sum()),
            "By_over_B0y": {"first": float(by[0]), "window_start": float(by[sel][0]), "last": float(by[-1]),
                            "min": float(by[sel].min()), "max": float(by[sel].max())},
            "slope": s, "stderr": e, "slope_plus_2se": s + 2.0 * e, "limit": limit,
            "pass_slope": bool(s + 2.0 * e <= limit),
            "halves": [{"slope": s1, "stderr": e1}, {"slope": s2, "stderr": e2}],
            "halves_difference": abs(s1 - s2), "halves_tolerance": tol, "pass_halves": bool(abs(s1 - s2) <= tol)}


def slab_numbers(prof, info, t, slab, t_min: float, limit: float = 0.002) -> dict:
    d_i = info["d_i"]
    sel = t >= t_min
    out = {"slab_di": list(slab), "from_omega_ci_t": t_min, "outputs": int(sel.sum())}
    zn = (prof.x("n:e") >= slab[0] * d_i) & (prof.x("n:e") <= slab[1] * d_i)
    n0 = float(prof["n:i"][0][zn].mean())
    drift = {}
    for c in "xyz":
        d = (prof[f"flux{c}:e"][:, zn].mean(axis=1) - prof[f"flux{c}:i"][:, zn].mean(axis=1)) / (n0 * info["v0"])
        drift[c] = {"mean": float(d[sel].mean()), "stderr": float(d[sel].std(ddof=1) / np.sqrt(sel.sum())),
                    "first_output": float(d[0])}
    out["drift_over_n0_v0"] = drift
    out["drift_limit"] = limit
    out["pass_drift"] = bool(all(abs(v["mean"]) <= limit for v in drift.values()))
    zb = (prof.x("By") >= slab[0] * d_i) & (prof.x("By") <= slab[1] * d_i)
    by, bz = prof["By"][:, zb], prof["Bz"][:, zb]
    rms_y = (by - by.mean(axis=1, keepdims=True)).std(axis=1) / info["b"]
    rms_z = np.sqrt((bz * bz).mean(axis=1)) / info["b"]
    out["dB_perp_over_B0"] = {"By": float(rms_y[sel].mean()), "Bz": float(rms_z[sel].mean()),
                              "both": float(np.sqrt(rms_y ** 2 + rms_z ** 2)[sel].mean()),
                              "first_output": float(np.sqrt(rms_y[0] ** 2 + rms_z[0] ** 2))}
    return out


def velocity(run, info, w_ci: float, t_min: float) -> dict:
    h = run.series("velocity_e")
    cen = [0.5 * (e[1:] + e[:-1]) for e in h.edges]
    (_, ny, nz), dx = run.info["n_cells"], run.info["dx"]
    per_weight = info["per_column"] / (2 * ny * nz) / (dx[0] * dx[1] * dx[2])  # particles per unit of weight
    cnt = np.asarray(h["e"], dtype=np.float64) * per_weight  # particles per bin
    t = h.t * w_ci
    sel = t >= t_min
    c, s = info["bx"] / info["b"], info["by"] / info["b"]
    vth = np.sqrt(info["theta_e"])
    tot = cnt.sum(axis=(1, 2, 3))
    par_, perp, mean = [], [], []
    for w, n in zip(cnt, tot):
        m = [np.tensordot(w.sum(axis=tuple(a for a in range(3) if a != k)), cen[k], 1) / n for k in range(3)]
        d = [cen[k] - m[k] for k in range(3)]
        var = [np.tensordot(w.sum(axis=tuple(a for a in range(3) if a != k)), d[k] ** 2, 1) / n for k in range(3)]
        cxy = float(d[0] @ w.sum(axis=2) @ d[1]) / n
        tp = c * c * var[0] + s * s * var[1] + 2.0 * c * s * cxy
        par_.append(tp)
        perp.append(0.5 * (sum(var) - tp))
        mean.append(m)
    par_, perp, mean = np.array(par_), np.array(perp), np.array(mean)
    out = {"omega_ci_t": [float(t[0]), float(t[-1])], "outputs": int(sel.sum()), "from_omega_ci_t": t_min,
           "particles": {"first": float(tot[0]), "window_mean": float(tot[sel].mean())},
           "loaded": {"T_par_over_theta": float(par_[0] / info["theta_e"]),
                      "T_perp_over_theta": float(perp[0] / info["theta_e"]),
                      "mean_u": [float(v) for v in mean[0]]},
           "mean_u_window": [float(v) for v in mean[sel].mean(axis=0)]}
    for name, v in (("T_par", par_), ("T_perp", perp)):
        r = v[sel] / v[0]
        out[name] = {"mean": float(r.mean()), "min": float(r.min()), "max": float(r.max())}
    out["pass_temperature"] = bool(all(abs(out[k][q] - 1.0) <= 0.02 for k in ("T_par", "T_perp")
                                       for q in ("min", "max")))
    r = np.sqrt((cen[0][:, None, None] + info["u0"]) ** 2 + cen[1][None, :, None] ** 2
                + cen[2][None, None, :] ** 2) / vth
    tails = {}
    for k in (2, 3, 4):
        above = cnt[:, r > k].sum(axis=1)
        f0, n0 = above[0] / tot[0], tot[0]
        fp, n_p = above[sel].sum() / tot[sel].sum(), tot[sel].sum()
        sigma = float(np.sqrt(f0 * (1.0 - f0) / n0 + fp * (1.0 - fp) / n_p))
        tails[str(k)] = {"first": float(f0), "window": float(fp), "sigma": sigma, "z": float((fp - f0) / sigma)}
    out["share_above_k_vth"] = tails
    out["pass_tails"] = bool(all(abs(v["z"]) <= 3.0 for v in tails.values()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--json")
    ap.add_argument("--zone", type=float, nargs=2, default=(30.0, 37.0), help="d_i (control)")
    ap.add_argument("--window", type=float, nargs=2, default=(4.0, 40.0), help="1 / Omega_ci (control)")
    ap.add_argument("--slab", type=float, nargs=2, default=(20.0, 30.0), help="d_i (drift, dB_perp)")
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
            if (t >= a.window[0]).sum() >= 8 and x[-1] >= a.zone[1] * d_i:
                out["control"] = control(prof, info, t, a.zone, (a.window[0], min(a.window[1], t[-1] + 1e-9)))
        if info["by"] and (t >= a.window[0]).sum() >= 8 and x[-1] >= a.slab[1] * d_i:
            out["slab"] = slab_numbers(prof, info, t, a.slab, a.window[0])
        if "velocity_e" in run.names("histogram"):
            out["velocity"] = velocity(run, info, w_ci, a.window[0])
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
