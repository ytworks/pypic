"""Checks asked for by review R2 of the exploratory results (reviews/2026-11-13-results-r2.md).

They were written after the runs had been read, so none of them is a registered test: they answer the
reviewer's questions, add the reports that protocol Sect. 4 asks for, and give the counts by which the results
document was corrected.
  r2null <null run>           retention of the magnetic moment without a driver   -> analysis/r2null.json
  r2 <run> --null <null run>  the checks of a run with a front                    -> analysis/r2.json
"""
from __future__ import annotations

import json
import warnings
from math import erfc, sqrt
from pathlib import Path

import census
import numpy as np
import shockkin as sk

SPEED_EDGES = (0.0, 0.1, 0.2, 0.3, 1.0)  # the speed bins of the command null
THERMAL_BIN = 1  # the bin that holds the bulk at T = 0.01 (v_te = 0.1)
V_HT_SIGNAL = sqrt(0.2)  # a cold electron leaves a reflection with gamma = 1.5 at this V_HT
RAMP_HALF_DI = 3.0  # the ramp is looked for within this distance of the fitted mid-ramp surface [d_i]


def _pct(a, qs=(10, 50, 90)) -> dict:
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"n": 0}
    return {"n": int(len(a)), **{f"p{q:g}": float(np.percentile(a, q)) for q in qs}}


def _strict(obj):
    """obj with numpy scalars as Python numbers and every float that is not finite as None (strict JSON)."""
    if isinstance(obj, dict):
        return {str(k): _strict(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_strict(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _strict(obj.tolist())
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def _write(path: Path, out: dict) -> None:
    path.write_text(json.dumps(_strict(out), indent=1, allow_nan=False) + "\n")
    print(f"wrote {path}")


# ------------------------------------------------------------------ the magnetic moment without a driver
def window_ratio(mu, n_w: int, n_gap: int, starts):
    """mu after over mu before for orbits mu (n_t, n): each is the median over n_w samples (one gyro-period, as
    orbits.measure takes mu_in and mu_out), the two windows n_gap samples apart. One row per start sample; NaN
    where the second window does not fit."""
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # orbits that left the box are NaN
        for s in starts:
            b = s + n_w + n_gap
            if b + n_w > mu.shape[0]:
                rows.append(np.full(mu.shape[1], np.nan))
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                rows.append(np.nanmedian(mu[b:b + n_w], axis=0) / np.nanmedian(mu[s:s + n_w], axis=0))
    return np.array(rows)


def kept(ratio, tol: float = 0.2) -> dict:
    """Share of the ratios within tol of 1 (the criterion of orbits.summarize), and their percentiles."""
    r = np.asarray(ratio, float)
    r = r[np.isfinite(r)]
    if not len(r):
        return {"n": 0}
    return {"n": int(len(r)), "within 20 %": float(np.mean(np.abs(r - 1.0) <= tol)),
            **{f"p{q}": float(np.percentile(r, q)) for q in (10, 50, 90)}}


def cmd_r2null(args) -> int:
    """R2 question 1: how many electrons keep mu within 20 % over the time a reflection takes, in a uniform
    field without a front. The answer is the numerical loss alone."""
    run_dir = Path(args.run)
    run = census._open(run_dir)
    b0 = float(run.meta["planner"]["b0"])
    tr = run.tracks("tracks_e", "e")
    t = np.asarray(tr.t, float)
    dts = float(t[1] - t[0])
    gyro = 2.0 * np.pi / b0
    n_w = max(4, int(round(gyro / dts)))
    gaps = [0.0, 3.0 * gyro] + [float(g) for g in args.gaps]
    starts = [int(round(s / dts)) for s in args.starts]
    n = len(tr.ids)
    ratios: dict = {g: [] for g in gaps}
    which = []
    for a in range(0, n, args.chunk):
        sel = slice(a, min(n, a + args.chunk))
        u = np.stack([np.asarray(tr.f[k][:, sel], dtype=float) for k in ("ux", "uy", "uz")])
        b = np.stack([np.asarray(tr.f[k][:, sel], dtype=float) for k in ("Bx", "By", "Bz")])
        gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
        stay = np.isfinite(gam).all(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            speed = np.nanmedian(np.sqrt(np.sum(u ** 2, axis=0)) / gam, axis=0)
        mu = sk.magnetic_moment(u, b)
        mu[:, ~stay] = np.nan
        which.append(np.where(stay, np.digitize(speed, SPEED_EDGES) - 1, -1))
        for g in gaps:
            ratios[g].append(window_ratio(mu, n_w, int(round(g / dts)), starts))
    which = np.concatenate(which)
    rows = []
    for g in gaps:
        r = np.concatenate(ratios[g], axis=1)
        row = {"gap": g, "gap in gyro-periods": g / gyro, "starts used": int(np.isfinite(r).any(axis=1).sum()),
               "all": kept(r[:, which >= 0]), "by speed bin": []}
        for ib in range(len(SPEED_EDGES) - 1):
            row["by speed bin"].append({"v": [SPEED_EDGES[ib], SPEED_EDGES[ib + 1]], **kept(r[:, which == ib])})
        rows.append(row)
    out = {"run": str(run_dir), "commit": run.info.get("git_commit"), "gyro_period": gyro, "window": n_w * dts,
           "sample_step": dts, "starts": [float(s) for s in args.starts], "n_stay": int((which >= 0).sum()),
           "retention": rows,
           "note": "mu = u_perp^2 / 2|B| with the local B, median over one gyro-period; the gap is the time "
                   "between the end of the first window and the start of the second"}
    print(f"mu after a gap over mu before, {out['n_stay']} electrons that stay, window {out['window']:.1f}")
    print("    gap  gyro-p.   all: within 20 %, p10 / p50 / p90        thermal bin (v 0.1-0.2): within 20 %")
    for row in rows:
        a, th = row["all"], row["by speed bin"][THERMAL_BIN]
        print(f"  {row['gap']:5.0f}  {row['gap in gyro-periods']:5.1f}   {a['within 20 %']:.3f}   {a['p10']:.2f} / "
              f"{a['p50']:.2f} / {a['p90']:.2f}     {th['within 20 %']:.3f} (n {th['n']})")
    dest = run_dir / "analysis"
    dest.mkdir(exist_ok=True)
    _write(dest / "r2null.json", out)
    return 0


# ------------------------------------------------------------------ the front: bins before the valid domain
def prevalid(z, d_i: float, omega_ci: float) -> dict:
    """The bins with a front before the first valid time of their polar angle, by the criteria they fail
    (R2 change 5). Only the polar angles that reach the valid domain are counted."""
    det, valid = np.isfinite(z["r"]), z["valid"].astype(bool)
    first = np.where(valid.any(axis=0), valid.argmax(axis=0), 0)
    pre = det & (np.arange(det.shape[0])[:, None] < first[None, :])
    fb, fc, fe = (pre & ~z[k].astype(bool) for k in ("b", "c", "e"))
    with np.errstate(invalid="ignore", divide="ignore"):
        mach = z["v_n"] / z["v_f"]
        need = 0.7 * z["v_n"] / omega_ci / d_i
        lead = (z["r"] - z["r_debris"]) / d_i
    pc = pre & z["c"].astype(bool)
    never = det & ~valid.any(axis=0)[None, :]
    return {"bins": int(pre.sum()), "fail b": int(fb.sum()), "fail c": int(fc.sum()), "fail e": int(fe.sum()),
            "fail b alone": int((fb & ~fc & ~fe).sum()), "fail b and c": int((fb & fc).sum()),
            "fail b and e": int((fb & fe).sum()), "fail none of b, c, e (too short a row)":
                int((pre & ~fb & ~fc & ~fe).sum()),
            "pass c": int(pc.sum()), "fast Mach of those that pass c": _pct(mach[pc], (0, 50, 100)),
            "lead required by b [d_i]": _pct(need[pre], (0, 50, 100)), "lead measured [d_i]": _pct(lead[pre]),
            "bins at polar angles that never reach the valid domain": int(never.sum())}


def sect4_reports(z, fj: dict, d_i: float, omega_ci: float) -> dict:
    """The reports protocol Sect. 4 asks for (R2 change 4): the density front against the mid-ramp surface,
    the debris edge at its three thresholds, the equatorial plane against the B_z definition, the usable bins,
    and criterion (b) at the bin next to the equator."""
    valid, usable = z["valid"].astype(bool), z["usable"].astype(bool)
    deb = z["r_debris_thr"]
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)  # outputs without a front
        d_eq = np.nanmedian(z["eq_n"] - z["eq_b"], axis=1) if "eq_n" in z.files else np.array([])
        need = 0.7 * z["v_n"] / omega_ci
    seen = np.isfinite(d_eq)
    rows = []
    for i in range(len(z["t"])):
        if not np.isfinite(z["r"][i, -1]):
            continue
        rows.append({"t": float(z["t"][i]), "R": float(z["r"][i, -1]),
                     "debris edge at 0.03, 0.1, 0.3 n0": deb[i, :, -1],
                     "lead required": float(need[i, -1]), **{k: bool(z[k][i, -1]) for k in ("b", "c", "e", "valid")}})
    return {"valid bins": int(valid.sum()), "usable bins": int(usable.sum()),
            "density front minus mid-ramp surface, valid bins [d_e]": _pct(z["gap_b"][valid]),
            "debris edge, 0.3 n0 minus 0.03 n0, valid bins [d_e]": _pct((deb[:, 2] - deb[:, 0])[valid]),
            "debris edge, 0.1 n0 minus 0.03 n0, valid bins [d_e]": _pct((deb[:, 1] - deb[:, 0])[valid]),
            "valid bins without a debris edge at 0.03, 0.1, 0.3 n0":
                [int((valid & ~np.isfinite(deb[:, k])).sum()) for k in range(3)],
            "front minus debris edge (0.1 n0), valid bins [d_i]": _pct(((z["r"] - z["r_debris"]) / d_i)[valid]),
            "equatorial plane, density front minus B_z >= 1.25 B0 [d_e]":
                {"outputs": int(seen.sum()), **({"median": float(np.median(d_eq[seen])),
                                                 "min": float(d_eq[seen].min()), "max": float(d_eq[seen].max())}
                                                if seen.any() else {})},
            "fit rms, density front and mid-ramp surface [d_e]": [fj.get("fit_rms"), fj.get("fit_rms_B")],
            "bin next to the equator": rows}


# ------------------------------------------------------------------ the regime of the valid shocks
def perpendicular_compression(m_a: float, beta: float, gamma: float = 5.0 / 3.0) -> float:
    """Density and field compression of a perpendicular MHD shock with Alfven Mach number m_a and upstream
    beta = 2 p / B^2: the positive root of
    2 (2 - g) X^2 + g [2 beta + (g - 1) beta M_s^2 + 2] X - g (g + 1) beta M_s^2 = 0, M_s^2 = 2 M_A^2 / (g beta)."""
    m_s2 = 2.0 * m_a ** 2 / (gamma * beta)
    a = 2.0 * (2.0 - gamma)
    b = gamma * (2.0 * beta + (gamma - 1.0) * beta * m_s2 + 2.0)
    c = -gamma * (gamma + 1.0) * beta * m_s2
    return -c / b if a == 0 else (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a)


def regime(z, fj: dict, pl: dict) -> dict:
    """Dimensionless numbers of the valid shocks (R2 changes 6 and 8): Mach numbers, beta, V_n over the electron
    thermal speed, the mass-ratio measure, the thresholds in the natural energies, and the first valid front in
    three readings."""
    d_i, b0, v_a = float(pl["d_i"]), float(pl["b0"]), float(pl["v_a"])
    m_i = (b0 / v_a) ** 2
    omega_ci = b0 / m_i
    rho_d = float(pl["debris_gyroradius"])
    valid = z["valid"].astype(bool)
    it, ith = np.nonzero(valid)
    up = fj["upstream"]
    t_e = np.interp(z["t"][it], up["t"], up["T_e"])
    t_i = np.interp(z["t"][it], up["t"], up["T_i_used"])
    v_n = z["v_n"][valid]
    v_te = np.sqrt(t_e)
    vm, tm, tim = float(np.median(v_n)), float(np.median(t_e)), float(np.median(t_i))
    beta_e, beta = 2.0 * tm / b0 ** 2, 2.0 * (tm + tim) / b0 ** 2
    out = {"m_i": m_i, "omega_ci_inv": 1.0 / omega_ci, "R_m": float(pl["r_m"]), "rho_d": rho_d,
           "R_m / rho_d": float(pl["r_m"]) / rho_d, "R_B": float(pl["r_b"]),
           "V_n": _pct(v_n), "M_A": _pct(v_n / v_a), "fast Mach": _pct(v_n / z["v_f"][valid]),
           "theta_Bn [deg]": _pct(np.degrees(z["tbn"][valid]), (0, 50, 100)),
           "T_e upstream, median": tm, "T_i upstream, median": tim, "beta_e": beta_e, "beta": beta,
           "V_n / v_te": _pct(v_n / v_te), "(m_i / m_e) / M_A at the median V_n": m_i / (vm / v_a),
           "gamma - 1 = 0.5 over m_e V_n^2 / 2": 0.5 / (0.5 * vm ** 2),
           "gamma - 1 = 0.5 over m_i V_n^2 / 2": 0.5 / (0.5 * m_i * vm ** 2),
           "gamma - 1 = 0.5 over T_e": 0.5 / tm,
           "share of a Maxwellian faster than V_n along the normal": 0.5 * erfc(vm / (sqrt(2.0 * tm))),
           "compression b, valid bins": _pct(z["b_mean"][valid])}
    # near-perpendicular bins against the MHD jump (the two adiabatic indices bracket the kinetic answer)
    perp = valid & (np.degrees(z["tbn"]) >= 80.0)
    if perp.any():
        m_a = float(np.median(z["v_n"][perp])) / v_a
        out["theta_Bn >= 80 deg"] = {
            "bins": int(perp.sum()), "M_A": m_a, "b measured": _pct(z["b_mean"][perp]),
            "b of a perpendicular MHD shock, gamma 5/3 and 2": [perpendicular_compression(m_a, beta, g)
                                                                for g in (5.0 / 3.0, 2.0)]}
    first = np.where(valid.any(axis=0), valid.argmax(axis=0), -1)
    cols = np.flatnonzero(first >= 0)
    if len(cols):
        r0 = z["r"][first[cols], cols]
        t0 = z["t"][first[cols]]
        j = cols[-1]  # the bin closest to the equator
        r_eq, t_eq = float(z["r"][first[j], j]), float(z["t"][first[j]])
        out["first valid front"] = {
            "polar angle of the bin closest to the equator [deg]": float(np.degrees(z["theta"][j])),
            "t": t_eq, "t Omega_ci": t_eq * omega_ci, "R [d_e]": r_eq, "R / d_i": r_eq / d_i,
            "R / R_B": r_eq / float(pl["r_b"]), "R / rho_d": r_eq / rho_d, "R / R_m": r_eq / float(pl["r_m"]),
            "all polar angles: t": [float(t0.min()), float(t0.max())],
            "all polar angles: R [d_e]": [float(r0.min()), float(r0.max())]}
        sel = valid[:, j] & (z["t"] > 0)
        if sel.sum() >= 4:  # R ~ t^m in the valid domain next to the equator (1 free, 2/5 Sedov, 1/4 snowplough)
            m, _ = np.polyfit(np.log(z["t"][sel]), np.log(z["r"][sel, j]), 1)
            out["first valid front"]["exponent m of R ~ t^m next to the equator"] = float(m)
            out["first valid front"]["V_n next to the equator, first and last valid output"] = [
                float(z["v_n"][sel, j][0]), float(z["v_n"][sel, j][-1])]
    return out


# ------------------------------------------------------------------ the curvature parameter per band
def epsilon_bands(z, taus: dict) -> dict:
    """epsilon = V_n tan^2(theta_Bn) tau / R (protocol Sect. 3.1) in the valid bins, per 10-degree band of
    theta_Bn and for each dwell time of taus (R2 change 3)."""
    valid = z["valid"].astype(bool)
    tbn = z["tbn"][valid]
    with np.errstate(invalid="ignore", divide="ignore"):
        per_tau = z["v_n"][valid] * np.tan(tbn) ** 2 / z["r"][valid]
        v_ht = z["v_n"][valid] / np.cos(tbn)
    band = np.clip(np.degrees(tbn) // 10, 0, 8).astype(int)
    rows = []
    for b in range(9):
        m = band == b
        if not m.any():
            continue
        row = {"band": f"{10 * b}-{10 * b + 10}", "bins": int(m.sum()), "V_HT": _pct(v_ht[m], (50, 100))}
        for name, tau in taus.items():
            row[f"epsilon, {name}"] = float(np.median(per_tau[m]) * tau)
        rows.append(row)
    reach = v_ht >= V_HT_SIGNAL
    out = {"taus": taus, "bands": rows, "V_HT for gamma = 1.5 after one reflection of a cold electron": V_HT_SIGNAL,
           "valid bins that reach it": int(reach.sum()),
           "their theta_Bn [deg]": _pct(np.degrees(tbn[reach]), (0, 50, 100))}
    for name, tau in taus.items():
        out[f"valid bins with epsilon <= 0.3, {name}"] = int((per_tau * tau <= 0.3).sum())
        out[f"epsilon in the bins that reach it, {name}"] = _pct(per_tau[reach] * tau, (0, 50, 100))
    return out


# ------------------------------------------------------------------ the width of the ramp
def ramp_width(xi, p, lo: float = 0.1, hi: float = 0.9, search: float = 2.0):
    """Width of a normalised ramp p(xi), 0 upstream and 1 at the top, xi growing outwards: from the sample
    nearest to p = 0.5 within `search` of xi = 0, outwards to where p first falls below lo and inwards to where
    it first exceeds hi (linear interpolation). Returns (width, centre); NaN if a side is not reached."""
    xi, p = np.asarray(xi, float), np.asarray(p, float)
    near = np.flatnonzero((np.abs(xi) <= search) & np.isfinite(p))
    if not len(near):
        return np.nan, np.nan
    k0 = near[np.argmin(np.abs(p[near] - 0.5))]
    if not lo < p[k0] < hi:
        return np.nan, np.nan
    k = k0
    while k + 1 < len(p) and p[k + 1] >= lo:
        k += 1
    if k + 1 >= len(p) or not np.isfinite(p[k + 1]):
        return np.nan, np.nan
    x_lo = xi[k] + (lo - p[k]) * (xi[k + 1] - xi[k]) / (p[k + 1] - p[k])
    k = k0
    while k - 1 >= 0 and p[k - 1] <= hi:
        k -= 1
    if k - 1 < 0 or not np.isfinite(p[k - 1]):
        return np.nan, np.nan
    x_hi = xi[k] + (hi - p[k]) * (xi[k - 1] - xi[k]) / (p[k - 1] - p[k])
    return float(x_lo - x_hi), float(0.5 * (x_lo + x_hi))


def ramp_profiles(run, front_b, t_range, z, b0: float, centre, half: float, d_xi: float = 0.5):
    """|B| / B0 against xi = r - R_B(theta, t), R_B the fitted mid-ramp surface, in the bins of the polar angle:
    the cells of both meridional planes and of their four quadrants together (eight half-planes).
    Returns (xi, profile (n_t, n_xi, n_theta)); NaN at outputs outside the fit or without a valid bin."""
    n_th = len(z["theta"])
    step = float(z["theta"][1] - z["theta"][0])
    n_xi = int(round(2.0 * half / d_xi))
    xi = -half + (np.arange(n_xi) + 0.5) * d_xi
    tot = np.zeros((len(z["t"]), n_xi * n_th))
    cnt = np.zeros_like(tot)
    use = z["valid"].astype(bool).any(axis=1) & (z["t"] >= t_range[0]) & (z["t"] <= t_range[1])
    for name, axis in (("slice_y", 0), ("slice_x", 1)):
        try:
            ser = run.series(name)
        except (KeyError, FileNotFoundError, OSError):
            continue
        if not np.allclose(np.asarray(ser.t, float), z["t"]):
            raise SystemExit(f"{name} and front.npz are not at the same times")
        a, zz = ser.plane_axes("n:i")
        da = np.asarray(a, float)[:, None] - centre[axis]
        db = np.asarray(zz, float)[None, :] - centre[2]
        r, theta = np.hypot(da, db), np.arctan2(np.abs(da), np.abs(db))
        ith = np.clip(np.floor(theta / step).astype(int), 0, n_th - 1)
        for i in np.flatnonzero(use):
            bmag = np.sqrt(sum(ser.read(f"B{c}", i).astype(float) ** 2 for c in "xyz")) / b0
            k = np.floor((r - front_b(theta, float(z["t"][i]))[0] + half) / d_xi).astype(int)
            ok = (k >= 0) & (k < n_xi)
            idx = k[ok] * n_th + ith[ok]
            tot[i] += np.bincount(idx, weights=bmag[ok], minlength=n_xi * n_th)
            cnt[i] += np.bincount(idx, minlength=n_xi * n_th)
    with np.errstate(invalid="ignore", divide="ignore"):
        prof = np.where(cnt > 0, tot / cnt, np.nan).reshape(len(z["t"]), n_xi, n_th)
    prof[~use] = np.nan
    return xi, prof


def ramp_bands(xi, prof, z, layer: float) -> dict:
    """Ramp width (10-90 % and 20-80 % of the rise to the bin-average compression, along the normal) and its
    crossing time width / V_n in the valid bins, per 10-degree band of theta_Bn; `layer` is the width of the
    layer of the census (2 d_i), for the second reading of the registered tau of N3."""
    valid = z["valid"].astype(bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosn = z["r"] / np.hypot(z["r"], z["r_theta"])  # a radial length times this is the length along the normal
    rec: dict = {k: [] for k in ("band", "w", "w28", "tau", "tau_layer", "v_n")}
    shifted: dict = {}
    for i, j in zip(*np.nonzero(valid)):
        top = z["b_mean"][i, j]
        if not np.isfinite(top) or top <= 1.05 or not np.isfinite(prof[i, :, j]).all():
            continue
        p = (prof[i, :, j] - 1.0) / (top - 1.0)
        w, mid = ramp_width(xi, p)
        w28, _ = ramp_width(xi, p, 0.2, 0.8)
        band = int(np.clip(np.degrees(z["tbn"][i, j]) // 10, 0, 8))
        rec["band"].append(band)
        rec["w"].append(w * cosn[i, j])
        rec["w28"].append(w28 * cosn[i, j])
        rec["v_n"].append(z["v_n"][i, j])
        rec["tau"].append(w * cosn[i, j] / z["v_n"][i, j])
        rec["tau_layer"].append(layer / z["v_n"][i, j])
        if np.isfinite(mid):
            shifted.setdefault(band, []).append((np.interp(xi + mid, xi, p), cosn[i, j]))
    rec = {k: np.array(v, float) for k, v in rec.items()}
    rows = []
    for b in sorted(set(rec["band"].astype(int).tolist())):
        m = rec["band"] == b
        row = {"band": f"{10 * b}-{10 * b + 10}", "bins": int(m.sum()),
               "bins with a ramp found": int(np.isfinite(rec["w"][m]).sum()),
               "width 10-90 % [d_e]": _pct(rec["w"][m]), "width 20-80 % [d_e]": _pct(rec["w28"][m]),
               "V_n": float(np.median(rec["v_n"][m])), "crossing time, ramp 10-90 %": _pct(rec["tau"][m]),
               "crossing time, layer of the census": _pct(rec["tau_layer"][m])}
        if b in shifted:
            mean = np.mean([s for s, _ in shifted[b]], axis=0)
            row["width 10-90 % of the mean profile [d_e]"] = ramp_width(xi, mean)[0] * float(
                np.median([c for _, c in shifted[b]]))
        rows.append(row)
    return {"bands": rows, "all valid bins": {
        "bins": int(len(rec["w"])), "bins with a ramp found": int(np.isfinite(rec["w"]).sum()),
        "width 10-90 % [d_e]": _pct(rec["w"]), "width 20-80 % [d_e]": _pct(rec["w28"]),
        "crossing time, ramp 10-90 %": _pct(rec["tau"]), "crossing time, layer of the census": _pct(rec["tau_layer"])},
        "layer of the census [d_e]": layer,
        "note": "xi is measured from the fitted mid-ramp surface, so the fit error inside a bin and the output "
                "average over eight half-planes widen the ramp: the widths are upper bounds"}


def thermal_nu_tau(nu_null, ramp: dict, dwell: dict) -> dict:
    """N3 as registered (protocol Sect. 5): nu of the thermal electrons of the null run times the crossing time
    of the layer; and the same rate times the dwell times of the reflections. nu is given at the three lags."""
    nu = [float(v) for v in nu_null["nu_local"][THERMAL_BIN]]
    out = {"nu thermal (v 0.1-0.2) at lags 4, 16, 64": nu, "bands": []}
    for row in ramp["bands"]:
        ent = {"band": row["band"], "bins": row["bins"]}
        for name in ("crossing time, ramp 10-90 %", "crossing time, layer of the census"):
            tau = row[name].get("p50")
            ent[name] = tau
            ent["nu tau, " + name.split(", ")[1]] = None if tau is None else [min(nu) * tau, max(nu) * tau]
        out["bands"].append(ent)
    for name in ("crossing time, ramp 10-90 %", "crossing time, layer of the census"):
        tau = ramp["all valid bins"][name].get("p50")
        out["all valid bins, nu tau, " + name.split(", ")[1]] = None if tau is None else [min(nu) * tau,
                                                                                         max(nu) * tau]
    out["nu tau at the dwell times"] = {name: [min(nu) * tau, max(nu) * tau] for name, tau in dwell.items()}
    return out


# ------------------------------------------------------------------ the energetic electrons by id
def id_history(cz, t_valid0: float) -> dict:
    """The ambient electrons with gamma >= 1.5 of the census per dump and class, what becomes of those of the
    largest dump, and how many of the later ones were already energetic before the first valid shock
    (R2 changes 1 and 6, question 4). The census holds an electron in a dump only with gamma >= 1.25."""
    ids, gam, dump, cls = cz["id"], cz["gamma"], cz["dump"], cz["cls"]
    times = np.asarray(cz["times"], float)
    amb = cz["origin"] <= 2
    hot = amb & (gam >= 1.5)
    names = ("upstream", "layer", "downstream", "no front")
    which = np.where(cls <= 1, 0, cls - 1)
    series = []
    for k, t in enumerate(times):
        m = hot & (dump == k)
        series.append({"t": float(t), "n": int(m.sum()),
                       **{nm: int((m & (which == c)).sum()) for c, nm in enumerate(names)}})
    peak = int(np.argmax([s["n"] for s in series]))
    ids_peak = ids[hot & (dump == peak)]
    follow = []
    for k in range(peak, len(times)):
        m = amb & (dump == k) & np.isin(ids, ids_peak)
        g = gam[m]
        n15 = int((g >= 1.5).sum())
        if k > peak and not len(g) and follow and follow[-1]["gamma >= 1.25"] == 0:
            continue
        follow.append({"t": float(times[k]), "gamma >= 1.5": n15, "gamma >= 1.25": int(len(g)),
                       "median gamma of those above 1.25": float(np.median(g)) if len(g) else None,
                       "median distance to the front [d_e]": float(np.nanmedian(cz["delta"][m])) if len(g) else None,
                       **{nm: int((which[m] == c).sum()) for c, nm in enumerate(names)}})
    early = dump < int(np.searchsorted(times, t_valid0 - 1e-9))
    ids15, ids125 = np.unique(ids[hot & early]), np.unique(ids[amb & early])
    late = hot & ~early
    late_ids = np.unique(ids[late])
    up = late & (cls <= 1)
    out = {"per dump": [s for s in series if s["n"]], "largest dump": series[peak],
           "the electrons of the largest dump later": follow, "first valid time": float(t_valid0),
           "from the first valid time on": {
               "electron-dumps with gamma >= 1.5": int(late.sum()), "different electrons": int(len(late_ids)),
               "of these, already above 1.5 before": int(np.isin(late_ids, ids15).sum()),
               "of these, already above 1.25 before": int(np.isin(late_ids, ids125).sum()),
               "upstream electron-dumps": int(up.sum()), "different upstream electrons": int(len(np.unique(ids[up]))),
               "upstream, already above 1.25 before": int(np.isin(np.unique(ids[up]), ids125).sum())}}
    # which upstream electrons the comparison with the count model sees (baseline.measured_counts)
    upv = hot & (cls <= 1)
    out["upstream electron-dumps with gamma >= 1.5"] = {
        "all": int(upv.sum()), "in the valid domain": int((upv & (cz["flag"] >= 2)).sum()),
        "in the valid domain, with a reflection-point band": int((upv & (cz["flag"] >= 2)
                                                                   & (cz["band"] < census.NO_FOOT)).sum()),
        "by result": {census.RESULTS[r]: int((upv & (cz["result"] == r)).sum()) for r in range(len(census.RESULTS))
                      if (upv & (cz["result"] == r)).any()}}
    k1 = cz["k1"].astype(bool)
    out["K1"] = {"electron-dumps": int(k1.sum()), "electrons": [
        {"id": int(i), "t": [float(v) for v in cz["t"][k1 & (ids == i)]],
         "gamma": [float(v) for v in gam[k1 & (ids == i)]],
         "in the valid domain": int((k1 & (ids == i) & (cz["flag"] >= 2)).sum()),
         "beyond the debris": int((k1 & (ids == i) & cz["beyond_debris"].astype(bool)).sum())}
        for i in np.unique(ids[k1])[:20]]}
    return out


def downstream_branch(cj: dict, t_valid0: float) -> dict:
    """The downstream branch of C1 (R2 change 4): the weight above the energy cut behind the front against three
    times the Maxwell-Juttner tail at the upper fit temperature, in the outputs from the first valid time on."""
    out = {}
    for name in ("behind the front", "behind the front, outside the debris"):
        rows = [r[name] | {"t": r["t"]} for r in cj.get("downstream", []) if r["t"] >= t_valid0 - 1e-9
                and r[name].get("cells")]
        if not rows:
            out[name] = {"outputs": 0}
            continue
        worst = max(rows, key=lambda r: r["measured"])
        out[name] = {"outputs": len(rows), "measured, sum over the outputs": sum(r["measured"] for r in rows),
                     "three times the tail at 1.1 T, sum over the outputs": sum(3.0 * r["expected"][2] for r in rows),
                     "largest output": {"t": worst["t"], "measured": worst["measured"],
                                        "three times the tail at 1.1 T": 3.0 * worst["expected"][2],
                                        "cells": worst["cells"]},
                     "outputs with measured >= 100": sum(r["measured"] >= 100 for r in rows)}
    out["e_cut"] = cj["downstream"][0]["e_cut"] if cj.get("downstream") else None
    out["note"] = "weights of macro-electrons in the slab |y| < 8 d_e; all polar angles, not only the valid ones"
    return out


def tau_rule(dest: Path) -> dict:
    """The replacement rule of protocol Sect. 3.2 (R2 change 4): the results of the upstream electrons of the
    census with tau_dwell = 3 gyro-periods and with the dwell of the tracked reflections (census --tau)."""
    if not (dest / "census-tau.npz").exists():
        return {"applied": False}
    out = {"applied": True}
    for stem in ("census", "census-tau"):
        cz = np.load(dest / f"{stem}.npz")
        cj = json.loads((dest / f"{stem}.json").read_text())
        up = (cz["origin"] <= 2) & (cz["cls"] <= 1)
        ent = {"tau_dwell": cj["tau_dwell"]}
        for name, m in (("upstream electron-dumps, gamma >= 1.25", up),
                        ("upstream electron-dumps, gamma >= 1.5", up & (cz["gamma"] >= 1.5)),
                        ("upstream, valid domain, gamma >= 1.25", up & (cz["flag"] >= 2))):
            ent[name] = {"n": int(m.sum()), **{census.RESULTS[r]: int((m & (cz["result"] == r)).sum())
                                               for r in range(len(census.RESULTS)) if (m & (cz["result"] == r)).any()}}
        out[stem] = ent
    return out


# ------------------------------------------------------------------ the command
def cmd_r2(args) -> int:
    run_dir, null_dir = Path(args.run), Path(args.null)
    dest = run_dir / "analysis"
    run = census._open(run_dir)
    pl = run.meta["planner"]
    d_i, b0, v_a = float(pl["d_i"]), float(pl["b0"]), float(pl["v_a"])
    omega_ci = b0 / (b0 / v_a) ** 2
    z = np.load(dest / "front.npz")
    fj = json.loads((dest / "front.json").read_text())
    cj = json.loads((dest / "census.json").read_text())
    tj = json.loads((dest / "tracks.json").read_text())
    nj = json.loads((null_dir / "analysis" / "null.json").read_text())
    valid = z["valid"].astype(bool)
    if not valid.any():
        print("no valid bin: nothing to check")
        return 1
    t_valid0 = float(z["t"][valid.any(axis=1)].min())
    gyro = 2.0 * np.pi / b0
    refl = tj["reflections, all energies"]["dwell in the layer"]
    dwell = {"protocol (3 gyro-periods)": 3.0 * gyro, "tracked returns, median": refl["p50"],
             "tracked returns, 90 %": refl["p90"]}
    hot = tj["reflections"]["dwell in the layer"]
    if hot.get("n"):
        dwell["tracked reflections with gamma_out >= 1.25, 90 %"] = hot["p90"]
    model = census.FrontModel.from_run(run, run_dir)
    out = {"run": str(run_dir), "null": str(null_dir), "commit": run.info.get("git_commit"),
           "gyro_period": gyro, "dwell": dwell,
           "pre-valid bins": prevalid(z, d_i, omega_ci), "protocol Sect. 4": sect4_reports(z, fj, d_i, omega_ci),
           "regime": regime(z, fj, pl), "epsilon": epsilon_bands(z, dwell)}
    if model.front_b is not None:
        xi, prof = ramp_profiles(run, model.front_b, model.t_range_b, z, b0, pl["center"], RAMP_HALF_DI * d_i)
        out["ramp"] = ramp_bands(xi, prof, z, 2.0 * d_i)
        out["N3 for thermal electrons"] = thermal_nu_tau(nj["tracks"], out["ramp"], dwell)
    out["id history"] = id_history(np.load(dest / "census.npz"), t_valid0)
    out["downstream branch of C1"] = downstream_branch(cj, t_valid0)
    out["tau_dwell rule"] = tau_rule(dest)
    report(out)
    _write(dest / "r2.json", out)
    return 0


def report(out: dict) -> None:
    pv, s4, rg = out["pre-valid bins"], out["protocol Sect. 4"], out["regime"]
    print(f"R2 checks of {out['run']}")
    print(f"\npre-valid bins: {pv['bins']}; fail b {pv['fail b']}, c {pv['fail c']}, e {pv['fail e']}; b alone "
          f"{pv['fail b alone']}; pass c {pv['pass c']} with fast Mach {pv['fast Mach of those that pass c']}")
    print(f"  lead required {pv['lead required by b [d_i]']}, measured {pv['lead measured [d_i]']}")
    print(f"\nSect. 4: valid {s4['valid bins']}, usable {s4['usable bins']}; density front minus mid-ramp "
          f"{s4['density front minus mid-ramp surface, valid bins [d_e]']}")
    print(f"  equatorial plane {s4['equatorial plane, density front minus B_z >= 1.25 B0 [d_e]']}")
    print(f"  valid bins without a debris edge {s4['valid bins without a debris edge at 0.03, 0.1, 0.3 n0']}")
    print("\nregime of the valid bins")
    for k, v in rg.items():
        print(f"  {k}: {v}")
    print("\nepsilon per band of theta_Bn (bins; V_HT median, max; epsilon at each dwell time)")
    for row in out["epsilon"]["bands"]:
        eps = "  ".join(f"{v:7.2f}" for k, v in row.items() if k.startswith("epsilon"))
        print(f"  {row['band']:>6s} {row['bins']:4d}  {row['V_HT']['p50']:.3f} {row['V_HT']['p100']:.3f}  {eps}")
    print("  " + "; ".join(f"{k}: {v}" for k, v in out["epsilon"].items() if k not in ("bands", "taus")))
    if "ramp" in out:
        print("\nramp per band: bins, found, width 10-90 % (p10/p50/p90), mean profile, crossing time p50, "
              "thermal nu tau (ramp; layer of the census)")
        for row, n3 in zip(out["ramp"]["bands"], out["N3 for thermal electrons"]["bands"]):
            w = row["width 10-90 % [d_e]"]
            ws = f"{w['p10']:.1f}/{w['p50']:.1f}/{w['p90']:.1f}" if w.get("n") else "-"
            print(f"  {row['band']:>6s} {row['bins']:4d} {row['bins with a ramp found']:4d}  {ws}  "
                  f"{row.get('width 10-90 % of the mean profile [d_e]')}  {n3['crossing time, ramp 10-90 %']}  "
                  f"{n3['nu tau, ramp 10-90 %']}; {n3['nu tau, layer of the census']}")
        print(f"  all valid bins: {out['ramp']['all valid bins']}")
        print(f"  nu tau at the dwell times: {out['N3 for thermal electrons']['nu tau at the dwell times']}")
    h = out["id history"]
    print("\nambient electrons with gamma >= 1.5 per dump (t: n = upstream + layer + downstream + no front)")
    print("  " + "; ".join(f"{s['t']:.1f}: {s['n']} = {s['upstream']}+{s['layer']}+{s['downstream']}+{s['no front']}"
                           for s in h["per dump"]))
    print("the electrons of the largest dump later (t: gamma >= 1.5 / >= 1.25, median gamma, distance)")
    print("  " + "; ".join(f"{f['t']:.1f}: {f['gamma >= 1.5']}/{f['gamma >= 1.25']}, "
                           f"{f['median gamma of those above 1.25']}, {f['median distance to the front [d_e]']}"
                           for f in h["the electrons of the largest dump later"]))
    for k in ("from the first valid time on", "upstream electron-dumps with gamma >= 1.5", "K1"):
        print(f"{k}: {h[k]}")
    print(f"\ndownstream branch of C1: {out['downstream branch of C1']}")
    print(f"\ntau_dwell rule: {out['tau_dwell rule']}")


# ------------------------------------------------------------------ checks on synthetic data
def selfcheck(check) -> None:
    # a tanh ramp of scale 1.5: its 10-90 % width is 2 * 1.5 * atanh(0.8), centred where it was put
    xi = -12.0 + (np.arange(48) + 0.5) * 0.5
    p = 0.5 * (1.0 - np.tanh((xi - 0.7) / 1.5))
    w, mid = ramp_width(xi, p)
    w28, _ = ramp_width(xi, p, 0.2, 0.8)
    exact, exact28 = 3.0 * np.arctanh(0.8), 3.0 * np.arctanh(0.6)
    none = ramp_width(xi, np.full(len(xi), 0.02))[0]
    check("R2, width of a ramp", abs(w - exact) < 0.1 and abs(w28 - exact28) < 0.1 and abs(mid - 0.7) < 0.1
          and np.isnan(none),
          f"10-90 % {w:.2f} ({exact:.2f}), 20-80 % {w28:.2f} ({exact28:.2f}), centre {mid:.2f} (0.7)")
    # the perpendicular jump: 1 at a fast Mach number of 1, (g + 1) / (g - 1) for a strong shock
    g, beta = 5.0 / 3.0, 0.49
    weak = perpendicular_compression(sqrt(1.0 + 0.5 * g * beta), beta, g)
    strong = perpendicular_compression(1e4, beta, g)
    check("R2, perpendicular MHD jump", abs(weak - 1.0) < 1e-9 and abs(strong - 4.0) < 1e-3
          and abs(perpendicular_compression(1e4, beta, 2.0) - 3.0) < 1e-3,
          f"compression {weak:.6f} at fast Mach 1 (1), {strong:.4f} for a strong shock (4)")
    # mu kept exactly, and mu lost by a random walk
    rng = np.random.default_rng(3)
    mu = np.ones((400, 500))
    walk = np.exp(np.cumsum(rng.normal(0.0, 0.05, (400, 500)), axis=0))
    r0 = kept(window_ratio(mu, 20, 100, [0, 50]))
    r1 = kept(window_ratio(walk, 20, 100, [0, 50]))
    short = window_ratio(mu, 20, 370, [0])
    check("R2, retention of mu", r0["within 20 %"] == 1.0 and 0.1 < r1["within 20 %"] < 0.5
          and np.isnan(short).all(),
          f"constant mu kept {r0['within 20 %']:.2f} (1); random walk kept {r1['within 20 %']:.2f}")
