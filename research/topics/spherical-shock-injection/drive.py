"""How long the debris drives the front of a run (shock-forming-run.md, Sect. 5): the observables O1-O4.

  drive <run>   after `front`; reads analysis/front.npz and the planner block, writes analysis/drive.json

The front and its valid domain are those of protocol Sect. 4, unchanged. A bin is *driven* while the debris edge
(0.1 n0) still moves outwards at DRIVEN times the radial speed of the front or faster; both speeds are slopes
over the same five outputs. A bin is *strong* at a fast Mach number of MACH_STRONG or more. The duration of
k outputs in a row is (k - 1) output intervals.
"""
from __future__ import annotations

from pathlib import Path

import census
import numpy as np
import shockkin as sk
from r2checks import _pct, _write

DRIVEN = 0.5
MACH_STRONG = 2.0


def longest_run(mask, t):
    """Per column of mask (n_t, n_theta): the longest duration over which it is true in a row (0 for one
    output, NaN for none), in the units of t."""
    m = np.asarray(mask, bool)
    best = np.full(m.shape[1], np.nan)
    start = np.full(m.shape[1], -1)
    for i in range(m.shape[0]):
        start = np.where(m[i] & (start < 0), i, np.where(m[i], start, -1))
        span = np.where(start >= 0, t[i] - t[np.maximum(start, 0)], np.nan)
        best = np.fmax(best, span)
    return best


def _span(mask, t, theta, omega_ci: float) -> dict:
    """Bins of a mask and its longest duration at one polar angle, in 1 / Omega_ci."""
    runs = longest_run(mask, t) * omega_ci
    if not np.isfinite(runs).any():
        return {"bins": 0}
    j = int(np.nanargmax(runs))
    return {"bins": int(mask.sum()), "polar angles": int(np.isfinite(runs).sum()),
            "longest [1/Omega_ci]": float(runs[j]), "at polar angle [deg]": float(np.degrees(theta[j])),
            "median over the polar angles [1/Omega_ci]": float(np.nanmedian(runs))}


def drive(z, pl: dict) -> dict:
    """O1-O4 from the arrays of front.npz (z) and the planner block (pl)."""
    d_i, b0, v_a = float(pl["d_i"]), float(pl["b0"]), float(pl["v_a"])
    omega_ci = b0 / (b0 / v_a) ** 2
    r_m, r_b, rho_d = float(pl["r_m"]), float(pl["r_b"]), float(pl["debris_gyroradius"])
    t, theta, r, r_deb = z["t"], z["theta"], z["r"], z["r_debris"]
    valid = np.asarray(z["valid"], bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        mach = z["v_n"] / z["v_f"]
        ratio = sk.front_kinematics(theta, t, r_deb)["r_t"] / z["r_t"]
        driven = valid & (ratio >= DRIVEN)
        strong = valid & (mach >= MACH_STRONG)
    out = {"R_m": r_m, "rho_d": rho_d, "R_m / rho_d": r_m / rho_d, "R_B": r_b, "omega_ci_inv": 1.0 / omega_ci,
           "driven: debris-edge speed over front speed at least": DRIVEN, "strong: fast Mach at least": MACH_STRONG}
    if not valid.any():
        out["O1 first valid front"] = None
        return out
    i0 = int(np.flatnonzero(valid.any(axis=1))[0])
    sel = valid[i0]
    r0, e0 = float(np.median(r[i0, sel])), float(np.nanmedian(r_deb[i0, sel]))
    out["O1 first valid front"] = {
        "t": float(t[i0]), "t Omega_ci": float(t[i0] * omega_ci), "polar angles": int(sel.sum()),
        "R, median [d_e]": r0, "R / R_m": r0 / r_m, "R / rho_d": r0 / rho_d, "R / d_i": r0 / d_i, "R / R_B": r0 / r_b,
        "debris edge / R_m": e0 / r_m, "debris-edge speed over front speed": _pct(ratio[i0, sel], (50,)).get("p50"),
        "fast Mach, median": float(np.median(mach[i0, sel]))}
    out["O2 valid"] = _span(valid, t, theta, omega_ci)
    out["O2 valid and driven"] = _span(driven, t, theta, omega_ci)
    out["O2 debris-edge speed over front speed, valid bins"] = _pct(ratio[valid])
    out["O3 fast Mach, valid bins"] = _pct(mach[valid])
    out["O3 fast Mach, driven bins"] = _pct(mach[driven])
    out["O3 valid and strong"] = _span(strong, t, theta, omega_ci)
    out["O3 valid, driven and strong"] = _span(driven & strong, t, theta, omega_ci)
    expo = []
    for j in range(valid.shape[1]):  # R ~ t^m over the valid outputs of a polar angle (1 free, 2/5 Sedov)
        s = valid[:, j] & (t > 0)
        if s.sum() >= 4:
            expo.append(np.polyfit(np.log(t[s]), np.log(r[s, j]), 1)[0])
    out["O4 exponent m of R ~ t^m per polar angle"] = _pct(expo, (0, 50, 100))
    rows = []
    for i in np.flatnonzero(np.isfinite(r).any(axis=1)):  # the history, one row per output with a front
        a = np.isfinite(r[i])
        with np.errstate(invalid="ignore"):
            rows.append({"t": float(t[i]), "t Omega_ci": float(t[i] * omega_ci), "front": int(a.sum()),
                         "valid": int(valid[i].sum()), "driven": int(driven[i].sum()),
                         "strong": int(strong[i].sum()), "R / R_m": float(np.median(r[i, a])) / r_m,
                         "debris edge / R_m": float(np.nanmedian(r_deb[i, a])) / r_m,
                         "speed ratio": _pct(ratio[i, a], (50,)).get("p50"),
                         "fast Mach": _pct(mach[i, a], (50,)).get("p50")})
    out["history"] = rows
    return out


def report(out: dict) -> None:
    first = out["O1 first valid front"]
    print(f"R_m / rho_d = {out['R_m / rho_d']:.2f} (R_m {out['R_m']:.1f}, rho_d {out['rho_d']:.1f} d_e), "
          f"1 / Omega_ci = {out['omega_ci_inv']:.1f}")
    if first is None:
        print("O1: no valid front")
        return
    print(f"O1 first valid front: t = {first['t']:.1f} = {first['t Omega_ci']:.2f} / Omega_ci at "
          f"{first['polar angles']} polar angles; R = {first['R, median [d_e]']:.1f} d_e = "
          f"{first['R / R_m']:.2f} R_m = {first['R / rho_d']:.1f} rho_d = {first['R / d_i']:.1f} d_i = "
          f"{first['R / R_B']:.2f} R_B; debris edge at {first['debris edge / R_m']:.2f} R_m, moving at "
          f"{first['debris-edge speed over front speed']:.2f} of the front speed; fast Mach "
          f"{first['fast Mach, median']:.2f}")
    for key in ("O2 valid", "O2 valid and driven", "O3 valid and strong", "O3 valid, driven and strong"):
        s = out[key]
        if s["bins"]:
            print(f"{key}: {s['bins']} bins at {s['polar angles']} polar angles; longest "
                  f"{s['longest [1/Omega_ci]']:.2f} / Omega_ci at {s['at polar angle [deg]']:.1f} deg, median "
                  f"{s['median over the polar angles [1/Omega_ci]']:.2f}")
        else:
            print(f"{key}: no bin")
    for key in ("O2 debris-edge speed over front speed, valid bins", "O3 fast Mach, valid bins",
                "O3 fast Mach, driven bins", "O4 exponent m of R ~ t^m per polar angle"):
        print(f"{key}: {out[key]}")
    print("history: t, t Omega_ci, bins with a front / valid / driven / strong, R / R_m, debris edge / R_m, "
          "speed ratio, fast Mach (medians over the bins with a front)")

    def show(x):
        return "    -" if x is None else f"{x:5.2f}"

    for row in out["history"]:
        print(f"  {row['t']:6.1f} {row['t Omega_ci']:5.2f}  {row['front']:2d} {row['valid']:2d} {row['driven']:2d} "
              f"{row['strong']:2d}  {row['R / R_m']:5.2f} {show(row['debris edge / R_m'])} "
              f"{show(row['speed ratio'])} {show(row['fast Mach'])}")


def cmd_drive(args) -> int:
    run_dir = Path(args.run)
    dest = run_dir / "analysis"
    run = census._open(run_dir)
    out = {"run": str(run_dir), "commit": run.info.get("git_commit"),
           **drive(np.load(dest / "front.npz"), run.meta["planner"])}
    report(out)
    _write(dest / "drive.json", out)
    return 0


# ------------------------------------------------------------------ check on synthetic data
def selfcheck(check) -> None:
    # a front at 0.2 c; the debris edge follows at 0.15 c and stops at t = 200; valid from t = 100;
    # fast Mach 2.5 before t = 300 and 1.8 after it. With slopes over five outputs the edge keeps half the
    # front speed up to t = 190: 10 driven outputs (90 time units), 20 strong ones (190), 31 valid (300).
    t = np.arange(41) * 10.0
    theta = np.radians(np.arange(4) * 22.5 + 11.25)
    r = np.repeat((50.0 + 0.2 * t)[:, None], 4, axis=1)
    edge = np.repeat(np.minimum(45.0 + 0.15 * t, 75.0)[:, None], 4, axis=1)
    kin = sk.front_kinematics(theta, t, r)
    v_f = np.repeat(np.where(t < 300, 0.08, 0.2 / 1.8)[:, None], 4, axis=1)
    z = {"t": t, "theta": theta, "r": r, "r_debris": edge, "valid": np.repeat((t >= 100)[:, None], 4, axis=1),
         "v_n": np.full(r.shape, 0.2), "v_f": v_f, "r_t": kin["r_t"]}
    pl = {"d_i": 4.0, "b0": 0.4, "v_a": 0.1, "r_m": 80.0, "r_b": 147.0, "debris_gyroradius": 10.0}
    out = drive(z, pl)  # 1 / Omega_ci = 40
    got = [out[k].get("longest [1/Omega_ci]") for k in ("O2 valid", "O2 valid and driven", "O3 valid and strong")]
    first = out["O1 first valid front"]
    check("drive, driven and strong durations", np.allclose(got, [7.5, 2.25, 4.75])
          and out["O2 valid and driven"]["bins"] == 40 and first["t"] == 100.0
          and abs(first["R / R_m"] - 70.0 / 80.0) < 1e-12
          and abs(first["debris-edge speed over front speed"] - 0.75) < 1e-9,
          f"valid, driven, strong: {got} / Omega_ci (7.5, 2.25, 4.75); first valid at t = {first['t']:.0f}")
