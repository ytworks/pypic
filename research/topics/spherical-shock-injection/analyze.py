"""Analysis of the exploratory runs of the topic spherical-shock-injection (CPU only).

Usage (from the repository root):
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py selftest
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py predict --out <json>
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py front <run directory>
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py cadence <run directory>

Results of a command on a run go to <run directory>/analysis/ (not in git).

The definitions and criteria are fixed in exploratory-protocol.md; the kinematics is in shockkin.py.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import shockkin as sk  # noqa: E402

GAMMA_EDGES = [1.5, 2.0, 3.0, 5.0, np.inf]  # one threshold for all runs: gamma >= 1.5
BAND_EDGES_DEG = list(range(0, 91, 10))  # foot-point theta_Bn bands


# ------------------------------------------------------------------ self-test
def _sample_mj(rng, temp, n):
    """Momenta of a Maxwell-Juttner distribution by rejection from a wider Gaussian (weights returned)."""
    s = 1.15 * np.sqrt(temp)
    u = rng.normal(0.0, s, size=(3, n))
    gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    w = np.exp(-(gam - 1.0) / temp + np.sum(u ** 2, axis=0) / (2.0 * s ** 2))
    return u, w / w.sum()


def _boris_orbits(rng, e0, b0, n, steps, dt_step, every):
    """Electron orbits of the Boris scheme in uniform fields, recorded as the engine records tracks:
    the momentum half a step before the fields. Returns u, E, B as (3, n_t, n) and the times."""
    u = rng.normal(0.0, 0.1, (3, n))
    us, ts = [], []
    for step in range(steps + 1):
        if step % every == 0:
            us.append(u.copy())
            ts.append(step * dt_step)
        um = u - 0.5 * dt_step * e0[:, None]
        tv = -0.5 * dt_step * b0[:, None] / np.sqrt(1 + np.sum(um ** 2, axis=0))
        up = um + np.cross(um, tv, axis=0)
        up = um + np.cross(up, 2 * tv / (1 + np.sum(tv ** 2, axis=0)), axis=0)
        u = up - 0.5 * dt_step * e0[:, None]
    us = np.stack(us, 1)
    shape = (1, us.shape[1], n)
    return us, np.tile(e0[:, None, None], shape), np.tile(b0[:, None, None], shape), np.array(ts)


def selftest() -> int:
    rng = np.random.default_rng(7)
    results = []

    def check(name, ok, detail):
        results.append(bool(ok))
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {detail}")

    # 1 the reflection map is its own inverse and keeps the energy in the frame of the mirror
    n = 20000
    u_par, u_perp = rng.normal(0, 0.3, n), np.abs(rng.normal(0, 0.3, n))
    v_ht = rng.uniform(0.05, 0.95, n)
    enc, _, u_r, g_r = sk.ht_reflect(u_par, u_perp, v_ht, 3.0)
    u_back, g_back, _ = sk.before_reflection(u_r, u_perp, v_ht)
    err = np.max(np.abs(u_back - u_par)[enc])
    check("involution", err < 1e-12, f"max |u_par - map(map(u_par))| = {err:.1e} over {int(enc.sum())} encounters")
    err = np.max(np.abs(g_r - np.sqrt(1 + u_r ** 2 + u_perp ** 2)))
    check("energy after reflection", err < 1e-12, f"max |gamma_r - sqrt(1 + u_r^2 + u_perp^2)| = {err:.1e}")

    # 2 a cold electron leaves with gamma = (1 + V^2) / (1 - V^2)
    v = np.array([0.1, 0.447214, 0.8])
    _, _, _, g_r = sk.ht_reflect(0.0, 0.0, v, 3.0)
    err = np.max(np.abs(g_r - (1 + v ** 2) / (1 - v ** 2)))
    check("cold electron", err < 1e-12, f"max error {err:.1e}")

    # 3 flux on a slow mirror with a closed loss cone against the non-relativistic formula
    temp, vm = 1e-4, 0.015
    sig = np.sqrt(temp)
    flux, enc_flux = sk.reflected_flux(vm, temp, 1e12, 0.0, [1.0, np.inf], [-1.0, 1.0])
    from math import erf, exp, pi, sqrt
    x = vm / sig
    exact = vm * 0.5 * (1 + erf(x / sqrt(2))) + sig * exp(-0.5 * x * x) / sqrt(2 * pi)
    rel = abs(flux.sum() - exact) / exact
    check("planar flux", rel < 0.01 and abs(enc_flux - exact) / exact < 0.01,
          f"model {flux.sum():.6e}, encounter {enc_flux:.6e}, formula {exact:.6e} (rel {rel:.1e})")

    # 4 the quadrature against a Monte Carlo sum, finite loss cone, relativistic mirror
    temp, vm, b = 0.01, 0.47, 3.0
    edges = np.array(GAMMA_EDGES)
    flux, _ = sk.reflected_flux(vm, temp, b, 0.0, edges, [-1.0, 1.0])
    mc = np.zeros(len(edges) - 1)
    for _ in range(8):
        u, w = _sample_mj(rng, temp, 2_000_000)
        up, uperp = u[2], np.hypot(u[0], u[1])
        gam = np.sqrt(1 + up ** 2 + uperp ** 2)
        _, refl, _, g_r = sk.ht_reflect(up, uperp, vm, b)
        ww = w * np.maximum(vm - up / gam, 0) * refl
        idx = np.searchsorted(edges, g_r, side="right") - 1
        ok = idx >= 0
        mc += np.bincount(idx[ok], weights=ww[ok], minlength=len(mc))[:len(mc)] / 8
    rel = abs(flux.sum() - mc.sum()) / mc.sum()
    rel0 = abs(flux[0, 0] - mc[0]) / mc[0]
    check("quadrature vs Monte Carlo", rel < 0.03 and rel0 < 0.03,
          f"gamma >= 1.5: {flux.sum():.4e} vs {mc.sum():.4e} (rel {rel:.1e}); first bin rel {rel0:.1e}")

    # 5 a fitted sphere gives V_HT = V / cos(theta) and theta_Bn = theta
    r0, vs = 20.0, 0.2
    th = np.repeat(np.radians(np.arange(2.5, 90, 5.0)), 30)
    tt = np.tile(np.linspace(100, 600, 30), 18)
    rad = r0 + vs * tt + rng.normal(0, 1.0, th.size)
    front, rms, _ = sk.fit_front(th, tt, rad)
    tq = np.radians(np.array([10.0, 40.0, 70.0]))
    f = sk.foot(front, tq, 350.0)
    e_v = np.max(np.abs(f["v_ht"] * np.cos(tq) / vs - 1))
    e_t = np.max(np.abs(np.degrees(f["tbn"] - tq)))
    check("front fit", e_v < 0.03 and e_t < 1.0 and abs(rms - 1.0) < 0.1,
          f"V_HT cos(theta)/V - 1 = {e_v:.1e}, theta_Bn - theta = {e_t:.2f} deg, rms {rms:.2f} (noise 1.0)")

    # 6 tracing back recovers synthetic reflections off a sphere
    front = sk.sphere(r0, vs)
    table = sk.FootTable(front, np.arange(0.0, 601.0, 1.0), rho_max=150.0)
    n = 5000
    t_r = rng.uniform(150, 500, n)
    th_r = np.radians(rng.uniform(15, 75, n))
    f = sk.foot(front, th_r, t_r)
    u_perp = np.abs(rng.normal(0, 0.3, n)) + 0.05
    enc, _, u_r, g_r = sk.ht_reflect(rng.normal(0, 0.1, n), u_perp, f["v_ht"], 3.0)
    v_r = u_r / g_r
    t_now = 600.0
    z_now = f["z"] + v_r * (t_now - t_r)
    t_est, v_est, tbn_est = table.trace_back(f["rho"], z_now, v_r, t_now)
    found = np.isfinite(t_est[enc])  # only an electron the foot point caught can have been reflected
    e_t = np.nanmax(np.abs(t_est - t_r)[enc])
    e_a = np.nanmax(np.abs(np.degrees(tbn_est - th_r))[enc])
    e_v = np.nanmax(np.abs(v_est / f["v_ht"] - 1)[enc])
    check("trace back", found.all() and e_t < 0.5 and e_a < 0.5 and e_v < 0.02,
          f"found {int(found.sum())}/{int(enc.sum())}, max errors: time {e_t:.2f}, theta_Bn {e_a:.2f} deg, "
          f"V_HT {e_v:.1e}")
    # an electron that never met the front (ahead of it at all times) has no root
    t_none, _, _ = table.trace_back(np.array([50.0]), np.array([190.0]), np.array([0.01]), t_now)
    check("trace back, no encounter", np.isnan(t_none).all(), f"t_r = {t_none}")

    # 7 the tail of the loaded distribution (measured in the 16-step check of L: 0.019628 at gamma >= 1.05)
    frac = sk.mj_tail_fraction(0.01, 1.05)
    check("Maxwell-Juttner tail", abs(frac / 0.019628 - 1) < 0.005, f"model {frac:.6f}, loaded 0.019628")

    # 8 count model on a sphere: every reflection is counted once (time step and angle step converge)
    kw = dict(temp=0.01, b=3.0, phi=0.0, density=1.0, z_half=1e9, gamma_edges=[1.0, np.inf],
              band_edges_deg=[0, 90])
    c1 = sk.reflection_counts(sk.sphere(20.0, 0.2), np.arange(0.5, 200, 1.0), np.radians(np.arange(0.5, 90, 1.0)),
                              [200.0], **kw)
    c2 = sk.reflection_counts(sk.sphere(20.0, 0.2), np.arange(0.25, 200, 0.5),
                              np.radians(np.arange(0.25, 90, 0.5)), [200.0], **kw)
    rel = abs(c1["reflected"][0] / c2["reflected"][0] - 1)
    inbox = abs(c2["counts"][0, 0, 0] / c2["reflected"][0] - 1)
    check("count model", rel < 0.02 and inbox < 1e-6,
          f"reflected {c1['reflected'][0]:.4e} vs {c2['reflected'][0]:.4e} at half the steps (rel {rel:.1e}); "
          f"in an infinite box all are kept (rel {inbox:.1e})")

    # 9 the front of a synthetic shell with the noise of 2 particles per cell (measured: 0.29 per cell)
    r_of = lambda th, t: (21.0 + 0.2 * t) * (1.0 + 0.08 * np.cos(2 * th))  # noqa: E731
    t_out = 100.0 + 18.5 * np.arange(20)
    ax = 0.5 * np.arange(768)
    res = measure_front(_synthetic_frames(rng, ax, r_of, t_out), ax, ax, (192.0, 192.0), **_SYNTH)
    th = res["theta"]
    true = r_of(th[None, :], t_out[:, None])
    err = res["r_fold"] - true
    spread = np.nanmax(np.abs(res["r_quad"] - true[:, None, :]))
    check("front position", np.isfinite(err).all() and np.abs(err).max() < 1.5 and spread < 2.0,
          f"detected {int(np.isfinite(err).sum())}/{err.size}, error {err.mean():+.2f} +- {err.std():.2f} d_e "
          f"(max {np.abs(err).max():.2f}; the mean over d_i moves a sharp jump out by 0.67), "
          f"single quadrants max {spread:.2f}")
    f_true = sk.foot(lambda a, b: (r_of(a, b), 0.2 * (1 + 0.08 * np.cos(2 * a)) + 0 * b,
                                   -(21.0 + 0.2 * b) * 0.16 * np.sin(2 * a)), th[None, :], t_out[:, None])
    mid = slice(2, -2)
    e_v = np.nanmax(np.abs(res["v_n"][mid] / f_true["v_n"][mid] - 1))
    e_a = np.nanmax(np.abs(np.degrees(res["tbn"][mid] - f_true["tbn"][mid])))
    fit_a = np.max(np.abs(np.degrees(sk.foot(res["front"], th[None, :], t_out[:, None])["tbn"] - f_true["tbn"])))
    check("front speed and normal", e_v < 0.08 and e_a < 6.0 and res["rms"] < 1.0 and fit_a < 2.0,
          f"V_n within {e_v:.1%}, theta_Bn within {e_a:.1f} deg (bins) and {fit_a:.1f} deg (fit, rms "
          f"{res['rms']:.2f} d_e)")
    inside = (true >= _SYNTH["r_lo"]) & (true <= _SYNTH["r_hi"])
    inside[:4] = inside[-4:] = False  # the speed needs five outputs, the domain five in a row
    lost = int((inside & ~res["valid"]).sum())
    e_b = np.nanmax(np.abs(res["b_mean"] - 2.5))
    gap = np.nanmean(res["r_fold"] - res["r_debris"])
    check("valid domain, compression, debris edge", lost == 0 and e_b < 0.15 and abs(gap - 12.0) < 1.5,
          f"{lost} bins lost inside the radial range, b = 2.5 within {e_b:.2f}, front ahead of the debris by "
          f"{gap:.1f} d_e (12)")
    flat = measure_front(_synthetic_frames(rng, ax, None, t_out[:6]), ax, ax, (192.0, 192.0), **_SYNTH)
    n_false = int(np.isfinite(flat["r_quad"]).sum())
    check("no front in noise", n_false == 0, f"{n_false} detections in {flat['r_quad'].size} single-quadrant bins")

    # 10 Boris orbits in uniform fields: work of the recorded field against the recorded energy; the moment
    b0 = np.array([0.0, 0.0, 0.2857])
    us, es, bs, ts = _boris_orbits(rng, np.array([0.0, 0.004, -0.004]), b0, 400, 800, 0.25, 4)
    wc = sk.work_closure(us, es, ts, 0.25, b=bs)
    par = np.nanmedian(np.abs(wc["work_par"] / wc["dgamma"] - 1)[np.abs(wc["dgamma"]) >= 0.1])
    check("work closure", wc["n_used"] > 300 and wc["median_rel_error"] < 0.01 and par < 0.05,
          f"{wc['n_used']} orbits with |d gamma| >= 0.1, median error {wc['median_rel_error']:.1e}; "
          f"the work along B is the gain within {par:.1%} (the other component only drifts)")
    us, es, bs, ts = _boris_orbits(rng, np.zeros(3), b0, 400, 800, 0.25, 4)
    mu = sk.magnetic_moment(us, bs)
    drift = np.max(np.abs(mu[-1] / mu[0] - 1))
    nu = sk.pitch_scatter_rate(us, bs, ts, lags=(4, 16))
    check("magnetic moment", drift < 1e-9 and nu[4][0] < 1e-12,
          f"largest change in a uniform field {drift:.1e}; pitch-angle scattering rate {nu[4][0]:.1e}")

    # 11 temperature of the bulk from a histogram, and the thermal bound of the test B1
    u, w = _sample_mj(rng, 0.012, 4_000_000)
    edges = np.logspace(-6, 1, 225)
    hist, _ = np.histogram(np.sqrt(1 + np.sum(u ** 2, axis=0)) - 1, bins=edges, weights=w * 5.8e6)
    temp, _ = sk.fit_temperature(hist, edges, 0.01)
    e_th = sk.thermal_bound(0.01, 9.06e8)
    check("temperature fit, thermal bound", abs(temp / 0.012 - 1) < 0.01 and abs(e_th - 0.20) < 0.005,
          f"T = {temp:.5f} (0.012); E_th = {e_th:.3f} for T = 0.01 and 9.06e8 electrons (0.20)")

    # 12 the per-electron test: reflected thermal electrons are explained, a faster one is not
    u, w = _sample_mj(rng, 0.01, 400_000)
    u_perp = np.hypot(u[0], u[1])
    _, refl, u_r, g_r = sk.ht_reflect(u[2], u_perp, 0.5, 3.0)
    keep = refl & (np.sqrt(1 + np.sum(u ** 2, axis=0)) - 1 <= e_th)
    expl, e_before, _ = sk.s0_test(u_r[keep], u_perp[keep], 0.5, 3.0, e_th)
    twice = sk.ht_reflect(-u_r[keep], u_perp[keep], 0.5, 3.0)[2]  # reflected a second time
    expl2, _, _ = sk.s0_test(twice, u_perp[keep], 0.5, 3.0, e_th)
    check("S0 test", expl.all() and not expl2.any(),
          f"{int(expl.sum())}/{int(keep.sum())} once-reflected thermal electrons explained, "
          f"{int(expl2.sum())} of the twice-reflected")

    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


# ------------------------------------------------------------------ the measured front
_SYNTH = dict(d_i=4.0, v_a=0.0714286, c_s=0.0456, r_lo=44.0, r_hi=172.0)


def _synthetic_frames(rng, ax, r_of, times, noise=0.29):
    """Meridional frames of a shell behind a front r_of(theta, t) (None: undisturbed plasma)."""
    dx, dz = ax[:, None] - 192.0, ax[None, :] - 192.0
    rad = np.hypot(dx, dz)
    theta = np.arccos(np.clip(np.abs(dz) / np.maximum(rad, 1e-9), 0, 1))
    for t in times:
        if r_of is None:
            n = np.ones(rad.shape)
            shell = cavity = np.zeros(rad.shape, bool)
        else:
            front = r_of(theta, t)
            shell, cavity = (rad <= front) & (rad > front - 8.0), rad <= front - 8.0
            n = np.where(shell, 2.5, np.where(cavity, 0.3, 1.0))
        yield {"t": t, "ni": n + noise * np.sqrt(n) * rng.normal(size=n.shape),
               "ndi": (rad <= (r_of(theta, t) - 12.0 if r_of else 0.0)).astype(float),
               "bmag": np.where(shell, 2.5, np.where(cavity, 0.2, 1.0)) + 0.03 * rng.normal(size=n.shape)}


def measure_front(frames, across, along, centre, d_i, v_a, c_s, r_lo, r_hi, n_front=1.5, n_debris=0.1,
                  behind_di=3.0, ahead_di=1.0) -> dict:
    """Front, debris edge, compression, kinematics and valid domain from meridional frames (protocol Sect. 4).

    frames: dicts with t, ni (ambient ion density over n0), ndi (debris ions, or None), bmag (|B| / B0), on
    the plane (across, along) with `along` the direction of B0. Returns arrays on (time, theta bin)."""
    pm = sk.PolarMap(across, along, centre)
    rows: dict = {k: [] for k in ("t", "r_fold", "r_quad", "r_debris", "r_debris_quad", "b_mean", "b_cell")}
    for fr in frames:
        rows["t"].append(fr["t"])
        r_f = sk.outermost_crossing(pm.fold_mean(fr["ni"]), pm.r, n_front)
        rows["r_fold"].append(r_f)
        rows["r_quad"].append(sk.outermost_crossing(np.moveaxis(pm.quadrants(pm.mean(fr["ni"])), 1, 0), pm.r, n_front))
        if fr.get("ndi") is None:
            rows["r_debris"].append(np.full(pm.n_q, np.nan))
            rows["r_debris_quad"].append(np.full((4, pm.n_q), np.nan))
        else:
            rows["r_debris"].append(sk.outermost_crossing(pm.fold_mean(fr["ndi"]), pm.r, n_debris))
            rows["r_debris_quad"].append(
                sk.outermost_crossing(np.moveaxis(pm.quadrants(pm.mean(fr["ndi"])), 1, 0), pm.r, n_debris))
        layer = (pm.r[:, None] >= r_f - behind_di * d_i) & (pm.r[:, None] <= r_f + ahead_di * d_i)
        for key, binned in (("b_mean", pm.fold_mean(fr["bmag"])), ("b_cell", pm.quadrants(pm.max(fr["bmag"])).max(0))):
            top = np.where(layer & np.isfinite(binned), binned, -np.inf).max(axis=0)
            rows[key].append(np.where(np.isfinite(top), top, np.nan))
    out = {k: np.array(v) for k, v in rows.items()}
    out["theta"] = pm.theta
    out.update(sk.front_kinematics(pm.theta, out["t"], out["r_fold"]))
    out.update(sk.valid_domain(out["r_fold"], out["r_debris"], out["v_n"], out["tbn"], v_a, c_s, d_i, r_lo, r_hi))
    found = np.isfinite(out["r_fold"])
    out["front"], out["rms"] = None, float("nan")
    if found.sum() >= 40:
        tt, th = np.broadcast_arrays(out["t"][:, None], pm.theta[None, :])
        out["front"], out["rms"], _ = sk.fit_front(th[found], tt[found], out["r_fold"][found])
    return out


def _open(run_dir):
    import hdf5plugin  # noqa: F401  (the files are compressed with zstd)

    from pypic.analysis import open_run
    return open_run(str(run_dir))


def _slice_frames(series, b0, with_debris):
    for i, t in enumerate(series.t):
        bmag = np.sqrt(sum(series.read(f"B{c}", i).astype(float) ** 2 for c in "xyz")) / b0
        yield {"t": float(t), "ni": series.read("n:i", i), "bmag": bmag,
               "ndi": series.read("n:di", i) if with_debris else None}


def cmd_front(args) -> int:
    run = _open(args.run)
    pl = run.meta["planner"]
    d_i, b0, v_a, centre = pl["d_i"], pl["b0"], pl["v_a"], pl["center"]
    c_s = sk.sound_speed(v_a, pl["beta"])
    half = 0.5 * float(run.L[2])
    r_lo, r_hi = 0.3 * pl["r_b"], half - 5.0 * d_i
    kw = dict(d_i=d_i, v_a=v_a, c_s=c_s, r_lo=r_lo, r_hi=r_hi)
    mer = run.series("slice_y")
    x, z = mer.plane_axes("n:i")
    res = measure_front(_slice_frames(mer, b0, "n:di" in mer.keys()), x, z, (centre[0], centre[2]), **kw)
    t, th = res["t"], res["theta"]
    deg = np.degrees(th)

    # the equatorial plane: one front per azimuth, and the example's definition (B_z >= 1.25 B0)
    eq = run.series("slice_z")
    xe, ye = eq.plane_axes("n:i")
    pe = sk.PolarMap(xe, ye, (centre[0], centre[1]))
    eq_n, eq_b = [], []
    for i in range(len(eq.t)):
        eq_n.append(sk.outermost_crossing(pe.mean(eq.read("n:i", i)), pe.r, 1.5))
        eq_b.append(sk.outermost_crossing(pe.mean(eq.read("Bz", i)) / b0, pe.r, 1.25))
    eq_n, eq_b = np.array(eq_n), np.array(eq_b)

    spread = np.nanstd(np.where(np.isfinite(res["r_quad"]).sum(axis=1, keepdims=True) >= 3, res["r_quad"], np.nan),
                       axis=1)
    summary = {"run": str(args.run), "commit": run.info.get("git_commit"), "d_i": d_i, "v_a": v_a, "c_s": c_s,
               "r_lo": r_lo, "r_hi": r_hi, "n_outputs": len(t), "fit_rms": res["rms"],
               "fit": res["front"].to_dict() if res["front"] is not None else None,
               "detected_bins": int(np.isfinite(res["r_fold"]).sum()), "valid_bins": int(res["valid"].sum()),
               "bins": int(res["valid"].size),
               "quadrant_spread_median": float(np.nanmedian(spread)) if np.isfinite(spread).any() else None,
               "criteria": {k: int(res[k].sum()) for k in "abcde"}}
    print(f"front of {args.run}: {summary['detected_bins']}/{summary['bins']} (time, theta) bins with a front, "
          f"{summary['valid_bins']} in the valid domain; criteria a-e: {summary['criteria']}")
    print(f"fit rms {res['rms']:.2f} d_e (needs < d_i = {d_i}); spread of the four quadrants, median "
          f"{summary['quadrant_spread_median']} d_e")
    with np.errstate(invalid="ignore"):
        mach = res["v_n"] / sk.fast_speed(v_a, c_s, res["tbn"])
    print("\nper output: time, then per theta bin R [d_e] (* = valid); '.' = no front")
    print("  t      " + " ".join(f"{d:5.1f}" for d in deg))
    for i in range(len(t)):
        cells = [f"{res['r_fold'][i, j]:4.0f}{'*' if res['valid'][i, j] else ' '}"
                 if np.isfinite(res["r_fold"][i, j]) else "   . " for j in range(len(th))]
        print(f"  {t[i]:6.1f} " + " ".join(cells))
    table = []
    print("\nper theta bin, averaged over the valid outputs: n, V_n/c, fast Mach, theta_Bn, V_HT, b (bin, cell), "
          "front - debris [d_e]")
    for j in range(len(th)):
        v = res["valid"][:, j]
        row = {"theta_deg": float(deg[j]), "n_valid": int(v.sum()), "t_first": None, "t_last": None}
        if v.any():
            sub = np.isfinite(res["v_ht"][v, j]) & (res["v_ht"][v, j] < 1)
            row.update(t_first=float(t[v][0]), t_last=float(t[v][-1]), v_n=float(np.mean(res["v_n"][v, j])),
                       mach_fast=float(np.mean(mach[v, j])), tbn_deg=float(np.degrees(np.mean(res["tbn"][v, j]))),
                       v_ht=float(np.median(res["v_ht"][v, j])), subluminal=int(sub.sum()),
                       b_mean=float(np.nanmean(res["b_mean"][v, j])), b_cell=float(np.nanmean(res["b_cell"][v, j])),
                       debris_gap=float(np.nanmean(res["r_fold"][v, j] - res["r_debris"][v, j]))
                       if np.isfinite(res["r_debris"][v, j]).any() else None)
            print(f"  {deg[j]:5.1f}  {row['n_valid']:3d}  {row['v_n']:.4f}  {row['mach_fast']:5.2f}  "
                  f"{row['tbn_deg']:5.1f}  {row['v_ht']:7.3f}  {row['b_mean']:5.2f} {row['b_cell']:5.2f}  "
                  f"{row['debris_gap'] if row['debris_gap'] is None else round(row['debris_gap'], 1)}")
        else:
            print(f"  {deg[j]:5.1f}    0  (no valid output)")
        table.append(row)
    summary["per_theta"] = table
    # equator: the meridional bin next to the equator against the equatorial plane and the B_z definition
    te = eq.t
    eq_row = {"t": [float(a) for a in te], "n_i_mean": np.nanmean(eq_n, axis=1).tolist(),
              "n_i_std": np.nanstd(eq_n, axis=1).tolist(), "bz_mean": np.nanmean(eq_b, axis=1).tolist()}
    summary["equator"] = eq_row
    seen = np.flatnonzero(np.isfinite(res["r_fold"][:, -1]))
    last = int(seen[-1]) if len(seen) else None
    if last is not None:
        print(f"\nequator at t = {t[last]:.0f}: meridional bin {res['r_fold'][last, -1]:.1f}, equatorial plane "
              f"{np.nanmean(eq_n[last]):.1f} +- {np.nanstd(eq_n[last]):.1f} over the azimuth, B_z >= 1.25 B0 at "
              f"{np.nanmean(eq_b[last]):.1f} d_e")
    if res["front"] is not None:
        lum = []
        fine = np.radians(np.arange(0.25, 90, 0.5))
        for ti in t:
            f = sk.foot(res["front"], fine, ti)
            sup = f["v_ht"] >= 1.0
            lum.append(float(np.degrees(f["tbn"][sup].min())) if sup.any() else None)
        summary["luminal_tbn_deg"] = lum
        shown = [a for a in lum if a is not None]
        if shown:
            print(f"luminal angle (foot points with V_HT >= c) from the fit: theta_Bn >= {min(shown):.1f}-"
                  f"{max(shown):.1f} deg over the run")
    out = Path(args.run) / "analysis"
    out.mkdir(exist_ok=True)
    arrays = {k: v for k, v in res.items() if isinstance(v, np.ndarray)}
    np.savez(out / "front.npz", eq_t=te, eq_n=eq_n, eq_b=eq_b, **arrays)
    (out / "front.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {out / 'front.npz'} and front.json")
    return 0


# ------------------------------------------------------------------ cadence of the tracks
def _read_tracks(tr, keys, sel=None):
    """Tracked quantities as float64 arrays (len(keys), n_t, n)."""
    if sel is None:
        return np.stack([np.asarray(tr[k], dtype=float) for k in keys])
    return np.stack([np.asarray(tr.f[k][:, sel], dtype=float) for k in keys])


def cmd_cadence(args) -> int:
    """Work closure of tracks recorded every step, read at every 1, 2, 3 and 4 steps (review R1, change 5).
    Fixes the cadence of tracks_e: the coarsest one at which the criterion N4 (10 %) holds per orbit."""
    run = _open(args.run)
    tr = run.tracks("tracks_e", "e")
    t = np.asarray(tr.t, float)
    dt_step = float(run.dt)
    step = int(round((t[1] - t[0]) / dt_step))
    if step != 1:
        raise SystemExit(f"the tracks are recorded every {step} steps; this check needs every step")
    u = _read_tracks(tr, ("ux", "uy", "uz"))
    e = _read_tracks(tr, ("Ex", "Ey", "Ez"))
    b = _read_tracks(tr, ("Bx", "By", "Bz"))
    b0 = float(run.meta["planner"]["b0"])
    bmax = np.nanmax(np.sqrt(np.sum(b ** 2, axis=0)), axis=0)
    gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    stay = np.isfinite(gam).all(axis=0)
    top = np.sort(bmax[stay])[-max(1, int(stay.sum()) // 100):]
    print(f"{u.shape[2]} tracked electrons, {u.shape[1]} samples to t = {t[-1]:.1f}; {int(stay.sum())} stay in the "
          f"box; largest |B| on a track {np.nanmax(bmax) / b0:.2f} B0 (median of the top 1 %: "
          f"{np.median(top) / b0:.2f})")
    summary = {"run": str(args.run), "n_tracked": int(u.shape[2]), "n_stay": int(stay.sum()), "cadences": {}}
    strong = bmax >= 2.0 * b0  # orbits that met the compressed field
    ref = None
    print("  every  orbits            |d gamma| >=  number  |work/d gamma - 1|: median  90 %  within 10 %  "
          "parallel work against every 1 (of d gamma): median  90 %")
    cadences = (1, 2, 3, 4)
    for c in cadences:
        wc = sk.work_closure(u[:, ::c], e[:, ::c], t[::c], dt_step, b=b[:, ::c], dgamma_min=0.0)
        if ref is None:
            ref = wc
        rows = {}
        for name, subset in (("all", stay), ("|B| >= 2 B0 met", stay & strong)):
            for dg in (0.02, 0.05, 0.1, 0.2):
                ok = subset & np.isfinite(wc["work"]) & (np.abs(wc["dgamma"]) >= dg)
                n = int(ok.sum())
                key = f"{name}, {dg}"
                if n == 0:
                    rows[key] = {"n": 0}
                    continue
                rel = np.abs(wc["work"][ok] / wc["dgamma"][ok] - 1.0)
                par = np.abs(wc["work_par"][ok] - ref["work_par"][ok]) / np.abs(ref["dgamma"][ok])
                r = rows[key] = {
                    "n": n, "median": float(np.median(rel)), "p90": float(np.percentile(rel, 90)),
                    "within_10pct": float(np.mean(rel <= 0.1)), "par_median": float(np.median(par)),
                    "par_p90": float(np.percentile(par, 90)),
                    "par_share": float(np.median(wc["work_par"][ok] / wc["dgamma"][ok]))}
                print(f"  {c:5d}  {name:16s}  {dg:12.2f}  {n:6d}  {r['median']:26.3f}  {r['p90']:.3f}  "
                      f"{r['within_10pct']:10.1%}  {r['par_median']:44.3f}  {r['par_p90']:.3f}")
        # samples per gyration in the upstream field and where the field is strongest
        rows["samples_per_gyration_b0"] = 2 * np.pi / (b0 * c * dt_step)
        rows["samples_per_gyration_bmax"] = 2 * np.pi / (float(np.nanmax(bmax)) * c * dt_step)
        summary["cadences"][str(c)] = rows
    # the criterion N4 is applied per orbit: a cadence passes if 95 % of the orbits with |d gamma| >= 0.1
    # close within 10 %, among all orbits and among those that met the compressed field
    passing = []
    for c in cadences:
        rows = [summary["cadences"][str(c)][f"{name}, 0.1"] for name in ("all", "|B| >= 2 B0 met")]
        if all(r["n"] >= 30 and r["within_10pct"] >= 0.95 for r in rows):
            passing.append(c)
    summary["coarsest_passing"] = max(passing) if passing else None
    print(f"coarsest cadence at which 95 % of the orbits with |d gamma| >= 0.1 close within 10 % (all, and those "
          f"that met |B| >= 2 B0; >= 30 orbits each): {summary['coarsest_passing']}")
    out = Path(args.run) / "analysis"
    out.mkdir(exist_ok=True)
    (out / "cadence.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {out / 'cadence.json'}")
    return 0


# ------------------------------------------------------------------ prediction before the runs
def predict(args) -> int:
    """Count model B2 for a nominal sphere. The numbers are the same for L, M and H: the three runs share
    V_d, v_te, the shell and the box, and the nominal front does not depend on the field strength."""
    v_sh = args.front_fraction * args.v_d
    t_end = args.steps * args.dt
    t_grid = np.arange(0.5 * args.t_step, t_end, args.t_step)
    theta = np.radians(np.arange(0.25, 90.0, 0.5))
    t_dumps = np.arange(args.every, args.steps + 1, args.every) * args.dt
    out = {"model": "sphere", "r0": args.r0, "v_sh": v_sh, "temp": args.temp, "density": args.density,
           "z_half": args.z_half, "gamma_edges": [g if np.isfinite(g) else None for g in GAMMA_EDGES],
           "band_edges_deg": BAND_EDGES_DEG, "t_dumps": t_dumps.tolist(),
           "luminal_angle_deg": float(np.degrees(np.arccos(v_sh))), "cases": []}
    print(f"nominal sphere: R = {args.r0} + {v_sh:.4f} t, luminal angle {out['luminal_angle_deg']:.1f} deg, "
          f"T = {args.temp}, {args.density} macro-electrons per d_e^3")
    for b in args.b:
        res = sk.reflection_counts(sk.sphere(args.r0, v_sh), t_grid, theta, t_dumps, args.temp, b, 0.0,
                                   args.density, args.z_half, GAMMA_EDGES, BAND_EDGES_DEG)
        out["cases"].append({"b": b, "phi": 0.0, "counts": res["counts"].tolist(),
                             "reflected": res["reflected"].tolist(), "encountered": res["encountered"].tolist()})
        peak = int(np.argmax(res["counts"].sum(axis=(1, 2))))
        print(f"\nb = {b}: macro-electrons with gamma >= 1.5 in the box, by foot-point theta_Bn band")
        print("  band      reflected   in box at the last dump, gamma bins 1.5-2 / 2-3 / 3-5 / >5"
              f"      in box at the peak (t = {t_dumps[peak]:.0f})")
        for ib in range(len(BAND_EDGES_DEG) - 1):
            if res["reflected"][ib] < 0.5:
                continue
            last = " / ".join(f"{x:9.0f}" for x in res["counts"][-1, ib])
            print(f"  {BAND_EDGES_DEG[ib]:2d}-{BAND_EDGES_DEG[ib + 1]:2d}  {res['reflected'][ib]:10.0f}   {last}"
                  f"      {res['counts'][peak, ib].sum():10.0f}")
        tot = res["counts"].sum(axis=(1, 2))
        print("  total in box per dump:", " ".join(f"{x:.0f}" for x in tot))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
        print(f"\nwrote {args.out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest", help="checks of the kinematics and the count model on synthetic data")
    p = sub.add_parser("predict", help="reflection baseline for a nominal sphere (before the runs)")
    p.add_argument("--v-d", type=float, default=0.25)
    p.add_argument("--front-fraction", type=float, default=0.87, help="front speed over V_d")
    p.add_argument("--r0", type=float, default=21.12, help="radius of the front at t = 0 (the shell)")
    p.add_argument("--temp", type=float, default=0.01, help="upstream electron temperature, m_e c^2")
    p.add_argument("--b", type=float, nargs="+", default=[2.0, 3.0, 4.0], help="compression B_max / B_0")
    p.add_argument("--density", type=float, default=16.0, help="macro-electrons per d_e^3")
    p.add_argument("--z-half", type=float, default=192.0)
    p.add_argument("--steps", type=int, default=2960)
    p.add_argument("--dt", type=float, default=0.25)
    p.add_argument("--every", type=int, default=148, help="cadence of the dump hot_e")
    p.add_argument("--t-step", type=float, default=1.0)
    p.add_argument("--out", default=None)
    p = sub.add_parser("front", help="front, valid shock domain and compression of a run (protocol Sect. 4)")
    p.add_argument("run", type=Path)
    p = sub.add_parser("cadence", help="work closure of tracks read every 1-4 steps (fixes the track cadence)")
    p.add_argument("run", type=Path)
    args = ap.parse_args(argv)
    if args.cmd == "selftest":
        return selftest()
    if args.cmd == "front":
        return cmd_front(args)
    if args.cmd == "cadence":
        return cmd_cadence(args)
    if args.cmd == "predict":
        return predict(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
