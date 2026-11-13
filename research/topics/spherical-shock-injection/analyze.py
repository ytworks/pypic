"""Analysis of the exploratory runs of the topic spherical-shock-injection (CPU only).

Usage (from the repository root):
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py selftest
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py predict --out <json>
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py front <run directory>
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py cadence <run directory>
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py null <run directory>

Results of a command on a run go to <run directory>/analysis/ (not in git).

The definitions and criteria are fixed in exploratory-protocol.md; the kinematics is in shockkin.py.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import baseline  # noqa: E402
import campaign  # noqa: E402
import census  # noqa: E402
import disk  # noqa: E402
import drive  # noqa: E402
import orbits  # noqa: E402
import potential  # noqa: E402
import r2checks  # noqa: E402
import record  # noqa: E402
import routes  # noqa: E402
import shockkin as sk  # noqa: E402
import sites  # noqa: E402
import upstream  # noqa: E402

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

    def plane(times=t_out, front=r_of, **kw):
        return measure_plane(_synthetic_frames(rng, ax, front, times, **kw), ax, ax, (192.0, 192.0), _SYNTH["d_i"])

    first = plane()
    res = measure_front([first, plane()], **_SYNTH)
    th = res["theta"]
    true = r_of(th[None, :], t_out[:, None])
    err = res["r"] - true
    spread = np.nanmax(np.abs(res["r_quad"] - true[None, :, None, :]))
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
    pull = float(np.sqrt(np.mean((res["v_n"][mid] - f_true["v_n"][mid]) ** 2)) / np.mean(res["sigma_v"][mid]))
    check("front speed and normal", e_v < 0.08 and e_a < 6.0 and res["rms"] < 1.0 and fit_a < 2.0
          and 0.5 < pull < 2.0,
          f"V_n within {e_v:.1%}, theta_Bn within {e_a:.1f} deg (bins) and {fit_a:.1f} deg (fit, rms "
          f"{res['rms']:.2f} d_e); rms error of V_n over the mean slope error {pull:.2f}")
    inside = np.zeros(true.shape, bool)
    inside[2:-2] = True  # the speed needs five outputs
    inside &= (t_out >= 112.0)[:, None]  # the age of the shock, 2 / Omega_ci
    lost = int((inside & ~res["valid"]).sum())
    e_b = np.nanmax(np.abs(res["b_mean"] - 2.5))
    gap = np.nanmean(res["r"] - res["r_debris"])
    check("valid domain, compression, debris edge", lost == 0 and res["n_row"] == 6 and e_b < 0.15
          and abs(gap - 16.0) < 1.5,
          f"{lost} bins lost of {int(inside.sum())} that can be valid ({res['n_row']} outputs in a row), b = 2.5 "
          f"within {e_b:.2f}, front ahead of the debris by {gap:.1f} d_e (16)")
    gap_b, flow = float(np.nanmedian(res["gap_b"])), float(np.nanmedian(res["flow"][inside]))
    check("mid-ramp surface, ion flow", abs(gap_b - 3.0) < 1.0 and res["rms_B"] < 1.0 and abs(flow - 1.0) < 0.1,
          f"density front ahead of the mid-ramp surface by {gap_b:.2f} d_e (3; fit rms {res['rms_B']:.2f}); ion "
          f"velocity over V_n (1 - n0/n): median {flow:.3f}")
    sub = inside & (f_true["v_ht"] < 0.7)
    shares = {k: float((res[k] & sub).sum() / sub.sum()) for k in ("ok_planes", "ok_thresholds", "ok_local",
                                                                    "ok_quadrants", "usable")}
    check("agreement flags on a clean front", shares["ok_planes"] == 1.0 and shares["ok_thresholds"] == 1.0
          and shares["usable"] > 0.9,
          "share of the valid bins with V_HT < 0.7 that pass: " + ", ".join(f"{k} {v:.2f}" for k, v in shares.items()))
    other = measure_front([first, plane(scale=1.15)], **_SYNTH)
    soft = measure_front([plane(ramp=20.0)], **_SYNTH)
    near = measure_front([plane(debris_gap=5.0)], **_SYNTH)
    bad = (int(other["ok_planes"].sum()), int(soft["ok_thresholds"].sum()), int(near["b"].sum()),
           int(other["usable"].sum() + soft["usable"].sum() + near["usable"].sum()))
    check("fronts that must not be used", bad == (0, 0, 0, 0),
          f"bins passed: {bad[0]} with planes 15 % apart, {bad[1]} with a ramp of 20 d_e (thresholds), {bad[2]} "
          f"with the debris 5 d_e behind the front; usable in the three cases {bad[3]}")
    flat = plane(times=t_out[:6], front=None)
    n_false = int(np.isfinite(flat["r_quad"]).sum() + np.isfinite(flat["r_thr"]).sum())
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

    # 13 canonical radius: constant around a gyro-orbit, and it follows a displaced guiding centre
    n = 20000
    b0 = 0.2857
    gx, gy = rng.uniform(-150, 150, (2, n))
    gx[np.hypot(gx, gy) < 10] += 20.0
    u_perp = np.abs(rng.normal(0.0, 0.5, n))
    r_l = u_perp / b0
    ph = rng.uniform(0, 2 * np.pi, n)  # an electron in B along +z: r = r_gc + r_L e, u = u_perp (-sin, cos)
    x, y = gx + r_l * np.cos(ph), gy + r_l * np.sin(ph)
    rho_p = sk.canonical_radius(x, y, -u_perp * np.sin(ph), u_perp * np.cos(ph), b0)
    want = np.sqrt(gx ** 2 + gy ** 2 - r_l ** 2)
    gc = sk.guiding_centre(x, y, -u_perp * np.sin(ph), u_perp * np.cos(ph), b0)
    err_gc = float(np.max(np.hypot(gc[0] - gx, gc[1] - gy)))
    err = float(np.nanmax(np.abs(rho_p - want)))
    spread = float(np.max(np.abs(np.hypot(x, y) - np.hypot(gx, gy))))
    check("canonical radius", err < 1e-9 and err_gc < 1e-9 and spread > 3.0,
          f"rho_P = sqrt(rho_gc^2 - rho_L^2) within {err:.1e} d_e while the position moves by up to {spread:.1f}")

    # 14 birth cell from the id, and the origin classes
    n_cells, ppc = (768, 768, 768), 2
    cells = rng.integers(0, 768, (3, 1000))
    ids = ppc * ((cells[0] * 768 + cells[1]) * 768 + cells[2]) + rng.integers(0, ppc, 1000)
    back = sk.birth_cell(ids.astype(np.uint64), n_cells, ppc)
    rad = sk.birth_radius(ids, n_cells, (0.5, 0.5, 0.5), ppc, (192.0, 192.0))
    want = np.hypot((cells[0] + 0.5) * 0.5 - 192.0, (cells[1] + 0.5) * 0.5 - 192.0)
    org = _origin(np.array([0, 2 * 768 ** 3 - 1, 2 * 768 ** 3, 20 * 768 ** 3], dtype=np.uint64), n_cells, ppc)
    check("birth cell", all((a == b).all() for a, b in zip(back, cells)) and np.allclose(rad, want)
          and org.tolist() == [0, 0, 2, 1], "1000 ids decoded; origins ambient / face-injected / other")

    census.selfcheck(check)
    orbits.selfcheck(check)
    baseline.selfcheck(check)
    upstream.selfcheck(check)
    potential.selfcheck(check)
    record.selfcheck(check)
    r2checks.selfcheck(check)
    drive.selfcheck(check)
    sites.selfcheck(check)
    th = np.radians(np.array([90.0, 0.0, 45.0]))
    cos_d, tbn_d = sk.dipole_tbn(th, np.full(3, 100.0), np.zeros(3))
    check("dipole theta_Bn on a sphere", np.allclose(np.degrees(tbn_d), [90.0, 0.0, np.degrees(np.arctan(0.5))]),
          f"{np.degrees(tbn_d).round(3).tolist()} at polar angles 90, 0, 45 (90, 0, 26.565)")
    bx_, by_, bz_ = sk.dipole_field(np.array([100.0, 0.0]), np.zeros(2), np.array([0.0, 100.0]), 1e6)
    check("dipole field, equator and pole", np.allclose([bz_[0], bz_[1]], [-1.0, 2.0]) and np.allclose(bx_, 0.0),
          f"B_z {bz_.round(6).tolist()} at r = 100 on the equator and the pole for m = 1e6 (-1, 2)")
    campaign.selfcheck(check)
    routes.selfcheck(check)
    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


# ------------------------------------------------------------------ the measured front
N_FRONT = (1.3, 1.5, 1.7)  # thresholds of the ambient ion density; the front is the middle one
N_DEBRIS = (0.03, 0.1, 0.3)  # thresholds of the debris density; the edge is the middle one
_SYNTH = dict(d_i=4.0, v_a=0.0714286, c_s=0.0456, omega_ci=1.0 / 56.0, r_outer=184.0)


def _synthetic_frames(rng, ax, r_of, times, noise=0.29, scale=1.0, ramp=0.0, debris_gap=16.0):
    """Frames on a plane through the centre of a shell behind a front scale * r_of(theta, t) (None: undisturbed).

    The density jumps to 2.5 at the front (or rises over `ramp` d_e behind it), the field 3 d_e behind it; the
    shell is 11 d_e thick, with a cavity behind; the debris ends debris_gap d_e behind the front; the ions of
    the shell move outwards with V (1 - 1/2.5), V = 0.2 (1 + 0.08 cos 2 theta)."""
    dx, dz = ax[:, None] - 192.0, ax[None, :] - 192.0
    rad = np.maximum(np.hypot(dx, dz), 1e-9)
    theta = np.arccos(np.clip(np.abs(dz) / rad, 0, 1))
    for t in times:
        if r_of is None:
            n, bmag, deb, u_r = np.ones(rad.shape), np.ones(rad.shape), np.zeros(rad.shape), np.zeros(rad.shape)
        else:
            front = scale * r_of(theta, t)
            shell, cavity = (rad <= front) & (rad > front - 11.0), rad <= front - 11.0
            rise = np.clip((front - rad) / ramp, 0.0, 1.0) if ramp else 1.0
            n = np.where(shell, 1.0 + 1.5 * rise, np.where(cavity, 0.3, 1.0))
            bmag = np.where(shell & (rad <= front - 3.0), 2.5, np.where(cavity, 0.2, 1.0))
            deb = (rad <= front - debris_gap).astype(float)
            u_r = np.where(shell, 0.2 * (1 + 0.08 * np.cos(2 * theta)) * 0.6, 0.0)
        ni = n + noise * np.sqrt(n) * rng.normal(size=n.shape)
        yield {"t": t, "ni": ni, "ndi": deb, "nde": None, "bmag": bmag + 0.03 * rng.normal(size=n.shape),
               "fa": ni * u_r * dx / rad, "fb": ni * u_r * dz / rad}


def measure_plane(frames, across, along, centre, d_i, behind_di=3.0, ahead_di=1.0) -> dict:
    """Front, debris edge, compression, mid-ramp surface and ion flow on one plane through the centre.

    frames: dicts with t, ni (ambient ion density over n0), ndi and nde (debris ions and electrons, or None),
    bmag (|B| / B0), fa and fb (ambient ion flux along the two axes), on the plane (across, along) with `along`
    the direction of B0. Returns arrays with the time as the first axis and the theta bin as the last:
    r_thr (per threshold of N_FRONT), r_quad (four quadrants, at N_FRONT[1]), r_debris (per threshold of
    N_DEBRIS), b_mean and r_bmax (largest bin average of |B|/B0 in the layer, and where), b_cell (largest
    single cell), r_b (outermost radius within 3 d_i of the front where |B|/B0 averaged over d_i reaches
    (1 + b_mean) / 2), n_peak and u_r (largest bin density within 3 d_i behind the front; the radial ion
    velocity in that bin)."""
    pm = sk.PolarMap(across, along, centre)
    da = np.asarray(across, float)[:, None] - centre[0]
    db = np.asarray(along, float)[None, :] - centre[1]
    rad = np.maximum(np.hypot(da, db), 1e-9)
    rows: dict = {k: [] for k in ("t", "r_thr", "r_quad", "r_debris", "b_mean", "r_bmax", "b_cell", "r_b", "n_peak",
                                  "u_r")}
    for fr in frames:
        rows["t"].append(fr["t"])
        n_fold = pm.fold_mean(fr["ni"])
        r_thr = np.array([sk.outermost_crossing(n_fold, pm.r, thr) for thr in N_FRONT])
        r_f = r_thr[1]
        rows["r_thr"].append(r_thr)
        rows["r_quad"].append(
            sk.outermost_crossing(np.moveaxis(pm.quadrants(pm.mean(fr["ni"])), 1, 0), pm.r, N_FRONT[1]))
        edges = np.full((len(N_DEBRIS), pm.n_q), np.nan)
        for key in ("ndi", "nde"):
            if fr.get(key) is not None:
                folded = pm.fold_mean(fr[key])
                edges = np.fmax(edges, np.array([sk.outermost_crossing(folded, pm.r, thr) for thr in N_DEBRIS]))
        rows["r_debris"].append(edges)
        b_fold = pm.fold_mean(fr["bmag"])
        rr = pm.r[:, None]
        layer = (rr >= r_f - behind_di * d_i) & (rr <= r_f + ahead_di * d_i)
        masked = np.where(layer & np.isfinite(b_fold), b_fold, -np.inf)
        top = masked.max(axis=0)
        b_mean = np.where(np.isfinite(top), top, np.nan)
        rows["b_mean"].append(b_mean)
        rows["r_bmax"].append(np.where(np.isfinite(top), pm.r[masked.argmax(axis=0)], np.nan))
        cell = np.where(layer, pm.quadrants(pm.max(fr["bmag"])).max(axis=0), -np.inf).max(axis=0)
        rows["b_cell"].append(np.where(np.isfinite(cell), cell, np.nan))
        ramp = np.where(rr <= r_f + 3.0 * d_i, b_fold - 0.5 * (1.0 + b_mean), -np.inf)
        r_b = sk.outermost_crossing(ramp, pm.r, 0.0)
        rows["r_b"].append(np.where(r_b >= r_f - 3.0 * d_i, r_b, np.nan))
        f_fold = pm.fold_mean((fr["fa"] * da + fr["fb"] * db) / rad)
        dens = np.where((rr >= r_f - behind_di * d_i) & (rr <= r_f) & np.isfinite(n_fold), n_fold, -np.inf)
        k = dens.argmax(axis=0)[None, :]
        peak = np.take_along_axis(n_fold, k, axis=0)[0]
        seen = np.isfinite(dens.max(axis=0))
        rows["n_peak"].append(np.where(seen, peak, np.nan))
        with np.errstate(invalid="ignore", divide="ignore"):
            rows["u_r"].append(np.where(seen, np.take_along_axis(f_fold, k, axis=0)[0] / peak, np.nan))
    out = {k: np.array(v) for k, v in rows.items()}
    out["theta"] = pm.theta
    return out


def _fit(theta, t, radius, mask=None, n_min=40):
    """Fit of the finite (and masked) bins of radius (n_t, n_theta); (None, nan) with fewer than n_min bins."""
    use = np.isfinite(radius) if mask is None else np.isfinite(radius) & mask
    if use.sum() < n_min:
        return None, float("nan")
    tt, th = np.broadcast_arrays(t[:, None], theta[None, :])
    front, rms, _ = sk.fit_front(th[use], tt[use], radius[use])
    return front, rms


def _fit_vht(front, theta, t):
    """V_HT of a fit on the bins; NaN without a fit and outside its time range."""
    if front is None:
        return np.full((len(t), len(theta)), np.nan)
    t0, t1 = front.t_range
    v = sk.foot(front, theta[None, :], t[:, None])["v_ht"]
    return np.where(((t >= t0) & (t <= t1))[:, None], v, np.nan)


def measure_front(planes, d_i, v_a, c_s, omega_ci, r_outer, dipole=None) -> dict:
    """The front of protocol Sect. 4 from one or two meridional planes (results of measure_plane, same outputs).

    The front is the mean of the planes where all of them show one. Returns arrays on (time, theta bin):
    r, the kinematics from differences (sk.front_kinematics), the valid domain (sk.valid_domain, with c_s a
    number or one value per output), and the agreement flags
      ok_planes      the planes agree within d_i (true with one plane);
      ok_thresholds  R moves by less than d_i and V_HT of the fits by less than 10 % between 1.3 and 1.7 n0;
      ok_local       V_HT of the fit and from differences of neighbouring bins agree within 15 %;
      ok_quadrants   V_HT of the fits of the four quadrants of each plane agree within 15 %;
      usable         valid and all flags and the rms residual of the fit below d_i: the bins used for B1'.
    Also the fits front (density, 1.5 n0), front_B (mid-ramp surface over the bins with b >= 1.5) and their rms,
    gap_b (density front minus mid-ramp surface) and flow (ion velocity along the normal behind the front over
    V_n (1 - n0 / n))."""
    t, th = planes[0]["t"], planes[0]["theta"]
    stack = {k: np.array([p[k] for p in planes]) for k in planes[0] if k not in ("t", "theta")}
    out: dict = {"t": t, "theta": th, "n_planes": len(planes)}
    r_all = stack["r_thr"]  # (plane, time, threshold, theta)
    r_thr = np.where(np.isfinite(r_all).all(axis=0), r_all.mean(axis=0), np.nan)
    out["r_thr"], out["r"] = r_thr, r_thr[:, 1]
    out["r_plane"] = r_all[:, :, 1]
    out["r_quad"] = stack["r_quad"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # bins without a front in any plane stay NaN
        for k in ("b_mean", "r_bmax", "r_b", "n_peak", "u_r"):
            out[k] = np.nanmean(stack[k], axis=0)
        out["b_cell"] = np.nanmax(stack["b_cell"], axis=0)
        out["r_debris_thr"] = np.nanmax(stack["r_debris"], axis=0)
    out["r_debris"] = out["r_debris_thr"][:, 1]
    out.update(sk.front_kinematics(th, t, out["r"]))
    omega_row = None
    if dipole is not None:
        # a dipole along z: theta_Bn from the local field direction, v_A and Omega_ci from |B| at the front
        out["cos_tbn"], out["tbn"] = sk.dipole_tbn(th[None, :], out["r"], out["r_theta"])
        with np.errstate(invalid="ignore", divide="ignore"):
            out["v_ht"] = out["v_n"] / out["cos_tbn"]
            bmag = dipole["moment"] / out["r"] ** 3 * np.sqrt(1.0 + 3.0 * np.cos(th)[None, :] ** 2)
        out["b_local"] = bmag
        v_a = bmag / np.sqrt(dipole["m_i"])
        omega_ci = bmag / dipole["m_i"]
        omega_row = dipole["moment"] / dipole["r_ref"] ** 3 / dipole["m_i"]
    fits = [_fit(th, t, r_thr[:, k]) for k in range(len(N_FRONT))]
    out["front"], out["rms"] = fits[1]
    out["fits_thr"] = [f for f, _ in fits]
    vht = np.array([_fit_vht(f, th, t) for f, _ in fits])
    out["v_ht_fit"] = vht[1]
    with np.errstate(invalid="ignore"):
        out["ok_planes"] = np.isfinite(out["r"]) & (np.ptp(r_all[:, :, 1], axis=0) < d_i)
        out["ok_thresholds"] = (np.ptp(r_thr, axis=1) < d_i) & sk.vht_agree(vht[1], vht[[0, 2]], 0.10)
    out["ok_local"] = sk.vht_agree(vht[1], out["v_ht"][None], 0.15)
    ok_q = np.ones(out["r"].shape, bool)
    for p in planes:
        quad = np.array([_fit_vht(_fit(th, t, p["r_quad"][:, q])[0], th, t) for q in range(4)])
        ok_q &= sk.vht_agree(vht[1], quad, 0.15)
    out["ok_quadrants"] = ok_q
    out["front_B"], out["rms_B"] = _fit(th, t, out["r_b"], out["b_mean"] >= 1.5)
    out["gap_b"] = out["r"] - out["r_b"]
    out.update(sk.valid_domain(out["r"], out["r_debris"], out["r_bmax"], out["v_n"], out["sigma_v"], out["tbn"], t,
                               v_a, c_s, omega_ci, r_outer, omega_row=omega_row))
    with np.errstate(invalid="ignore", divide="ignore"):
        # a flow along the normal has the radial part u_n R / sqrt(R^2 + R_theta^2)
        u_n = out["u_r"] * np.sqrt(out["r"] ** 2 + out["r_theta"] ** 2) / out["r"]
        out["flow"] = u_n / (out["v_n"] * (1.0 - 1.0 / out["n_peak"]))
    out["usable"] = (out["valid"] & out["ok_planes"] & out["ok_thresholds"] & out["ok_local"] & out["ok_quadrants"]
                     & bool(out["rms"] < d_i))
    return out


def _open(run_dir):
    import hdf5plugin  # noqa: F401  (the files are compressed with zstd)

    from pypic.analysis import open_run
    return open_run(str(run_dir))


def _slice_frames(series, b0, comps):
    """Frames of a slice series for measure_plane; comps: the two in-plane axes, e.g. "xz". b0 is a number (the
    uniform field) or an array of the plane (the static field of a dipole run): |B| is divided by it."""
    keys = series.keys()
    for i, t in enumerate(series.t):
        bmag = np.sqrt(sum(series.read(f"B{c}", i).astype(float) ** 2 for c in "xyz")) / b0
        yield {"t": float(t), "ni": series.read("n:i", i), "bmag": bmag,
               "ndi": series.read("n:di", i) if "n:di" in keys else None,
               "nde": series.read("n:de", i) if "n:de" in keys else None,
               "fa": series.read(f"flux{comps[0]}:i", i), "fb": series.read(f"flux{comps[1]}:i", i)}


def upstream_temperature(run, species=("e", "i")) -> dict:
    """Temperatures of the far upstream plasma against time, from the spectra of the corner column (spec_up).

    The ion spectrum is in the same energy unit as the electron one if its first fit gives the initial
    temperature; that is checked by the caller on the first output (mean energy 1.5 T)."""
    s = run.series("spec_up")
    edges = np.asarray(s.edges[0], float)
    cen = np.sqrt(edges[1:] * edges[:-1])
    out: dict = {"t": np.asarray(s.t, float)}
    for sp in species:
        temps, means = [], []
        for i in range(len(s.t)):
            h = np.asarray(s.read(sp, i), float)
            means.append(float((h * cen).sum() / h.sum()))
            temp = means[-1] / 1.5  # first guess from the mean energy
            for _ in range(3):
                temp = float(sk.fit_temperature(h, edges, temp)[0])
            temps.append(temp)
        out[f"T_{sp}"], out[f"mean_{sp}"] = np.array(temps), np.array(means)
    return out


def _count(mask) -> int:
    return int(np.asarray(mask).sum())


def cmd_front(args) -> int:
    run = _open(args.run)
    pl = run.meta["planner"]
    d_i, b0, v_a, centre = pl["d_i"], pl["b0"], pl["v_a"], pl["center"]
    m_i = (b0 / v_a) ** 2  # n0 = 1: v_a = B0 / sqrt(m_i)
    omega_ci = b0 / m_i
    half = 0.5 * float(run.L[2])
    r_outer = half - 8.0
    names = [n for n in ("slice_y", "slice_x") if _has_series(run, n)]
    dips = run.config.get("external_fields", {}).get("dipoles") or []
    dipole = None
    if dips:
        dipole = {"moment": float(dips[0]["moment"][2]), "m_i": m_i, "r_ref": float(args.r_ref)}
    planes = []
    for name in names:
        ser = run.series(name)
        comps = "xz" if name == "slice_y" else "yz"
        a, z = ser.plane_axes("n:i")
        norm = b0
        if dipole is not None:
            aa, zz = np.meshgrid(np.asarray(a, float) - centre[0 if name == "slice_y" else 1],
                                 np.asarray(z, float) - centre[2], indexing="ij")
            bx, by, bz = sk.dipole_field(aa, 0.0 * aa, zz, dipole["moment"], b0)
            norm = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
        planes.append(measure_plane(_slice_frames(ser, norm, comps), a, z,
                                    (centre[0 if name == "slice_y" else 1], centre[2]), d_i))
    t = planes[0]["t"]
    c_s0 = sk.sound_speed(v_a, pl["beta"])
    up = None
    if _has_series(run, "spec_up"):
        up = upstream_temperature(run)
        # the ion spectrum may be in units of the ion rest energy: bring it to m_e c^2 by the first output
        t_i = up["T_i"] * (m_i if up["T_i"][0] < 0.3 * up["T_e"][0] else 1.0)
        c_s = np.interp(t, up["t"], sk.sound_speed_of(up["T_e"], t_i, m_i))
        up["T_i_used"] = t_i
    else:
        c_s = np.full(len(t), c_s0)
    res = measure_front(planes, d_i, v_a, c_s, omega_ci, r_outer, dipole=dipole)
    th = res["theta"]
    deg = np.degrees(th)

    # the equatorial plane: one front per azimuth, and the example's definition (B_z >= 1.25 B0)
    eq = run.series("slice_z")
    xe, ye = eq.plane_axes("n:i")
    pe = sk.PolarMap(xe, ye, (centre[0], centre[1]))
    eq_n, eq_b = [], []
    for i in range(len(eq.t)):
        eq_n.append(sk.outermost_crossing(pe.mean(eq.read("n:i", i)), pe.r, N_FRONT[1]))
        if dipole is None:
            eq_b.append(sk.outermost_crossing(pe.mean(eq.read("Bz", i)) / b0, pe.r, 1.25))
        else:  # |B_z| over the static |B_z| of the dipole (plus B0) in the equatorial plane
            static = np.abs(-dipole["moment"] / np.maximum(pe.r, 1e-9) ** 3 + b0)
            ring = pe.mean(np.abs(eq.read("Bz", i).astype(float)))
            eq_b.append(sk.outermost_crossing(ring / static.reshape((-1,) + (1,) * (ring.ndim - 1)), pe.r, 1.25))
    eq_n, eq_b = np.array(eq_n), np.array(eq_b)

    flags = ("ok_planes", "ok_thresholds", "ok_local", "ok_quadrants")
    v = res["valid"]
    summary = {"run": str(args.run), "commit": run.info.get("git_commit"), "planes": names, "d_i": d_i, "v_a": v_a,
               "omega_ci_inv": 1.0 / omega_ci, "r_outer": r_outer, "n_row": res["n_row"],
               "c_s_initial": c_s0, "c_s_first": float(c_s[0]), "c_s_last": float(c_s[-1]),
               "upstream": None if up is None else {k: np.asarray(a).tolist() for k, a in up.items()},
               "n_outputs": len(t), "fit_rms": res["rms"], "fit_rms_B": res["rms_B"],
               "fit": res["front"].to_dict() if res["front"] is not None else None,
               "fit_B": res["front_B"].to_dict() if res["front_B"] is not None else None,
               "fits_thresholds": {str(n): f.to_dict() if f is not None else None
                                   for n, f in zip(N_FRONT, res["fits_thr"])},
               "bins": int(v.size), "detected_bins": _count(np.isfinite(res["r"])), "valid_bins": _count(v),
               "usable_bins": _count(res["usable"]),
               "criteria": {k: _count(res[k]) for k in "abcde"},
               "flags_in_valid": {k: _count(res[k] & v) for k in flags}}
    print(f"front of {args.run} from {' and '.join(names)}: {summary['detected_bins']}/{summary['bins']} (time, "
          f"theta) bins with a front; criteria a-e {summary['criteria']} ({res['n_row']} outputs in a row); "
          f"{summary['valid_bins']} valid, {summary['usable_bins']} usable for B1'")
    print(f"flags inside the valid domain: {summary['flags_in_valid']}; fit rms {res['rms']:.2f} d_e (density) and "
          f"{res['rms_B']:.2f} d_e (mid-ramp), needed < d_i = {d_i}")
    if up is not None:
        print(f"upstream (corner column): T_e {up['T_e'][0]:.4f} -> {up['T_e'][-1]:.4f}, T_i {up['T_i_used'][0]:.4f}"
              f" -> {up['T_i_used'][-1]:.4f}; c_s {c_s[0]:.4f} -> {c_s[-1]:.4f} (initial {c_s0:.4f})")
    with np.errstate(invalid="ignore"):
        mach = res["v_n"] / res["v_f"]
    print("\nper output: time, then per theta bin R [d_e] (* = usable, + = valid only); '.' = no front")
    print("  t      " + " ".join(f"{d:5.1f}" for d in deg))
    for i in range(len(t)):
        cells = [f"{res['r'][i, j]:4.0f}{'*' if res['usable'][i, j] else '+' if v[i, j] else ' '}"
                 if np.isfinite(res["r"][i, j]) else "   . " for j in range(len(th))]
        print(f"  {t[i]:6.1f} " + " ".join(cells))
    table = []
    print("\nper theta bin over the valid outputs: n valid, n usable, V_n/c, sigma_V, fast Mach, theta_Bn, V_HT "
          "(differences; fit), b (bin, cell), front - debris, front - mid-ramp [d_e], ion flow ratio")

    def med(key, sel, j):
        vals = res[key][sel, j]
        return float(np.nanmedian(vals)) if np.isfinite(vals).any() else None

    def show(x, fmt):
        return "    -" if x is None else format(x, fmt)

    for j in range(len(th)):
        sel = v[:, j]
        row = {"theta_deg": float(deg[j]), "n_valid": _count(sel), "n_usable": _count(res["usable"][:, j]),
               **{k: _count(res[k][sel, j]) for k in flags}}
        if sel.any():
            row.update(t_first=float(t[sel][0]), t_last=float(t[sel][-1]), v_n=med("v_n", sel, j),
                       sigma_v=med("sigma_v", sel, j), mach_fast=float(np.nanmedian(mach[sel, j])),
                       tbn_deg=float(np.degrees(np.nanmedian(res["tbn"][sel, j]))), v_ht=med("v_ht", sel, j),
                       v_ht_fit=med("v_ht_fit", sel, j), b_mean=med("b_mean", sel, j), b_cell=med("b_cell", sel, j),
                       debris_gap=float(np.nanmedian(res["r"][sel, j] - res["r_debris"][sel, j]))
                       if np.isfinite(res["r_debris"][sel, j]).any() else None,
                       gap_b=med("gap_b", sel, j), flow=med("flow", sel, j))
            print(f"  {deg[j]:5.1f} {row['n_valid']:3d} {row['n_usable']:3d}  {row['v_n']:.4f} {row['sigma_v']:.4f} "
                  f"{row['mach_fast']:5.2f}  {row['tbn_deg']:5.1f}  {row['v_ht']:7.3f} {show(row['v_ht_fit'], '7.3f')}"
                  f"  {row['b_mean']:5.2f} {row['b_cell']:5.2f}  {show(row['debris_gap'], '5.1f')} "
                  f"{show(row['gap_b'], '5.1f')}  {show(row['flow'], '5.2f')}")
        else:
            print(f"  {deg[j]:5.1f}   0   0  (no valid output)")
        table.append(row)
    summary["per_theta"] = table
    with np.errstate(invalid="ignore"):
        ok_flow = np.abs(res["flow"] - 1.0) <= 0.3
    summary["flow_within_30pct"] = {"valid": _count(ok_flow & v), "of": _count(v)}
    edge = res["r_debris_thr"]
    if np.isfinite(edge).any():
        spread = float(np.nanmedian(edge[:, 2] - edge[:, 0]))
        summary["debris_edge_0.3_minus_0.03_median"] = spread
        print(f"\ndebris edge: 0.3 n0 against 0.03 n0, median difference {spread:+.1f} d_e; ion flow within 30 % of "
              f"V_n (1 - n0/n) in {summary['flow_within_30pct']['valid']}/{_count(v)} valid bins")
    # equator: the meridional bin next to the equator against the equatorial plane and the B_z definition
    te = eq.t
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # outputs without a front stay NaN
        summary["equator"] = {"t": [float(a) for a in te], "n_i_mean": np.nanmean(eq_n, axis=1).tolist(),
                              "n_i_std": np.nanstd(eq_n, axis=1).tolist(),
                              "bz_mean": np.nanmean(eq_b, axis=1).tolist()}
    seen = np.flatnonzero(np.isfinite(res["r"][:, -1]))
    if len(seen):
        last = int(seen[-1])
        print(f"\nequator at t = {t[last]:.0f}: meridional bin {res['r'][last, -1]:.1f}, equatorial plane "
              f"{np.nanmean(eq_n[last]):.1f} +- {np.nanstd(eq_n[last]):.1f} over the azimuth, B_z >= 1.25 B0 at "
              f"{np.nanmean(eq_b[last]):.1f} d_e")
    if res["front"] is not None:
        lum = []
        fine = np.radians(np.arange(0.25, 90, 0.5))
        for ti in t:
            f = sk.foot(res["front"], fine, ti)
            sup = f["v_ht"] >= 1.0
            inside = res["front"].t_range[0] <= ti <= res["front"].t_range[1]
            lum.append(float(np.degrees(f["tbn"][sup].min())) if sup.any() and inside else None)
        summary["luminal_tbn_deg"] = lum
        shown = [a for a in lum if a is not None]
        if shown:
            print(f"luminal angle (foot points with V_HT >= c) from the fit: theta_Bn >= {min(shown):.1f}-"
                  f"{max(shown):.1f} deg over the run")
    out = Path(args.run) / "analysis"
    out.mkdir(exist_ok=True)
    arrays = {k: a for k, a in res.items() if isinstance(a, np.ndarray)}
    np.savez(out / "front.npz", eq_t=te, eq_n=eq_n, eq_b=eq_b, c_s=c_s, **arrays)
    (out / "front.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {out / 'front.npz'} and front.json")
    return 0


def _has_series(run, name) -> bool:
    try:
        run.series(name)
    except (KeyError, FileNotFoundError, OSError):
        return False
    return True


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


# ------------------------------------------------------------------ numerical floor (null run)
SPEED_EDGES = [0.0, 0.1, 0.2, 0.3, 1.0]  # speed bins of the scattering rate (protocol N3)
LAGS = (4.0, 16.0, 64.0)  # lags in time units (protocol N3)


def _origin(ids, n_cells, ppc):
    """0 ambient (loaded at the start), 1 injected at a face (ids from 20 N), 2 anything else."""
    n = int(np.prod(n_cells))
    ids = np.asarray(ids, dtype=np.uint64)
    return np.where(ids < np.uint64(ppc * n), 0, np.where(ids >= np.uint64(20 * n), 1, 2))


def _thin(rows, every=4):
    """Every `every`-th row and the last one."""
    return rows[::every] + ([] if (len(rows) - 1) % every == 0 else rows[-1:])


def _percentiles(a, qs=(50, 90, 99, 99.9)):
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return {"n": 0}
    out = {"n": int(len(a)), "max": float(a.max())}
    out.update({f"p{q}": float(v) for q, v in zip(qs, np.percentile(a, qs))})
    return out


def floor_counts(run, thresholds=(1.25, 1.5, 2.0)):
    """N1: macro-electrons of the species e in the dump hot_e above each threshold, by origin and dump."""
    cfg = run.config
    half = np.asarray(run.L, float) / 2
    p = run.particles("hot_e")
    rows = []
    for it, tt in zip(p.iterations, p.times):
        d = p.get(it, "e")
        gam = np.sqrt(1.0 + d["ux"].astype(float) ** 2 + d["uy"].astype(float) ** 2 + d["uz"].astype(float) ** 2)
        org = _origin(d["id"], cfg["n_cells"], cfg["ppc"])
        to_face = np.min([half[i] - np.abs(d[c] - half[i]) for i, c in enumerate("xyz")], axis=0)
        row = {"t": float(tt), "n_dumped": int(len(gam)), "gamma_max": float(gam.max()) if len(gam) else None}
        for thr in thresholds:
            for name, code in (("ambient", 0), ("injected", 1), ("other", 2)):
                row[f"{name}>={thr}"] = int(np.sum((gam >= thr) & (org == code)))
        sel = gam >= 1.5
        row["to_face_median>=1.5"] = float(np.median(to_face[sel])) if sel.any() else None
        rows.append(row)
    return rows


def floor_temperature(run, t0):
    """N2: electron temperature against time in the corner column (spec_up), and its map in the meridional
    slab (spec_meridional) by distance from the centre."""
    out = {}
    s = run.series("spec_up")
    edges = np.asarray(s.edges[0], float)
    cen = np.sqrt(edges[1:] * edges[:-1])
    rows = []
    for i, tt in enumerate(s.t):
        h = np.asarray(s.read("e", i), float)
        temp = t0
        for _ in range(3):
            temp = float(sk.fit_temperature(h, edges, temp)[0])
        mean_e = float((h * cen).sum() / h.sum())
        n_tot = float(h.sum())
        tail = {f"E>={e}": float(h[cen >= e].sum()) for e in (0.1, 0.2, 0.5)}
        expect = {f"mj E>={e}": float(sk.mj_tail_fraction(temp, 1.0 + e) * n_tot) for e in (0.1, 0.2, 0.5)}
        rows.append({"t": float(tt), "T_fit": temp, "mean_E": mean_e, "weight": n_tot, **tail, **expect})
    out["corner"] = rows
    s = run.series("spec_meridional")
    ex, ez, ee = (np.asarray(a, float) for a in s.edges)
    cx, cz = 0.5 * (ex[1:] + ex[:-1]), 0.5 * (ez[1:] + ez[:-1])
    centre = run.meta["planner"]["center"]
    rr = np.hypot(cx[:, None] - centre[0], cz[None, :] - centre[2])
    inside = (np.abs(cx[:, None] - centre[0]) < centre[0] - 16) & (np.abs(cz[None, :] - centre[2]) < centre[2] - 16)
    rings = [(0, 50), (50, 100), (100, 150), (150, 176)]
    rows = []
    for i, tt in enumerate(s.t):
        h = np.asarray(s.read("e", i), float)
        guess = t0  # one fitting range for the whole map: the temperature of the summed spectrum
        for _ in range(3):
            guess = float(sk.fit_temperature(h[inside].sum(axis=0), ee, guess)[0])
        temp = sk.fit_temperature(h, ee, guess)[0]
        row = {"t": float(tt), "T_all": guess}
        for lo, hi in rings:
            sel = inside & (rr >= lo) & (rr < hi) & np.isfinite(temp)
            row[f"T {lo}-{hi}"] = {"median": float(np.median(temp[sel])), "p5": float(np.percentile(temp[sel], 5)),
                                   "p95": float(np.percentile(temp[sel], 95))}
        rows.append(row)
    out["meridional"] = rows
    return out


def floor_tracks(run, chunk=40000):
    """N3 and N5 from the tracked electrons of a run without a driver: pitch-angle scattering rate and
    energy diffusion per speed bin, and the canonical radius against the birth radius and against time."""
    cfg, pl = run.config, run.meta["planner"]
    b0, centre = float(pl["b0"]), pl["center"]
    tr = run.tracks("tracks_e", "e")
    t = np.asarray(tr.t, float)
    dts = float(t[1] - t[0])
    lags = [int(round(lag / dts)) for lag in LAGS]
    ids = np.asarray(tr.ids)
    n = len(ids)
    nb = len(SPEED_EDGES) - 1
    acc = {k: np.zeros((nb, len(lags))) for k in ("mu_local", "mu_axis", "pairs_local", "pairs_axis", "dg2", "dg",
                                                   "pairs_g")}
    n_bin = np.zeros(nb, int)
    every = max(1, int(round(18.5 / dts)))  # the times of the particle dumps
    rows_t = np.arange(0, len(t), every)
    rho_p = np.full((len(rows_t), n), np.nan)
    stay = np.zeros(n, bool)
    dgam_run = np.full(n, np.nan)
    for a in range(0, n, chunk):
        sel = slice(a, min(n, a + chunk))
        u = _read_tracks(tr, ("ux", "uy", "uz"), sel)
        b = _read_tracks(tr, ("Bx", "By", "Bz"), sel)
        gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
        stay[sel] = np.isfinite(gam).all(axis=0)
        dgam_run[sel] = gam[-1] - gam[0]
        speed = np.nanmedian(np.sqrt(np.sum(u ** 2, axis=0)) / gam, axis=0)
        which = np.digitize(speed, SPEED_EDGES) - 1
        axis = np.zeros_like(b)
        axis[2] = 1.0
        for ib in range(nb):
            m = which == ib
            n_bin[ib] += int(m.sum())
            if not m.any():
                continue
            for name, field in (("local", b[:, :, m]), ("axis", axis[:, :, m])):
                res = sk.pitch_scatter_rate(u[:, :, m], field, t, lags=lags)
                for il, lag in enumerate(lags):
                    nu, pairs = res.get(lag, (np.nan, 0))
                    if pairs:
                        acc[f"mu_{name}"][ib, il] += nu * pairs
                        acc[f"pairs_{name}"][ib, il] += pairs
            for il, lag in enumerate(lags):
                dg = gam[lag:, m] - gam[:-lag, m]
                ok = np.isfinite(dg)
                acc["dg2"][ib, il] += float(np.sum(dg[ok] ** 2))
                acc["dg"][ib, il] += float(np.sum(dg[ok]))
                acc["pairs_g"][ib, il] += int(ok.sum())
        xy = _read_tracks(tr, ("x", "y"), sel)[:, rows_t]
        rho_p[:, sel] = sk.canonical_radius(xy[0] - centre[0], xy[1] - centre[1], u[0, rows_t], u[1, rows_t], b0)
        if a == 0:  # the birth cell from the id against the position at t = 0
            bx, by, _ = sk.birth_cell(ids[sel], cfg["n_cells"], cfg["ppc"])
            off = max(float(np.nanmax(np.abs(xy[0, 0] - (bx + 0.5) * cfg["dx"][0]))),
                      float(np.nanmax(np.abs(xy[1, 0] - (by + 0.5) * cfg["dx"][1]))))
    out = {"n_tracked": int(n), "n_stay": int(stay.sum()), "sample_step": dts, "lags_time": list(LAGS),
           "birth_cell_offset_max": off, "speed_edges": SPEED_EDGES, "n_per_speed_bin": n_bin.tolist()}
    with np.errstate(divide="ignore", invalid="ignore"):
        out["nu_local"] = (acc["mu_local"] / acc["pairs_local"]).tolist()
        out["nu_axis"] = (acc["mu_axis"] / acc["pairs_axis"]).tolist()
        out["energy_diffusion"] = (acc["dg2"] / acc["pairs_g"] / np.asarray(LAGS)[None, :]).tolist()
        out["heating_rate"] = (acc["dg"] / acc["pairs_g"] / np.asarray(LAGS)[None, :]).tolist()
    out["dgamma_run"] = _percentiles(np.abs(dgam_run))
    out["n_dgamma>=0.1"] = int(np.sum(np.abs(dgam_run) >= 0.1))
    rho_b = sk.birth_radius(ids, cfg["n_cells"], cfg["dx"], cfg["ppc"], centre)
    ranges = {"all": rho_b >= 4.0, "4-30": (rho_b >= 4) & (rho_b < 30), "30-100": (rho_b >= 30) & (rho_b < 100),
              ">=100": rho_b >= 100}
    rows = []
    for j, it in enumerate(rows_t):
        row = {"t": float(t[it])}
        for name, m in ranges.items():
            row[f"birth {name}"] = _percentiles(np.abs(rho_p[j, m] - rho_b[m]))
        row["start all"] = _percentiles(np.abs(rho_p[j] - rho_p[0])[ranges["all"]])
        row["mean signed"] = float(np.nanmean((rho_p[j] - rho_b)[ranges["all"]]))
        rows.append(row)
    out["delta_rho"] = rows
    return out


def cmd_null(args) -> int:
    """Numerical floor of a run without a driver (protocol Sect. 5): N1 counts of energetic electrons,
    N2 heating, N3 pitch-angle scattering, N5 the canonical radius against the birth radius."""
    run = _open(args.run)
    t0 = float(run.config["species"][0]["beams"][0]["theta"])
    summary = {"run": str(args.run), "T0": t0}
    rows = floor_counts(run)
    summary["N1"] = rows
    print("N1 macro-electrons of the species e in hot_e, by origin (ambient / injected at a face / other)")
    print("      t   dumped   gamma>=1.25         gamma>=1.5          gamma>=2     gamma max   to the face (>=1.5)")
    for r in rows:
        def tri(thr, r=r):
            return " / ".join(f"{r[f'{k}>={thr}']:5d}" for k in ("ambient", "injected", "other"))
        face = "-" if r["to_face_median>=1.5"] is None else f"{r['to_face_median>=1.5']:.1f}"
        gmax = "-" if r["gamma_max"] is None else f"{r['gamma_max']:.3f}"
        print(f"  {r['t']:5.0f}  {r['n_dumped']:7d}  {tri(1.25)}  {tri(1.5)}  {tri(2.0)}  {gmax:>8s}  {face:>8s}")
    worst = {k: max(r[f"{k}>=1.5"] for r in rows) for k in ("ambient", "injected", "other")}
    summary["N1_max_per_dump>=1.5"] = worst
    print(f"largest count per dump at gamma >= 1.5: {worst}")

    temp = floor_temperature(run, t0)
    summary["N2"] = temp
    print(f"\nN2 electron temperature (T0 = {t0:.4f}); corner column spec_up")
    print("      t    T fit   2/3 <E>   electrons above E = 0.1 / 0.2 / 0.5 (Maxwell-Juttner at the fitted T)")
    for r in _thin(temp["corner"]):
        tails = " / ".join(f"{r[f'E>={e}']:.3g} ({r[f'mj E>={e}']:.3g})" for e in (0.1, 0.2, 0.5))
        print(f"  {r['t']:5.0f}  {r['T_fit']:.5f}  {2 / 3 * r['mean_E']:.5f}   {tails}")
    print("  meridional slab spec_meridional: median T (5 %-95 % of the map) by distance from the centre")
    for r in temp["meridional"]:
        cells = "  ".join(f"{k[2:]}: {v['median']:.4f} ({v['p5']:.4f}-{v['p95']:.4f})"
                          for k, v in r.items() if k.startswith("T "))
        print(f"  {r['t']:5.0f}  {cells}")

    trk = floor_tracks(run)
    summary["tracks"] = trk
    print(f"\nN3 {trk['n_tracked']} tracked electrons ({trk['n_stay']} stay in the box), samples every "
          f"{trk['sample_step']} time units; birth cell from the id against the position at t = 0: largest "
          f"offset {trk['birth_cell_offset_max']:.3f} d_e")
    print(f"  pitch-angle scattering rate nu at lags {LAGS} (against the local field | against the axis), "
          "energy diffusion <(d gamma)^2>/lag and heating <d gamma>/lag")
    for ib in range(len(SPEED_EDGES) - 1):
        loc = " ".join(f"{v:.2e}" for v in trk["nu_local"][ib])
        ax = " ".join(f"{v:.2e}" for v in trk["nu_axis"][ib])
        de = " ".join(f"{v:.2e}" for v in trk["energy_diffusion"][ib])
        he = " ".join(f"{v:+.2e}" for v in trk["heating_rate"][ib])
        print(f"  v {SPEED_EDGES[ib]:.1f}-{SPEED_EDGES[ib + 1]:.1f} ({trk['n_per_speed_bin'][ib]:6d}): nu {loc} | {ax}"
              f"   diffusion {de}   heating {he}")
    print(f"  |gamma(end) - gamma(0)|: {trk['dgamma_run']}; orbits with 0.1 or more: {trk['n_dgamma>=0.1']}")
    print("\nN5 |rho_P - rho_birth| of the tracked electrons born at rho >= 4 d_e: median / 99 % / 99.9 % / max, "
          "then |rho_P - rho_P(0)|: 99.9 % / max, and the mean of rho_P - rho_birth")
    for r in _thin(trk["delta_rho"]):
        a, s0 = r["birth all"], r["start all"]
        print(f"  {r['t']:5.0f}  {a['p50']:.3f} / {a['p99']:.3f} / {a['p99.9']:.3f} / {a['max']:.3f}    "
              f"{s0['p99.9']:.3f} / {s0['max']:.3f}   {r['mean signed']:+.4f}   (n = {a['n']})")
    out = Path(args.run) / "analysis"
    out.mkdir(exist_ok=True)
    (out / "null.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {out / 'null.json'}")
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
    p.add_argument("--r-ref", type=float, default=70.0, help="dipole runs: radius [d_e] of the equatorial ion "
                   "gyro-frequency that sets the persistence of 2 / Omega_ci (dipole-design.md)")
    p = sub.add_parser("cadence", help="work closure of tracks read every 1-4 steps (fixes the track cadence)")
    p.add_argument("run", type=Path)
    p = sub.add_parser("census", help="energetic electrons by origin, class, band and result (protocol Sect. 6)")
    p.add_argument("run", type=Path)
    p.add_argument("--null", type=Path, help="null.json of the run without a driver (floors N1, N3, N5)")
    p.add_argument("--tau", type=float, help="tau_dwell that replaces 3 gyro-periods (the rule of protocol "
                   "Sect. 3.2); the result goes to analysis/census-tau.json and .npz, census.json is kept")
    p = sub.add_parser("campaign", help="campaign A: one pass over the tracks (entries, gamma_max, frame, sites)")
    p.add_argument("run", type=Path)
    p.add_argument("--chunk", type=int, default=60000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--summary-only", action="store_true", help="rewrite campaign.json from the saved tables")
    p = sub.add_parser("disk", help="dipole runs: histories and mechanism of the equatorial disk (disk_e dump)")
    p.add_argument("run", type=Path)
    p = sub.add_parser("routes", help="M6 of global-plan.md: routes of the top 1 % electrons (after campaign)")
    p.add_argument("run", type=Path)
    p = sub.add_parser("prune", help="delete tracks and particle dumps of an analysed run (campaign A)")
    p.add_argument("run", type=Path)
    p.add_argument("--null", action="store_true", help="a null run: counts.json is not required")
    p = sub.add_parser("sites", help="site boundaries along each ray from the |B| profile (sites.md); after front")
    p.add_argument("run", type=Path)
    p = sub.add_parser("sitegain", help="net change of gamma of the tracked electrons per site and theta_Bn")
    p.add_argument("run", type=Path)
    p.add_argument("--chunk", type=int, default=100000)
    p = sub.add_parser("counts", help="electrons per dump by class and energy, without back-tracing (large runs)")
    p.add_argument("run", type=Path)
    p = sub.add_parser("n3c", help="fluctuation power upstream of the front against the null run (N3c)")
    p.add_argument("run")
    p.add_argument("null")
    p = sub.add_parser("phi", help="rise of the potential along the field lines through the layer (Sect. 3.1)")
    p.add_argument("run")
    p.add_argument("--null", help="directory of the null run: its E along the same paths gives the noise")
    p = sub.add_parser("baseline", help="count model B2 on the measured front against the census (Sect. 3.4)")
    p.add_argument("run")
    p.add_argument("--workers", type=int, default=4)
    p = sub.add_parser("tracks", help="measures on the tracked orbits (protocol Sect. 7, N4, f_FN); after census")
    p.add_argument("run")
    p.add_argument("--chunk", type=int, default=20000, help="orbits measured per block")
    p = sub.add_parser("null", help="numerical floor N1, N2, N3, N5 of a run without a driver")
    p.add_argument("run", type=Path)
    p = sub.add_parser("r2null", help="review R2: retention of the magnetic moment in a run without a driver")
    p.add_argument("run")
    p.add_argument("--gaps", type=float, nargs="+", default=[150.0, 250.0, 420.0],
                   help="times between the two windows, besides 0 and 3 gyro-periods")
    p.add_argument("--starts", type=float, nargs="+", default=[0.0, 100.0, 200.0], help="start of the first window")
    p.add_argument("--chunk", type=int, default=40000)
    p = sub.add_parser("r2", help="review R2: pre-valid bins, regime, epsilon, ramp width, thermal nu tau, ids")
    p.add_argument("run")
    p.add_argument("--null", required=True, help="directory of the null run (its analysis/null.json)")
    p = sub.add_parser("drive", help="how long the debris drives the front (shock-forming-run.md); after front")
    p.add_argument("run")
    p = sub.add_parser("record", help="write runs/<date>-<commit>-<name>.json and runs/index.json (clean tree)")
    p.add_argument("run")
    p.add_argument("name")
    args = ap.parse_args(argv)
    if args.cmd == "record":
        return record.cmd_record(args)
    if args.cmd == "drive":
        return drive.cmd_drive(args)
    if args.cmd == "r2null":
        return r2checks.cmd_r2null(args)
    if args.cmd == "r2":
        return r2checks.cmd_r2(args)
    if args.cmd == "n3c":
        return upstream.cmd_n3c(args)
    if args.cmd == "phi":
        return potential.cmd_phi(args)
    if args.cmd == "baseline":
        return baseline.cmd_baseline(args)
    if args.cmd == "tracks":
        return orbits.cmd_tracks(args)
    if args.cmd == "census":
        return census.cmd_census(args)
    if args.cmd == "campaign":
        return campaign.write_summary(args.run) if args.summary_only else campaign.cmd_campaign(args)
    if args.cmd == "disk":
        return disk.cmd_disk(args)
    if args.cmd == "routes":
        return routes.cmd_routes(args)
    if args.cmd == "prune":
        return campaign.cmd_prune(args)
    if args.cmd == "sites":
        return sites.cmd_sites(args)
    if args.cmd == "sitegain":
        return sites.cmd_sitegain(args)
    if args.cmd == "counts":
        return census.cmd_counts(args)
    if args.cmd == "selftest":
        return selftest()
    if args.cmd == "front":
        return cmd_front(args)
    if args.cmd == "cadence":
        return cmd_cadence(args)
    if args.cmd == "null":
        return cmd_null(args)
    if args.cmd == "predict":
        return predict(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
