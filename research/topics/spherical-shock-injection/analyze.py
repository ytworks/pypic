"""Analysis of the exploratory runs of the topic spherical-shock-injection (CPU only).

Usage (from the repository root):
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py selftest
    .venv/bin/python research/topics/spherical-shock-injection/analyze.py predict --out <json>

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

    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


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
    args = ap.parse_args(argv)
    if args.cmd == "selftest":
        return selftest()
    if args.cmd == "predict":
        return predict(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
