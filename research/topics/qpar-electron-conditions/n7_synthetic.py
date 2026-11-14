"""Null test N7: the diffusion measurement (diffusion.py) on synthetic fields.

    python research/topics/qpar-electron-conditions/n7_synthetic.py --json <file> [--n 20000] [--seed 7]

CPU only (NumPy). Three parts; the first is N7 as the README fixes it, the other two cover what it leaves out.

1. ``known``: a uniform field B0 at the angle theta to x, plus a magnetic field without memory (drawn anew for
   every particle and step). The direction of a particle then diffuses on the sphere at a known rate nu, and

       D / D_B = (2/3) [cos^2(theta) Omega / nu + sin^2(theta) (nu / Omega) / (1 + (nu / Omega)^2)]

   in the conventions of diffusion.py (derivation: the velocity correlation is exp(-nu t) along B0 and
   exp(-nu t) cos(Omega t) across it). The time step keeps nu dt <= 0.01, the order of the error of the rate.
   Fixed before the run (README, N7): the measured D is within 10 % of this for D between 0.1 and 10 D_B, and it is
   taken from an interval of a factor 3 in time with a log-slope within 1 +- 0.1.
2. ``streaming``: the same uniform field without the random part. The motion is ballistic; the rule of the
   diffusive interval must return no value.
3. ``fields``: the path that part 1 does not use, fields that vary along x. (a) A field linear in x, cut and
   mirrored: the interpolated value against the exact one, for a component at the nodes and one at the half cells.
   (b) A static field (0, By(x), Bz(x)): u_y + q a_z(x) and u_z - q a_y(x) with a_z' = Bz, a_y' = By are constants
   of the motion. Their largest change over five gyrations is given for a set of time steps, and once more for a
   push that reads the two components half a cell away from where they are (``wrong_stagger``), which shows what
   the number is sensitive to. No threshold is fixed for part 3.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import diffusion as df  # noqa: E402

U0, B0 = 10.0, 0.0667  # p = m_i v_pt and the upstream field of V1 (m_e c, m_e c omega_pe / e)
RATES = (1 / 15, 2 / 9, 2 / 3, 20 / 9, 20 / 3)  # nu / Omega: D = 10, 3, 1, 0.3, 0.1 D_B at theta = 0
ANGLES = (0.0, 30.0)


def exact(ratio: float, theta: float) -> float:
    c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    return (2.0 / 3.0) * (c * c / ratio + s * s * ratio / (1.0 + ratio * ratio))


def last_slope(t, msd, factor: float = 3.0) -> float:
    sel = t >= t[-1] / factor
    return float(np.polyfit(np.log(t[sel]), np.log(msd[sel]), 1)[0])


def known(n: int, seed: int) -> list:
    g = np.sqrt(1.0 + U0 * U0)
    omega, d_b = B0 / g, df.bohm(U0, B0)
    rows = []
    for theta in ANGLES:
        c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
        for k, ratio in enumerate(RATES):
            rng = np.random.default_rng([seed, int(theta), k])
            nu = ratio * omega
            dt = min(0.05, 0.01 / ratio) / omega
            steps = int(round(100.0 / (nu * dt)))
            win = df.uniform((B0 * c, B0 * s, 0.0))
            x, u = rng.random(n) * win.length, df.isotropic(n, U0, rng)
            t, mean, msd = df.push(win, x, u, dt, steps, df.sample_steps(steps), noise=g * np.sqrt(nu / dt), rng=rng)
            iv = df.diffusive_interval(t, msd)
            row = {"theta_deg": theta, "nu_over_omega": ratio, "nu_dt": nu * dt, "omega_dt": omega * dt,
                   "steps": steps, "exact": exact(ratio, theta), "interval": iv,
                   "last_log_slope": last_slope(t, msd), "mean_over_rms": float(mean[-1] / np.sqrt(msd[-1])),
                   "u_change_max": float(np.abs(np.sqrt((u * u).sum(axis=0)) - U0).max() / U0)}
            if iv:
                row["measured"] = iv["d"] / d_b
                row["measured_line"] = iv["d_line"] / d_b
                row["deviation"] = row["measured"] / row["exact"] - 1.0
            row["in_scope"] = bool(0.1 * (1 - 1e-9) <= row["exact"] <= 10.0 * (1 + 1e-9))
            row["pass"] = bool(iv and abs(row["deviation"]) <= 0.10)
            rows.append(row)
    return rows


def streaming(n: int, seed: int) -> dict:
    rng = np.random.default_rng([seed, 99])
    g = np.sqrt(1.0 + U0 * U0)
    dt, steps = 0.05 * g / B0, 20000
    win = df.uniform((B0 * np.cos(np.radians(30.0)), B0 * np.sin(np.radians(30.0)), 0.0))
    x, u = rng.random(n) * win.length, df.isotropic(n, U0, rng)
    t, _, msd = df.push(win, x, u, dt, steps, df.sample_steps(steps))
    iv = df.diffusive_interval(t, msd)
    return {"steps": steps, "last_log_slope": last_slope(t, msd), "interval": iv, "pass": iv is None}


def integral(win: df.Window, name: str, x):
    """Integral from 0 to x of the interpolated component ``name`` (piecewise linear, periodic)."""
    f = win.f[name]
    nxt = np.roll(f, -1)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f + nxt))]) * win.dx
    s = x / win.dx - win.stagger[name]
    i = np.floor(s)
    w = s - i
    i = i.astype(np.int64)
    k = i % win.n
    # the knots sit at (i + stagger) dx: the integral from 0 to the first knot is added by the caller's difference
    return (i // win.n) * cum[-1] + cum[k] + win.dx * (f[k] * w + 0.5 * (nxt[k] - f[k]) * w * w)


def fields(seed: int) -> dict:
    dx, cells = 0.1, 24000
    xs = np.arange(cells + 1) * dx
    # (a) a field linear in x, cut to [i0, i1) and mirrored
    i0, i1 = 3000, 21000
    lin = {"Bx": 0.3 + 0.002 * xs, "By": 0.3 + 0.002 * (xs + 0.5 * dx)}
    win = df.mirrored(lin, i0, i1, dx)
    half = (i1 - i0) * dx
    rng = np.random.default_rng([seed, 5])
    x = rng.random(4000) * 2.0 * half
    lab = np.where(x < half, x, 2.0 * half - x) + i0 * dx  # position in the lineout that the mirror maps to
    err = {}
    for name in lin:
        inner = (np.abs(x - half) > dx) & (x > dx) & (x < 2.0 * half - dx)  # away from the two kinks
        err[name] = float(np.abs(win.at(name, x) - (0.3 + 0.002 * lab))[inner].max())
    # (b) constants of the motion in a static field (0, By(x), Bz(x))
    length = cells * dx
    xh = (np.arange(cells + 1) + 0.5) * dx
    by = B0 * (0.6 + 0.5 * np.sin(2 * np.pi * 3 * xh / length + 1.0))
    bz = B0 * (0.8 * np.cos(2 * np.pi * 5 * xh / length) + 0.3 * np.sin(2 * np.pi * 11 * xh / length))
    win = df.mirrored({"By": by, "Bz": bz}, i0, i1, dx)
    g = np.sqrt(1.0 + U0 * U0)
    period = 2.0 * np.pi * g / B0
    wrong = df.Window(win.f, dx, stagger=df.STAGGER_X | {"By": 0.0, "Bz": 0.0})
    out = {"linear_field_max_error": err, "omega_dt_per_dt": B0 / g, "invariants": []}
    for dt, pushed in ((3.2, win), (1.6, win), (0.8, win), (0.4, win), (0.2, win), (0.1, win), (0.4, wrong)):
        rng = np.random.default_rng([seed, 6])
        n = 2000
        x, u = rng.random(n) * win.length, df.isotropic(n, U0, rng)

        def consts(x, u, dt=dt):
            xm = x - 0.5 * dt * u[0] / g  # the position at the time of u (half a step before x)
            return u[1] - integral(win, "Bz", xm), u[2] + integral(win, "By", xm)  # q = -1

        before = consts(x, u)
        df.push(pushed, x, u, dt, int(round(5.0 * period / dt)))
        after = consts(x, u)
        out["invariants"].append({"dt": dt, "wrong_stagger": pushed is wrong,
                                  "max_change_over_u": [float(np.abs(a - b).max() / U0)
                                                        for a, b in zip(after, before)]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json")
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    def git(*c):
        return subprocess.run(("git", "-C", here) + c, capture_output=True, text=True).stdout.strip()

    t0 = time.time()
    out = {"commit": git("rev-parse", "HEAD"), "clean_tree": git("status", "--short") == "", "n": a.n,
           "seed": a.seed, "u0": U0, "b0": B0, "statistical_error": float(np.sqrt(2.0 / a.n)),
           "known": known(a.n, a.seed), "streaming": streaming(a.n, a.seed), "fields": fields(a.seed)}
    scope = [r for r in out["known"] if r["in_scope"]]
    out["n7_pass"] = bool(all(r["pass"] for r in scope) and out["streaming"]["pass"])
    out["seconds"] = time.time() - t0
    text = json.dumps(out, indent=1)
    print(text)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
