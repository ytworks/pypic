"""Analyse a Langmuir-oscillation run (Landau damping, Dawson oscillation, wave breaking).

python examples/landau_damping/analyze.py output/landau_damping/landau [--movie]
    landau:       omega and damping rate of the k mode vs the kinetic root (gate: gamma 10 %, omega 2 %)
    nonlinear:    damping, then trapping oscillations (bounce period reported)
    dawson:       per-sheet frequency and amplitude vs the exact cold solution (gate: 1 %, 2 %)
    wavebreaking: first sheet crossing time vs arcsin(omega_pe / k v0) (gate: 5 %)
    relativistic: sheet period vs the exact relativistic oscillator (gate: 0.5 %)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import theory as T  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402
from pypic.analysis.spectral import fit_damped  # noqa: E402


def _fit_sine(t, y, w0):
    """Least-squares y = c + a sin(w t) + b cos(w t) with w refined around w0 -> (w, amplitude, c)."""
    from scipy.optimize import minimize_scalar

    def resid(w):
        A = np.stack([np.ones_like(t), np.sin(w * t), np.cos(w * t)], 1)
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        return float(np.sum((A @ coef - y) ** 2)), coef

    r = minimize_scalar(lambda w: resid(w)[0], bounds=(0.8 * w0, 1.2 * w0), method="bounded",
                        options={"xatol": 1e-10})
    _, coef = resid(r.x)
    return float(r.x), float(math.hypot(coef[1], coef[2])), float(coef[0])


def _sheet_frequencies(tr, info, t_max, n_max=600):
    """(x0, frequency, amplitude) of tracked electrons spread evenly over x0 (cold cases, started at
    equilibrium, so x0 = x(t=0)), fitted over t <= t_max."""
    t = tr.t
    sel = t <= t_max
    x = tr["x"][sel]
    L = info["L"]
    x = x[0] + np.unwrap((x - x[0]) * 2 * np.pi / L, axis=0) * L / (2 * np.pi)
    ok = np.nonzero(np.all(np.isfinite(x), axis=0))[0]
    ok = ok[np.argsort(x[0, ok])]
    cols = ok[np.linspace(0, len(ok) - 1, min(n_max, len(ok))).astype(int)]
    out = []
    for j in cols:
        w, a, _ = _fit_sine(t[sel], x[:, j], 1.0)
        out.append((x[0, j], w, a))
    return np.array(out)


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    case = info["case"]
    res = {"case": case}
    md = run.series("modes")
    z = md["Ex"][:, 0].astype(complex)
    t = md.t
    figs = []
    if case in ("landau", "nonlinear"):
        w_th, g_th = info["omega_theory"], info["gamma_theory"]
        if case == "landau":
            f = fit_damped(t, z, omega0=w_th, gamma0=g_th)
            res.update(omega_fit=f["omega"], gamma_fit=f["gamma"], omega_theory=w_th, gamma_theory=g_th,
                       omega_rel_error=(f["omega"] - w_th) / w_th, gamma_rel_error=(f["gamma"] - g_th) / g_th,
                       fit_window=f["t_window"])
            res["passed"] = bool(abs(res["gamma_rel_error"]) <= 0.10 and abs(res["omega_rel_error"]) <= 0.02)
        else:
            f = fit_damped(t, z, omega0=w_th, gamma0=g_th, t_max=15.0)
            amp = np.abs(z)
            mins = [i for i in range(1, len(amp) - 1) if amp[i] < amp[i - 1] and amp[i] < amp[i + 1]]
            env_t = t[mins]
            res.update(initial_gamma_fit=f["gamma"], linear_gamma_theory=g_th,
                       bounce_period_estimate=float(np.median(np.diff(env_t[env_t > 15]))) if len(env_t) > 3 else None,
                       bounce_period_linear_theory=2 * math.pi / math.sqrt(info["alpha"]))
        fig, ax = P.plt.subplots(figsize=(6.5, 4))
        P.style()
        ax.semilogy(t, np.abs(z), color=P.color(0), label=r"$|E_x(k)|$")
        tt = np.linspace(*res.get("fit_window", (0, 15)), 50)
        a0 = np.abs(z[np.searchsorted(t, tt[0])])
        ax.semilogy(tt, a0 * np.exp(g_th * (tt - tt[0])), **P.THEORY, label=rf"theory $\gamma$={g_th:.4f}")
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.legend()
        figs.append(("mode", fig))
    else:
        tr = run.tracks("tracks", "electron")
        if case == "wavebreaking":
            t_th = T.crossing_time(info["k_v0"])
            x = tr["x"]
            L = info["L"]
            ux = tr["ux"]
            ok = np.all(np.isfinite(x), axis=0)
            x = x[:, ok]
            x = x[0] + np.unwrap((x - x[0]) * 2 * np.pi / L, axis=0) * L / (2 * np.pi)
            # sheets cross when the map x0 -> x(t) stops being monotonic; bin-averaged sheets (64 bins in
            # x0, ~60 particles each) remove the particle-noise swaps of nearly coincident particles
            nb = 64
            b = np.minimum((x[0] / L * nb).astype(int), nb - 1)
            cnt = np.bincount(b, minlength=nb)
            xm = np.stack([np.bincount(b, weights=row, minlength=nb) for row in x]) / np.maximum(cnt, 1)
            cross = np.any(np.diff(xm[:, cnt > 0], axis=1) < 0, axis=1)
            idx = np.nonzero(cross)[0]
            t_cross = float(tr.t[idx[0]]) if idx.size else float("nan")
            res.update(t_cross=t_cross, t_cross_theory=t_th, rel_error=(t_cross - t_th) / t_th,
                       max_gamma_end=float(np.nanmax(np.sqrt(1 + ux[-1] ** 2))))
            res["passed"] = bool(abs(res["rel_error"]) <= 0.05)
        else:
            k, v0 = info["k"], info["v0"]
            t_max = 3 * 2 * math.pi
            sh = _sheet_frequencies(tr, info, t_max)
            big = np.abs(np.sin(k * sh[:, 0])) > 0.7  # sheets with near-maximal amplitude
            if case == "dawson":
                w_err = np.abs(sh[big, 1] - 1.0)
                a_th = v0 * np.abs(np.sin(k * sh[big, 0]))
                a_err = np.abs(sh[big, 2] - a_th) / a_th
                res.update(sheets=int(big.sum()), freq_rel_error_median=float(np.median(w_err)),
                           amp_rel_error_median=float(np.median(a_err)))
                res["passed"] = bool(res["freq_rel_error_median"] <= 0.01 and res["amp_rel_error_median"] <= 0.02)
            else:
                u0 = info["u0"]
                # amplitude of each sheet from its own initial proper velocity u0 sin(k x0)
                us = u0 * np.abs(np.sin(k * sh[big, 0]))
                T_th = np.array([T.relativistic_period(T.amplitude_from_u0(u)) for u in us])
                w_th = 2 * math.pi / T_th
                err = np.abs(sh[big, 1] - w_th) / w_th
                res.update(sheets=int(big.sum()), beta_max=T.beta_max(T.amplitude_from_u0(u0)),
                           omega_theory_max_amplitude=float(2 * math.pi / T.relativistic_period(
                               T.amplitude_from_u0(u0))), freq_rel_error_median=float(np.median(err)),
                           nonrelativistic_shift=float(np.median(1 - w_th)))
                res["passed"] = bool(res["freq_rel_error_median"] <= 0.005)
        fig, ax = P.plt.subplots(figsize=(6.5, 4))
        P.style()
        x = tr["x"]
        for j in range(0, min(x.shape[1], 40), 4):
            ax.plot(tr.t, x[:, j], color=P.color(j // 4 % 8), lw=0.8)
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.set_ylabel("x [c/ωpe]")
        ax.set_title("electron sheets")
        figs.append(("sheets", fig))
    s = run.scalars
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    if show_plot:
        res["plots"] = [P.save(f, os.path.join(run_dir, f"{name}.png")) for name, f in figs]
    else:
        for _, f in figs:
            P.plt.close(f)
    return res


def movie(run_dir: str) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title=f"Langmuir oscillation: {run.meta['planner']['case']}")
    m.phase("phase", "electron").xt("lineout", "Ex").spectrum("spectrum", "electron").history()
    return m.save(os.path.join(run_dir, "langmuir.mp4"), fps=20, n_frames=150)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
