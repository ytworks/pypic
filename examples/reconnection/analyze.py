"""Analyse a relativistic reconnection run.

python examples/reconnection/analyze.py output/reconnection/pair_sigma10 [--movie]
gate: before the tearing onset (t <= 10 / omega_pe) the y-profile of B_x stays within 5 % of B0 of the
      analytic double Harris sheet (equilibrium), and Gauss's law holds.
reported (empirical comparisons): peak reconnection rate d(flux)/dt / (B0 v_A) vs ~0.1, magnetic energy
      converted, power-law index of dN/dgamma vs Werner et al. 2016, E_parallel share of the work done on
      the most energetic tracked electrons.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402
from pypic.analysis.fields import flux_function  # noqa: E402
from pypic.analysis.particles import fit_power_law, parallel_perpendicular_work, spectrum  # noqa: E402
from pypic.profiles import evaluate  # noqa: E402


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    B0 = info["B0"]
    y1, y2 = info["y_sheets"]
    pr = run.series("profile_y")
    y = pr.x("Bx")
    bx_th = B0 * evaluate({"shape": "double_tanh", "centers": [y1, y2], "width": info["delta"]}, y)
    early = np.nonzero(pr.t <= 10.0)[0]
    dev = max(float(np.max(np.abs(pr.read("Bx", i) - bx_th))) for i in early) / B0 if early.size else float("nan")
    # reconnected flux of the lower sheet: along y = y1, flux function maximum minus minimum
    sl = run.series("slice_xy")
    xa, ya = sl.plane_axes("Bx")
    j1 = int(np.argmin(np.abs(ya - y1)))
    dxs, dys = xa[1] - xa[0], ya[1] - ya[0]
    psi = []
    for i in range(len(sl.t)):
        A = flux_function(sl.read("Bx", i), sl.read("By", i), dxs, dys)
        line = A[:, j1]
        psi.append(float(line.max() - line.min()))
    psi = np.asarray(psi)
    rate = np.gradient(psi, sl.t) / (B0 * info["v_A"])
    s = run.scalars
    eb = s["magnetic"]
    sp = run.series("spectrum")
    G, f = spectrum(sp.edges[0], sum(sp.read(k, -1) for k in sp.keys()))
    sig = info["sigma"]
    pl = fit_power_law(G - 1.0, f, lo=0.5 * sig, hi=5.0 * sig)  # the non-thermal range of sigma-scaling studies
    active = np.isfinite(rate) & (rate > 0.02)  # times when flux is being reconnected
    res = {"equilibrium_max_dev_over_B0": dev, "reconnection_rate_peak": float(np.nanmax(rate)),
           "reconnection_rate_mean_active": float(np.mean(rate[active])) if active.any() else float("nan"),
           "reconnection_rate_literature": 0.1, "magnetic_energy_converted": float(1 - eb[-1] / eb[0]),
           "power_law_index": pl["p"], "power_law_err": pl["p_err"], "power_law_range": [pl["lo"], pl["hi"]],
           "werner2016_index": info["werner_index"], "sigma": sig}
    tr = run.tracks("tracks", "electron")
    if len(tr.t) > 2:
        w = parallel_perpendicular_work(tr, -1.0, 1.0)
        gend = np.nan_to_num(w["gamma"][-1])
        top = gend >= np.quantile(gend, 0.99)
        wp, wq = np.nansum(w["W_par"][-1][top]), np.nansum(w["W_perp"][-1][top])
        res["E_parallel_share_top1pct"] = float(wp / (wp + wq)) if (wp + wq) != 0 else float("nan")
    res["passed"] = bool(dev <= 0.05)
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(1, 3, figsize=(15, 4.2))
        i_mid = int(np.argmax(rate)) if np.isfinite(rate).any() else len(sl.t) // 2
        P.field_map(ax[0], sl.read("Jz", i_mid), sl.extent("Jz"), signed=True, label="J_z",
                    xlabel="x [c/ωpe]", ylabel="y [c/ωpe]")
        ax[0].set_title(f"J_z at t={sl.t[i_mid]:.0f} (peak rate)")
        ax[1].plot(sl.t, rate, color=P.color(0), label="PIC")
        ax[1].axhline(0.1, **P.THEORY, label="~0.1 (published)")
        ax[1].set_xlabel(r"$\omega_{pe} t$")
        ax[1].set_ylabel(r"$d\Psi/dt / (B_0 v_A)$")
        ax[1].legend()
        P.spectrum_plot(ax[2], G - 1.0, f, label=f"t={sp.t[-1]:.0f}", fit=pl, xlabel=r"$\gamma - 1$")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "reconnection.png"))
    return res


def summary(results: dict[str, dict], path: str) -> str:
    """Power-law index p of dN/dgamma vs sigma: the Werner et al. (2016) fit 1.9 + 0.7 / sqrt(sigma) and
    the Sironi & Spitkovsky (2014) values, with the measured index of every run at its sigma."""
    from theory import SIRONI_SPITKOVSKY_2014, werner_index

    P.style()
    fig, ax = P.plt.subplots(figsize=(6.4, 4.2))
    s = np.geomspace(1, 60, 200)
    ax.plot(s, [werner_index(v) for v in s], **P.THEORY, label="Werner et al. 2016: 1.9 + 0.7/√σ")
    ss = np.array(sorted(SIRONI_SPITKOVSKY_2014))
    ax.plot(ss, [SIRONI_SPITKOVSKY_2014[k] for k in ss], "s", ms=7, mfc="none", color=P.INK["secondary"],
            label="Sironi & Spitkovsky 2014")
    names = list(results)
    for i, n in enumerate(names):
        r = results[n]
        sig = r.get("sigma", 10.0) * (1 + 0.04 * (i - (len(names) - 1) / 2))  # small offsets to separate points
        ax.errorbar(sig, r["power_law_index"], yerr=r["power_law_err"], fmt="o", ms=7, color=P.color(i),
                    capsize=3, label=f"PIC {n}: p = {r['power_law_index']:.2f}")
    ax.set_xscale("log")
    ax.set_xlabel("magnetization σ")
    ax.set_ylabel(r"power-law index $p$  ($dN/d\gamma \propto \gamma^{-p}$)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def movie(run_dir: str, out: str | None = None, fps: int = 15, n_frames: int = 80,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="relativistic pair reconnection")
    m.slice("slice_xy", "Jz").slice("slice_xy", "n:electron", signed=False).spectrum("spectrum", "electron")
    m.history()
    return m.save(out or os.path.join(run_dir, "reconnection.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
