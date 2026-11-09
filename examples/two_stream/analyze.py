"""Analyse a two-stream run: growth rate vs linear theory, energy/conservation history, phase space.

python examples/two_stream/analyze.py output/two_stream/ei [--movie]
    -> summary.png (+ two_stream.mp4) in the run directory, verdict JSON on stdout
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from pypic.analysis import open_run
from pypic.analysis import plot as P
from pypic.analysis.spectral import fit_growth


def analyze(run_dir: str, show_plot: bool = True):
    run = open_run(run_dir)
    lin = run.meta.get("planner", {})
    nwl = run.meta.get("params", {}).get("n_wavelengths", 4)
    s = run.scalars
    md = run.series("modes")
    modes, mt = md["Ex"], md.t
    res = {}
    line = None
    if modes.size:
        m = min(nwl, modes.shape[1]) - 1  # mode index nwl corresponds to k_max
        g_fit, line = fit_growth(mt, modes[:, m], omega_r=lin.get("omega_r", 0.0))
        gth = lin.get("gamma_max", float("nan"))
        res.update(gamma_fit_kmax=g_fit, gamma_theory=gth, relative_error=(g_fit - gth) / gth)
    b = s["budget"]
    res["energy_error_max"] = float(np.max(np.abs(b - b[0])) / abs(b[0]))
    g = s["gauss"][np.isfinite(s["gauss"])]
    d = s["divB"][np.isfinite(s["divB"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    res["divB_max"] = float(d.max()) if d.size else float("nan")

    if show_plot:
        P.style()
        fig = P.plt.figure(figsize=(14, 8.5))
        P.energy_history(fig.add_subplot(2, 3, 1), run)
        P.conservation(fig.add_subplot(2, 3, 2), run)
        ax = fig.add_subplot(2, 3, 3)
        for j in range(min(8, modes.shape[1])):
            ax.semilogy(mt, np.abs(modes[:, j]), color=P.color(j), lw=2.2 if j == m else 0.8, label=f"m={j + 1}")
        if line is not None:
            ax.semilogy(line[0], line[1] * 2, **P.THEORY,
                        label=f"fit $\\gamma$={g_fit:.3f} (theory {res['gamma_theory']:.3f})")
        ax.legend(ncol=2)
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.set_title(r"$|E_x(k_m)|$, $k_m = 2\pi m/L_x$")
        ph = run.series("phase")
        sp = "electron" if "electron" in ph.keys() else ph.keys()[0]
        nt = len(ph.t)
        for i, idx in enumerate([nt // 2, nt - 1]):
            ax = fig.add_subplot(2, 3, 4 + i)
            P.phase_space(ax, ph.read(sp, idx), ph.edges, xlabel="x [c/ωpe]", ylabel="u_x / c")
            ax.set_title(f"{sp} (x, u_x) at t={ph.t[idx]:.1f}")
        lo = run.series("lineout")
        ax = fig.add_subplot(2, 3, 6)
        x = lo.x("Ex")
        P.field_map(ax, lo["Ex"].T, [x[0], x[-1], lo.t[0], lo.t[-1]], signed=True, label="E_x",
                    xlabel="x [c/ωpe]", ylabel=r"$\omega_{pe} t$")
        ax.set_title(r"$E_x(x, t)$ (y-z average)")
        fig.suptitle(os.path.basename(os.path.normpath(run_dir)))
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "summary.png"), dpi=110)
    return res


def summary(runs: dict[str, str], results: dict[str, dict], out_dir: str) -> list[str]:
    """Cross-configuration comparison with linear theory (one chart per question):
    growth.png  |E_x(k_max)|(t) of every configuration with the fitted rate and the theory slope;
    growth_ratio.png  gamma_fit / gamma_theory with the 10 % gate band;
    energy.png  relative total-energy error, explicit vs implicit, against the implicit 1e-8 contract."""
    P.style()
    paths = []
    names = list(runs)
    nc = min(3, len(names))
    nr = -(-len(names) // nc)
    fig, axes = P.plt.subplots(nr, nc, figsize=(4.2 * nc, 3.3 * nr), squeeze=False)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    for ax, name in zip(axes.flat, names):
        run = open_run(runs[name])
        lin = run.meta.get("planner", {})
        nwl = run.meta.get("params", {}).get("n_wavelengths", 4)
        md = run.series("modes")
        a = np.abs(md["Ex"][:, min(nwl, md["Ex"].shape[1]) - 1])
        g_fit, line = fit_growth(md.t, md["Ex"][:, min(nwl, md["Ex"].shape[1]) - 1], omega_r=lin.get("omega_r", 0.0))
        g_th = results[name]["gamma_theory"]
        ax.semilogy(md.t, a, color=P.color(0), lw=1.2, label=r"PIC $|E_x(k_{max})|$")
        tl, env = line
        tm = 0.5 * (tl[0] + tl[-1])
        em = float(np.exp(np.interp(tm, tl, np.log(env))))
        tt = np.linspace(tl[0] - 0.25 * (tl[-1] - tl[0]), tl[-1] + 0.25 * (tl[-1] - tl[0]), 50)
        ax.semilogy(tt, 2 * em * np.exp(g_th * (tt - tm)), **P.THEORY, label=rf"theory $\gamma$={g_th:.3f}")
        ax.axvspan(tl[0], tl[-1], color=P.INK["grid"], alpha=0.6, lw=0, label="fit window")
        err = results[name]["relative_error"]
        ax.set_title(f"{name}: fit {g_fit:.3f} ({err:+.1%})")
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.legend(loc="lower right")
    fig.tight_layout()
    paths.append(P.save(fig, os.path.join(out_dir, "growth.png"), dpi=110))

    fig, ax = P.plt.subplots(figsize=(6.0, 0.55 * len(names) + 1.2))
    y = np.arange(len(names))[::-1]
    ax.axvspan(0.9, 1.1, color=P.INK["grid"], alpha=0.7, lw=0, label="gate ±10 %")
    ax.axvline(1.0, **P.THEORY, label="linear theory")
    ratio = [results[n]["gamma_fit_kmax"] / results[n]["gamma_theory"] for n in names]
    ax.plot(ratio, y, "o", ms=8, color=P.color(0), label="PIC fit")
    for yi, r in zip(y, ratio):
        ax.annotate(f"{r - 1:+.1%}", (r, yi), xytext=(0, 9), textcoords="offset points", ha="center",
                    fontsize=8, color=P.INK["secondary"], bbox={"fc": P.INK["surface"], "ec": "none", "pad": 0.5})
    ax.set_yticks(y, names)
    ax.set_xlim(0.85, 1.15)
    ax.set_ylim(y.min() - 0.5, y.max() + 0.6)
    ax.set_xlabel(r"$\gamma_{fit} / \gamma_{theory}$ (k_max mode)")
    ax.legend(loc="upper left", ncol=3, bbox_to_anchor=(0, 1.18))
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    paths.append(P.save(fig, os.path.join(out_dir, "growth_ratio.png"), dpi=110))

    pair = [n for n in ("ei", "ei_implicit") if n in runs]
    if len(pair) == 2:
        fig, ax = P.plt.subplots(figsize=(6.5, 3.6))
        for i, n in enumerate(pair):
            s = open_run(runs[n]).scalars
            b = s["budget"]
            ax.semilogy(s.t, np.abs(b - b[0]) / abs(b[0]) + 1e-17, color=P.color(i),
                        label=("explicit (ei)" if n == "ei" else "implicit (ei_implicit)"))
        ax.axhline(1e-8, **P.THEORY, label="implicit contract 1e-8")
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.set_ylabel("|ΔE_total| / E_total(0)")
        ax.set_ylim(1e-16, 1)
        ax.legend(loc="lower right")
        fig.tight_layout()
        paths.append(P.save(fig, os.path.join(out_dir, "energy.png"), dpi=110))
    return paths


def movie(run_dir: str, out: str | None = None, fps: int = 15, n_frames: int = 300, dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="two-stream " + os.path.basename(os.path.normpath(run_dir)))
    m.phase("phase", "electron").xt("lineout", "Ex").slice("slice_xy", "Ex").history()
    return m.save(out or os.path.join(run_dir, "two_stream.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
