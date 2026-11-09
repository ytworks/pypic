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


def movie(run_dir: str, out: str | None = None) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="two-stream " + os.path.basename(os.path.normpath(run_dir)))
    m.phase("phase", "electron").xt("lineout", "Ex").slice("slice_xy", "Ex").history()
    return m.save(out or os.path.join(run_dir, "two_stream.mp4"), fps=15)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
