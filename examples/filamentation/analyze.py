"""Analyse a filamentation run: growth rate of each k_perp mode of B vs cold linear theory,
saturation level vs the magnetic-trapping estimate, filament structure.

Many modes grow at once, and all of them leave the linear regime when the fastest ones saturate, so
the growth rates are fitted over one common linear window taken from the total magnetic energy: from
two e-folds (in amplitude) above its initial noise to two e-folds below its peak.

python examples/filamentation/analyze.py output/filamentation/pair [--movie]
gate: median |Gamma_fit - Gamma_theory| / Gamma_theory <= 10 % over the modes m = 2..10
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from pypic.analysis import open_run
from pypic.analysis import plot as P


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    md = run.series("modes")
    t = md.t
    n = info["params"]["n_modes"]
    gth = np.asarray(info["gamma_theory"])
    ks = np.asarray(info["k_modes"])
    bz, by = md["Bz"], md["By"]
    s = run.scalars
    eb_t, eb = s.t, s["magnetic"]
    noise = np.median(eb[1: max(3, len(eb) // 20)])
    i_pk = int(np.argmax(eb))
    above = np.nonzero(eb[:i_pk] > noise * np.exp(4.0))[0]
    t_a = float(eb_t[above[0]]) if above.size else float(t[1])
    below = np.nonzero(eb[:i_pk] > eb[i_pk] * np.exp(-4.0))[0]
    t_b = float(eb_t[below[0]]) if below.size else float(eb_t[i_pk])
    win = (t >= t_a) & (t <= t_b)
    fits, amps = [], []
    for m in range(n):
        a_m = np.sqrt(np.abs(bz[:, m]) ** 2 + np.abs(by[:, n + m]) ** 2) / np.sqrt(2)  # k along y: B_z; along z: B_y
        fits.append(float(np.polyfit(t[win], np.log(np.maximum(a_m[win], 1e-300)), 1)[0]) if win.sum() > 4
                    else float("nan"))
        amps.append(a_m)
    fits = np.asarray(fits)
    sel = slice(1, min(10, n))
    rel = (fits[sel] - gth[sel]) / gth[sel]
    ek0 = s["kinetic"][0]
    res = {"k": ks.tolist(), "gamma_fit": fits.tolist(), "gamma_theory": gth.tolist(), "linear_window": [t_a, t_b],
           "median_rel_error_m2_10": float(np.nanmedian(np.abs(rel))),
           "max_rel_error_m2_10": float(np.nanmax(np.abs(rel))),
           "eps_B_peak": float(eb.max() / ek0), "t_peak": float(s.t[int(np.argmax(eb))]),
           "B_sat_trapping_estimate": info["B_sat_trapping"]}
    res["passed"] = bool(res["median_rel_error_m2_10"] <= 0.10)
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(1, 3, figsize=(15, 4.2))
        kk = np.linspace(ks[0] * 0.5, ks[-1], 200)
        from pypic.theory.streaming import filamentation_growth

        wp2 = info["wp_tot"] ** 2
        u0 = info["u0"]
        if info["mixture"] == "e-p":
            comps = [(wp2 / 4, u0), (wp2 / 4, -u0), (wp2 / 4, u0), (wp2 / 4, -u0)]
        else:
            comps = [(0.5, u0), (0.5, -u0), (1.0 / info["params"]["mass_ratio"], 0.0)]
        ax[0].plot(kk, filamentation_growth(kk, comps), **P.THEORY, label="cold theory")
        ax[0].plot(ks, fits, "o", color=P.color(0), label="PIC (B modes)")
        ax[0].set_xlabel(r"$k_\perp c/\omega_{pe}$")
        ax[0].set_ylabel(r"$\Gamma/\omega_{pe}$")
        ax[0].legend()
        for i, m in enumerate((1, 3, 7)):
            if m < n:
                ax[1].semilogy(t, amps[m], color=P.color(i), label=f"m={m + 1}")
        ax[1].set_xlabel(r"$\omega_{pe} t$")
        ax[1].set_title(r"$|B(k_\perp)|$")
        ax[1].legend()
        sl = run.series("slice_yz")
        i_pk = int(np.argmin(np.abs(sl.t - res["t_peak"])))
        P.field_map(ax[2], sl.read("Bz", i_pk), sl.extent("Bz"), signed=True, label="B_z",
                    xlabel="y [c/ωpe]", ylabel="z [c/ωpe]")
        ax[2].set_title(f"B_z at t={sl.t[i_pk]:.1f} (magnetic peak)")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "filamentation.png"))
    return res


def movie(run_dir: str) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="current filamentation (Weibel)")
    m.slice("slice_yz", "Bz").slice("slice_yz", "Jx").phase("phase", "electron").history()
    return m.save(os.path.join(run_dir, "filamentation.mp4"), fps=15, n_frames=120)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
