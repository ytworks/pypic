"""Analyse a filamentation run: growth rate of each k_perp mode of B vs cold linear theory,
saturation level vs the magnetic-trapping estimate, filament structure.

Many modes grow at once, and all of them leave the linear regime when the fastest ones saturate, so
the growth rates are fitted over one common linear window taken from the total magnetic energy: from
two e-folds (in amplitude) above its initial noise to two e-folds below its peak. With a magnetic
field along the drift the total magnetic energy contains that field, and the window is taken from
the energy of the analysed modes instead.

python examples/filamentation/analyze.py output/filamentation/pair [--movie]
gate: median |Gamma_fit - Gamma_theory| / Gamma_theory <= 10 % over the modes m = 2..10 (from
gate_m_min with a field); above the field of marginal stability: the energy of the modes m = 2..10
at the end over that at 10 % of the run <= e^2
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from pypic.analysis import open_run
from pypic.analysis import plot as P


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    md = run.series("modes")
    t = md.t
    n = info["params"]["n_modes"]
    gth = np.asarray(info["gamma_theory"])
    ks = np.asarray(info["k_modes"])
    bz, by = md["Bz"], md["By"]
    s = run.scalars
    # amplitude of mode m: k along y shows in B_z, k along z in B_y
    amps = [np.sqrt(np.abs(bz[:, m]) ** 2 + np.abs(by[:, n + m]) ** 2) / np.sqrt(2) for m in range(n)]
    g = s["gauss"][np.isfinite(s["gauss"])]
    gauss = float(g.max()) if g.size else float("nan")
    if info.get("stable"):
        return _stable(run_dir, t, amps, info, float(s.t[-1]), gauss, show_plot, fig_dir)
    field = bool(info.get("B0"))
    eb_t, eb = (t, np.sum(np.square(amps), axis=0)) if field else (s.t, s["magnetic"])
    noise = np.median(eb[1: max(3, len(eb) // 20)])
    i_pk = int(np.argmax(eb))
    above = np.nonzero(eb[:i_pk] > noise * np.exp(4.0))[0]
    t_a = float(eb_t[above[0]]) if above.size else float(t[1])
    below = np.nonzero(eb[:i_pk] > eb[i_pk] * np.exp(-4.0))[0]
    t_b = float(eb_t[below[0]]) if below.size else float(eb_t[i_pk])
    win = (t >= t_a) & (t <= t_b)
    fits = np.asarray([float(np.polyfit(t[win], np.log(np.maximum(a_m[win], 1e-300)), 1)[0]) if win.sum() > 4
                       else float("nan") for a_m in amps])
    m_lo = info.get("gate_modes", [2, 10])[0]
    sel = slice(m_lo - 1, min(10, n))
    rel = (fits[sel] - gth[sel]) / gth[sel]
    ek0 = s["kinetic"][0]
    # the energy of the field along the drift is the magnetic energy of step 0 (its flux is conserved)
    eb_all = s["magnetic"] - s["magnetic"][0] if field else s["magnetic"]
    key = f"rel_error_m{m_lo}_10"
    res = {"k": ks.tolist(), "gamma_fit": fits.tolist(), "gamma_theory": gth.tolist(), "linear_window": [t_a, t_b],
           "median_" + key: float(np.nanmedian(np.abs(rel))), "max_" + key: float(np.nanmax(np.abs(rel))),
           "eps_B_peak": float(eb_all.max() / ek0), "t_peak": float(s.t[int(np.argmax(eb_all))]),
           "B_sat_trapping_estimate": info["B_sat_trapping"]}
    res["passed"] = bool(res["median_" + key] <= 0.10)
    res["gauss_max"] = gauss
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
        ax[0].plot(kk, filamentation_growth(kk, comps, info.get("omega_c", 0.0)), **P.THEORY, label="cold theory")
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
        res["plot"] = P.save(fig, os.path.join(fig_dir or run_dir, "filamentation.png"))
    return res


STABLE_LIMIT = float(np.exp(2.0))  # growth of the mode energy allowed above the field of marginal stability


def stable_ratio(t, amps, t_end: float) -> float:
    """Energy of the modes m = 2..10 at the last output over that at the first output from 10 % of the run."""
    w = np.sum(np.square(amps[1:10]), axis=0)
    return float(w[-1] / w[int(np.argmax(t >= 0.1 * t_end))])


def _stable(run_dir, t, amps, info, t_end, gauss, show_plot, fig_dir) -> dict:
    ratio = stable_ratio(t, amps, t_end)
    res = {"mode_energy_ratio_m2_10": ratio, "limit": STABLE_LIMIT, "omega_c": info["omega_c"],
           "omega_c_critical": info["omega_c_critical"], "passed": bool(ratio <= STABLE_LIMIT), "gauss_max": gauss}
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(figsize=(6, 4.2))
        for i, m in enumerate((1, 3, 7)):
            ax.semilogy(t, amps[m], color=P.color(i), label=f"m={m + 1}")
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.set_title(r"$|B(k_\perp)|$ above the field of marginal stability")
        ax.legend()
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(fig_dir or run_dir, "filamentation.png"))
    return res


def movie(run_dir: str, out: str | None = None, fps: int = 15, n_frames: int = 120,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="current filamentation (Weibel)")
    m.slice("slice_yz", "Bz").slice("slice_yz", "Jx").phase("phase", "electron").history()
    return m.save(out or os.path.join(run_dir, "filamentation.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
