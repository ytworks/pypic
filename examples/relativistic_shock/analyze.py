"""Analyse a relativistic shock run: jump conditions and particle acceleration.

python examples/relativistic_shock/analyze.py output/relativistic_shock/pair_2d [--movie]
gate: shock speed and density compression within 10 % of the jump conditions of the run's
      effective dimensionality (Synge gas, theory.py); reported: downstream temperature, the
      non-thermal tail (power-law index of dN/dgamma) and gamma_max(t), magnetic energy fraction.
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
from pypic.analysis.fields import fit_velocity, front_position  # noqa: E402
from pypic.analysis.particles import energy_cutoff, fit_power_law, fit_temperature, spectrum  # noqa: E402


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    skin = info["skin"]
    jc = info["jump"]
    lo = run.series("density")
    t = lo.t
    x = lo.x("n:electron")
    ne = lo["n:electron"]
    xs = np.array([front_position(x, ne[i], 2.0) for i in range(len(t))])
    t_valid = info["t_valid"]
    sel = (t >= 0.4 * t[-1]) & (t <= 0.95 * t_valid) & np.isfinite(xs)
    v_sh, _ = fit_velocity(t[sel], xs[sel])
    i_last = int(np.nonzero(sel)[0][-1])
    xl = xs[i_last]
    down = (x > 0.3 * xl) & (x < 0.8 * xl)
    # upstream: between the shock and the receding plasma edge L - beta0 t (vacuum gap beyond it)
    edge = info["Lx"] - np.sqrt(1 - 1 / info["gamma0"] ** 2) * t[i_last]
    up = (x > xl + 5 * skin) & (x < edge - 5 * skin)
    n2 = float(np.median(ne[i_last][down]))
    # the jump conditions refer to the incoming flow, whose lab density is 1 by construction; the density
    # measured just ahead of the shock also contains particles reflected into the upstream (precursor)
    n1_ahead = float(np.median(ne[i_last][up])) if up.sum() > 5 else float("nan")
    n1 = 1.0
    comp = n2 / n1
    res = {"dims_eff": info["dims_eff"], "beta_sh": v_sh, "beta_sh_theory": jc["beta_sh"],
           "beta_sh_rel_error": (v_sh - jc["beta_sh"]) / jc["beta_sh"], "compression": comp,
           "compression_theory": jc["compression"],
           "compression_rel_error": (comp - jc["compression"]) / jc["compression"],
           "compression_gamma_3_2": info["jump_gamma_3_2"]["compression"], "density_ahead_of_shock": n1_ahead,
           "beta_sh_gamma_3_2": info["jump_gamma_3_2"]["beta_sh"],
           "t_measured": float(t[i_last])}
    # downstream spectrum from the (x, gamma) histogram at the last valid time
    xg = run.series("x_gamma")
    ex, eg = xg.edges
    ig = int(np.argmin(np.abs(xg.t - t[i_last])))
    xc = 0.5 * (ex[1:] + ex[:-1])
    xs_g = float(np.interp(xg.t[ig], t, xs))
    dmask = (xc > 0.3 * xs_g) & (xc < 0.8 * xs_g)
    h = sum(xg.read(sp, ig)[dmask].sum(axis=0) for sp in xg.keys())
    G, f = spectrum(eg, h)
    theta, nth = fit_temperature(G, f)
    pl = fit_power_law(G - 1.0, f)
    gmax = []
    for i in range(len(xg.t)):
        xs_i = float(np.interp(xg.t[i], t, xs))
        m = (xc > 0.3 * xs_i) & (xc < 0.8 * xs_i) if np.isfinite(xs_i) and xs_i > 0 else np.zeros_like(xc, bool)
        hh = sum(xg.read(sp, i)[m].sum(axis=0) for sp in xg.keys()) if m.any() else np.zeros(len(eg) - 1)
        Gi, fi = spectrum(eg, hh)
        gmax.append(energy_cutoff(Gi, fi, 1e-4) if hh.sum() > 0 else np.nan)
    # tail fitted over an explicit window above the thermal peak: from 4x the peak of gamma dN/dgamma
    # to half the cut-off (where N(>gamma) falls to 1e-4); needs at least half a decade
    g_pk = float(G[int(np.argmax(G * f))])
    g_cut = energy_cutoff(G, f, 1e-4)
    plw = {"p": float("nan"), "p_err": float("nan"), "lo": 4 * g_pk, "hi": 0.5 * g_cut}
    if np.isfinite(g_cut) and plw["hi"] > 3.0 * plw["lo"]:
        plw = fit_power_law(G, f, lo=plw["lo"], hi=plw["hi"])
    # maximum energy growth gamma_max ~ t^alpha (Sironi et al. 2013: alpha ~ 1/2) over the second half
    gm = np.asarray(gmax, dtype=float)
    tg = xg.t
    okg = np.isfinite(gm) & (tg >= 0.5 * tg[-1]) & (tg <= t_valid)
    alpha = float(np.polyfit(np.log(tg[okg]), np.log(gm[okg]), 1)[0]) if okg.sum() >= 4 else float("nan")
    s = run.scalars
    res.update(theta_downstream_fit=theta, theta2_theory=jc["theta2"], nonthermal_fraction=nth,
               tail_index_p=pl["p"], tail_index_err=pl["p_err"], tail_range=[pl["lo"], pl["hi"]],
               tail_index_window_p=plw["p"], tail_index_window_err=plw["p_err"],
               tail_window=[plw["lo"], plw["hi"]], tail_index_literature_2d=2.4,
               gamma_max_t=list(map(float, gmax)), gamma_max_growth_exponent=alpha,
               gamma_max_growth_literature=0.5, eps_B_end=float(s["magnetic"][-1] / s["kinetic"][0]))
    res["passed"] = bool(abs(res["beta_sh_rel_error"]) <= 0.10 and abs(res["compression_rel_error"]) <= 0.10)
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(1, 3, figsize=(15, 4.2))
        P.field_map(ax[0], ne.T, [x[0], x[-1], t[0], t[-1]], signed=False, label="n_e / n0",
                    xlabel="x [c/ωpe]", ylabel=r"$\omega_{pe} t$")
        tt = np.linspace(0, t[-1], 50)
        ax[0].plot(jc["beta_sh"] * tt, tt, **P.THEORY, label=f"jump cond. β_sh={jc['beta_sh']:.3f}")
        ax[0].legend()
        ax[1].plot(x, ne[i_last], color=P.color(0), label=f"n_e, t={t[i_last]:.0f}")
        ax[1].axhline(jc["compression"] * n1, **P.THEORY, label=f"n2/n1={jc['compression']:.2f}")
        ax[1].set_xlabel("x [c/ωpe]")
        ax[1].legend()
        P.spectrum_plot(ax[2], G, f, label="downstream", fit=plw if np.isfinite(plw["p"]) else None,
                        xlabel=r"$\gamma$")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "shock.png"))
    return res


def movie(run_dir: str) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="relativistic pair shock")
    m.slice("slice_xy", "Bz").slice("slice_xy", "n:electron", signed=False).phase("phase", "electron")
    m.xt("density", "n:electron", signed=False)
    return m.save(os.path.join(run_dir, "shock.mp4"), fps=15, n_frames=100)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
