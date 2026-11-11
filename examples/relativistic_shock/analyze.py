"""Analyse a relativistic shock run: jump conditions and particle acceleration.

python examples/relativistic_shock/analyze.py output/relativistic_shock/pair_2d [--movie]
gate: shock speed and density compression within 10 % of the jump conditions of the run's
      effective dimensionality (Synge gas, theory.py); with a magnetised flow (sigma > 0) also the
      mean downstream field over B1 within 10 % of the compression (flux freezing).
reported: downstream temperature, the non-thermal tail (power-law index of dN/dgamma) and
      gamma_max(t), magnetic energy fraction; sigma > 0: the field ahead of the shock (precursor).
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


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
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
    # formed shock, still far from the receding plasma edge (the reflected-particle precursor runs ahead)
    sel = (t >= 0.3 * t_valid) & (t <= 0.7 * t_valid) & np.isfinite(xs)
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
    sigma = info.get("sigma", 0.0)
    if sigma:
        # flux freezing: the field along B1 is compressed with the plasma, B2 / B1 = n2 / n1
        name, B1 = "B" + info["field_component"], info["B1"]
        xb, b = lo.x(name), lo[name][i_last] / B1
        b2 = float(np.mean(b[(xb > 0.3 * xl) & (xb < 0.8 * xl)]))
        ahead = b[(xb > xl + 5 * skin) & (xb < edge - 5 * skin)]
        res.update(sigma=sigma, field_component=name, field_compression=b2,
                   field_compression_rel_error=(b2 - jc["compression"]) / jc["compression"],
                   field_ahead_mean=float(np.mean(ahead)), field_ahead_rms=float(np.std(ahead)),
                   field_ahead_max_deviation=float(np.max(np.abs(ahead - 1.0))))
    # downstream spectrum from the (x, gamma) histogram at the latest time the run is still valid
    xg = run.series("x_gamma")
    ex, eg = xg.edges
    okt = np.nonzero(xg.t <= 0.95 * t_valid)[0]
    ig = int(okt[-1]) if okt.size else len(xg.t) - 1
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
    # tail fitted over an explicit window above the thermal peak: from 3x the peak of gamma dN/dgamma
    # to 0.7x the cut-off (where N(>gamma) falls to 1e-4); needs at least a factor 2 in gamma
    g_pk = float(G[int(np.argmax(G * f))])
    g_cut = energy_cutoff(G, f, 1e-4)
    plw = {"p": float("nan"), "p_err": float("nan"), "lo": 3 * g_pk, "hi": 0.7 * g_cut}
    if np.isfinite(g_cut) and plw["hi"] > 2.0 * plw["lo"]:
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
               tail_window=[plw["lo"], plw["hi"]], tail_index_literature_2d=2.4, t_spectrum=float(xg.t[ig]),
               gamma_max_t=list(map(float, gmax)), gamma_max_growth_exponent=alpha,
               gamma_max_growth_literature=0.5, eps_B_end=float(s["magnetic"][-1] / s["kinetic"][0]))
    res["passed"] = bool(abs(res["beta_sh_rel_error"]) <= 0.10 and abs(res["compression_rel_error"]) <= 0.10
                         and abs(res.get("field_compression_rel_error", 0.0)) <= 0.10)
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
        if sigma:
            ax[1].plot(xb, b, color=P.color(1), lw=0.8, label=f"{name} / B1")
        ax[1].axhline(jc["compression"] * n1, **P.THEORY, label=f"n2/n1={jc['compression']:.2f}")
        ax[1].set_xlabel("x [c/ωpe]")
        ax[1].legend()
        P.spectrum_plot(ax[2], G, f, label="downstream", fit=plw if np.isfinite(plw["p"]) else None,
                        xlabel=r"$\gamma$")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(fig_dir or run_dir, "shock.png"))
    return res


def summary(results: dict[str, dict], path: str) -> str:
    """Shock speed and compression of every run divided by the Synge-gas jump conditions of the run's
    effective dimensionality, with the 10 % gate band (two panels, one per quantity)."""
    P.style()
    names = list(results)
    y = np.arange(len(names))[::-1]
    fig, axes = P.plt.subplots(1, 2, figsize=(10, 0.5 * len(names) + 1.6), sharey=True)
    for ax, (key, th, lab) in zip(axes, (("beta_sh", "beta_sh_theory", r"$\beta_{sh}$"),
                                         ("compression", "compression_theory", r"$n_2/n_1$"))):
        ax.axvspan(0.9, 1.1, color=P.INK["grid"], alpha=0.7, lw=0, label="gate ±10 %")
        ax.axvline(1.0, **P.THEORY, label="jump conditions (Synge gas)")
        for yi, n in zip(y, names):
            r = results[n]
            ratio = r[key] / r[th]
            ax.plot(ratio, yi, "o", ms=8, color=P.color(0), label=("PIC" if yi == y[0] else None))
            ax.annotate(f"{r[key]:.3f} vs {r[th]:.3f} ({ratio - 1:+.1%}), {r['dims_eff']}D", (ratio, yi),
                        xytext=(0, 10), textcoords="offset points", ha="center", fontsize=7,
                        color=P.INK["secondary"], bbox={"fc": P.INK["surface"], "ec": "none", "pad": 0.5})
        ax.set_xlim(0.8, 1.2)
        ax.set_xlabel(f"{lab}: PIC / theory")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y, names)
    axes[0].set_ylim(y.min() - 0.6, y.max() + 0.7)
    axes[0].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def movie(run_dir: str, out: str | None = None, fps: int = 15, n_frames: int = 100,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="relativistic pair shock")
    m.slice("slice_xy", "Bz").slice("slice_xy", "n:electron", signed=False).phase("phase", "electron")
    m.xt("density", "n:electron", signed=False)
    return m.save(out or os.path.join(run_dir, "shock.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
