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


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
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
            # strong (alpha = 0.5) Landau damping benchmark (Cheng & Knorr 1976 and later Vlasov codes,
            # k = 0.5, alpha = 0.5): the envelope of |E_k| first decays at gamma_1 ~ -0.29, then grows
            # by particle trapping at gamma_2 ~ +0.08 (t ~ 20-40)
            amp = np.abs(z)
            nwin = max(1, int(np.pi / w_th / np.median(np.diff(t))))
            env = np.array([amp[max(0, i - nwin): i + nwin + 1].max() for i in range(len(amp))])
            g1 = float(np.polyfit(t[(t > 2) & (t < 15)], np.log(env[(t > 2) & (t < 15)]), 1)[0])
            w2 = (t > 20) & (t < 40)
            g2 = float(np.polyfit(t[w2], np.log(env[w2]), 1)[0]) if w2.sum() > 4 else float("nan")
            res.update(gamma_1_fit=g1, gamma_1_reference=-0.29, gamma_2_fit=g2, gamma_2_reference=0.08,
                       linear_gamma_theory=g_th, note="nonlinear benchmark: reported, not gated")
        fig, ax = P.plt.subplots(figsize=(6.5, 4))
        P.style()
        ax.semilogy(t, np.abs(z), color=P.color(0), label=r"$|E_x(k)|$")
        tt = np.linspace(*res.get("fit_window", (0, 15)), 50)
        a0 = np.abs(z[np.searchsorted(t, tt[0])])
        ax.semilogy(tt, a0 * np.exp(g_th * (tt - tt[0])), **P.THEORY, label=rf"linear theory $\gamma$={g_th:.4f}")
        if case == "nonlinear":
            # reference envelope slopes (Cheng & Knorr 1976), anchored on the measured envelope
            for (ta, tb), gref, ls in (((2.0, 15.0), -0.29, ":"), ((20.0, 40.0), 0.08, "-.")):
                tr_ = np.linspace(ta, tb, 30)
                e0 = float(np.interp(ta, t, env)) * 1.6
                ax.semilogy(tr_, e0 * np.exp(gref * (tr_ - ta)), color=P.INK["primary"], ls=ls, lw=1.2,
                            label=rf"Cheng & Knorr $\gamma$={gref:+.2f}")
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
        res["plots"] = [P.save(f, os.path.join(fig_dir or run_dir, f"{name}.png")) for name, f in figs]
    else:
        for _, f in figs:
            P.plt.close(f)
    return res


def summary(results: dict[str, dict], path: str) -> str:
    """|measured - theory| / theory of every gated quantity against its tolerance (log scale);
    the nonlinear benchmark rates are reported without a tolerance (hollow markers)."""
    rows = []  # (label, error, tolerance or None)
    r = results.get("landau")
    if r:
        rows += [("landau: damping rate γ", r["gamma_rel_error"], 0.10),
                 ("landau: frequency ω", r["omega_rel_error"], 0.02)]
    r = results.get("nonlinear")
    if r:
        rows += [("nonlinear: γ₁ vs −0.29", (r["gamma_1_fit"] + 0.29) / -0.29, None),
                 ("nonlinear: γ₂ vs +0.08", (r["gamma_2_fit"] - 0.08) / 0.08, None)]
    r = results.get("dawson")
    if r:
        rows += [("dawson: frequency (median)", r["freq_rel_error_median"], 0.01),
                 ("dawson: amplitude (median)", r["amp_rel_error_median"], 0.02)]
    r = results.get("wavebreaking")
    if r:
        rows.append(("wavebreaking: crossing time", r["rel_error"], 0.05))
    r = results.get("relativistic")
    if r:
        rows.append(("relativistic: period (median)", r["freq_rel_error_median"], 0.005))
    P.style()
    fig, ax = P.plt.subplots(figsize=(7.2, 0.42 * len(rows) + 1.2))
    y = np.arange(len(rows))[::-1]
    for i, (yi, (_lab, err, tol)) in enumerate(zip(y, rows)):
        e = max(abs(err), 1e-5)
        ax.plot(e, yi, "o", ms=8, color=P.color(0), mfc=P.color(0) if tol else "none",
                label=("measured |error|" if i == 0 else None))
        if tol:
            ax.plot(tol, yi, "|", ms=16, mew=2, color=P.INK["primary"], label=("gate tolerance" if i == 0 else None))
        ax.annotate(f"{err:+.2%}", (1.0, yi), xycoords=("axes fraction", "data"), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=7.5, color=P.INK["secondary"])
    ax.plot([], [], "o", color=P.color(0), mfc="none", label="reported (not gated)")
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 1.0)
    ax.set_yticks(y, [r_[0] for r_ in rows], fontsize=8)
    ax.set_xlabel("|PIC − theory| / theory   (right: signed error)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def sheets_movie(run_dir: str, out: str, fps: int = 10, n_frames: int = 90, dpi: int = 70,
                 n_sheets: int = 1500) -> str:
    """Cold-plasma sheets from the tracked electrons (cold cases): the phase space (x, u_x) coloured by
    the initial position x0, the Lagrangian map x0 -> x (monotonic until the sheets cross, i.e. until the
    wave breaks) and the electron density from the sheet positions (it spikes where they cross)."""
    from matplotlib import animation

    run = open_run(run_dir)
    info = run.meta["planner"]
    tr = run.tracks("tracks", "electron")
    t, L = tr.t, info["L"]
    x, ux = tr["x"], tr["ux"]
    ok = np.nonzero(np.all(np.isfinite(x), axis=0))[0]
    ok = ok[np.argsort(x[0, ok])]
    cols = ok[np.linspace(0, len(ok) - 1, min(n_sheets, len(ok))).astype(int)]
    x = x[:, cols]
    x = x[0] + np.unwrap((x - x[0]) * 2 * np.pi / L, axis=0) * L / (2 * np.pi)
    ux = ux[:, cols]
    x0 = x[0]
    t_cross = T.crossing_time(info["k_v0"]) if info.get("k_v0", 0) > 1 else None
    P.style()
    fig, ax = P.plt.subplots(1, 3, figsize=(13, 3.6))
    umax = 1.15 * float(np.nanmax(np.abs(ux)))
    sc = ax[0].scatter(np.mod(x[0], L), ux[0], c=x0, s=3, cmap=P.SEQUENTIAL, vmin=0, vmax=L)
    ax[0].set_xlim(0, L)
    ax[0].set_ylim(-umax, umax)
    ax[0].set_xlabel("x [c/ωpe]")
    ax[0].set_ylabel("u_x / c")
    ax[0].set_title("phase space (colour: initial position x₀)")
    ln, = ax[1].plot(x0, x[0], color=P.color(0), lw=1.4, label="PIC sheets")
    ax[1].plot(x0, x0, **{**P.THEORY, "lw": 0.8}, label="x = x₀")
    span = float(np.nanmax(np.abs(x - x0)))
    ax[1].set_xlim(0, L)
    ax[1].set_ylim(-span, L + span)
    ax[1].set_xlabel("x₀ [c/ωpe]")
    ax[1].set_ylabel("x [c/ωpe]")
    ax[1].set_title("Lagrangian map: sheets cross where it folds")
    ax[1].legend(loc="upper left", fontsize=7)
    edges = np.linspace(0, L, 129)
    xc = 0.5 * (edges[1:] + edges[:-1])

    def density(i):
        h, _ = np.histogram(np.mod(x[i], L), bins=edges)
        return h / (len(x0) / (len(edges) - 1))

    nl, = ax[2].plot(xc, density(0), color=P.color(1), lw=1.4)
    ax[2].set_xlim(0, L)
    ax[2].set_ylim(0, max(4.0, 1.1 * max(float(density(i).max()) for i in range(0, len(t), max(1, len(t) // 60)))))
    ax[2].set_xlabel("x [c/ωpe]")
    ax[2].set_ylabel("n_e / n₀")
    ax[2].set_title("electron density")
    title = fig.suptitle("")
    fig.tight_layout()
    idx = np.linspace(0, len(t) - 1, n_frames).astype(int)

    def update(f):
        i = idx[f]
        sc.set_offsets(np.column_stack([np.mod(x[i], L), ux[i]]))
        ln.set_ydata(x[i])
        nl.set_ydata(density(i))
        s = f"{info['case']}: ω_pe t = {t[i]:.2f}"
        if t_cross is not None:
            s += f"   (theory: sheets first cross at {t_cross:.3f}" + (", crossed)" if t[i] >= t_cross else ")")
        title.set_text(s)
        return sc, ln, nl, title

    ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=False)
    ani.save(out, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    P.plt.close(fig)
    return out


def movie(run_dir: str, out: str | None = None, fps: int = 20, n_frames: int = 150,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title=f"Langmuir oscillation: {run.meta['planner']['case']}")
    m.phase("phase", "electron").xt("lineout", "Ex").spectrum("spectrum", "electron").history()
    return m.save(out or os.path.join(run_dir, "langmuir.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
