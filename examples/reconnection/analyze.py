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


def flux_maps(run_dir: str):
    """(run, slice series, x, y, A(t) list): the flux function A_z of every slice output; its contours
    are the in-plane magnetic field lines (B = curl(A_z z))."""
    run = open_run(run_dir)
    sl = run.series("slice_xy")
    xa, ya = sl.plane_axes("Bx")
    dxs, dys = xa[1] - xa[0], ya[1] - ya[0]
    A = [flux_function(sl.read("Bx", i), sl.read("By", i), dxs, dys) for i in range(len(sl.t))]
    return run, sl, xa, ya, A


def null_line(A: np.ndarray, ya: np.ndarray, y_sheet: float, half_width: float):
    """(a, y, valid, across_max) on the magnetic null line B_x = dA/dy = 0 of one current sheet: the extremum
    of A across the band |y - y_sheet| < half_width in every column. ``valid`` is False where that extremum
    lies on the band edge: the null line has left the band and the column has no stationary point."""
    band = np.nonzero(np.abs(ya - y_sheet) < half_width)[0]
    sub = A[:, band]
    mid = sub[:, len(band) // 2]
    # across the sheet A has a maximum or a minimum on the null line, the same kind everywhere
    across_max = bool(np.median(sub.max(axis=1) - mid) < np.median(mid - sub.min(axis=1)))
    j = np.argmax(sub, axis=1) if across_max else np.argmin(sub, axis=1)
    a = sub[np.arange(sub.shape[0]), j]
    valid = (j > 0) & (j < len(band) - 1)
    return a, ya[band][j], valid, across_max


def null_points(A: np.ndarray, xa: np.ndarray, ya: np.ndarray, y_sheet: float, half_width: float,
                prominence: float) -> tuple[list, list]:
    """X and O points of one current sheet: the extrema of A along its null line (``null_line``) are
    X points where A is a saddle and O points (island centres) where A is an extremum in both directions.
    Candidates on the edge of the search band are not stationary points and are dropped."""
    from scipy.signal import find_peaks

    a, yl, valid, across_max = null_line(A, ya, y_sheet, half_width)
    n = len(a)
    ext = np.concatenate([a, a, a])  # periodic in x
    hi, _ = find_peaks(ext, prominence=prominence)
    lo, _ = find_peaks(-ext, prominence=prominence)
    hi = [i - n for i in hi if n <= i < 2 * n and valid[i - n]]
    lo = [i - n for i in lo if n <= i < 2 * n and valid[i - n]]
    # an extremum of the same kind as across the sheet is an O point; the opposite kind is an X point
    o_idx, x_idx = (hi, lo) if across_max else (lo, hi)
    return [(xa[i], yl[i]) for i in x_idx], [(xa[i], yl[i]) for i in o_idx]


def xo_flux(A: np.ndarray, ya: np.ndarray, y_sheet: float, half_width: float) -> float:
    """Reconnected flux of one sheet measured between its X and O points: the range of A along the null
    line, i.e. A at the main O point minus A at the deepest X point (0 for the unperturbed sheet)."""
    a, _, valid, _ = null_line(A, ya, y_sheet, half_width)
    return float(a[valid].max() - a[valid].min()) if valid.any() else float("nan")


def _half_width(info: dict) -> float:
    return 4 * info["delta"] + 6.0


def _levels(A0: np.ndarray, n: int = 26) -> np.ndarray:
    return np.linspace(float(A0.min()), float(A0.max()), n + 2)[1:-1]


def _draw_field(ax, data, A, xa, ya, levels, info, vlim, label, signed=True, step=2):
    ext = [xa[0], xa[-1], ya[0], ya[-1]]
    im = P.field_map(ax, data, ext, signed=signed, vlim=vlim, label=label, xlabel="x [c/ωpe]",
                     ylabel="y [c/ωpe]")
    cs = ax.contour(xa[::step], ya[::step], A[::step, ::step].T, levels=levels, colors=P.INK["primary"],
                    linewidths=0.45, alpha=0.75, linestyles="solid")
    marks = []
    prom = 0.02 * float(levels[-1] - levels[0])
    for ys in info["y_sheets"]:
        xs, os_ = null_points(A, xa, ya, ys, _half_width(info), prom)
        if xs:
            marks += ax.plot(*zip(*xs), "x", ms=7, mew=2, color="#eda100", label=None)
        if os_:
            marks += ax.plot(*zip(*os_), "o", ms=5, mfc="none", mew=1.5, color="#008300", label=None)
    return im, cs, marks


def fieldlines(run_dir: str, path: str, n_snap: int = 6) -> str:
    """Key chart of reconnection: J_z with the magnetic field lines (contours of the flux function, the
    same levels in every panel) at n_snap times, X points (×) and O points (○) of both sheets, and the
    reconnected flux of the lower sheet against time with the snapshot times marked."""
    run, sl, xa, ya, A = flux_maps(run_dir)
    info = run.meta["planner"]
    t = sl.t
    psi = np.array([xo_flux(a, ya, info["y_sheets"][0], _half_width(info)) for a in A])
    i_pk = int(np.argmax(np.gradient(psi, t)))
    # snapshots where the reconnected flux reaches given fractions of its final value
    fr = np.linspace(0.15, 0.9, n_snap - 2)
    picks = [0] + [int(np.argmax(psi >= f * psi[-1])) for f in fr] + [len(t) - 1]
    picks = sorted(set(picks))
    levels = _levels(A[0])
    vlim = P.robust_limits([sl.read("Jz", i) for i in picks], signed=True)
    P.style()
    fig = P.plt.figure(figsize=(16, 9.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.8])
    for k, i in enumerate(picks):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        _draw_field(ax, sl.read("Jz", i), A[i], xa, ya, levels, info, vlim, "J_z")
        ax.set_title(f"t = {t[i]:.0f} / ω_pe,  Ψ = {psi[i] / psi[-1]:.0%} of final"
                     + ("  (peak rate)" if i == i_pk else ""))
    ax = fig.add_subplot(gs[2, :])
    ax.plot(t, psi / (info["B0"] * info["Lx"]), color=P.color(0), label="reconnected flux Ψ = A(O) − A(X), lower sheet")
    for i in picks:
        ax.axvline(t[i], color=P.INK["muted"], lw=0.8, ls=":")
    ax.set_xlabel(r"$\omega_{pe} t$")
    ax.set_ylabel(r"$\Psi / (B_0 L_x)$")
    ax.plot([], [], "x", color="#eda100", mew=2, label="X point (field lines reconnect)")
    ax.plot([], [], "o", mfc="none", color="#008300", mew=1.5, label="O point (plasmoid centre)")
    ax.plot([], [], color=P.INK["primary"], lw=0.8, label="field lines (contours of A_z, fixed levels)")
    ax.legend(loc="lower right")
    fig.tight_layout()
    return P.save(fig, path, dpi=100)


class FieldLinePanel:
    """Movie panel: a slice component with the magnetic field lines (fixed contour levels of A_z) and the
    X / O points redrawn at every frame (pypic.analysis.animate.Panel interface)."""

    def __init__(self, run_dir, comp="Jz", signed=True, flux=None):
        self.run, self.sl, self.xa, self.ya, self.A = flux or flux_maps(run_dir)
        self.info = self.run.meta["planner"]
        self.comp, self.signed = comp, signed
        self.t = self.sl.t
        self.levels = _levels(self.A[0])
        self.vlim = P.robust_limits([self.sl.read(comp, i) for i in range(0, len(self.t), max(1, len(self.t) // 12))],
                                    signed=signed)
        self.art = []

    def times(self):
        return self.t

    def prepare(self, frame_times):
        pass

    def init(self, ax):
        self.ax = ax

    def update(self, t):
        i = max(0, int(np.searchsorted(self.t, t + 1e-9, side="right")) - 1)
        self.ax.clear()
        lab = {"Jz": "J_z", "n:electron": "n_e"}.get(self.comp, self.comp)
        if not hasattr(self, "_cbar"):
            im, cs, marks = _draw_field(self.ax, self.sl.read(self.comp, i), self.A[i], self.xa, self.ya, self.levels,
                                        self.info, self.vlim, lab, self.signed)
            self._cbar = True
        else:
            ext = [self.xa[0], self.xa[-1], self.ya[0], self.ya[-1]]
            cmap = P.DIVERGING if self.signed else P.SEQUENTIAL
            self.ax.imshow(np.asarray(self.sl.read(self.comp, i)).T, origin="lower", extent=ext, aspect="auto",
                           cmap=cmap, vmin=self.vlim[0], vmax=self.vlim[1], interpolation="nearest")
            self.ax.contour(self.xa[::2], self.ya[::2], self.A[i][::2, ::2].T, levels=self.levels,
                            colors=P.INK["primary"], linewidths=0.45, alpha=0.75, linestyles="solid")
            prom = 0.02 * float(self.levels[-1] - self.levels[0])
            for ys in self.info["y_sheets"]:
                xs, os_ = null_points(self.A[i], self.xa, self.ya, ys, _half_width(self.info), prom)
                if xs:
                    self.ax.plot(*zip(*xs), "x", ms=7, mew=2, color="#eda100")
                if os_:
                    self.ax.plot(*zip(*os_), "o", ms=5, mfc="none", mew=1.5, color="#008300")
            self.ax.set_xlabel("x [c/ωpe]")
            self.ax.set_ylabel("y [c/ωpe]")
        self.ax.set_title(f"{lab} + field lines (× X point, ○ O point)  t={self.t[i]:.0f}")


class FluxPanel:
    """Movie panel: reconnected flux of the lower sheet, revealed up to the frame time."""

    def __init__(self, flux):
        run, sl, xa, ya, A = flux
        info = run.meta["planner"]
        self.t = sl.t
        hw = _half_width(info)
        self.psi = np.array([xo_flux(a, ya, info["y_sheets"][0], hw) for a in A]) / (info["B0"] * info["Lx"])

    def times(self):
        return self.t

    def prepare(self, frame_times):
        pass

    def init(self, ax):
        ax.plot(self.t, self.psi, color=P.INK["grid"], lw=1.0)
        (self.line,) = ax.plot([], [], color=P.color(0), lw=1.6)
        ax.set_xlim(self.t[0], self.t[-1])
        ax.set_ylim(0, 1.1 * float(self.psi.max()))
        ax.set_xlabel(r"$\omega_{pe} t$")
        ax.set_ylabel(r"$\Psi / (B_0 L_x)$")
        ax.set_title("reconnected flux A(O) − A(X), lower sheet")

    def update(self, t):
        i = max(0, int(np.searchsorted(self.t, t + 1e-9, side="right")) - 1)
        self.line.set_data(self.t[: i + 1], self.psi[: i + 1])


def movie(run_dir: str, out: str | None = None, fps: int = 15, n_frames: int = 80,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    flux = flux_maps(run_dir)
    m = Movie(flux[0], title="relativistic pair reconnection", panel_size=(6.4, 3.6))
    m.add(FieldLinePanel(run_dir, "Jz", True, flux)).add(FieldLinePanel(run_dir, "n:electron", False, flux))
    m.add(FluxPanel(flux)).spectrum("spectrum", "electron")
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
