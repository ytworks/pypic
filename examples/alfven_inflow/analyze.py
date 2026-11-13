"""Analyse an Alfvénic-inflow run.

python examples/alfven_inflow/analyze.py output/alfven_inflow/gate [--movie]

Without a wall (the gate): the wave that the flow carries through the box against the incident one
of the fluctuation table (theory.incident), at the last output. Per x row (the plane z = 0 of the
slice output) the complex amplitude of e^{i k y} of B_x and B_z, a(x), and of the incident field, b.
- gate (README, docs/plans/2026-11-13-upstream-state-restart.md, fixed before the first run):
  |a - b| / mean|b| <= 12 % on every row with x > 10 c/omega_pe (the outflow face x- carries no
  incident field) and <= 4 % on average over the rows 10 < x < lx - 4; the densities of electrons and
  ions in bins of 2 c/omega_pe from x = 10 within 10 % of 1; Gauss's law <= 1e-10; div B <= 1e-12.
- reported: the mean amplitude ratio |a| / |b| and phase a / b of each component.

With a wall (the piston): reported quantities, not a gate (the verdict is the conservation contracts
Gauss <= 1e-10 and div B <= 1e-12).
- shock position x_sh(t): the outermost x where the ion density, averaged over d_i / 2, exceeds
  (1 + r) / 2 (r: the compression of the jump conditions); its speed from a fit over the second half
  of the run; the downstream density and B_y between 0.2 and 0.8 x_sh at the last profile.
- the wave ahead of and behind the shock: mean |a| / |b| of B_x and B_z over the rows from x_sh + 4
  to x_sh + 20 c/omega_pe (ahead; the foot of the reflected ions lies there), from x_sh + 20 to
  lx - 4 (far upstream) and from 0.2 x_sh to x_sh - 4 (behind), at every slice output, their means
  over the outputs with both ahead and behind measured, and behind / ahead (transmission). A
  frozen-in field gives B_x continuous (normal to the shock) and B_z compressed by r (tangential).
  b is the incident wave of the planner's table also in the null run without the wave (``wave:
  false``), whose levels are then the shock's own k-mode in units of the wave.
- the downstream density (the incoming flow has density 1; the density far upstream is reported),
  B_y / B0 there, the downstream temperatures T_e, T_i (2/3 of the mean kinetic energy per particle
  in the frame of the ions) against the total T_e + T_i of the jump conditions, and the electron
  energy spectra downstream and far upstream.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import theory as T  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402
from pypic.analysis.fields import fit_velocity, front_position  # noqa: E402

GAUSS_MAX, DIVB_MAX = 1e-10, 1e-12
# gate
X_OUT = 10.0  # rows beyond this distance from the outflow face (c/omega_pe)
FACE_MARGIN = 4.0  # the mean leaves out the rows this close to the inflow face
ROW_MAX, MEAN_MAX, DENSITY_TOL = 0.12, 0.04, 0.10
BIN = 2.0  # density bins (c/omega_pe)
COMPS = ("Bx", "Bz")
# piston
GAP, AHEAD = 4.0, 16.0  # rows from x_sh + GAP to x_sh + GAP + AHEAD (ahead), 0.2 x_sh to x_sh - GAP (behind)
DOWN = (0.2, 0.8)  # downstream part of [0, x_sh] for the density, B_y and the temperatures


def _finite_max(a) -> float:
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    return float(a.max()) if a.size else float("nan")


def row_modes(F, y, k) -> np.ndarray:
    """Complex amplitude of e^{i k y} on every row of F[x, y] (y spans whole wavelengths)."""
    return 2.0 * (np.asarray(F, dtype=float) * np.exp(-1j * k * np.asarray(y))[None, :]).mean(axis=1)


def plane_modes(run, info, j, comp):
    """x, the row modes of the slice output j of ``comp`` and those of the incident field at its time."""
    s = run.series("plane")
    x, y = s.plane_axes(comp)
    F = s.read(comp, j)
    cfg = run.config
    inc = T.incident(info["table"], comp, x[:, None], y[None, :], 0.0, float(s.t[j]), cfg["E0"], cfg["B0"])
    return x, row_modes(F, y, info["k"]), row_modes(inc, y, info["k"])


def binned_density(x, n, lo, hi, width=BIN):
    """Mean of a lineout over bins of ``width`` from lo to hi: (bin centres, means)."""
    edges = np.arange(lo, hi + 1e-9, width)
    idx = np.digitize(x, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    cnt = np.bincount(idx[ok], minlength=len(edges) - 1)
    tot = np.bincount(idx[ok], weights=np.asarray(n, dtype=float)[ok], minlength=len(edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 0.5 * (edges[1:] + edges[:-1]), tot / cnt


def gate_measures(x, a, b, lx) -> dict:
    """Row errors |a - b| / mean|b| beyond X_OUT: the largest, and the mean of the rows short of
    lx - FACE_MARGIN; nan without rows."""
    sel = x > X_OUT
    if not sel.any():
        return {"row_max": float("nan"), "mean": float("nan"), "ratio": float("nan"), "phase": float("nan")}
    ref = np.abs(b[sel]).mean()
    err = np.abs(a - b)[sel] / ref
    inner = x[sel] < lx - FACE_MARGIN
    return {"row_max": float(err.max()), "mean": float(err[inner].mean()) if inner.any() else float("nan"),
            "ratio": float(np.abs(a[sel]).mean() / ref), "phase": float(np.angle(np.mean(a[sel] / b[sel]))),
            "row_max_x": float(x[sel][int(np.argmax(err))])}


def gate_checks(res: dict) -> dict:
    """The criteria on the measured values; a nan fails."""
    return {
        "gauss": bool(res["gauss_max"] <= GAUSS_MAX),
        "divB": bool(res["divB_max"] <= DIVB_MAX),
        "row_error": bool(all(res[f"{c}_row_max"] <= ROW_MAX for c in COMPS)),
        "mean_error": bool(all(res[f"{c}_mean"] <= MEAN_MAX for c in COMPS)),
        "density": bool(res["density_deviation"] <= DENSITY_TOL),
    }


def _stats(run) -> dict:
    s = run.scalars
    out = {"gauss_max": _finite_max(s["gauss"]), "divB_max": _finite_max(s["divB"]),
           "particles_start": float(s["n_particles"][0]), "particles_end": float(s["n_particles"][-1]),
           "fill_max": float(np.max(s["fill"]))}
    stats = run.meta.get("run")
    if stats:
        out.update(ms_per_step=1e3 * stats["wall_time_s"] / stats["steps"],
                   ns_per_particle_step=stats["ns_per_particle_step"])
    return out


def analyze_open(run, info, run_dir, show_plot, fig_dir) -> dict:
    s = run.series("plane")
    j = len(s.t) - 1
    res = {"t": float(s.t[j]), "k_di": info["k_di"], "omega": info["omega"], "omega_over_kva": info["omega_over_kva"],
           "v0": info["v0"], "v_a": info["v_a"], "flow_mach": info["params"]["flow_mach"],
           "transit_fraction": float(s.t[j]) / info["transit_time"]}
    modes = {}
    for comp in COMPS:
        x, a, b = plane_modes(run, info, j, comp)
        modes[comp] = (x, a, b)
        for k, v in gate_measures(x, a, b, info["lx"]).items():
            res[f"{comp}_{k}"] = v
    lo = run.series("profile")
    dens = {}
    for sp in ("e", "i"):
        xc, n = binned_density(lo.x(f"n:{sp}"), lo[f"n:{sp}"][-1], X_OUT, info["lx"])
        dens[sp] = (xc, n)
        # np.max keeps a nan: an empty bin is no measurement
        res[f"density_{sp}_deviation"] = float(np.max(np.abs(n - 1.0)))
        res[f"density_{sp}_min"], res[f"density_{sp}_max"] = float(np.min(n)), float(np.max(n))
    res["density_deviation"] = max(res["density_e_deviation"], res["density_i_deviation"])
    res.update(_stats(run))
    checks = gate_checks(res)
    res.update(checks=checks, passed=all(checks.values()))
    if show_plot:
        res["plot"] = _figure_open(run, info, modes, dens, res, os.path.join(fig_dir or run_dir, "alfven_inflow.png"))
    return res


def shock_history(run, info):
    """(t, x_sh) from the ion-density lineouts, the level (1 + r) / 2 of the jump conditions."""
    lo = run.series("profile")
    x = lo.x("n:i")
    w = max(1, int(round(0.5 * info["d_i"] / run.dx[0])))
    level = 0.5 * (1.0 + info["shock"]["compression"])
    n = lo["n:i"]
    xs = np.array([front_position(x, np.convolve(n[i], np.ones(w) / w, "same"), level) for i in range(len(lo.t))])
    return lo.t, xs


def wave_across(run, info, t_sh, x_sh):
    """Mean |a| / |b| of B_x and B_z ahead of the shock, far upstream and behind it at every slice
    output (nan where a window holds no row)."""
    s = run.series("plane")
    out = {c: {"ahead": [], "far": [], "behind": []} for c in COMPS}
    for j, t in enumerate(s.t):
        xs = float(np.interp(t, t_sh, x_sh)) if np.isfinite(x_sh).any() else np.nan
        for comp in COMPS:
            x, a, b = plane_modes(run, info, j, comp)
            ref = np.abs(b).mean()
            ahead = (x > xs + GAP) & (x < xs + GAP + AHEAD)
            far = (x > xs + GAP + AHEAD) & (x < run.L[0] - GAP)
            behind = (x > DOWN[0] * xs) & (x < xs - GAP)
            for name, m in (("ahead", ahead), ("far", far), ("behind", behind)):
                out[comp][name].append(float(np.abs(a[m]).mean() / ref) if m.any() else np.nan)
    return s.t, out


def spectra(run, x_sh):
    """Electron energy spectra dN/dE (per bin width) downstream (DOWN of [0, x_sh]) and far upstream
    (x_sh + GAP + AHEAD ... lx - GAP) at the last histogram output: (energy centres, down, up)."""
    h = run.series("spectrum")
    ex, ee = h.edges
    xc = 0.5 * (ex[1:] + ex[:-1])
    H = np.asarray(h.read("e", len(h.t) - 1), dtype=float)
    lx = run.L[0]
    down = (xc > DOWN[0] * x_sh) & (xc < DOWN[1] * x_sh)
    up = (xc > x_sh + GAP + AHEAD) & (xc < lx - GAP)
    de = np.diff(ee)
    ec = np.sqrt(ee[1:] * ee[:-1])

    def norm(m):  # normalised to an integral of 1 over the energy
        f = H[m].sum(axis=0) / de
        return f / max(float((f * de).sum()), 1e-300) if m.any() else np.full(len(ec), np.nan)

    return ec, norm(down), norm(up), float(h.t[-1])


def temperatures(run, x_sh) -> dict:
    """Downstream T_e and T_i: 2/3 of the mean kinetic energy per particle; the ions' mean flow is
    taken away from both (non-relativistic, the electrons follow the ions)."""
    lo = run.series("profile")
    x = lo.x("n:i")
    m = (x > DOWN[0] * x_sh) & (x < DOWN[1] * x_sh)
    out = {}
    if not m.any():
        return {"Te_down": np.nan, "Ti_down": np.nan, "flow_down": np.nan}
    ni, ne = lo["n:i"][-1][m].sum(), lo["n:e"][-1][m].sum()
    v = lo["fluxx:i"][-1][m].sum() / ni
    mi = run.species_info["i"]["mass"]
    out["flow_down"] = float(v)
    out["Te_down"] = float(2.0 / 3.0 * (lo["ekin:e"][-1][m].sum() / ne - 0.5 * v * v))
    out["Ti_down"] = float(2.0 / 3.0 * (lo["ekin:i"][-1][m].sum() / ni - 0.5 * mi * v * v))
    return out


def analyze_piston(run, info, run_dir, show_plot, fig_dir) -> dict:
    sh = info["shock"]
    t, xs = shock_history(run, info)
    ok = np.isfinite(xs) & (t >= 0.5 * t[-1])
    v_sh = fit_velocity(t[ok], xs[ok])[0] if ok.sum() >= 3 else float("nan")
    x_end = float(xs[-1])
    lo = run.series("profile")
    x = lo.x("n:i")
    down = (x > DOWN[0] * x_end) & (x < DOWN[1] * x_end)
    up = (x > x_end + GAP + AHEAD) & (x < info["lx"] - GAP)
    xb = lo.x("By")
    downb = (xb > DOWN[0] * x_end) & (xb < DOWN[1] * x_end)
    n2 = float(np.mean(lo["n:i"][-1][down])) if down.any() else float("nan")
    n1 = float(np.mean(lo["n:i"][-1][up])) if up.any() else float("nan")
    b2 = float(np.mean(lo["By"][-1][downb]) / info["b0"]) if downb.any() else float("nan")
    tw, wave = wave_across(run, info, t, xs)
    res = {"t_end": float(t[-1]), "omega_ci_t_end": info["omega_ci_t_end"], "v0": info["v0"], "v_a": info["v_a"],
           "flow_mach": info["params"]["flow_mach"], "k_di": info["k_di"], "wave": info.get("wave", True),
           "shock_speed": v_sh, "shock_speed_theory": sh["v_shock"], "shock_speed_ratio": v_sh / sh["v_shock"],
           "shock_position_end": x_end, "mach_alfven_shock": sh["mach_alfven"],
           "mach_alfven_measured": (info["v0"] + v_sh) / info["v_a"],
           # the incoming flow has the density 1 (the load and the injection); n1 is the measured one
           "compression": n2, "density_upstream": n1, "compression_theory": sh["compression"],
           "field_compression": b2, "temperature_theory": sh["temperature"]}
    res.update(temperatures(run, x_end))
    res["temperature_ratio"] = (res["Te_down"] + res["Ti_down"]) / sh["temperature"]
    for comp in COMPS:
        w = {k: np.asarray(v) for k, v in wave[comp].items()}
        both = np.isfinite(w["ahead"]) & np.isfinite(w["behind"])
        for k, v in w.items():
            res[f"{comp}_{k}_end"] = float(v[-1])
            res[f"{comp}_{k}_mean"] = float(np.nanmean(v[both])) if np.isfinite(v[both]).any() else float("nan")
        res[f"{comp}_transmission"] = res[f"{comp}_behind_mean"] / res[f"{comp}_ahead_mean"]
        res[f"{comp}_samples"] = int(both.sum())
    ec, fd, fu, t_sp = spectra(run, x_end)
    res.update(spectrum_t=t_sp, Te_up=float(info["params"]["theta_e"]))
    res.update(_stats(run))
    checks = {"gauss": bool(res["gauss_max"] <= GAUSS_MAX), "divB": bool(res["divB_max"] <= DIVB_MAX)}
    res.update(checks=checks, passed=all(checks.values()))
    if show_plot:
        res["plot"] = _figure_piston(run, info, t, xs, tw, wave, (ec, fd, fu, t_sp), res,
                                     os.path.join(fig_dir or run_dir, "alfven_piston.png"))
    return res


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    if info["wall"]:
        return analyze_piston(run, info, run_dir, show_plot, fig_dir)
    return analyze_open(run, info, run_dir, show_plot, fig_dir)


# ------------------------------------------------------------------------------------------ figures

def _bz_map(ax, run, info, j, title):
    s = run.series("plane")
    P.field_map(ax, np.asarray(s.read("Bz", j), dtype=float) / info["b0"], s.extent("Bz"), signed=True,
                label="B_z / B0", xlabel="x [c/ωpe]", ylabel="y [c/ωpe]")
    ax.set_title(f"{title}, t = {s.t[j]:.0f}")
    ax.set_aspect("equal")


def _figure_open(run, info, modes, dens, res, path):
    P.style()
    fig, ax = P.plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw={"width_ratios": [1.5, 1, 1]})
    _bz_map(ax[0], run, info, len(run.series("plane").t) - 1, "B_z of the wave carried in from x+")
    g = ax[1]
    for k, comp in enumerate(COMPS):
        x, a, b = modes[comp]
        ref = np.abs(b).mean()
        g.plot(x, np.abs(a) / ref, color=P.color(k), label=f"|{comp} mode| / incident")
        g.plot(x, np.abs(a - b) / ref, color=P.color(k), ls=":", lw=1.0, label=f"|{comp} − incident| / incident")
    g.axhline(1.0, **P.THEORY, label="incident")
    g.axhspan(0, ROW_MAX, color=P.INK["grid"], alpha=0.6, lw=0, label=f"gate: row error ≤ {ROW_MAX:.0%}")
    g.axvline(X_OUT, color=P.INK["muted"], lw=0.8)
    g.set_ylim(0, 1.3)
    g.set_xlabel("x [c/ωpe]")
    g.set_title("k-mode of each row against the incident wave")
    g.legend(loc="lower right", fontsize=7)
    g = ax[2]
    for k, sp in enumerate(("e", "i")):
        xc, n = dens[sp]
        g.plot(xc, n, "o-", ms=3, color=P.color(k), label=f"n_{sp} (bins of {BIN:g} c/ωpe)")
    g.axhspan(1 - DENSITY_TOL, 1 + DENSITY_TOL, color=P.INK["grid"], alpha=0.6, lw=0, label="gate: within 10 %")
    g.set_ylim(0.8, 1.2)
    g.set_xlabel("x [c/ωpe]")
    g.set_title("densities")
    g.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    return P.save(fig, path)


def _figure_piston(run, info, t, xs, tw, wave, spec, res, path):
    P.style()
    fig, ax = P.plt.subplots(2, 2, figsize=(14, 9))
    s = run.series("plane")
    j = len(s.t) - 1
    _bz_map(ax[0, 0], run, info, j, "B_z: the wave ahead of and behind the shock")
    ax[0, 0].axvline(float(np.interp(s.t[j], t, xs)), color=P.INK["primary"], ls="--", lw=1.0)
    lo = run.series("profile")
    g = ax[0, 1]
    n = lo["n:i"]
    P.field_map(g, n.T, [lo.x("n:i")[0], lo.x("n:i")[-1], lo.t[0], lo.t[-1]], signed=False, label="n_i / n0",
                xlabel="x [c/ωpe]", ylabel="t [1/ωpe]")
    g.set_aspect("auto")
    g.plot(xs, t, color=P.INK["primary"], lw=1.0, label="shock (ion density)")
    ok = np.isfinite(xs) & (t >= 0.5 * t[-1])
    if ok.any():
        t0, x0 = t[ok][0], xs[ok][0]
        g.plot(x0 + res["shock_speed_theory"] * (t - t0), t, **P.THEORY, label="jump conditions")
    g.set_xlim(lo.x("n:i")[0], lo.x("n:i")[-1])
    g.set_title(f"ion density; shock speed {res['shock_speed']:.4f} c (theory {res['shock_speed_theory']:.4f} c)")
    g.legend(loc="upper right", fontsize=7)
    g = ax[1, 0]
    for k, comp in enumerate(COMPS):
        g.plot(tw, wave[comp]["ahead"], color=P.color(k), label=f"|{comp}| ahead of the shock")
        g.plot(tw, wave[comp]["far"], color=P.color(k), ls=":", lw=1.0, label=f"|{comp}| far upstream")
        g.plot(tw, wave[comp]["behind"], color=P.color(k), ls="--", label=f"|{comp}| behind the shock")
    g.axhline(1.0, color=P.INK["muted"], lw=0.8)
    g.axhline(res["compression_theory"], color=P.INK["muted"], lw=0.8, ls="-.", label="frozen-in B_z: r")
    g.set_xlabel("t [1/ωpe]")
    g.set_ylabel("mean |k-mode| / incident wave")
    g.set_title("the wave across the shock" + ("" if info.get("wave", True) else " (null run: no wave)"))
    g.legend(loc="upper left", fontsize=7)
    g = ax[1, 1]
    ec, fd, fu, t_sp = spec
    g.loglog(ec, fd, color=P.color(0), label="downstream")
    g.loglog(ec, fu, color=P.color(1), label="far upstream")
    T_ = res["Te_down"]
    if np.isfinite(T_) and T_ > 0:
        g.loglog(ec, 2.0 * np.sqrt(ec / np.pi) * T_ ** -1.5 * np.exp(-ec / T_), color=P.color(0), ls=":", lw=1.0,
                 label=f"Maxwellian of the downstream T_e = {T_:.2e}")
    top = np.nanmax(np.concatenate([fd, fu]))
    g.set_ylim(1e-4 * top, 3 * top)
    g.set_xlim(ec[0], ec[-1])
    g.set_xlabel("electron kinetic energy [m_e c²]")
    g.set_ylabel("dN/dE (normalised)")
    g.set_title(f"electron spectra, t = {t_sp:.0f}")
    g.legend(fontsize=7)
    fig.tight_layout()
    return P.save(fig, path)


def movie(run_dir: str, out: str | None = None, fps: int = 12, n_frames: int = 100, dpi: int = 110) -> str:
    """B_z and the ion density in the plane z = 0, and the density lineouts."""
    from pypic.analysis import Movie
    from pypic.analysis.animate import LineoutPanel, SlicePanel

    run = open_run(run_dir)
    info = run.meta["planner"]
    m = Movie(run, title="Alfvénic wave carried into a piston shock" if info["wall"] else
              "Alfvénic wave carried in by a flow", ncols=3, panel_size=(5.4, 3.6))
    a = 1.5 * info["amplitude"]
    m.add(SlicePanel(run, "plane", "Bz", scale=1.0 / info["b0"], vlim=(-a, a), label="B_z / B0"))
    top = 1.15 * info["shock"]["compression"] + 0.5 if info["wall"] else 1.5
    m.add(SlicePanel(run, "plane", "n:i", signed=False, vlim=(0.0, top), label="n_i / n0"))
    m.add(LineoutPanel(run, "profile", ["n:e", "n:i"], ylim=(0.0, top)))
    return m.save(out or os.path.join(run_dir, "alfven_inflow.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
