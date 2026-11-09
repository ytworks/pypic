"""Analyse a plasma-wave run: omega-k spectra of the lineouts vs linear theory.

python examples/plasma_waves/analyze.py output/plasma_waves/unmagnetized [--movie]
    -> omega_k.png (+ movie) in the run directory; verdict JSON on stdout

Every theory branch is located in the measured spectrum by a peak search within +-15 % of the
theoretical frequency at each k; the verdict compares the peak with theory where the branch is
resolved, with tolerance max(3 %, 2 d_omega / omega), d_omega = 2 pi / t_run.
"""
from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np

from pypic.analysis import open_run
from pypic.analysis import plot as P
from pypic.analysis.spectral import branch_ridge, omega_k
from pypic.theory import waves as W


def branches(info: dict, k: np.ndarray, dx: float, dt: float) -> dict:
    """{component: [(label, k, omega_theory, checked[, (w_lo, w_hi)]), ...]} for the run's geometry
    (the optional bounds keep the peak search inside a Bernstein band)."""
    geo = info["geometry"]
    vt = math.sqrt(info["theta_e"])
    kpos = k[k > 0]
    out = {}
    if geo == "ion_acoustic":
        kia = kpos[kpos * vt <= 1.0]
        out["Ex"] = [("ion acoustic (kinetic)", kia, W.ion_acoustic(kia, info["theta_e"], info["Ti"],
                                                                     info["mass_ratio"]).real, True)]
        return out
    kl = kpos[kpos * vt <= 0.45]  # weakly damped Langmuir waves only
    lang = ("Langmuir (kinetic)", kl, W.langmuir(kl, vt).real, True)
    em = ("EM (Yee)", kpos, W.em_wave_yee(kpos, dx, dt), True)
    if geo == "unmagnetized":
        out["Ex"] = [lang]
        out["Ey"] = [em]
        out["Ez"] = [em]
    elif geo == "perpendicular":
        wce = info["wce"]
        cp_ = W.cold_perpendicular(kpos, 1.0, wce)
        out["Ez"] = [("O mode", kpos, cp_["O"], True)]
        out["Ey"] = [("X upper", kpos, cp_["X_upper"], True), ("X lower", kpos, cp_["X_lower"], False)]
        rho = vt / wce  # thermal Larmor radius: Bernstein structure for k rho ~ 0.3 .. 3
        kb = kpos[(kpos * rho >= 0.3) & (kpos >= 2.0) & (kpos <= info.get("k_bernstein_max", 40.0))]
        bern = W.bernstein(kb, vt, 1.0, wce, n_bands=4)
        # search strictly inside each band: the cyclotron harmonics n omega_ce carry strong lines
        # bands 1-3 are gated; band 4 is weakly excited thermally at these k (reported only)
        out["Ex"] = [(f"Bernstein {b + 1}", kb, bern[b], b < 3, ((b + 1.03) * wce, (b + 1.97) * wce))
                     for b in range(bern.shape[0])]
    elif geo == "parallel":
        wce = info["wce"]
        cpar = W.cold_parallel(kpos, 1.0, wce)
        ksep = kpos <= 2.0  # R and L windows disjoint (beyond, the +-15 % L window reaches the R ridge)
        wh = cpar["whistler"]
        brs = [("whistler", kpos, np.where(wh < 0.8 * wce, wh, np.nan), True),
               ("R", kpos[ksep], cpar["R_upper"][ksep], True), ("L", kpos[ksep], cpar["L"][ksep], True)]
        out["Ey"] = brs
        out["Ez"] = brs
        out["Ex"] = [lang]
    return out


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    lo = run.series("lineout")
    t = lo.t
    dt_out = float(np.median(np.diff(t)))
    t_run = t[-1] - t[0]
    dw = 2 * math.pi / t_run
    kwin = info["k_window"]
    res = {"geometry": info["geometry"], "d_omega": dw, "branches": {}}
    spectra = {}
    ok_all = True
    for comp in lo.components:
        k, w, Pw = omega_k(lo[comp], dt_out, run.dx[0])
        spectra[comp] = (k, w, Pw)
    br = branches(info, np.unique(np.abs(spectra["Ex"][0])), run.dx[0], run.dt)
    for comp, lst in br.items():
        k, w, Pw = spectra[comp]
        for label, kt, wt, checked, *bounds in lst:
            sel = ((kt <= kwin) | bool(bounds)) & np.isfinite(wt)  # Bernstein branches carry their own k range
            km, wm, wth = branch_ridge(k, w, Pw, kt[sel], wt[sel], bounds=bounds[0] if bounds else None)
            found = np.isfinite(wm)
            err = np.abs(wm[found] - wth[found]) / wth[found]
            tol = np.maximum(0.03, 2 * dw / wth[found])
            frac_ok = float(np.mean(err <= tol)) if found.any() else 0.0
            entry = {"component": comp, "points": int(len(km)), "found": int(found.sum()),
                     "median_rel_error": float(np.median(err)) if found.any() else float("nan"),
                     "max_rel_error": float(err.max()) if found.any() else float("nan"),
                     "fraction_within_tolerance": frac_ok, "checked": checked}
            if checked:
                entry["passed"] = bool(found.sum() >= max(3, 0.5 * len(km)) and frac_ok >= 0.9)
                ok_all &= entry["passed"]
            res["branches"][f"{comp}: {label}"] = entry
    res["passed"] = bool(ok_all)
    s = run.scalars
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    d = s["divB"][np.isfinite(s["divB"])]
    res["divB_max"] = float(d.max()) if d.size else float("nan")

    if show_plot:
        P.style()
        comps = list(br)
        fig, axes = P.plt.subplots(1, len(comps), figsize=(5.2 * len(comps), 4.2), squeeze=False)
        wmax = info["w_max"]
        for ax, comp in zip(axes[0], comps):
            k, w, Pw = spectra[comp]
            curves = [(b[0], np.concatenate([-b[1][::-1], b[1]]), np.concatenate([b[2][::-1], b[2]]))
                      for b in br[comp]]
            kmax = info.get("k_bernstein_max", kwin) if any(len(b) > 4 for b in br[comp]) else kwin
            P.omega_k_map(ax, k, w, Pw, curves=curves, k_max=kmax, w_max=min(wmax, 3.0) if kmax > kwin else wmax)
            ax.set_title(f"{comp}: {info['geometry']}")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "omega_k.png"))
    return res


def summary(results: dict[str, dict], path: str) -> str:
    """Median frequency error of every dispersion branch (all configurations) against the 3 % tolerance:
    filled markers are gated branches, hollow ones reported only."""
    P.style()
    rows = []
    for i, (name, res) in enumerate(results.items()):
        for label, b in res["branches"].items():
            rows.append((i, f"{name} · {label}", b["median_rel_error"], b["checked"], b["found"], b["points"]))
    fig, ax = P.plt.subplots(figsize=(9.5, 0.28 * len(rows) + 1.3))
    y = np.arange(len(rows))[::-1]
    ax.axvline(0.03, **P.THEORY, label="tolerance 3 %")
    seen = set()
    for yi, (i, _, err, checked, found, pts) in zip(y, rows):
        name = list(results)[i]
        lab = name if name not in seen else None
        seen.add(name)
        e = max(err, 1e-5) if np.isfinite(err) else np.nan
        ax.plot(e, yi, "o", ms=7, color=P.color(i), mfc=P.color(i) if checked else "none", label=lab)
        ax.annotate(f"{found}/{pts}", (e, yi), xytext=(8, 0), textcoords="offset points", va="center",
                    fontsize=7, color=P.INK["muted"])
    ax.plot([], [], "o", color=P.INK["secondary"], mfc="none", label="reported (not gated)")
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 0.3)
    ax.set_yticks(y, [r[1] for r in rows], fontsize=7.5)
    ax.set_xlabel("median |ω_PIC − ω_theory| / ω_theory   (label: k points found / searched)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def movie(run_dir: str, out: str | None = None, fps: int = 20, n_frames: int = 200,
          dpi: int = 110) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="thermal plasma waves")
    m.xt("lineout", "Ex").xt("lineout", "Ey").lineout("lineout", ["Ex", "Ey"]).history()
    return m.save(out or os.path.join(run_dir, "plasma_waves.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
