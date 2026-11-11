"""Analyse a bow-shock run: upstream state, sheath, magnetopause and shock on the subsolar line.

python examples/bow_shock/analyze.py output/bow_shock/gate [--movie]
gate (definitions in docs/plans/2026-11-10-bow-shock.md, lengths in units of the planned magnetopause
distance R_mp): Gauss law <= 1e-10; upstream density and u_x of both species within 5 % of the wind
over the last crossing; ion density compression >= 2 in the sheath; magnetopause (where the density
falls back to the wind density) within 25 % of R_mp; shock upstream of the magnetopause.
reported: MHD compression and the stand-off of Spreiter et al. (1966) against the measured ones, the
magnetopause by magnetic pressure, the charge of the body, absorbed particles per step, largest fill.
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
from pypic.analysis.fields import front_position  # noqa: E402

GAUSS_MAX, UPSTREAM_TOL, COMPRESSION_MIN, MAGNETOPAUSE_TOL = 1e-10, 0.05, 2.0, 0.25
BAND, SMOOTH = 0.1, 0.1  # half width across and running mean along the subsolar line, in R_mp


def _running_mean(a, w):
    return np.convolve(np.asarray(a, dtype=float), np.ones(w) / w, "valid")


def _fall_to(x, p, i0, level):
    """First point after index i0 where p has fallen to level (linear interpolation); nan if it never does."""
    below = np.nonzero(p[i0:] <= level)[0]
    if below.size == 0:
        return float("nan")
    i = i0 + int(below[0])
    if i == i0:
        return float(x[i])
    return float(x[i - 1] + (level - p[i - 1]) * (x[i] - x[i - 1]) / (p[i] - p[i - 1]))


def _mean(a) -> float:
    """Mean of the blocks of a region; nan for a region without a block (nothing was measured)."""
    return float(a.astype(float).mean()) if a.size else float("nan")


def upstream_state(run, info) -> dict:
    """Largest deviation of the density (against t = 0) and of u_x (against the wind) of each species over the
    outputs of the last crossing, in the blocks of the coarse output that lie inside the upstream region.
    nan where there is no measurement: no block lies inside the region (``blocks`` is 0: the blocks are
    too large for it) or a value is not finite. A nan deviation fails the gate."""
    f = run.fields("coarse")
    its = f.iterations
    c = next(d for d in run.config["output"]["diagnostics"] if d["name"] == "coarse")["coarsen"]
    r_mp, v = info["r_mp"], info["v_sw"]
    sel = T.upstream_blocks(run.n_cells, c, run.dx, r_mp)
    region = np.ix_(*sel)
    last = [it for it in its if it >= its[-1] - info["crossing_steps"]]
    out = {"blocks": int(np.prod([len(s) for s in sel])), "outputs": len(last)}
    for sp in "ei":
        n0 = _mean(f.get(its[0], f"n:{sp}").data[region])
        with np.errstate(divide="ignore", invalid="ignore"):  # an empty region or a zero density gives nan
            n = np.array([_mean(f.get(it, f"n:{sp}").data[region]) for it in last]) / n0
            ux = np.array([_mean(f.get(it, f"fluxx:{sp}").data[region]) for it in last]) / (n * n0 * v)
        # np.max keeps a nan (the builtin max drops it, which reported a missing measurement as no deviation)
        out[sp] = {"density": float(np.max(np.abs(n - 1.0))), "u_x": float(np.max(np.abs(ux - 1.0))),
                   "n_end": float(n[-1]), "u_x_end": float(ux[-1])}
    out["deviation"] = float(np.max([[out[sp]["density"], out[sp]["u_x"]] for sp in "ei"]))
    return out


def subsolar(run, info) -> dict:
    """Profiles on the subsolar line from the slice z = z_c: band average across, running mean along x,
    mean over the outputs of the last half crossing. Density in units of the wind density n_0 (the
    slice at t = 0 upstream of x_c - R_mp), magnetic pressure of the total field in units of the ram
    pressure n_0 m_i v_sw^2."""
    sz = run.series("slice_z")
    r_mp, (xc, yc, _) = info["r_mp"], info["center"]
    x, y = sz.plane_axes("n:i")
    band = np.abs(y - yc) <= BAND * r_mp + 1e-9
    w = max(1, int(round(SMOOTH * r_mp / (x[1] - x[0]))))
    step = sz.step
    late = np.nonzero(step >= step[-1] - info["crossing_steps"] / 2)[0]
    sel = slice(int(late[0]), int(late[-1]) + 1)
    n0 = float(sz.read("n:i", 0)[x <= xc - r_mp].astype(float).mean())
    xs = _running_mean(x, w)
    ni = _running_mean(sz.read("n:i", sel)[:, :, band].astype(float).mean(axis=(0, 2)), w) / n0
    b2 = np.zeros_like(xs)
    for comp in ("Bx", "By", "Bz"):  # staggered by half a cell: brought to the positions of the density
        xa, ya = sz.plane_axes(comp)
        bb = np.abs(ya - yc) <= BAND * r_mp + 1e-9
        prof = _running_mean(sz.read(comp, sel)[:, :, bb].astype(float).mean(axis=(0, 2)), w)
        b2 += np.interp(xs, _running_mean(xa, w), prof) ** 2
    return {"x": xs, "n_i": ni, "pb": 0.5 * b2 / info["ram_pressure"], "n0": n0, "outputs": len(late),
            "band_samples": int(band.sum()), "smooth_samples": w, "sel": sel}


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    r_mp, xc, r_body = info["r_mp"], info["center"][0], info["body_radius"]
    s = run.scalars
    g = s["gauss"][np.isfinite(s["gauss"])]
    res = {"gauss_max": float(g.max()) if g.size else float("nan"), "r_mp": r_mp}
    up = upstream_state(run, info)
    res.update(upstream_deviation=up["deviation"], upstream={sp: up[sp] for sp in "ei"})
    # subsolar line, upstream of the surface of the body
    line = subsolar(run, info)
    front = line["x"] < xc - r_body
    x, p, pb = line["x"][front], line["n_i"][front], line["pb"][front]
    i_max = int(np.argmax(p))
    x_mp = _fall_to(x, p, i_max, 1.0)
    x_sh = front_position(x, p, COMPRESSION_MIN, from_right=False)
    x_mp_b = front_position(x, pb, 1.0, from_right=False)
    x_th = info["compression_mhd"]
    res.update(
        compression=float(p[i_max]), x_compression=float(x[i_max]), compression_mhd=x_th,
        x_magnetopause=x_mp, magnetopause=(xc - x_mp) / r_mp, x_shock=x_sh, shock=(xc - x_sh) / r_mp,
        shock_spreiter=1.0 + 1.1 / x_th, sheath=(x_mp - x_sh) / (xc - x_mp), sheath_spreiter=1.1 / x_th,
        magnetopause_magnetic=(xc - x_mp_b) / r_mp, wind_density=line["n0"])
    step = s["step"]
    n_abs = s["n_absorbed"].astype(float)
    cross = info["crossing_steps"]
    n_cross = max(1, int(round(step[-1] / cross)))  # the last record lies up to one interval before the end
    at = [min(int(np.searchsorted(step, k * cross)), len(step) - 1) for k in range(n_cross + 1)]
    res.update(body_charge=float(s["absorbed_charge"][-1]),
               absorbed_per_step=[float((n_abs[b] - n_abs[a]) / (step[b] - step[a])) for a, b in zip(at, at[1:])],
               fill_max=float(np.max(s["fill"])), particles_start=float(s["n_particles"][0]),
               particles_end=float(s["n_particles"][-1]))
    checks = {
        "gauss": bool(res["gauss_max"] <= GAUSS_MAX),
        "upstream": bool(up["deviation"] <= UPSTREAM_TOL),  # False for nan: no measurement, no pass
        "compression": bool(res["compression"] >= COMPRESSION_MIN),
        "magnetopause": bool(np.isfinite(x_mp) and abs(res["magnetopause"] - 1.0) <= MAGNETOPAUSE_TOL),
        "shock_upstream": bool(np.isfinite(x_sh) and np.isfinite(x_mp) and x_sh < x_mp),
    }
    res.update(checks=checks, passed=all(checks.values()))
    if show_plot:
        res["plot"] = _figure(run, info, line, res, os.path.join(fig_dir or run_dir, "bow_shock.png"))
    return res


def _mark_body(ax, info, a, b):
    """Body (filled) and the planned magnetopause distance (dashed circle) in the plane of axes a, b."""
    c = (info["center"][a], info["center"][b])
    ax.add_patch(P.plt.Circle(c, info["body_radius"], fc=P.INK["muted"], ec="none"))
    ax.add_patch(P.plt.Circle(c, info["r_mp"], fc="none", ec=P.INK["primary"], ls="--", lw=0.8))
    ax.set_aspect("equal")


def _figure(run, info, line, res, path):
    P.style()
    xc = info["center"][0]
    fig, ax = P.plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for k, (name, a, b) in enumerate((("slice_z", 0, 1), ("slice_y", 0, 2))):
        sl = run.series(name)
        d = sl.read("n:i", line["sel"]).astype(float).mean(axis=0) / line["n0"]
        P.field_map(ax[k], d, sl.extent("n:i"), signed=False, vlim=(0.0, max(4.0, 1.1 * res["compression"])),
                    label="n_i / n_0", xlabel="x [c/ωpe]", ylabel=f"{'xyz'[b]} [c/ωpe]")
        _mark_body(ax[k], info, a, b)
        ax[k].set_title(f"ion density, plane {'xyz'[3 - a - b]} = centre\n(mean of the last half crossing)")
    yc = info["center"][1]
    for xm, i in ((res["x_shock"], 1), (res["x_magnetopause"], 3)):
        if np.isfinite(xm):
            ax[0].plot([xm, xm], [yc - 0.5 * info["r_mp"], yc + 0.5 * info["r_mp"]], color=P.color(i), lw=1.2)
    x = (xc - line["x"]) / info["r_mp"]
    a = ax[2]
    a.plot(x, line["n_i"], color=P.color(0), label="n_i / n_0")
    a.plot(x, line["pb"], color=P.color(2), label="B² / 2 over the ram pressure")
    a.axhline(res["compression_mhd"], **P.THEORY,
              label=f"perpendicular MHD shock, n2/n1 = {res['compression_mhd']:.2f}")
    a.axvline(res["shock_spreiter"], color=P.INK["primary"], ls=":", lw=1.0, label="shock (Spreiter et al. 1966)")
    a.axvline(1.0, color=P.INK["primary"], ls="--", lw=0.8, label="magnetopause (pressure balance)")
    for xm, lab, i in ((res["shock"], "shock", 1), (res["magnetopause"], "magnetopause", 3)):
        if np.isfinite(xm):
            a.axvline(xm, color=P.color(i), lw=1.2, label=f"{lab}, measured: {xm:.2f} R_mp")
    a.axvspan(0.0, info["body_radius"] / info["r_mp"], color=P.INK["grid"], lw=0)
    a.set_xlim(min(2.4, xc / info["r_mp"]), 0.0)
    a.set_ylim(0.0, max(4.5, 1.2 * res["compression"]))
    a.set_xlabel("distance from the centre on the subsolar line [R_mp]")
    a.set_title("subsolar line")
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7)
    fig.tight_layout()
    return P.save(fig, path)


def summary(results: dict[str, dict], path: str) -> str:
    """Magnetopause distance, shock distance and largest density on the subsolar line of every run
    against the pressure balance, the stand-off of Spreiter et al. (1966) and the perpendicular MHD shock."""
    P.style()
    names = list(results)
    y = np.arange(len(names))[::-1]
    panels = (("magnetopause", None, "magnetopause distance / R_mp", "pressure balance"),
              ("shock", "shock_spreiter", "shock distance / Spreiter et al. (1966)", "Spreiter et al. (1966)"),
              ("compression", "compression_mhd", "largest density / perpendicular MHD shock", "MHD shock"))
    fig, axes = P.plt.subplots(1, 3, figsize=(14, 0.5 * len(names) + 1.8), sharey=True)
    for k, (ax, (key, th, lab, th_lab)) in enumerate(zip(axes, panels)):
        if k == 0:
            ax.axvspan(1 - MAGNETOPAUSE_TOL, 1 + MAGNETOPAUSE_TOL, color=P.INK["grid"], alpha=0.7, lw=0,
                       label="gate ±25 %")
        ax.axvline(1.0, **P.THEORY, label=th_lab)
        for yi, n in zip(y, names):
            r = results[n]
            ref = r[th] if th else 1.0
            ax.plot(r[key] / ref, yi, "o", ms=8, color=P.color(0), label=("PIC" if yi == y[0] else None))
            side = -1 if r[key] < ref else 1  # the label stays on the point's side of the reference line
            ax.annotate(f"{r[key]:.2f} vs {ref:.2f}", (r[key] / ref, yi), xytext=(8 * side, 8),
                        textcoords="offset points", ha="right" if side < 0 else "left", fontsize=7,
                        color=P.INK["secondary"])
        ax.set_xlim(0.5, 1.5)
        ax.set_xlabel(lab)
        ax.grid(axis="y", visible=False)
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=7)
    axes[0].set_yticks(y, names)
    axes[0].set_ylim(y.min() - 0.6, y.max() + 0.7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def movie(run_dir: str, out: str | None = None, fps: int = 12, n_frames: int = 100, dpi: int = 110) -> str:
    from pypic.analysis import Movie
    from pypic.analysis.animate import SlicePanel

    run = open_run(run_dir)
    info = run.meta["planner"]
    m = Movie(run, title="bow shock of a magnetised planet")
    for diag, comp, scale, label, log in (("slice_z", "n:i", 1.0, "n_i / n_0 (z = centre)", False),
                                          ("slice_y", "n:i", 1.0, "n_i / n_0 (y = centre)", False),
                                          ("slice_z", "Bz", 1.0 / info["b_imf"], "B_z / B_imf (z = centre)", True)):
        m.add(SlicePanel(run, diag, comp, signed=False, log=log, label=label, scale=scale,
                         blank=(info["center"], info["body_radius"])))
    return m.save(out or os.path.join(run_dir, "bow_shock.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
