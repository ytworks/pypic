"""Analyse a wakefield run: the wake behind the driver against the linear theory.

python examples/wakefield/analyze.py output/wakefield/laser_linear
gate (laser): on the axis E_x = c A_th(x) cos(k x + phi) + b fitted behind the pulse, with A_th the
      theory at the place (local a0 and waist of the beam): |c - 1| <= 5 %, |k v_g - 1| <= 2 %; on
      the lines at r = waist / 2, with the same k and phi, E_x = c_x A_x cos + b and
      E_r = -c_r A_r sin + b': |c_r / c_x - 1| <= 10 %.
gate (bunch): on the axis E_x = c A_th cos(k (x - x_b) + phi) + b over five wavelengths behind the
      bunch: |c - 1| <= 5 %, |k beta - 1| <= 2 %; the amplitude on the lines at r = sigma_r over
      the one on the axis within 10 % of R(sigma_r) / R(0).
null run (drive: false): max |E_x| on the axis <= 1 % of the predicted amplitude, mean electron
      density on the plateau within 1 % of n0.
reported: the fitted phase, the amplitude at the output before the last, the place of the driver.
The places are those of the laboratory (OUT-D5). No interval of a fit or of the density comes
closer than 1 c/omega_pe to the lower x face of the box, and the density is taken up to
1 c/omega_pe behind the front of the plasma of a moving window.

regime "blowout" (no gate; everything is reported next to the scalings of Lu et al. 2007):
cavity: the largest connected region of the plane with electron density below 0.5 n0 between the
      lower face + 1 and the centre of the pulse (the density averaged along x over one spacing
      across); its radius is half of its largest extent across.
field: max and min of E_x on the axis behind the pulse; the zero at which E_x rises inside the
      cavity (the nearest to the widest place) and the slope of a line fitted within 1 c/omega_pe
      of it (a sphere without electrons: 1 / 2).
electrons: from the histogram over u_x, those with u_x >= 2 omega0 / omega_pe (twice the gamma of
      the wake): their number in n0 (c/omega_pe)^3, mean and largest u_x, the most populated bin;
      the number of electrons and the charge for a laser of 0.8 um.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy import constants, ndimage
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import plan as PL  # noqa: E402
import theory as T  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402
from pypic.analysis.particles import spectrum  # noqa: E402


def line(plane: np.ndarray, coords: np.ndarray, pos: float) -> np.ndarray:
    """The row of ``plane`` (first axis x, second axis the coordinate ``coords``) at ``pos``,
    linear between the two neighbouring samples."""
    f = (pos - coords[0]) / (coords[1] - coords[0])
    i = min(max(int(np.floor(f)), 0), len(coords) - 2)
    w = f - i
    return (1.0 - w) * plane[:, i] + w * plane[:, i + 1]


def _design(x, amp, k, phi=None):
    one = np.ones_like(x)
    if phi is None:
        return np.stack([amp * np.cos(k * x), amp * np.sin(k * x), one], axis=1)
    return np.stack([amp * np.cos(k * x + phi), one], axis=1)


def fit_wave(x, y, amp) -> dict:
    """Least squares of y = c amp(x) cos(k x + phi) + b with free c > 0, k, phi, b. For a given k
    the model is linear in (c cos phi, -c sin phi, b); k is searched within one bin around the
    peak of the spectrum of y."""
    x, y, amp = (np.asarray(v, dtype=float) for v in (x, y, amp))
    length = x[-1] - x[0]
    kk = 2.0 * np.pi * np.fft.rfftfreq(16 * len(x), d=(x[1] - x[0]))
    spec = np.abs(np.fft.rfft((y - y.mean()) * np.hanning(len(y)), 16 * len(x)))
    k0 = float(kk[1 + int(np.argmax(spec[1:]))])

    def residual(k):
        m = _design(x, amp, k)
        return float(np.sum((m @ np.linalg.lstsq(m, y, rcond=None)[0] - y) ** 2))

    dk = 2.0 * np.pi / length
    k = float(minimize_scalar(residual, bounds=(max(k0 - dk, 0.1 * k0), k0 + dk), method="bounded",
                              options={"xatol": 1e-10 * k0}).x)
    a, b, off = np.linalg.lstsq(_design(x, amp, k), y, rcond=None)[0]
    c = float(np.hypot(a, b))
    return {"c": c, "k": k, "phi": float(np.arctan2(-b, a)), "offset": float(off), "k_spectrum": k0,
            "rms_residual": float(np.sqrt(residual(k) / len(x)) / (c * np.mean(amp)))}


def fit_amplitude(x, y, amp, k: float, phi: float) -> float:
    """c of y = c amp(x) cos(k x + phi) + b for a given k and phi (c has a sign)."""
    return float(np.linalg.lstsq(_design(np.asarray(x, dtype=float), np.asarray(amp, dtype=float), k, phi),
                                 np.asarray(y, dtype=float), rcond=None)[0][0])


def _window(p, x_driver: float, x_low: float = 0.0):
    """(x_lo, x_hi) of the measurement behind the driver in a box that begins at x_low."""
    if p.driver == "bunch":
        hi = x_driver - 5.0 * p.sigma_x
        return max(p.start + 30.0, hi - 5.0 * 2.0 * np.pi * PL.speed(p), x_low + 1.0), hi
    return max(p.ramp[1], x_low) + 1.0, x_driver - 3.0 * p.duration


def _measure(p, s, i: int) -> dict:
    """Fits of the output i of the plane through the axis."""
    xe, ze = s.plane_axes("Ex", i)
    zc, speed = 0.5 * p.width, PL.speed(p)
    r = 0.5 * p.waist if p.driver == "laser" else p.sigma_r
    x_drv = PL.driver_position(p, float(s.t[i]))
    lo, hi = _window(p, x_drv, xe[0] - 0.5 * (xe[1] - xe[0]))
    sel = (xe >= lo) & (xe <= hi)
    x = xe[sel]
    origin = x_drv if p.driver == "bunch" else 0.0
    ex = s.read("Ex", i)
    on = fit_wave(x - origin, line(ex, ze, zc)[sel], PL.wake_amplitude(p, x)[0])
    out = {"t": float(s.t[i]), "x_driver": x_drv, "window": [float(lo), float(hi)], "amplitude_ratio": on["c"],
           "k": on["k"], "k_times_speed": on["k"] * speed, "phase": on["phi"], "offset": on["offset"],
           "rms_residual": on["rms_residual"], "k_spectrum": on["k_spectrum"], "fit": on, "x": x, "origin": origin}
    ax_r, ar_r = PL.wake_amplitude(p, x, r)
    cx = 0.5 * sum(fit_amplitude(x - origin, line(ex, ze, zc + sg * r)[sel], ax_r, on["k"], on["phi"])
                   for sg in (1.0, -1.0))
    out["side_amplitude_ratio"] = cx
    if p.driver == "bunch":
        # the amplitude at r = sigma_r over the one on the axis (both in units of their theory)
        out["R_ratio"] = cx * ax_r[0] / (on["c"] * PL.wake_amplitude(p, x)[0][0])
    else:
        xz, zz = s.plane_axes("Ez", i)
        selz = (xz >= lo) & (xz <= hi)
        ez = s.read("Ez", i)
        # E_r = -c_r A_r sin(k x + phi) + b' = c_r A_r cos(k x + phi + pi / 2) + b'
        cr = 0.5 * sum(sg * fit_amplitude(xz[selz], line(ez, zz, zc + sg * r)[selz],
                                          PL.wake_amplitude(p, xz[selz], r)[1], on["k"], on["phi"] + 0.5 * np.pi)
                       for sg in (1.0, -1.0))
        out["radial_amplitude_ratio"] = cr
        out["Er_over_Ex_ratio"] = cr / cx
    return out


def cavity(density, x, z, x_lo: float, x_hi: float, level: float = 0.5):
    """The largest connected region with density below ``level`` among the samples x_lo <= x <= x_hi of
    a plane (first axis x). The density is first averaged along x over one spacing across. None if
    there is no such sample; else the radius (half of the largest extent across), the x of that
    place (the middle of the columns that have it) and the first and last x of the region."""
    n = max(1, int(round((z[1] - z[0]) / (x[1] - x[0]))))
    sel = (x >= x_lo) & (x <= x_hi)
    if sel.sum() < 2:
        return None
    lab, count = ndimage.label(ndimage.uniform_filter1d(np.asarray(density, dtype=float), n, axis=0,
                                                        mode="nearest")[sel] < level)
    if count == 0:
        return None
    across = (lab == 1 + int(np.argmax(np.bincount(lab.ravel())[1:]))).sum(axis=1)
    cols, widest = np.nonzero(across)[0], np.nonzero(across == across.max())[0]
    xs = x[sel]
    return {"radius": float(0.5 * across.max() * (z[1] - z[0])), "x_widest": float(np.median(xs[widest])),
            "x_back": float(xs[cols[0]]), "x_front": float(xs[cols[-1]])}


def rising_zero(x, y, x_near: float, half: float = 1.0):
    """(place, slope) of the zero of y at which it rises that is nearest to x_near, from a line fitted
    to the samples within ``half`` of the crossing. NaN without such a crossing."""
    up = np.nonzero((y[:-1] < 0.0) & (y[1:] >= 0.0))[0]
    if not up.size:
        return float("nan"), float("nan")
    x0 = x[up[int(np.argmin(np.abs(x[up] - x_near)))]]
    sel = np.abs(x - x0) <= half
    slope, cut = np.polyfit(x[sel], y[sel], 1)
    return float(-cut / slope), float(slope)


def trapped(edges, counts, u_min: float) -> dict:
    """The particles of the bins at and above u_min of a histogram over u_x: number (sum of the
    weights), mean u_x, the centre of the most populated bin, the upper edge of the last occupied
    bin, and whether that is the last bin of the histogram (the range then cuts the spectrum)."""
    edges = np.asarray(edges, dtype=float)
    lo = int(np.searchsorted(edges, u_min, side="left"))
    c = np.asarray(counts, dtype=float)[lo:]
    number = float(c.sum())
    if number <= 0.0:
        return {"u_min": float(u_min), "number": 0.0, "u_mean": float("nan"), "u_peak": float("nan"),
                "u_max": float("nan"), "at_range_limit": False}
    centres = 0.5 * (edges[lo:-1] + edges[lo + 1:])
    top = int(np.nonzero(c)[0][-1])
    return {"u_min": float(u_min), "number": number, "u_mean": float(np.sum(c * centres) / number),
            "u_peak": float(centres[int(np.argmax(c))]), "u_max": float(edges[lo + top + 1]),
            "at_range_limit": top == len(c) - 1}


def _bubble(p, s, i: int, x_low: float | None = None) -> dict:
    """The pulse, the cavity behind it and E_x on the axis at the output i of the plane through the axis."""
    zc = 0.5 * p.width
    xy, zy = s.plane_axes("Ey", i)
    ey = line(s.read("Ey", i), zy, zc)
    out = {"t": float(s.t[i]), "a0": float(np.abs(ey).max() / p.omega0), "pulse_center": float("nan"),
           "cavity": None, "Ex_max": float("nan"), "Ex_min": float("nan"), "Ex_zero": float("nan"),
           "Ex_slope": float("nan")}
    if not np.any(ey):
        return out
    out["pulse_center"] = x_p = float(np.sum(xy * ey ** 2) / np.sum(ey ** 2))
    xn, zn = s.plane_axes("n:electron", i)
    if x_low is None:
        x_low = float(xn[0])
    x_lo = max(p.ramp[1], x_low) + 1.0
    out["cavity"] = cav = cavity(s.read("n:electron", i), xn, zn, x_lo, x_p)
    xe, ze = s.plane_axes("Ex", i)
    ex = line(s.read("Ex", i), ze, zc)
    behind = (xe >= x_lo) & (xe <= x_p)
    if behind.any():
        out["Ex_max"], out["Ex_min"] = float(ex[behind].max()), float(ex[behind].min())
        xb, eb = xe[behind], ex[behind]
        out["x_Ex_max"], out["x_Ex_min"] = float(xb[eb.argmax()]), float(xb[eb.argmin()])
    if cav is not None:
        inside = (xe >= cav["x_back"]) & (xe <= cav["x_front"])
        out["Ex_zero"], out["Ex_slope"] = rising_zero(xe[inside], ex[inside], cav["x_widest"])
    return out


def _blowout(run, p, info, show_plot: bool, fig_dir: str | None) -> dict:
    th = info["theory"]["blowout"]
    s, sp = run.series("axis_plane"), run.series("spectrum")
    last = len(s.t) - 1
    u_min = 2.0 * p.omega0
    edges = sp.edges[0]
    each = [_bubble(p, s, i) for i in range(last + 1)]
    beam = [trapped(edges, sp.read("electron", i), u_min) for i in range(last + 1)]
    b, tr = each[last], beam[last]
    g = run.scalars["gauss"]
    res = {"driver": p.driver, "drive": p.drive, "regime": p.regime, "t": b["t"], "theory": th,
           "gauss_max": float(np.nanmax(g)) if np.isfinite(g).any() else float("nan"),
           "pulse_center_measured": b["pulse_center"], "pulse_center_planned": PL.driver_position(p, b["t"]),
           "a0_measured": b["a0"], "cavity": b["cavity"],
           "cavity_radius": b["cavity"]["radius"] if b["cavity"] else float("nan"),
           "Ex_axis_max": b["Ex_max"], "Ex_axis_min": b["Ex_min"], "x_Ex_axis_max": b.get("x_Ex_max"),
           "x_Ex_axis_min": b.get("x_Ex_min"), "Ex_zero": b["Ex_zero"], "Ex_slope": b["Ex_slope"], "trapped": tr,
           "electrons_in_range": float(np.sum(sp.read("electron", last)))}
    res["cavity_radius_over_theory"] = res["cavity_radius"] / th["radius"]
    per_unit = T.electrons_per_unit(p.omega0)
    res["trapped_electrons_800nm"] = tr["number"] * per_unit
    res["trapped_charge_pC_800nm"] = tr["number"] * per_unit * constants.e * 1e12
    res["history"] = {"t": [e["t"] for e in each], "pulse_center": [e["pulse_center"] for e in each],
                      "a0": [e["a0"] for e in each],
                      "cavity_radius": [e["cavity"]["radius"] if e["cavity"] else float("nan") for e in each],
                      "Ex_min": [e["Ex_min"] for e in each], "Ex_slope": [e["Ex_slope"] for e in each],
                      "trapped_number": [t["number"] for t in beam], "u_max": [t["u_max"] for t in beam],
                      "u_peak": [t["u_peak"] for t in beam]}
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(3, 1, figsize=(11, 10))
        xn, zn = s.plane_axes("n:electron", last)
        xe, ze = s.plane_axes("Ex", last)
        P.field_map(ax[0], s.read("n:electron", last), [xn[0], xn[-1], zn[0], zn[-1]], signed=False, vlim=(0.0, 4.0),
                    label="electron density [n0]", ylabel="z [c/ωpe]")
        tw = ax[0].twinx()
        tw.plot(xe, line(s.read("Ex", last), ze, 0.5 * p.width), color=P.color(1), lw=0.8)
        if np.isfinite(b["Ex_zero"]):  # the slope of Lu et al. through the measured zero
            half = 0.5 * th["radius"]
            tw.plot([b["Ex_zero"] - half, b["Ex_zero"] + half], [-th["slope"] * half, th["slope"] * half], **P.THEORY)
        tw.set_ylabel("E_x on the axis")
        tw.yaxis.set_label_position("left")  # the right side holds the colour bar of the map
        tw.yaxis.tick_left()
        tw.spines["left"].set_position(("outward", 42))
        tw.grid(False)
        ax[0].set_title(f"t = {b['t']:.1f}: cavity radius {res['cavity_radius']:.2f} (2√a0 = {th['radius']:.2f}), "
                        f"slope of E_x {b['Ex_slope']:.3f} (1/2)", fontsize=9)
        ph = run.series("phase")
        pe = ph.edges
        P.phase_space(ax[1], ph.read("electron", last), [pe[0] + float(ph.window_shift[last]) * run.dx[0], pe[1]],
                      xlabel="x [c/ωpe]", ylabel="u_x")
        ax[1].set_xlim(ax[0].get_xlim())
        cen, f = spectrum(edges, sp.read("electron", last))
        ok = (f > 0) & (cen >= 0.0)
        ax[2].semilogy(cen[ok], f[ok], color=P.color(0))
        ax[2].axvline(u_min, **P.THEORY)
        ax[2].set_xlabel("u_x")
        ax[2].set_ylabel("dN/du_x [n0 (c/ωpe)³]")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(fig_dir or run.path, "wakefield.png"))
    return res


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    p = PL.params_from_dict(info["params"])
    if p.regime == "blowout":
        return _blowout(run, p, info, show_plot, fig_dir)
    th = info["theory"]
    s = run.series("axis_plane")
    last = len(s.t) - 1
    xe, ze = s.plane_axes("Ex", last)
    zc = 0.5 * p.width
    # the box at the last output: its lower face and the front of the plasma of a moving window
    x_low = float(xe[0] - 0.5 * (xe[1] - xe[0]))
    w = run.window
    front = (w["front_cell"] + int(s.window_shift[last])) * run.dx[0] if w else np.inf
    predicted = th["Ex_axis_focus"] if p.driver == "laser" else th["Ex_axis"]
    res = {"driver": p.driver, "drive": p.drive, "t": float(s.t[last]), "predicted_amplitude": predicted,
           "speed": info["speed"], "x_low": x_low}
    # mean electron density on the plateau, in the plane, without the two nodes next to each face
    xn, zn = s.plane_axes("n:electron", last)
    seln = (xn >= max(p.ramp[1], x_low) + 1.0) & (xn <= min(p.ramp[2], front) - 1.0)
    res["density_plateau"] = float(s.read("n:electron", last)[seln][:, 2:-2].mean())
    ex_axis = np.array([line(s.read("Ex", i), ze, zc) for i in range(last + 1)])
    res["Ex_axis_max"] = float(np.abs(ex_axis).max())
    g = run.scalars["gauss"]
    res["gauss_max"] = float(np.nanmax(g)) if np.isfinite(g).any() else float("nan")
    m = None
    if not p.drive:
        res["Ex_axis_max_over_predicted"] = res["Ex_axis_max"] / predicted
        res["passed"] = bool(res["Ex_axis_max_over_predicted"] <= 0.01 and abs(res["density_plateau"] - 1.0) <= 0.01)
    else:
        m = _measure(p, s, last)
        res.update({k: v for k, v in m.items() if k not in ("fit", "x", "origin")})
        before = _measure(p, s, last - 1)
        res["amplitude_ratio_previous_output"] = before["amplitude_ratio"]
        res["t_previous_output"] = before["t"]
        ok = abs(m["amplitude_ratio"] - 1.0) <= 0.05 and abs(m["k_times_speed"] - 1.0) <= 0.02
        if p.driver == "laser":
            res["Er_over_Ex_theory"] = th["Er_over_Ex_half_waist"]
            xy, zy = s.plane_axes("Ey", last)
            ey = line(s.read("Ey", last), zy, zc) ** 2
            res["pulse_center_measured"] = float(np.sum(xy * ey) / np.sum(ey))
            ok = ok and abs(m["Er_over_Ex_ratio"] - 1.0) <= 0.10
        else:
            res["R_ratio_theory"] = th["R_ratio"]
            res["R_ratio_rel_error"] = m["R_ratio"] / th["R_ratio"] - 1.0
            xb, zb = s.plane_axes("n:bunch", last)
            nb = line(s.read("n:bunch", last), zb, zc)
            res["bunch_center_measured"] = float(np.sum(xb * nb) / np.sum(nb))
            ok = ok and abs(res["R_ratio_rel_error"]) <= 0.10
        res["passed"] = bool(ok)
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
        P.field_map(ax[0], s.read("Ex", last), [xe[0], xe[-1], ze[0], ze[-1]], label="E_x", ylabel="z [c/ωpe]")
        ax[0].set_aspect("auto")
        ax[1].plot(xe, ex_axis[last], color=P.color(0), lw=0.8, label=f"E_x on the axis, t={s.t[last]:.1f}")
        if m is not None:
            f, x = m["fit"], m["x"]
            amp = PL.wake_amplitude(p, x)[0]
            ax[1].plot(x, f["c"] * amp * np.cos(f["k"] * (x - m["origin"]) + f["phi"]) + f["offset"], color=P.color(1),
                       lw=0.8, label=f"fit: c={f['c']:.3f}, k·v={m['k_times_speed']:.4f}")
            ax[1].plot(x, amp, **P.THEORY, label="linear theory (amplitude)")
            ax[1].plot(x, -amp, **P.THEORY)
        ax[1].set_xlabel("x [c/ωpe]")
        ax[1].legend(loc="upper left", fontsize=7)
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(fig_dir or run_dir, "wakefield.png"))
    return res


def movie(run_dir: str, out: str | None = None, fps: int = 6, n_frames: int = 100, dpi: int = 110) -> str:
    """The plane through the axis at its laboratory place: E_x and the driver (E_y of the pulse, or the density
    of the bunch); the radial field of a linear wake; the electron density and the phase space (x, u_x) of a
    blowout run."""
    from pypic.analysis import Movie
    from pypic.analysis.animate import PhasePanel, SlicePanel

    run = open_run(run_dir)
    p = PL.params_from_dict(run.meta["planner"]["params"])
    blowout = p.regime == "blowout"
    m = Movie(run, title=f"wake of a {'laser pulse' if p.driver == 'laser' else 'bunch'}", ncols=1,
              panel_size=(9.0, 2.6))
    if blowout:
        m.add(SlicePanel(run, "axis_plane", "n:electron", signed=False, vlim=(0.0, 4.0), label="electron density [n0]"))
    m.add(SlicePanel(run, "axis_plane", "Ex", label="E_x"))
    if not blowout:
        m.add(SlicePanel(run, "axis_plane", "Ez", label="E_z"))
    if p.driver == "laser":
        m.add(SlicePanel(run, "axis_plane", "Ey", label="E_y (pulse)"))
    else:
        m.add(SlicePanel(run, "axis_plane", "n:bunch", signed=False, vlim=(0.0, p.bunch_density),
                         label="bunch density [n0]"))
    if blowout:
        m.add(PhasePanel(run, "phase", "electron"))
    return m.save(out or os.path.join(run_dir, "wakefield.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
