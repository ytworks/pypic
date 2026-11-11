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
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import plan as PL  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402


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


def _window(p, x_driver: float):
    """(x_lo, x_hi) of the measurement behind the driver."""
    if p.driver == "bunch":
        hi = x_driver - 5.0 * p.sigma_x
        return max(p.start + 30.0, hi - 5.0 * 2.0 * np.pi * PL.speed(p)), hi
    return p.ramp[1] + 1.0, x_driver - 3.0 * p.duration


def _measure(p, s, i: int) -> dict:
    """Fits of the output i of the plane through the axis."""
    xe, ze = s.plane_axes("Ex")
    zc, speed = 0.5 * p.width, PL.speed(p)
    r = 0.5 * p.waist if p.driver == "laser" else p.sigma_r
    x_drv = PL.driver_position(p, float(s.t[i]))
    lo, hi = _window(p, x_drv)
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
        xz, zz = s.plane_axes("Ez")
        selz = (xz >= lo) & (xz <= hi)
        ez = s.read("Ez", i)
        # E_r = -c_r A_r sin(k x + phi) + b' = c_r A_r cos(k x + phi + pi / 2) + b'
        cr = 0.5 * sum(sg * fit_amplitude(xz[selz], line(ez, zz, zc + sg * r)[selz],
                                          PL.wake_amplitude(p, xz[selz], r)[1], on["k"], on["phi"] + 0.5 * np.pi)
                       for sg in (1.0, -1.0))
        out["radial_amplitude_ratio"] = cr
        out["Er_over_Ex_ratio"] = cr / cx
    return out


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    p = PL.params_from_dict(info["params"])
    th = info["theory"]
    s = run.series("axis_plane")
    xe, ze = s.plane_axes("Ex")
    zc = 0.5 * p.width
    last = len(s.t) - 1
    predicted = th["Ex_axis_focus"] if p.driver == "laser" else th["Ex_axis"]
    res = {"driver": p.driver, "drive": p.drive, "t": float(s.t[last]), "predicted_amplitude": predicted,
           "speed": info["speed"]}
    # mean electron density on the plateau, in the plane, without the two nodes next to each face
    xn, zn = s.plane_axes("n:electron")
    seln = (xn >= p.ramp[1] + 1.0) & (xn <= p.ramp[2] - 1.0)
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
            ey = line(s.read("Ey", last), s.plane_axes("Ey")[1], zc) ** 2
            res["pulse_center_measured"] = float(np.sum(s.plane_axes("Ey")[0] * ey) / np.sum(ey))
            ok = ok and abs(m["Er_over_Ex_ratio"] - 1.0) <= 0.10
        else:
            res["R_ratio_theory"] = th["R_ratio"]
            res["R_ratio_rel_error"] = m["R_ratio"] / th["R_ratio"] - 1.0
            nb = line(s.read("n:bunch", last), s.plane_axes("n:bunch")[1], zc)
            res["bunch_center_measured"] = float(np.sum(s.plane_axes("n:bunch")[0] * nb) / np.sum(nb))
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


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    a = ap.parse_args()
    print(json.dumps(analyze(a.run_dir), indent=2, default=float))
