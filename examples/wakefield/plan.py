"""Planner for a plasma wake behind a laser pulse or a particle bunch (static box).

A cold plasma of electrons and ions fills a plateau along x between two linear ramps; the x faces
are open and no plasma touches them. The driver enters from the low-x side:

- ``driver: laser``: an antenna (PIC-D20) in the vacuum before the ramp emits a Gaussian beam with
  a Gaussian envelope, polarised along y, focused at ``focus``. The run ends when the centre of
  the pulse is at ``x_end``; its travel time is the integral of 1 / v_g over the local density,
  with the group velocity of the Yee grid.
- ``driver: bunch``: a rigid Gaussian bunch (very heavy particles) starts at ``start`` inside the
  plateau on top of a partner of the opposite charge at rest (PIC-D21). The run ends when the
  bunch is at ``x_end``.

``drive: false`` runs the same box for the same time without the driver (null test).

Theory (theory.py): the linear wake of both drivers; the analysis compares the on-axis E_x, its
wavelength and the radial dependence with it.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

import numpy as np

from pypic.antenna import diffraction_wavenumber, wavenumber
from pypic.config import Beam, Laser, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402


@dataclass
class WakeParams:
    driver: str = "laser"  # "laser" or "bunch"
    drive: bool = True  # False: the same box and run time without the driver
    nx: int = 3200
    n_perp: int = 128  # cells along y and along z
    width: float = 24.0  # size along y and z, c/omega_pe
    dx: float | None = None  # cell along x; None: lambda0 / cells_per_wavelength (laser)
    cells_per_wavelength: float = 32.0
    transverse: str = "periodic"  # y and z faces: "periodic" or "reflecting"
    ramp: tuple = (4.0, 8.0, 58.0, 61.0)  # density 0, n0, n0, 0 at these x (linear in between)
    x_end: float = 54.0  # place of the driver at the end of the run
    ppc: int = 8
    ion_mass: float = 1836.0
    n_out: int = 8  # outputs of the plane through the axis
    # laser
    omega0: float = 10.0  # omega_pe
    a0: float = 0.1  # in vacuum at the focus
    duration: float = 2.0  # field envelope exp(-(t / duration)^2), 1/omega_pe
    waist: float = 4.0  # field envelope exp(-r^2 / waist^2) at the focus
    focus: float | None = None  # None: the centre of the box
    plane_x: float = 2.0  # antenna
    t_peak: float | None = None  # None: 4 * duration
    # bunch
    bunch_density: float = 0.05  # peak, n0
    sigma_x: float = 1.0
    sigma_r: float = 1.0
    gamma: float = 10.0
    bunch_mass: float = 1.0e6
    start: float = 12.0
    seed: int = 12345
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> WakeParams:
    p = WakeParams(**d)
    p.ramp = tuple(float(v) for v in p.ramp)
    if p.driver not in ("laser", "bunch"):
        raise ValueError("driver must be 'laser' or 'bunch'")
    if len(p.ramp) != 4 or list(p.ramp) != sorted(p.ramp):
        raise ValueError("ramp: four increasing x positions")
    if p.driver == "bunch" and p.dx is None:
        raise ValueError("driver 'bunch' needs dx")
    return p


def cells(p: WakeParams):
    """(dx, dy = dz, dt)."""
    dx = p.dx if p.dx is not None else 2.0 * math.pi / (p.omega0 * p.cells_per_wavelength)
    return dx, p.width / p.n_perp, 0.5 * dx


def density(p: WakeParams, x):
    """Density of the plasma at x in n0."""
    return np.interp(x, p.ramp, [0.0, 1.0, 1.0, 0.0])


def laser_geometry(p: WakeParams) -> dict:
    """Plane of the antenna, focus, wavenumbers of the grid in vacuum and Rayleigh length."""
    dx, _, dt = cells(p)
    k = wavenumber(p.omega0, dx, dt)
    kd = diffraction_wavenumber(k, dx)
    return {"x_plane": (math.floor(p.plane_x / dx) + 0.5) * dx,
            "focus": p.focus if p.focus is not None else 0.5 * p.nx * dx,
            "t_peak": p.t_peak if p.t_peak is not None else 4.0 * p.duration,
            "k": k, "k_d": kd, "rayleigh": 0.5 * kd * p.waist ** 2,
            "v_g": T.yee_group_velocity(p.omega0, dx, dt, 1.0),
            "v_g_vacuum": T.yee_group_velocity(p.omega0, dx, dt, 0.0)}


def _pulse_path(p: WakeParams, n: int = 8001):
    """(x, t): the time at which the centre of the pulse is at x, from the antenna to the far face."""
    dx, _, dt = cells(p)
    g = laser_geometry(p)
    x = np.linspace(g["x_plane"], p.nx * dx, n)
    slow = np.array([1.0 / T.yee_group_velocity(p.omega0, dx, dt, math.sqrt(d)) for d in density(p, x)])
    t = g["t_peak"] + np.concatenate([[0.0], np.cumsum(0.5 * (slow[1:] + slow[:-1]) * np.diff(x))])
    return x, t


def speed(p: WakeParams) -> float:
    """Speed of the driver on the plateau: the wake has the wavenumber 1 / speed."""
    if p.driver == "bunch":
        return math.sqrt(1.0 - 1.0 / p.gamma ** 2)
    return laser_geometry(p)["v_g"]


def arrival_time(p: WakeParams, x: float) -> float:
    """Time at which the driver (centre of the pulse or of the bunch) is at x."""
    if p.driver == "bunch":
        return (x - p.start) / speed(p)
    xs, ts = _pulse_path(p)
    return float(np.interp(x, xs, ts))


def driver_position(p: WakeParams, t: float) -> float:
    """Place of the driver at the time t."""
    if p.driver == "bunch":
        return p.start + speed(p) * t
    xs, ts = _pulse_path(p)
    return float(np.interp(t, ts, xs))


def wake_amplitude(p: WakeParams, x, r: float = 0.0):
    """Amplitudes (E_x, E_r) of the linear wake at the places x on the plateau, at the radius r.
    The pulse leaves at x the wake of its local a0 and waist; the bunch is rigid (E_r of the bunch
    is not evaluated: NaN)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    if p.driver == "bunch":
        ex = float(T.beam_wake(p.bunch_density, p.sigma_x, p.sigma_r, p.gamma, r)[0])
        return np.full(x.shape, ex), np.full(x.shape, np.nan)
    g = laser_geometry(p)
    out = []
    for xx in x:
        a, w = T.beam_at(xx - g["focus"], p.a0, p.waist, g["rayleigh"])
        out.append(T.laser_wake(T.amplitude_in_plasma(a, p.omega0), p.duration, w, r, g["v_g"]))
    out = np.asarray(out)
    return out[:, 0], out[:, 1]


def plan_wakefield(p: WakeParams, n_gpus: int = 1):
    dx, dy, dt = cells(p)
    Lx, c_perp = p.nx * dx, 0.5 * p.width
    n_steps = int(math.ceil(arrival_time(p, p.x_end) / dt / p.n_out)) * p.n_out
    plasma = {"axis": "x", "shape": "table", "x": list(p.ramp), "values": [0.0, 1.0, 1.0, 0.0]}
    species = [SpeciesConfig("electron", -1.0, 1.0, [Beam(1.0, profile=plasma)]),
               SpeciesConfig("ion", 1.0, p.ion_mass, [Beam(1.0, profile=plasma)])]
    lasers, fields = [], ["Ex", "Ez", "n:electron"]
    if p.driver == "laser":
        g = laser_geometry(p)
        if p.drive:
            lasers = [Laser(plane_x=p.plane_x, a0=p.a0, omega=p.omega0, polarization="y", waist=p.waist,
                            focus_x=g["focus"],
                            envelope={"shape": "gaussian", "duration": p.duration, "t_peak": g["t_peak"]})]
            fields = fields + ["Ey"]
        ex0, er0 = wake_amplitude(p, g["focus"], 0.5 * p.waist)
        theory = {**g, "Ex_axis_focus": float(wake_amplitude(p, g["focus"])[0][0]),
                  "Er_over_Ex_half_waist": float(er0[0] / ex0[0])}
        head = (f"laser a0={p.a0:g}, omega0={p.omega0:g}, duration={p.duration:g}, waist={p.waist:g}; antenna at "
                f"x={g['x_plane']:.3f}, peak at t={g['t_peak']:g}, focus x={g['focus']:.2f}, "
                f"Rayleigh length {g['rayleigh']:.1f}")
        pred = (f"linear wake: v_g={g['v_g']:.5f} (vacuum {g['v_g_vacuum']:.5f}), E_x on the axis at the focus "
                f"{theory['Ex_axis_focus']:.4e}, E_r/E_x at r=waist/2 {theory['Er_over_Ex_half_waist']:.4f}")
    else:
        beta = speed(p)
        if p.drive:
            shape = [{"shape": "gauss", "width": math.sqrt(2.0) * s, "axis": a, "center": c}
                     for a, s, c in (("x", p.sigma_x, p.start), ("y", p.sigma_r, c_perp), ("z", p.sigma_r, c_perp))]
            species += [SpeciesConfig("bunch", -1.0, p.bunch_mass,
                                      [Beam(p.bunch_density, drift=(p.gamma * beta, 0.0, 0.0), profile=shape)]),
                        SpeciesConfig("partner", 1.0, p.bunch_mass, [Beam(p.bunch_density, profile=shape)])]
            fields = fields + ["n:bunch"]
        R = T.beam_wake_profile([0.0, p.sigma_r], p.sigma_r, 1.0 / beta)
        theory = {"beta": beta, "k_w": 1.0 / beta, "Ex_axis": float(wake_amplitude(p, p.x_end)[0][0]),
                  "R_ratio": float(R[1] / R[0])}
        head = (f"bunch n_b={p.bunch_density:g}, sigma_x={p.sigma_x:g}, sigma_r={p.sigma_r:g}, gamma={p.gamma:g} "
                f"(mass {p.bunch_mass:g}) on a partner at rest; starts at x={p.start:g}")
        pred = (f"linear wake: k_w=1/beta={theory['k_w']:.5f}, E_x on the axis {theory['Ex_axis']:.5f}, "
                f"R(sigma_r)/R(0)={theory['R_ratio']:.4f}")
    diags = [{"name": "axis_plane", "kind": "slice", "axis": "y", "fields": fields, "every": {"count": p.n_out}}]
    tr = (p.transverse, p.transverse)
    bc = (("open", "open"), tr, tr)
    n_cells = (p.nx, p.n_perp, p.n_perp)
    try:
        dims = auto_dims(n_gpus, n_cells)
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=n_cells, dx=(dx, dy, dy), dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims, boundaries=bc, lasers=lasers, diag_interval=max(1, n_steps // 500),
                    gauss_interval=max(1, n_steps // 25), sort_interval=50,
                    output=OutputConfig(diagnostics=diags, timescales={"plasma_period": 2.0 * math.pi,
                                                                       "t_end": n_steps * dt}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    filled = 0.5 * (p.ramp[2] + p.ramp[3] - p.ramp[0] - p.ramp[1]) / Lx
    n_plasma = int(2 * p.nx * p.n_perp ** 2 * p.ppc * filled)
    info = {"driver": p.driver, "drive": p.drive, "dx": dx, "dy": dy, "dt": dt, "Lx": Lx, "speed": speed(p),
            "t_end": n_steps * dt, "x_driver_end": driver_position(p, n_steps * dt), "theory": theory,
            "params": asdict(p)}
    lines = ["=== pypic wakefield planner (static box) ===",
             head if p.drive else "no driver (null run); box and run time of: " + head,
             f"grid {cfg.n_cells}, dx={dx:.5f}, dy=dz={dy:.5f}, L=({Lx:.2f}, {p.width:g}, {p.width:g}), "
             f"dt={dt:.6f}, {n_steps} steps to t={n_steps * dt:.2f}; x faces open, y and z {p.transverse}",
             f"plasma: ramp {p.ramp[0]:g}..{p.ramp[1]:g}, plateau to {p.ramp[2]:g}, down to {p.ramp[3]:g}; "
             f"electrons and ions (mass {p.ion_mass:g}), {p.ppc} per cell each, about {n_plasma:,} particles",
             f"driver at x={info['x_driver_end']:.2f} at the end",
             pred]
    return cfg, "\n".join(lines), info
