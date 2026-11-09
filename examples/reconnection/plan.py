"""Planner for relativistic pair reconnection in a periodic double Harris sheet (quasi-2-D).

Two current sheets at y = L_y/4 and 3 L_y/4 (B_x = B0 [tanh((y-y1)/d) - tanh((y-y2)/d) - 1], periodic
in x and y) over a cold uniform pair background of magnetisation sigma. Reconnection starts from
noise (tearing); the plasmoid chain accelerates particles into a power law. A compensated current
filter suppresses numerical Cherenkov noise of the relativistic outflows.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

from pypic.config import Beam, CurrentFilter, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402


@dataclass
class ReconnectionParams:
    sigma: float = 10.0
    n_cs_over_nb: float = 5.0
    delta: float = 1.0  # sheet half-thickness in c/omega_pe (background)
    theta_b: float = 0.01  # background temperature
    Lx: float = 100.0
    Ly: float = 100.0
    cells_per_skin: float = 4.0
    nz: int = 4
    ppc: int = 4
    t_end: float = 400.0
    filter_passes: int = 2
    n_tracks: int = 20000
    seed: int = 17
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> ReconnectionParams:
    return ReconnectionParams(**d)


def plan_reconnection(p: ReconnectionParams, n_gpus: int = 1):
    h = T.harris(p.sigma, p.n_cs_over_nb, p.delta)
    dx = 1.0 / p.cells_per_skin
    nx = int(math.ceil(p.Lx / dx / 8)) * 8
    ny = int(math.ceil(p.Ly / dx / 8)) * 8
    Lx, Ly = nx * dx, ny * dx
    y1, y2 = 0.25 * Ly, 0.75 * Ly
    dt = 0.5 * dx
    n_steps = int(math.ceil(p.t_end / dt))
    ud = h["u_d"]
    s1 = {"axis": "y", "shape": "sech2", "center": y1, "width": p.delta}
    s2 = {"axis": "y", "shape": "sech2", "center": y2, "width": p.delta}

    def beams(sign):  # sign = -1 for electrons, +1 for positrons (J_z < 0 at y1, > 0 at y2)
        return [Beam(1.0, (0.0, 0.0, 0.0), p.theta_b),
                Beam(h["n_cs"], (0.0, 0.0, sign * -ud), h["theta_cs"], profile=s1),
                Beam(h["n_cs"], (0.0, 0.0, sign * ud), h["theta_cs"], profile=s2)]

    species = [SpeciesConfig("electron", -1.0, 1.0, beams(-1)), SpeciesConfig("positron", 1.0, 1.0, beams(+1))]
    field_profile = {"kind": "profile", "component": "Bx", "axis": "y", "shape": "double_tanh",
                     "centers": [y1, y2], "width": p.delta, "scale": h["B0"]}
    n_part = nx * ny * p.nz * p.ppc * 3
    g_hi = 50.0 * (1.0 + p.sigma)
    diags = [
        {"name": "slice_xy", "kind": "slice", "fields": ["Bx", "By", "Jz", "n:electron"], "axis": "z",
         "coarsen": [2, 2, 1], "every": {"count": 80}},
        {"name": "profile_y", "kind": "lineout", "fields": ["Bx", "Jz", "n:electron"], "axis": "y",
         "every": {"count": 200}},
        {"name": "spectrum", "kind": "histogram", "quantities": ["gamma"], "bins": [160], "log": [True],
         "ranges": [[1.0, g_hi]], "every": {"count": 80}},
        {"name": "tracks", "kind": "tracks", "species": ["electron"], "stride": max(1, n_part // p.n_tracks),
         "with_fields": True, "every": {"count": 800}},
    ]
    try:
        dims = auto_dims(n_gpus, (nx, ny, p.nz))
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=(nx, ny, p.nz), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims, initial_fields=[field_profile], diag_interval=max(1, n_steps // 500),
                    gauss_interval=max(1, n_steps // 25), sort_interval=50,
                    current_filter=CurrentFilter(p.filter_passes, True),
                    output=OutputConfig(diagnostics=diags, timescales={"plasma_period": 2 * math.pi,
                                                                       "t_end": p.t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info = {**h, "sigma": p.sigma, "delta": p.delta, "y_sheets": [y1, y2], "Lx": Lx, "Ly": Ly,
            "werner_index": T.werner_index(p.sigma), "params": asdict(p)}
    lines = ["=== pypic relativistic reconnection planner (double Harris sheet, pair plasma) ===",
             f"sigma={p.sigma:g}: B0={h['B0']:.3f}, v_A={h['v_A']:.3f} c; sheets n_cs={h['n_cs']:g} n_b, "
             f"delta={p.delta:g}, beta_d={h['beta_d']:.3f} (u_d={ud:.3f}), theta_cs={h['theta_cs']:.3f}",
             f"grid ({nx}, {ny}, {p.nz}), dx={dx:.3f}, L=({Lx:.0f}, {Ly:.0f}), dt={dt:.4f}, {n_steps} steps to "
             f"t={p.t_end:g}; {n_part:,} particle slots (sheet beams keep only weights > 1e-10)",
             f"published scalings: rate ~0.1 v_A, Werner+2016 p = 1.9 + 0.7/sqrt(sigma) = {info['werner_index']:.2f}"]
    return cfg, "\n".join(lines), info
