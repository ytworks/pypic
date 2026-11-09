"""Planner for test-particle orbits with exact relativistic solutions (theory.py).

cases (test electrons, no plasma, periodic box)
    gyration    uniform B_z; electrons with u_perp = 0.1, 1, 10 (non- to ultra-relativistic)
    exb         E_y and B_z with E/B = 0.5: guiding centre drifts at E x B / B^2
    hyperbolic  uniform E_x: hyperbolic motion from rest up to gamma ~ 10
    e_gt_b      E_y > B_z: unbounded acceleration (no frame without E)
    plane_wave  vacuum plane wave (exact Yee eigenmode) of normalised amplitude a0: figure-8 orbits,
                drift a0^2 / (4 + a0^2) for particles starting at a node of the vector potential
dt_scale < 1 refines the time step (and, for the plane wave, the cells with it) for the
second-order convergence check of the analysis.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

from pypic.config import OutputConfig, SimConfig, SpeciesConfig
from pypic.setup.planner import apply_overrides

CASES = ("gyration", "exb", "hyperbolic", "e_gt_b", "plane_wave")


@dataclass
class OrbitParams:
    case: str = "gyration"
    B: float = 1.0
    E: float = 0.5
    a0: float = 2.0  # plane wave: e E0 / (m omega c), omega = k c = 1
    n_particles: int = 8
    periods: float = 5.0  # run length in units of the characteristic period
    w_dt: float = 0.05  # characteristic frequency x dt at dt_scale = 1
    dt_scale: float = 1.0
    cells_per_wavelength: int = 64
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> OrbitParams:
    return OrbitParams(**d)


def plan_orbits(p: OrbitParams, n_gpus: int = 1):
    if p.case not in CASES:
        raise ValueError(f"case must be one of {CASES}")
    rng = np.random.default_rng(5)
    E0 = (0.0, 0.0, 0.0)
    B0 = (0.0, 0.0, 0.0)
    fields = []
    info = {"case": p.case, "params": asdict(p)}
    if p.case == "plane_wave":
        k = 1.0
        L = 2 * math.pi / k
        nx = int(round(p.cells_per_wavelength / p.dt_scale))
        dx = L / nx
        dt = 0.5 * dx  # c dt / dx = 0.5 (Yee CFL 1/sqrt(3) for 3-D cells)
        E_amp = p.a0  # omega = 1, m = e = 1
        fields = [{"kind": "plane_wave", "axis": "x", "mode": 1, "polarization": "y", "amplitude": E_amp}]
        n_cells = (nx, 4, 4)
        dxs = (dx, L / 4, L / 4)
        drift = p.a0 ** 2 / (4 + p.a0 ** 2)
        period = 2 * math.pi / (1 - drift)  # lab period of the figure-8 (phase advances at omega (1 - v_d))
        t_end = p.periods * period
        # start at rest at several phases; x0 = 0 is a node of A (E maximal): pure drift a0^2/(4+a0^2)
        xs = np.linspace(0.0, L, p.n_particles, endpoint=False)
        parts = [[float(x), 0.5 * L / 4, 0.5 * L / 4, 0.0, 0.0, 0.0] for x in xs]
        info.update(k=k, E_wave=E_amp, omega=1.0, drift_theory=drift)
    else:
        L = 16.0
        n_cells = (8, 8, 8)
        dxs = (L / 8,) * 3
        if p.case == "gyration":
            B0 = (0.0, 0.0, p.B)
            us = np.geomspace(0.1, 10.0, p.n_particles)
            parts = [[L / 2, L / 2, L / 2, float(u), 0.0, 0.05 * float(u)] for u in us]
            w = p.B  # non-relativistic gyro-frequency (relativistic ones are slower: longer period)
            period = 2 * math.pi / w
        elif p.case == "exb":
            B0 = (0.0, 0.0, p.B)
            E0 = (0.0, p.E * p.B, 0.0)
            parts = [[float(x), L / 2, L / 2, 0.0, 0.0, 0.0] for x in rng.uniform(0, L, p.n_particles)]
            period = 2 * math.pi / p.B
        elif p.case == "hyperbolic":
            E0 = (p.E, 0.0, 0.0)
            parts = [[float(x), L / 2, L / 2, 0.0, 0.0, 0.0] for x in rng.uniform(0, L, p.n_particles)]
            period = 2.0 / p.E  # time to reach u ~ 2
        else:  # e_gt_b
            E0 = (0.0, p.E, 0.0)
            B0 = (0.0, 0.0, p.B)
            if not p.E > p.B:
                raise ValueError("e_gt_b needs E > B")
            parts = [[float(x), L / 2, L / 2, 0.0, 0.0, 0.0] for x in rng.uniform(0, L, p.n_particles)]
            period = 2.0 / p.E
        dt = p.w_dt * period / (2 * math.pi) * p.dt_scale
        t_end = p.periods * period
    n_steps = int(math.ceil(t_end / dt))
    species = [SpeciesConfig("electron", -1.0, 1.0, [], test=True, particles=parts)]
    diags = [{"name": "tracks", "kind": "tracks", "species": ["electron"], "every": 1, "with_fields": True,
              "dtype": "float64"}]
    cfg = SimConfig(n_cells=n_cells, dx=dxs, dt=dt, n_steps=n_steps, species=species, ppc=1, B0=B0, E0=E0,
                    initial_fields=fields, diag_interval=0, gauss_interval=0, sort_interval=0,
                    output=OutputConfig(diagnostics=diags, timescales={"t_end": t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info.update(dt=dt, t_end=t_end, particles=parts, B0=B0, E0=E0, period=period)
    lines = ["=== pypic test-particle orbit planner ===",
             f"case={p.case}: {len(parts)} test electrons, dt={dt:.5g} (dt_scale={p.dt_scale:g}), "
             f"t_end={t_end:.4g} ({n_steps} steps), B0={B0}, E0={E0}" +
             (f", plane wave a0={p.a0:g}, drift {info.get('drift_theory', 0):.4f} c" if p.case == "plane_wave" else "")]
    return cfg, "\n".join(lines), info
