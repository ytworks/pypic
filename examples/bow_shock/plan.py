"""Planner for the bow shock of a magnetised planet in the solar wind (3-D, explicit scheme).

A uniform magnetised wind (electrons and ions, speed v_sw along +x, interplanetary field B_imf)
fills an open box and keeps entering through the upstream face and the four side faces; the
downstream face only removes particles. The planet is an absorbing sphere with a magnetic dipole
at its centre. The wind field (B_imf and the motional field E = -v x B) and the dipole are external
fields, so the grid holds the difference to them and the absorbing layers of the open faces relax
the fields to the wind values (PIC-D12 to PIC-D16).

The planner takes dimensionless numbers (mass ratio, Alfven Mach number, omega_pe / Omega_ce,
magnetopause distance in ion skin depths) and sets the dipole moment from the pressure balance at
the wanted magnetopause distance R_mp (theory.py). Lengths of the box, of the body and of the
diagnostics are given in units of R_mp.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

from pypic.config import (
    ID_LIMIT,
    SLOT_BYTES,
    SORT_BYTES,
    Beam,
    Body,
    CurrentFilter,
    Dipole,
    ExternalFields,
    Injection,
    OutputConfig,
    SimConfig,
    SpeciesConfig,
)
from pypic.injection import ids_per_step
from pypic.parallel.decomp import auto_dims
from pypic.setup import mhd
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402

WCE_DT_MAX = 2.0  # largest Omega_ce dt at the surface of the body the planner accepts
GHOST = 2  # ghost layers of the grid (pypic.grid.Grid)


@dataclass
class BowShockParams:
    # --- wind ---
    mass_ratio: float = 16.0  # m_i / m_e
    mach_alfven: float = 4.0  # v_sw / v_A, v_A the Alfven speed of the IMF with the ion mass
    wpe_wce: float = 4.0  # omega_pe / Omega_ce of the IMF (B_imf = 1 / wpe_wce)
    vte: float = 0.1  # sqrt(T_e / m_e c^2)
    ti_te: float = 1.0
    imf: tuple = (0.0, 0.0, 1.0)  # direction of the IMF
    # --- planet ---
    dp: float = 5.0  # magnetopause distance R_mp in ion skin depths
    body_radius: float = 1.0 / 1.5  # in R_mp
    core: float = 0.2  # regularisation length of the dipole in R_mp (inside the body)
    dipole: tuple = (0.0, 0.0, -1.0)  # direction of the moment
    # --- box and numerics ---
    box: tuple = (6.0, 8.0, 8.0)  # in R_mp
    upstream: float = 2.5  # distance of the planet's centre from the upstream face in R_mp
    cells_per_de: float = 2.0
    ppc: int = 8
    crossings: float = 3.0  # run length in crossing times L_x / v_sw
    filter_passes: int = 4  # binomial passes of the current filter (compensated)
    marder_zone: float = 0.125
    npml: int = 12
    capacity_factor: float = 2.0  # the sheath compresses the plasma: room for the fullest rank
    # --- output ---
    coarse: float = 0.2  # block size of the 3-D output in R_mp
    slice_coarsen: int = 1  # block size of the two slices in cells
    seed: int = 12345
    overrides: dict = field(default_factory=dict)  # any SimConfig field


def params_from_dict(d: dict) -> BowShockParams:
    d = dict(d)
    for k in ("imf", "dipole", "box"):
        if k in d:
            d[k] = tuple(float(v) for v in d[k])
    return BowShockParams(**d)


def _unit(v, name):
    n = math.sqrt(sum(c * c for c in v))
    if len(v) != 3 or n == 0.0:
        raise ValueError(f"{name} needs three components, not all zero")
    return tuple(c / n for c in v)


def _cells(length: float, dx: float) -> int:
    return max(8, int(round(length / dx / 8)) * 8)


def plan_bow_shock(p: BowShockParams, n_gpus: int = 1):
    b_imf = 1.0 / p.wpe_wce
    v = p.mach_alfven * mhd.alfven_speed(b_imf, p.mass_ratio)
    if not 0.0 < v < 1.0:
        raise ValueError(f"wind speed {v:g} c: lower mach_alfven or raise wpe_wce or mass_ratio")
    gamma = 1.0 / math.sqrt(1.0 - v * v)
    u = (gamma * v, 0.0, 0.0)
    b_dir = _unit(p.imf, "imf")
    B = tuple(b_imf * c for c in b_dir)
    E = (0.0, v * B[2], -v * B[1])  # -v x B for v along x
    r_mp = p.dp * math.sqrt(p.mass_ratio)  # d_i = sqrt(m_i / m_e) d_e
    dx = 1.0 / p.cells_per_de
    dt = 0.5 * dx
    n_cells = tuple(_cells(b * r_mp, dx) for b in p.box)
    L = tuple(n * dx for n in n_cells)
    center = (p.upstream * r_mp, 0.5 * L[1], 0.5 * L[2])
    radius, core = p.body_radius * r_mp, p.core * r_mp
    moment = T.dipole_moment(r_mp, v, p.mass_ratio)
    m_dir = _unit(p.dipole, "dipole")
    # the largest external field a particle can meet: the pole of the body's surface
    wce_dt = (2.0 * moment / radius ** 3 + b_imf) * dt
    if wce_dt > WCE_DT_MAX:
        raise ValueError(f"Omega_ce dt = {wce_dt:.2f} at the surface of the body (limit {WCE_DT_MAX:g}): "
                         f"raise body_radius or cells_per_de")
    t_cross = L[0] / v
    cross_steps = int(round(t_cross / dt))
    n_steps = int(round(p.crossings * cross_steps))
    every_3d, every_slice = max(1, cross_steps // 4), max(1, cross_steps // 16)
    theta_e = p.vte ** 2
    species = [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, u, theta_e)]),
               SpeciesConfig("i", 1.0, p.mass_ratio, [Beam(1.0, u, theta_e * p.ti_te / p.mass_ratio)])]
    c3 = max(1, int(round(p.coarse * r_mp / dx)))
    if not all(T.upstream_blocks(n_cells, (c3,) * 3, (dx,) * 3, r_mp)):
        raise ValueError(f"coarse = {p.coarse:g} R_mp leaves no block of the 3-D output inside the upstream region "
                         f"(x from {T.UPSTREAM_X[0]:g} to {T.UPSTREAM_X[1]:g} R_mp): lower coarse")
    cs = int(p.slice_coarsen)
    moments = ["n:e", "n:i", "flux:e", "flux:i", "ekin:e", "ekin:i", "E", "B", "rho"]
    in_plane = ["n:e", "n:i", "flux:i", "B"]
    diags = [
        {"name": "coarse", "kind": "fields", "fields": moments, "coarsen": [c3, c3, c3], "every": every_3d},
        {"name": "slice_z", "kind": "slice", "fields": in_plane, "axis": "z", "position": center[2],
         "coarsen": [cs, cs, 1], "every": every_slice},
        {"name": "slice_y", "kind": "slice", "fields": in_plane, "axis": "y", "position": center[1],
         "coarsen": [cs, 1, cs], "every": every_slice},
    ]
    try:
        dims = auto_dims(n_gpus, n_cells)
    except ValueError:
        dims = None
    omega_ci = b_imf / p.mass_ratio
    cfg = SimConfig(
        n_cells=n_cells, dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc, seed=p.seed, dims=dims,
        boundaries=(("open", "open"),) * 3, npml=p.npml, capacity_factor=p.capacity_factor,
        external_fields=ExternalFields(E=E, B=B, dipoles=[Dipole(tuple(moment * c for c in m_dir), center, core)]),
        body=Body(center=center, radius=radius),
        injection=Injection(faces={"x-": "paired", "y-": "flux", "y+": "flux", "z-": "flux", "z+": "flux"},
                            neutralizer="e"),
        current_filter=CurrentFilter(p.filter_passes, True), marder_zone=p.marder_zone,
        diag_interval=max(1, n_steps // 288), gauss_interval=every_3d,
        output=OutputConfig(diagnostics=diags, timescales={"crossing": t_cross, "ion_gyration": 2 * math.pi / omega_ci,
                                                           "plasma_period": 2 * math.pi}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    n_cell = math.prod(cfg.n_cells)
    n_part = 2 * n_cell * cfg.ppc
    ids = n_part + cfg.n_steps * ids_per_step(cfg)
    if ids >= ID_LIMIT:
        raise ValueError(f"the run needs {ids:.3g} particle ids, the limit is 2^53: fewer steps, cells or ppc")
    gpus = max(n_gpus, 1)
    # Device memory of a rank whose particle arrays are full (terms measured at the gate): the bytes
    # of a particle slot and of the sort of one species per particle, by precision (pypic.config);
    # 31 arrays of the grid and the arrays an output step holds (one per component of the 3-D
    # output), each with the ghost and absorbing layers; 1 GB of CUDA context, kernels and NCCL
    # buffers.
    slots = n_part * cfg.capacity_factor / gpus
    padded = math.prod(n / d + (2 if d == 1 else 1) * cfg.npml + 2 * GHOST
                       for n, d in zip(cfg.n_cells, cfg.dims or (1, 1, 1)))
    n_out = sum(3 if m.split(":")[0] in ("flux", "E", "B") else 1 for m in moments)
    mem = (slots * SLOT_BYTES[cfg.precision] + slots / 2 * SORT_BYTES[cfg.precision]
           + padded * 8 * (31 + n_out)) / 1e9 + 1.0
    beta = mhd.plasma_beta(p.vte, p.ti_te, b_imf)
    x_mhd = mhd.compression(p.mach_alfven, beta)
    lam_d = p.vte  # Debye length in d_e
    info = {
        "r_mp": r_mp, "d_i": math.sqrt(p.mass_ratio), "center": list(center), "body_radius": radius, "core": core,
        "v_sw": v, "b_imf": b_imf, "B": list(B), "E": list(E), "moment": moment,
        "b_standoff": T.standoff_field(v, p.mass_ratio), "ram_pressure": p.mass_ratio * v * v,
        "beta": beta, "mach_alfven": p.mach_alfven, "mach_fast": mhd.fast_mach(p.mach_alfven, beta),
        "compression_mhd": x_mhd, "sheath_spreiter": 1.1 / x_mhd,
        "shock_standoff_spreiter": T.shock_standoff(r_mp, x_mhd),
        "crossing_time": t_cross, "crossing_steps": cross_steps, "omega_ci_t_end": omega_ci * cfg.n_steps * cfg.dt,
        "wce_dt_surface": wce_dt, "debye_cells": lam_d / dx, "particles": n_part, "memory_per_gpu_GB": mem,
        "id_fraction": ids / ID_LIMIT, "coarse_cells": c3, "dims": cfg.dims, "params": asdict(p)}
    lines = [
        "=== pypic bow-shock planner (magnetised planet in the solar wind) ===",
        f"wind: m_i/m_e={p.mass_ratio:g}, M_A={p.mach_alfven:g} (fast Mach {info['mach_fast']:.2f}), "
        f"v_sw={v:.4f} c, B_imf={b_imf:g} along {tuple(round(c, 3) for c in b_dir)}, beta={beta:.2f}, "
        f"v_te={p.vte:g} c (Debye length {lam_d / dx:.2f} dx)",
        f"planet: R_mp={r_mp:g} d_e = {p.dp:g} d_i ({r_mp / dx:.0f} cells), body radius {radius:.2f}, "
        f"centre ({center[0]:g}, {center[1]:g}, {center[2]:g}); dipole moment {moment:.1f} along "
        f"{tuple(round(c, 3) for c in m_dir)}, core {core:g}; Omega_ce dt <= {wce_dt:.2f} at the surface",
        f"grid {cfg.n_cells}, dx={dx:g}, L=({L[0]:g}, {L[1]:g}, {L[2]:g}) d_e, dt={cfg.dt:g}, {cfg.n_steps} steps = "
        f"{p.crossings:g} crossings of {cross_steps} steps = {info['omega_ci_t_end']:.1f} / Omega_ci",
        f"{n_part:,} particles ({cfg.ppc} per cell and species), capacity factor {cfg.capacity_factor:g}; "
        f"memory up to {mem:.1f} GB/GPU on {gpus} GPU(s), rank grid {cfg.dims}; ids {info['id_fraction']:.1e} of 2^53",
        f"faces: x- paired inflow, y and z flux, x+ outflow; filter {p.filter_passes}+comp, Marder zone "
        f"{cfg.marder_zone:g}, npml {cfg.npml}",
        f"theory: perpendicular MHD shock n2/n1={x_mhd:.2f}; sheath (Spreiter) {info['sheath_spreiter']:.2f} R_mp, "
        f"shock at {info['shock_standoff_spreiter']:.1f} d_e from the centre",
        f"output: 3-D moments in blocks of {c3} cells every {every_3d} steps, slices z and y through the centre "
        f"every {every_slice} steps",
    ]
    return cfg, "\n".join(lines), info
