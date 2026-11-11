"""Planner for a blast wave around a magnetised body (3-D, explicit scheme).

A magnetised plasma at rest (electrons and ions, uniform field B0 along z) fills an open box and
keeps entering through its six faces. At the centre sits an absorbing sphere with a magnetic dipole.
A dense shell of debris around the sphere starts with a radial speed V_d = M_A v_A; it sweeps up
the ambient plasma, pushes the field out of a cavity and drives a shock outwards (Le et al. 2021;
Keenan et al. 2022 for a cloud in a uniform field). B0 and the dipole are external fields, so the
grid holds the difference to them (PIC-D12 to PIC-D16).

The planner takes dimensionless numbers (mass ratio, Alfven Mach number of the debris,
omega_pe / Omega_ce, equal-mass radius R_m in ion skin depths). Lengths of the shell and of the
dipole are given in units of R_m, the box in units of the stopping radius R_B = R_m M_A^(2/3), the
run in units of R_B / V_d (theory.py). ``debris: false`` and ``body: false`` leave out the shell and
the sphere with its dipole: the null tests of the faces and of the body.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

from pypic.config import (
    ID_LIMIT,
    MAX_CAPACITY,
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
from pypic.profiles import support
from pypic.setup import mhd
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402

WCE_DT_MAX = 2.0  # largest Omega_ce dt at the surface of the body the planner accepts
WPE_DT_MAX = 1.2  # largest omega_pe dt at the peak of the shell (ambient and debris electrons)
OVERLAP_MAX = 0.01  # largest fraction of the shell's particles inside the body
HALF_WIDTH_MIN = 1.3  # smallest half-width of the box in R_B: room for the front at the end of the run
BODY_CELLS_MIN = 6.0  # smallest radius of the body in cells
GHOST = 2  # ghost layers of the grid (pypic.grid.Grid)
FACES = ("x-", "x+", "y-", "y+", "z-", "z+")


@dataclass
class BlastWaveParams:
    # --- ambient plasma ---
    mass_ratio: float = 16.0  # m_i / m_e
    wpe_wce: float = 4.0  # omega_pe / Omega_ce of the uniform field (B0 = 1 / wpe_wce, along z)
    vte: float = 0.1  # sqrt(T_e / m_e c^2)
    ti_te: float = 1.0
    # --- debris ---
    debris: bool = True
    mach_alfven: float = 3.0  # V_d / v_A, v_A the Alfven speed of B0 with the ion mass
    rm: float = 8.0  # equal-mass radius R_m in ion skin depths
    shell_radius: float = 0.33  # in R_m
    shell_width: float = 0.1  # in R_m
    debris_beams: int = 4  # beams per debris species: the weight of a particle is 1 / debris_beams
    # --- body ---
    body: bool = True
    dipole_radius: float = 0.2  # distance r_dip at which the equatorial dipole field equals B0, in R_m
    body_radius: float = 1.0 / 1.5  # in r_dip
    core: float = 0.2  # regularisation length of the dipole in r_dip (inside the body)
    dipole: tuple = (0.0, 0.0, 1.0)  # direction of the moment
    # --- box and numerics ---
    half_width: float = 1.44  # of the cubic box, in R_B
    cells_per_de: float = 2.0
    ppc: int = 8
    run: float = 1.25  # run length in R_B / V_d
    filter_passes: int = 4  # binomial passes of the current filter (compensated)
    marder_zone: float = 0.125
    npml: int = 12
    capacity_factor: float = 1.5
    # --- output ---
    coarse: float = 0.1  # block size of the 3-D output in R_m
    slice_coarsen: int = 1  # block size of the two slices in cells
    seed: int = 12345
    overrides: dict = field(default_factory=dict)  # any SimConfig field


def params_from_dict(d: dict) -> BlastWaveParams:
    d = dict(d)
    if "dipole" in d:
        d["dipole"] = tuple(float(v) for v in d["dipole"])
    return BlastWaveParams(**d)


def _unit(v, name):
    n = math.sqrt(sum(c * c for c in v))
    if len(v) != 3 or n == 0.0:
        raise ValueError(f"{name} needs three components, not all zero")
    return tuple(c / n for c in v)


def _cells(length: float, dx: float) -> int:
    return max(8, int(round(length / dx / 8)) * 8)


def _block_cells(center, r_support, cfg) -> tuple:
    """Cells of the block a radial beam loads (``Simulation._beam_block``): in the box, and on the rank
    that holds the most of it (the box is split evenly along each axis)."""
    total, most = 1, 1
    for o, n, dx, d in zip(center, cfg.n_cells, cfg.dx, cfg.dims or (1, 1, 1)):
        lo = min(max(int(math.floor((o - r_support) / dx)), 0), n)
        hi = max(min(int(math.ceil((o + r_support) / dx)), n), lo)
        m = n // d
        total *= hi - lo
        most *= max(max(min(hi, (c + 1) * m) - max(lo, c * m), 0) for c in range(d))
    return total, most


def plan_blast_wave(p: BlastWaveParams, n_gpus: int = 1):
    b0 = 1.0 / p.wpe_wce
    v_a = mhd.alfven_speed(b0, p.mass_ratio)
    v_d = p.mach_alfven * v_a
    if not 0.0 < v_d < 1.0:
        raise ValueError(f"debris speed {v_d:g} c: lower mach_alfven or raise wpe_wce or mass_ratio")
    d_i = math.sqrt(p.mass_ratio)  # d_i = sqrt(m_i / m_e) d_e
    r_m = p.rm * d_i
    r_b = T.stopping_radius(r_m, p.mach_alfven)
    if p.half_width < HALF_WIDTH_MIN:
        raise ValueError(f"half_width = {p.half_width:g} R_B (limit {HALF_WIDTH_MIN:g}): the front leaves the box")
    dx = 1.0 / p.cells_per_de
    dt = 0.5 * dx
    n = _cells(2.0 * p.half_width * r_b, dx)
    n_cells = (n, n, n)
    L = n * dx
    center = (0.5 * L,) * 3
    B = (0.0, 0.0, b0)
    r_dip = p.dipole_radius * r_m
    radius, core = p.body_radius * r_dip, p.core * r_dip
    moment = b0 * r_dip ** 3
    m_dir = _unit(p.dipole, "dipole")
    wce_dt = b0 * dt
    if p.body:
        if radius < BODY_CELLS_MIN * dx:
            raise ValueError(f"body radius {radius / dx:.1f} cells (limit {BODY_CELLS_MIN:g}): raise dipole_radius, "
                             f"body_radius or cells_per_de")
        # the largest external field a particle can meet: the pole of the body's surface
        wce_dt = (2.0 * moment / radius ** 3 + b0) * dt
        if wce_dt > WCE_DT_MAX:
            raise ValueError(f"Omega_ce dt = {wce_dt:.2f} at the surface of the body (limit {WCE_DT_MAX:g}): "
                             f"raise body_radius or cells_per_de")
    t_b = r_b / v_d
    # ten 3-D outputs after the first; four slices and eight history records per 3-D output
    every_3d = max(8, int(round(p.run * t_b / dt / 80.0)) * 8)
    n_steps = 10 * every_3d
    every_slice, every_hist = every_3d // 4, every_3d // 8
    theta_e = p.vte ** 2
    theta_i = theta_e * p.ti_te / p.mass_ratio
    species = [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, (0.0, 0.0, 0.0), theta_e)]),
               SpeciesConfig("i", 1.0, p.mass_ratio, [Beam(1.0, (0.0, 0.0, 0.0), theta_i)])]
    ambient = ["e", "i"]
    r_s, w = p.shell_radius * r_m, p.shell_width * r_m
    n_pk = T.shell_peak_density(r_m, r_s, w)
    k_beams = int(p.debris_beams)
    shape = {"axis": "r", "shape": "gauss", "center": r_s, "width": w}
    wpe_dt = math.sqrt(1.0 + (n_pk if p.debris else 0.0)) * dt
    inside = T.shell_fraction_inside(r_s, w, radius) if p.debris and p.body else 0.0
    if p.debris:
        if k_beams < 1:
            raise ValueError("debris_beams must be 1 or more")
        if wpe_dt > WPE_DT_MAX:
            raise ValueError(f"omega_pe dt = {wpe_dt:.2f} at the peak of the shell (density {n_pk:.1f}, limit "
                             f"{WPE_DT_MAX:g}): widen the shell or raise cells_per_de")
        if inside > OVERLAP_MAX:
            raise ValueError(f"the shell overlaps the body: {100 * inside:.1f} % of its particles lie inside the "
                             f"body radius (limit {100 * OVERLAP_MAX:g} %): raise shell_radius or lower shell_width")
        # the same shape and positions for all beams of both species: neutral cell by cell; each beam
        # draws its own momenta, so debris_beams beams carry particles of 1 / debris_beams the weight
        u_d = v_d / math.sqrt(1.0 - v_d * v_d)
        for name, q, m, theta in (("de", -1.0, 1.0, theta_e), ("di", 1.0, p.mass_ratio, theta_i)):
            species.append(SpeciesConfig(name, q, m, [
                Beam(n_pk / k_beams, theta=theta, profile=dict(shape), origin=center, drift_radial=u_d)
                for _ in range(k_beams)]))
    names = [s.name for s in species]
    c3 = max(1, int(round(p.coarse * r_m / dx)))
    cs = int(p.slice_coarsen)
    moments = [f"{q}:{s}" for s in names for q in ("n", "flux", "ekin")] + ["E", "B", "rho"]
    in_plane = [f"n:{s}" for s in names] + ["flux:i", "E", "B"]
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
    omega_ci = b0 / p.mass_ratio
    dipoles = [Dipole(tuple(moment * c for c in m_dir), center, core)] if p.body else []
    cfg = SimConfig(
        n_cells=n_cells, dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc, seed=p.seed, dims=dims,
        boundaries=(("open", "open"),) * 3, npml=p.npml, capacity_factor=p.capacity_factor,
        external_fields=ExternalFields(B=B, dipoles=dipoles),
        body=Body(center=center, radius=radius) if p.body else None,
        injection=Injection(faces={f: "flux" for f in FACES}, species=ambient),
        current_filter=CurrentFilter(p.filter_passes, True), marder_zone=p.marder_zone,
        diag_interval=every_hist, gauss_interval=every_3d,
        output=OutputConfig(diagnostics=diags, timescales={"expansion": t_b, "ion_gyration": 2 * math.pi / omega_ci,
                                                           "plasma_period": 2 * math.pi}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    n_cell = math.prod(cfg.n_cells)
    n_part = len(ambient) * n_cell * cfg.ppc
    n_beams = sum(len(s.beams) for s in cfg.species)
    ids = n_beams * n_cell * cfg.ppc + cfg.n_steps * ids_per_step(cfg)
    if ids >= ID_LIMIT:
        raise ValueError(f"the run needs {ids:.3g} particle ids, the limit is 2^53: fewer steps, cells or ppc")
    gpus = max(n_gpus, 1)
    # Particle slots of one species on a rank, as the engine sizes them (PIC-D17): the cells the
    # species loads there times ppc and the capacity factor; of the block of the shell a rank counts
    # at least its even share. The debris numbers are those of the rank that holds the most of it.
    slots_a = cfg.capacity_factor * n_cell / gpus * cfg.ppc
    slots_d, n_debris, block = 0.0, 0, 0
    if p.debris:
        r_sup = support(shape, math.sqrt(sum(max(abs(o), abs(m * d - o)) ** 2
                                             for o, m, d in zip(center, cfg.n_cells, cfg.dx))))
        block, most = _block_cells(center, r_sup, cfg)
        cells = n_cell / gpus if block == n_cell else max(most, -(-block // gpus))
        slots_d = cfg.capacity_factor * k_beams * cells * cfg.ppc
        # particles kept after the cut of the shape: the cells of the ball of the support
        n_debris = int(2 * k_beams * cfg.ppc * min(block, 4.0 / 3.0 * math.pi * r_sup ** 3 / math.prod(cfg.dx)))
    if max(slots_a, slots_d) + 4096 > MAX_CAPACITY:
        raise ValueError(f"{max(slots_a, slots_d):.3g} particle slots per species and GPU exceed 2^31 - 1 "
                         f"(the kernels' particle index): more GPUs, fewer cells or a lower ppc")
    # Device memory of a rank whose particle arrays are full (terms measured at the bow-shock gate): 8
    # rows of float64 per particle slot; 40 bytes per particle of one species during its sort; 31
    # arrays of the grid and the arrays an output step holds (one per component of the 3-D output),
    # each with the ghost and absorbing layers; 1 GB of CUDA context, kernels and NCCL buffers.
    slots = len(ambient) * slots_a + 2 * slots_d
    padded = math.prod(m / d + (2 if d == 1 else 1) * cfg.npml + 2 * GHOST
                       for m, d in zip(cfg.n_cells, cfg.dims or (1, 1, 1)))
    n_out = sum(3 if m.split(":")[0] in ("flux", "E", "B") else 1 for m in moments)
    mem = (slots * 64 + max(slots_a, slots_d) * 40 + padded * 8 * (31 + n_out)) / 1e9 + 1.0
    beta = mhd.plasma_beta(p.vte, p.ti_te, b0)
    v_fast = mhd.fast_speed(v_a, beta)
    info = {
        "debris": p.debris, "body": p.body, "d_i": d_i, "r_m": r_m, "r_b": r_b, "center": list(center),
        "b0": b0, "v_a": v_a, "v_d": v_d, "v_fast": v_fast, "beta": beta, "mach_alfven": p.mach_alfven,
        "mach_fast": v_d / v_fast, "debris_gyroradius": p.mach_alfven * d_i,
        "shell_radius": r_s, "shell_width": w, "shell_peak_density": n_pk, "debris_beams": k_beams,
        "debris_extent_le": T.debris_extent(r_m, p.mach_alfven), "shell_inside_body": inside,
        "wpe_dt_peak": wpe_dt, "debris_particles": n_debris, "debris_block_cells": block,
        "debris_slots": slots_d, "ambient_slots": slots_a,
        "dipole_radius": r_dip, "body_radius": radius, "core": core, "moment": moment, "dipole": list(m_dir),
        "t_b": t_b, "t_b_steps": t_b / dt, "omega_ci_t_end": omega_ci * cfg.n_steps * cfg.dt,
        "wce_dt_surface": wce_dt, "debye_cells": p.vte / dx, "particles": n_part, "memory_per_gpu_GB": mem,
        "id_fraction": ids / ID_LIMIT, "coarse_cells": c3, "every_3d": every_3d, "every_slice": every_slice,
        "dims": cfg.dims, "params": asdict(p)}
    lines = [
        "=== pypic blast-wave planner (debris shell around a magnetised body) ===",
        f"ambient: m_i/m_e={p.mass_ratio:g}, B0={b0:g} along z, beta={beta:.2f}, v_te={p.vte:g} c "
        f"(Debye length {p.vte / dx:.2f} dx), v_A={v_a:.4f} c, fast speed {v_fast:.4f} c",
        f"debris ({'on' if p.debris else 'OFF'}): M_A={p.mach_alfven:g} (fast Mach {v_d / v_fast:.2f}), "
        f"V_d={v_d:.4f} c, R_m={r_m:g} d_e = {p.rm:g} d_i, R_B={r_b:.1f} d_e, gyroradius {p.mach_alfven * d_i:g} "
        f"({p.mach_alfven / p.rm:.3f} R_m); shell at {r_s:.2f}, width {w:.2f}, peak density {n_pk:.1f} in "
        f"{k_beams} beam(s) per species (omega_pe dt <= {wpe_dt:.2f}, {100 * inside:.2f} % inside the body); "
        f"extent after Le et al. 2021, Eq. 4 with D = 3: {info['debris_extent_le']:.1f} d_e",
        (f"body: radius {radius:.2f} ({radius / dx:.1f} cells) at the centre ({center[0]:g}); dipole moment "
         f"{moment:.1f} along {tuple(round(c, 3) for c in m_dir)}, field B0 at {r_dip:g}, core {core:g}; "
         f"Omega_ce dt <= {wce_dt:.2f} at the surface") if p.body else "body: OFF (no sphere, no dipole)",
        f"grid {cfg.n_cells}, dx={dx:g}, L={L:g} d_e (half-width {0.5 * L / r_b:.2f} R_B), dt={cfg.dt:g}, "
        f"{cfg.n_steps} steps = {cfg.n_steps * cfg.dt / t_b:.2f} R_B/V_d = {info['omega_ci_t_end']:.1f} / Omega_ci",
        f"{n_part:,} ambient particles ({cfg.ppc} per cell and species), about {n_debris:,} of the debris "
        f"(block of {block:,} cells), capacity factor {cfg.capacity_factor:g}; "
        f"memory up to {mem:.1f} GB/GPU on {gpus} GPU(s), rank grid {cfg.dims}; ids {info['id_fraction']:.1e} of 2^53",
        f"faces: six open faces with flux injection of the ambient species; filter {p.filter_passes}+comp, "
        f"Marder zone {cfg.marder_zone:g}, npml {cfg.npml}",
        f"output: 3-D moments in blocks of {c3} cells every {every_3d} steps, slices z and y through the centre "
        f"every {every_slice} steps, history every {every_hist} steps",
    ]
    return cfg, "\n".join(lines), info
