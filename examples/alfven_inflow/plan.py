"""Planner for an Alfvénic wave carried by a flow through an inflow face, optionally into a shock.

Electrons and ions flow at -V0 along x through B0 = B y (motional field E0 = V0 B z). The face x+
injects the flow (`flux` face, PIC-D13) and carries one cold two-fluid wave along B0 (k = k y, an R
or an L wave, |dB| / B0 = amplitude): the fluctuation table of PIC-D26 (theory.fluctuation). The
pattern is uniform along x, so the plasma carries it in unchanged.

- ``wall: false`` (the gate): x- is an open face where the plasma leaves. The box holds the same
  wave at t = 0 (``load: uniform``: profiles along y, the electrons displaced along z for E_z).
- ``wall: true`` (the piston): x- is a reflecting wall (E0 tangential to it is allowed because the
  open face of the axis injects the flow that carries it, PIC-D26). The flow piles up into a
  perpendicular shock that moves into the plasma that carries the wave. ``load: envelope`` loads the
  wave times an envelope h(x) that vanishes within ``wall_ramp`` of the wall and within
  ``face_ramp`` of the face (PIC-D24: vector potential, shaped velocities, displaced electrons; see
  theory.preload for why h vanishes at the face too). The face carries the wave from t = 0 on.

Quasi-2-D: two cells along z. The numbers of the plasma (theta_e = 1e-3, Debye length dx / 2, 32
particles per cell) are those of tests/gpu/test_inflow.py, whose measured errors the gate criteria
of the README follow.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

from pypic.config import Beam, CurrentFilter, Injection, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup import mhd
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402

LOADS = ("uniform", "envelope", "none")


@dataclass
class AlfvenInflowParams:
    # --- plasma ---
    mass_ratio: float = 16.0
    wpe_wce: float = 10.0  # omega_pe / Omega_ce: B0 = 1 / wpe_wce along y
    theta_e: float = 1e-3  # T_e / m_e c^2
    ti_te: float = 1.0
    flow_mach: float = 2.0  # V0 / v_A in the box (wall) frame
    # --- wave ---
    branch: str = "R"
    amplitude: float = 0.2  # |delta B| / B0
    wave: bool = True  # false: the null run, a table of zero amplitude (the analysis normalises by the wave)
    ly: float = 32.0  # box length along y = one wavelength, in c/omega_pe
    load: str = "uniform"  # the wave in the box at t = 0: uniform | envelope | none
    wall_ramp: tuple = (4.0, 12.0)  # envelope: 0 below the first, 1 above the second (c/omega_pe)
    face_ramp: float = 8.0  # envelope: falls to 0 over this distance before the face
    # --- box and numerics ---
    wall: bool = False  # reflecting x- wall (piston) instead of an open outflow face
    lx: float = 64.0
    nz: int = 2
    cells_per_de: float = 4.0
    cfl: float = 0.4  # c dt / dx
    ppc: int = 32
    t_end: float = 300.0
    filter_passes: int = 4
    marder_zone: float = 0.125
    npml: int = 12
    seed: int = 12345
    # --- output ---
    frames: int = 30  # slice outputs
    profiles: int = 100  # lineout outputs
    spectra: int = 20  # histogram outputs (piston)
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> AlfvenInflowParams:
    d = dict(d)
    if "wall_ramp" in d:
        d["wall_ramp"] = tuple(float(v) for v in d["wall_ramp"])
    return AlfvenInflowParams(**d)


def plan_alfven_inflow(p: AlfvenInflowParams, n_gpus: int = 1):
    if p.load not in LOADS:
        raise ValueError(f"load must be one of {LOADS}")
    if p.load == "envelope" and not p.wall:
        raise ValueError("load: envelope is for the piston (wall: true); the open box loads the uniform wave")
    if p.load == "uniform" and p.wall and p.wave:
        raise ValueError("a uniform wave has B_x on the reflecting wall, which zeroes it (div B): use the envelope")
    b = 1.0 / p.wpe_wce
    species = T.two_fluid(p.mass_ratio)
    v_a = T.alfven_speed(b, species)
    v0 = p.flow_mach * v_a
    if not 0.0 < v0 < 0.5:
        raise ValueError(f"flow speed {v0:g} c: give 0 < flow_mach v_A < 0.5 c")
    dx = 1.0 / p.cells_per_de
    dt = p.cfl * dx
    nx, ny = int(round(p.lx / dx)), int(round(p.ly / dx))
    if abs(nx * dx - p.lx) > 1e-9 or abs(ny * dx - p.ly) > 1e-9:
        raise ValueError("lx and ly must be whole numbers of cells")
    k = 2.0 * math.pi / p.ly
    wave = T.parallel_wave(k, b, species, p.branch, p.amplitude)
    table = T.fluctuation(wave, v0)
    # the null run keeps a table of zero amplitude: the face then holds the uniform motional field E0
    # (PIC-D26 allows E0 tangential to the wall only with a fluctuation on the open face of the axis)
    injected = table if p.wave else T.fluctuation(T.parallel_wave(k, b, species, p.branch, 0.0), v0)
    h = None
    if p.load == "envelope":
        w0, w1 = p.wall_ramp
        if not (0.0 <= w0 < w1 and w1 + p.face_ramp < p.lx and p.face_ramp > 0.0):
            raise ValueError("envelope: need 0 <= wall_ramp[0] < wall_ramp[1] < lx - face_ramp, face_ramp > 0")
        h = T.envelope(p.lx, (w0, w1), p.face_ramp)
    fields, pert = T.preload(table, p.ly, dx, h) if p.load != "none" and p.wave else ([], {"e": [], "i": []})
    u0 = v0 / math.sqrt(1.0 - v0 * v0)
    theta_i = p.theta_e * p.ti_te / p.mass_ratio
    sp = [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, (-u0, 0.0, 0.0), p.theta_e, perturbation=pert["e"])]),
          SpeciesConfig("i", 1.0, p.mass_ratio, [Beam(1.0, (-u0, 0.0, 0.0), theta_i, perturbation=pert["i"])])]
    n_steps = int(math.ceil(p.t_end / dt))
    beta = mhd.plasma_beta(math.sqrt(p.theta_e), p.ti_te, b)
    shock = T.piston_shock(v0, v_a, beta, mass=p.mass_ratio + 1.0) if p.wall else None
    plane = ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "n:e", "n:i"]
    diags = [
        {"name": "plane", "kind": "slice", "axis": "z", "position": 0.0, "fields": plane,
         "every": {"count": p.frames}},
        {"name": "profile", "kind": "lineout", "axis": "x", "average": False,
         "fields": ["n:e", "n:i", "By", "Ez", "flux:i", "ekin:e", "ekin:i"], "every": {"count": p.profiles}},
    ]
    if p.wall:
        diags.append({"name": "spectrum", "kind": "histogram", "quantities": ["x", "ekin"],
                      "bins": [nx // 4, 120], "ranges": [[0.0, p.lx], [1e-5, 1.0]], "log": [False, True],
                      "every": {"count": p.spectra}})
    bc = (("reflecting" if p.wall else "open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    try:
        dims = auto_dims(n_gpus, (nx, ny, p.nz))
    except ValueError:
        dims = None
    omega_ci = b / p.mass_ratio
    cfg = SimConfig(
        n_cells=(nx, ny, p.nz), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=sp, ppc=p.ppc, seed=p.seed, dims=dims,
        boundaries=bc, npml=p.npml, current_filter=CurrentFilter(p.filter_passes, True), marder_zone=p.marder_zone,
        B0=(0.0, b, 0.0), E0=(0.0, 0.0, v0 * b), initial_fields=fields,
        injection=Injection(faces={"x+": "flux"}, fluctuation=injected),
        # the downstream of a piston holds r times the density on most of a rank (the arrays grow if needed)
        capacity_factor=1.2 * shock["compression"] if shock else 1.3,
        diag_interval=max(1, n_steps // 500), gauss_interval=max(1, n_steps // 50), sort_interval=50,
        output=OutputConfig(diagnostics=diags, timescales={"ion_gyration": 2 * math.pi / omega_ci,
                                                           "wave_period": 2 * math.pi / wave["omega"]}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    d_i = T.ion_skin_depth(species)
    info = {
        "b0": b, "v_a": v_a, "v0": v0, "u0": u0, "beta": beta, "d_i": d_i, "k": k, "k_di": k * d_i,
        "omega": wave["omega"], "omega_over_kva": wave["omega"] / (k * v_a), "branch": p.branch,
        "amplitude": p.amplitude, "wave": p.wave, "table": table, "omega_ci": omega_ci,
        "omega_ci_t_end": omega_ci * n_steps * dt,
        "transit_time": p.lx / v0, "wall": p.wall, "load": p.load, "lx": p.lx, "ly": p.ly,
        "envelope": h, "shock": shock, "debye_cells": math.sqrt(p.theta_e) / dx,
        "particles": 2 * nx * ny * p.nz * p.ppc, "dims": cfg.dims, "params": asdict(p)}
    lines = [
        "=== pypic Alfvénic inflow planner (a wave carried in by the flow) ===",
        f"plasma: m_i/m_e={p.mass_ratio:g}, B0={b:g} along y, beta={beta:.2f}, theta_e={p.theta_e:g} (Debye length "
        f"{info['debye_cells']:.2f} dx), v_A={v_a:.4f} c, d_i={d_i:g}",
        f"flow: V0={v0:.4f} c = {p.flow_mach:g} v_A along -x; E0 = V0 B0 along z",
        f"wave: cold two-fluid {p.branch} wave, k d_i={k * d_i:.3f}, omega={wave['omega']:.5f} = "
        f"{info['omega_over_kva']:.3f} k v_A, |dB|/B0={p.amplitude:g}; "
        + (f"loaded: {p.load}; carried in through x+" if p.wave else "OFF (null run)"),
        f"grid {cfg.n_cells}, dx={dx:g}, dt={dt:g}, {n_steps} steps to t={n_steps * dt:g} "
        f"({info['omega_ci_t_end']:.2f} / Omega_ci, {n_steps * dt / info['transit_time']:.2f} transit times "
        f"lx / V0); {info['particles']:,} particles at t = 0",
        f"faces: x- {'reflecting wall' if p.wall else 'open'}, x+ flux inflow face with the fluctuation; "
        f"filter {p.filter_passes}+comp, Marder zone {p.marder_zone:g}, npml {p.npml}",
    ]
    if shock:
        lines.append(f"perpendicular MHD shock (gamma = 5/3): r={shock['compression']:.3f}, "
                     f"V_sh={shock['v_shock']:.4f} c, "
                     f"M_A={shock['mach_alfven']:.2f} (fast Mach {shock['mach_fast']:.2f}) in the shock frame")
    return cfg, "\n".join(lines), info
