"""Planner for a relativistic pair shock, unmagnetised or magnetised (reflecting-wall / piston method).

A cold pair plasma flowing at -beta0 (Lorentz factor gamma0) fills 0 <= x < L and hits the
reflecting wall at x = 0; the reflected and incoming streams are Weibel unstable and a shock forms
that moves back into the flow (the wall frame is the downstream frame). The right wall is
reflecting too: a vacuum gap opens there at ~c, so the run is useful until it meets the shock,
t < L / (1 + beta_sh). A compensated binomial current filter suppresses the numerical Cherenkov
instability of the relativistic flow (PIC-D10).

With sigma > 0 the flow carries a field B1 perpendicular to it and the motional field E = -v x B1,
which is tangential to the walls (PIC-D18): the left wall is a conductor at rest, the piston of a
perpendicular magnetised shock. field = "out": B1 along z, out of the simulation plane (the momenta
stay in the plane: 2-D gas); "in": B1 along y (the gyration takes them out of the plane: 3-D gas).

Theory (theory.py): jump conditions of a Synge gas (2-D for a box thinner than a few skin depths
along z, 3-D otherwise), with the magnetisation sigma = B1^2 / (gamma0 n1 m) (n1 m: lab mass
density of both species); Fermi acceleration tail p ~ 2.4 (2-D, Spitkovsky 2008b) is reported.
"""
from __future__ import annotations

import dataclasses
import math
import os
import sys
from dataclasses import asdict, dataclass

from pypic.config import Beam, CurrentFilter, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402


@dataclass
class ShockParams:
    gamma0: float = 15.0
    theta: float = 1e-4  # upstream temperature in its rest frame (T / m c^2)
    t_end: float = 400.0  # 1/omega_pe (omega_pe of the upstream lab-frame electron density)
    cells_per_skin: float = 10.0  # cells per relativistic skin depth sqrt(gamma0) c/omega_pe
    Ly_skin: float = 8.0  # transverse size in relativistic skin depths
    nz: int = 4  # thin along z (quasi-2-D); set nz_skin for 3-D
    nz_skin: float | None = None
    ppc: int = 4
    filter_passes: int = 4
    n_tracks: int = 20000
    sigma: float = 0.0  # magnetisation B1^2 / (gamma0 n1 m) of the flow; 0: no field
    field: str = "out"  # B1 out of the simulation plane (along z) or in it (along y); used with sigma > 0
    seed: int = 2024
    overrides: dict = dataclasses.field(default_factory=dict)


def params_from_dict(d: dict) -> ShockParams:
    return ShockParams(**d)


def plan_shock(p: ShockParams, n_gpus: int = 1):
    g0 = p.gamma0
    u0 = math.sqrt(g0 * g0 - 1)
    skin = math.sqrt(g0)  # relativistic skin depth in c/omega_pe
    dx = skin / p.cells_per_skin
    beta0 = u0 / g0
    if p.sigma < 0 or p.field not in ("out", "in"):
        raise ValueError("sigma must be >= 0 and field one of 'out', 'in'")
    thin = p.nz_skin is None or p.nz_skin < 2
    # dimensions of the downstream momenta: an in-plane field turns them out of the plane of a thin box
    dims_eff = 2 if thin and not (p.sigma and p.field == "in") else 3

    def jump(dims):
        return T.jump_conditions_magnetised(g0, p.sigma, dims) if p.sigma else T.jump_conditions(g0, dims)

    jc = jump(dims_eff)
    n_mass = 2.0  # lab mass density of the two species (n = 1, m = 1 each)
    B1 = math.sqrt(p.sigma * g0 * n_mass)
    b_axis = 2 if p.field == "out" else 1
    B0, E0 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
    B0[b_axis] = B1
    E0[3 - b_axis] = (-beta0 if b_axis == 2 else beta0) * B1  # E = -v x B for v = (-beta0, 0, 0)
    L = (1.0 + jc["beta_sh"]) * p.t_end * 1.05
    n_ranks = max(n_gpus, 1)
    nx = int(math.ceil(L / dx / (8 * n_ranks))) * 8 * n_ranks
    ny = int(math.ceil(p.Ly_skin * skin / dx / 8)) * 8
    nz = p.nz if p.nz_skin is None else int(math.ceil(p.nz_skin * skin / dx / 4)) * 4
    dt = 0.5 * dx  # c dt / dx = 0.5
    n_steps = int(math.ceil(p.t_end / dt))
    flow = [Beam(1.0, (-u0, 0.0, 0.0), p.theta)]
    species = [SpeciesConfig("electron", -1.0, 1.0, flow),
               SpeciesConfig("positron", 1.0, 1.0, [Beam(1.0, (-u0, 0.0, 0.0), p.theta)])]
    Lx = nx * dx
    n_e = nx * ny * nz * p.ppc
    diags = [
        {"name": "density", "kind": "lineout", "axis": "x",
         "fields": ["n:electron", "n:positron", "Bz", "Ex"] + (["By"] if B1 and b_axis == 1 else []),
         "every": {"count": 100}},
        {"name": "x_gamma", "kind": "histogram", "quantities": ["x", "gamma"], "bins": [256, 120],
         "ranges": [[0.0, Lx], [1.0, 100 * g0]], "log": [False, True], "every": {"count": 50}},
        {"name": "phase", "kind": "histogram", "species": ["electron"], "quantities": ["x", "ux"],
         "bins": [512, 256], "ranges": [[0.0, Lx], [-2 * u0, 2 * u0]], "every": {"count": 100}},
        {"name": "slice_xy", "kind": "slice", "fields": ["Bz", "n:electron"], "axis": "z", "coarsen": [2, 2, 1],
         "every": {"count": 60}},
        {"name": "tracks", "kind": "tracks", "species": ["electron"], "stride": max(1, n_e // p.n_tracks),
         "with_fields": True, "every": {"count": 400}},
    ]
    bc = (("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))
    try:
        dims = auto_dims(n_gpus, (nx, ny, nz))
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=(nx, ny, nz), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims, boundaries=bc, diag_interval=max(1, n_steps // 500),
                    gauss_interval=max(1, n_steps // 25), sort_interval=50,
                    current_filter=CurrentFilter(p.filter_passes, True), B0=tuple(B0), E0=tuple(E0),
                    output=OutputConfig(diagnostics=diags, timescales={"plasma_period": 2 * math.pi * skin,
                                                                       "t_end": p.t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info = {"gamma0": g0, "u0": u0, "skin": skin, "dims_eff": dims_eff, "jump": jc, "jump_2d": jump(2),
            "jump_3d": jump(3), "jump_gamma_3_2": T.jump_conditions_gamma(g0, 1.5),
            "t_valid": Lx / (1 + jc["beta_sh"]), "Lx": Lx, "sigma": p.sigma, "B1": B1, "field_component": "xyz"[b_axis],
            "params": asdict(p)}
    lines = ["=== pypic relativistic-shock planner (reflecting wall) ===",
             f"pair flow gamma0={g0:g} (u0={u0:.3f}), theta={p.theta:g}; "
             f"relativistic skin depth {skin:.2f} c/omega_pe",
             f"grid {cfg.n_cells}, dx={dx:.3f} ({p.cells_per_skin:g}/skin), "
             f"L=({Lx:.0f}, {ny * dx:.1f}, {nz * dx:.1f}), "
             f"dt={dt:.4f}, {n_steps} steps to t={p.t_end:g}; {2 * n_e:,} particles; filter {p.filter_passes}+comp",
             *([f"magnetised: sigma={p.sigma:g}, B1={B1:.3f} along {'xyz'[b_axis]} ({p.field.replace('out', 'out of')} "
                f"the plane), E{'xyz'[3 - b_axis]}={E0[3 - b_axis]:.3f} tangential to the walls (PIC-D18); "
                f"Larmor radius of the flow {u0 / B1 / dx:.1f} cells, Omega_c dt={B1 / g0 * dt:.3f} upstream"]
               if B1 else []),
             f"jump conditions ({dims_eff}-D Synge gas{f', sigma={p.sigma:g}' if B1 else ''}): "
             f"beta_sh={jc['beta_sh']:.3f}, n2/n1={jc['compression']:.3f}, theta2={jc['theta2']:.3f}"
             + ("; B2/B1 = n2/n1" if B1 else f"; Gamma=3/2: beta_sh={info['jump_gamma_3_2']['beta_sh']:.3f}, "
                                             f"n2/n1={info['jump_gamma_3_2']['compression']:.3f}"),
             f"valid until the vacuum gap reaches the shock: t < {info['t_valid']:.0f}"]
    return cfg, "\n".join(lines), info
