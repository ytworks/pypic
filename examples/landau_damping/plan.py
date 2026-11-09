"""Planner for Langmuir-oscillation experiments: Landau damping (warm) and the exact cold-plasma
oscillation / wave breaking (Dawson).

cases
    landau        warm electrons, density perturbation n = n0 (1 + alpha cos k x): linear Landau
                  damping for small alpha (gate: k lambda_De = 0.5, alpha = 0.02)
    nonlinear     the same with large alpha (trapping, O'Neil; e.g. alpha = 0.5)
    dawson        cold electrons started at equilibrium with v = v0 sin(k x0), k v0 < omega_pe:
                  every sheet oscillates at exactly omega_pe with amplitude v0 sin(k x0) / omega_pe
    wavebreaking  the same with k v0 > omega_pe: sheets cross at omega_pe t = arcsin(omega_pe / k v0)
    relativistic  cold, v0 relativistic: period of the relativistic sheet oscillator (beta_max ~ 0.3)
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

from pypic.config import Beam, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

CASES = ("landau", "nonlinear", "dawson", "wavebreaking", "relativistic")


@dataclass
class LangmuirParams:
    case: str = "landau"
    k_lambda_D: float = 0.5  # warm cases
    alpha: float = 0.02  # density perturbation (warm cases)
    theta_e: float = 1e-3  # T_e / m_e c^2 (warm cases)
    k_v0: float = 0.5  # cold cases: k v0 / omega_pe (crossing for >= 1)
    u0: float = 0.1  # cold cases: proper velocity amplitude (c)
    mass_ratio: float = 1836.15267
    cells_per_wavelength: int = 32
    n_transverse: int = 16
    ppc: int = 1024
    t_end: float | None = None
    n_tracks: int = 4000
    seed: int = 777
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> LangmuirParams:
    return LangmuirParams(**d)


def plan_langmuir(p: LangmuirParams, n_gpus: int = 1):
    if p.case not in CASES:
        raise ValueError(f"case must be one of {CASES}")
    warm = p.case in ("landau", "nonlinear")
    info = {"case": p.case, "params": asdict(p)}
    if warm:
        vt = math.sqrt(p.theta_e)
        k = p.k_lambda_D / vt
        from pypic.theory.waves import langmuir

        w = complex(langmuir([k], vt)[0])
        info.update(k=k, vt=vt, omega_theory=w.real, gamma_theory=w.imag, alpha=p.alpha)
        t_end = p.t_end or (min(8.0 / abs(w.imag), 60.0) if p.case == "landau" else 100.0)
        theta = p.theta_e
        pert = [{"axis": "x", "mode": 1, "density": p.alpha}]
        w_max = 1.3 * w.real
        timescales = {"plasma_period": 2 * math.pi, "shortest_period": 2 * math.pi / w_max,
                      "growth_time": 1.0 / abs(w.imag)}
    else:
        v0 = p.u0 / math.sqrt(1 + p.u0 * p.u0)
        k = p.k_v0 / v0
        info.update(k=k, v0=v0, u0=p.u0, k_v0=p.k_v0)
        t_end = p.t_end or {"dawson": 40.0, "wavebreaking": 15.0, "relativistic": 60.0}[p.case]
        theta = 0.0
        pert = [{"axis": "x", "mode": 1, "velocity": p.u0}]
        timescales = {"plasma_period": 2 * math.pi, "shortest_period": 2 * math.pi}
    L = 2 * math.pi / k
    nx = p.cells_per_wavelength
    dx = L / nx
    nt = p.n_transverse
    dt = 0.95 / math.sqrt(3.0) * dx
    dt = min(dt, 0.05)
    n_steps = int(math.ceil(t_end / dt))
    species = [SpeciesConfig("electron", -1.0, 1.0, [Beam(1.0, (0.0, 0.0, 0.0), theta, perturbation=pert)]),
               SpeciesConfig("ion", 1.0, p.mass_ratio, [Beam(1.0, (0.0, 0.0, 0.0), 0.0)])]
    n_e = nx * nt * nt * p.ppc
    diags = [
        {"name": "modes", "kind": "modes", "fields": ["Ex"], "modes": [[m, 0, 0] for m in range(1, 9)]},
        {"name": "lineout", "kind": "lineout", "fields": ["Ex"], "axis": "x"},
        {"name": "phase", "kind": "histogram", "species": ["electron"], "quantities": ["x", "ux"], "bins": [128, 128],
         "ranges": [[0.0, L], None], "every": {"count": 150}},
        {"name": "spectrum", "kind": "histogram", "species": ["electron"], "quantities": ["ekin"], "bins": [120],
         "ranges": [[1e-8, 1.0]], "log": [True], "every": {"count": 60}},
        {"name": "tracks", "kind": "tracks", "species": ["electron"], "stride": max(1, n_e // p.n_tracks)},
    ]
    if p.case == "wavebreaking":  # time the first crossing (t ~ 1) to 1 %: dense track output early on
        diags[-1].update(every={"t": 0.01}, t_end=3.0)
    try:
        dims = auto_dims(n_gpus, (nx, nt, nt))
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=(nx, nt, nt), dx=(dx, dx, dx), dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims, diag_interval=max(1, n_steps // 1000),
                    gauss_interval=max(1, n_steps // 50), sort_interval=50,
                    output=OutputConfig(diagnostics=diags, timescales={**timescales, "t_end": t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info.update(L=L, dx=dx, dt=cfg.dt, t_end=t_end)
    lines = ["=== pypic Langmuir-oscillation planner ===", f"case={p.case}: k={k:.4f} c^-1 omega_pe, L={L:.4f} "
             f"({nx} cells, dx={dx:.5f}), dt={cfg.dt:.5f}, t_end={t_end:g} ({n_steps} steps), "
             f"{n_e:,} electrons ({p.ppc}/cell)"]
    if warm:
        lines.append(f"Landau theory: k lambda_De={p.k_lambda_D:g}, omega={info['omega_theory']:.4f}, "
                     f"gamma={info['gamma_theory']:.4f} omega_pe; alpha={p.alpha:g}; dx/lambda_De={dx / vt:.3f}")
    else:
        lines.append(f"cold: u0={p.u0:g} (v0={v0:.4f} c), k v0={p.k_v0:g} omega_pe")
    return cfg, "\n".join(lines), info
