"""Planner for the plasma-wave dispersion experiment: a thermal plasma (optionally magnetised)
whose thermal fluctuations excite every normal mode; the omega-k spectrum of transverse-averaged
lineouts is compared with linear theory.

Geometry: k along x. ``geometry`` = unmagnetized | perpendicular (B0 along z: X mode and electron
Bernstein waves in Ex/Ey, O mode in Ez) | parallel (B0 along x: R/L/whistler in Ey/Ez, Langmuir in Ex)
| ion_acoustic (unmagnetised, T_e >> T_i, light ions, long run).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

from pypic.config import Beam, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

GEOMETRIES = ("unmagnetized", "perpendicular", "parallel", "ion_acoustic")


@dataclass
class WaveParams:
    geometry: str = "unmagnetized"
    theta_e: float = 0.01  # T_e / m_e c^2 (v_te = 0.1 c, lambda_De = 0.1 c/omega_pe)
    Ti_over_Te: float = 0.01
    mass_ratio: float = 1836.15267
    wce: float = 0.5  # omega_ce / omega_pe (magnetised geometries)
    dx_over_lambda_De: float = 1.0
    nx: int = 512
    n_transverse: int = 8  # cells along y and z (the lineouts average over them)
    ppc: int = 64
    t_end: float = 300.0  # run length in 1/omega_pe: frequency resolution 2 pi / t_end
    k_window: float = 3.0  # largest k c/omega_pe of the electromagnetic branches analysed
    k_bernstein_max: float = 40.0  # perpendicular: Bernstein waves analysed for 0.3 <= k rho_e and k <= this
    cfl: float = 0.95  # c dt / dx in units of 1/sqrt(3)
    seed: int = 4242
    output_budget_GB: float | None = None
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> WaveParams:
    return WaveParams(**d)


def frequency_window(p: WaveParams) -> float:
    """Highest frequency the analysis must resolve (sets the output cadence)."""
    k = p.k_window
    if p.geometry == "ion_acoustic":
        return 2.0 * math.sqrt(1.0 / p.mass_ratio)  # a few omega_pi
    w = math.sqrt(1.0 + k * k)  # EM / O branch
    if p.geometry in ("perpendicular", "parallel"):
        w = max(w, math.sqrt(1.0 + p.wce ** 2 + k * k) + p.wce, 5 * p.wce)
    return 1.1 * max(w, math.sqrt(1.0 + 3 * k * k * p.theta_e))


def plan_waves(p: WaveParams, n_gpus: int = 1):
    if p.geometry not in GEOMETRIES:
        raise ValueError(f"geometry must be one of {GEOMETRIES}")
    lam_De = math.sqrt(p.theta_e)
    dx = p.dx_over_lambda_De * lam_De
    nx, nt = p.nx, p.n_transverse
    dt = p.cfl / math.sqrt(3.0) * dx
    if p.geometry != "ion_acoustic":
        dt = min(dt, 0.2)
    n_steps = int(math.ceil(p.t_end / dt))
    B0 = (0.0, 0.0, 0.0)
    if p.geometry == "perpendicular":
        B0 = (0.0, 0.0, p.wce)
    elif p.geometry == "parallel":
        B0 = (p.wce, 0.0, 0.0)
    Ti = p.theta_e * p.Ti_over_Te
    species = [SpeciesConfig("electron", -1.0, 1.0, [Beam(1.0, (0.0, 0.0, 0.0), p.theta_e)]),
               SpeciesConfig("ion", 1.0, p.mass_ratio, [Beam(1.0, (0.0, 0.0, 0.0), Ti / p.mass_ratio)])]
    w_max = frequency_window(p)
    comps = ["Ex", "Ey", "Ez", "By", "Bz"]
    diags = [
        {"name": "lineout", "kind": "lineout", "fields": comps, "axis": "x"},  # auto: Nyquist of w_max, averaged
        {"name": "fields", "kind": "fields", "fields": ["E", "B"], "every": {"count": 3}},
        {"name": "spectrum", "kind": "histogram", "quantities": ["ekin"], "bins": [100], "log": [True],
         "ranges": [[1e-5, 1.0]], "every": {"count": 10}},
    ]
    timescales = {"plasma_period": 2 * math.pi, "shortest_period": 2 * math.pi / w_max, "t_end": p.t_end}
    try:
        dims = auto_dims(n_gpus, (nx, nt, nt))
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=(nx, nt, nt), dx=(dx, dx, dx), dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims, B0=B0, diag_interval=max(1, n_steps // 500),
                    gauss_interval=max(1, n_steps // 50), sort_interval=50,
                    output=OutputConfig(diagnostics=diags, timescales=timescales, budget_GB=p.output_budget_GB))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info = {"geometry": p.geometry, "theta_e": p.theta_e, "Ti": Ti, "mass_ratio": p.mass_ratio,
            "wce": p.wce if p.geometry in ("perpendicular", "parallel") else 0.0, "lambda_De": lam_De,
            "dx": dx, "dt": cfg.dt, "k_window": p.k_window, "k_bernstein_max": p.k_bernstein_max, "w_max": w_max,
            "frequency_resolution": 2 * math.pi / p.t_end, "params": asdict(p)}
    lines = [
        "=== pypic plasma-wave planner ===",
        f"geometry={p.geometry}: theta_e={p.theta_e:g} (lambda_De={lam_De:.3f} c/omega_pe), T_i/T_e={p.Ti_over_Te:g}, "
        f"m_i/m_e={p.mass_ratio:g}" + (f", omega_ce={p.wce:g} omega_pe, B0={B0}" if any(B0) else ""),
        f"grid {cfg.n_cells}, dx={dx:.4f} ({p.dx_over_lambda_De:g} lambda_De), L_x={nx * dx:.2f}: "
        f"dk={2 * math.pi / (nx * dx):.4f}, k_Nyquist={math.pi / dx:.1f}",
        f"dt={cfg.dt:.4f} (c dt/dx={cfg.dt / dx:.3f}), t_end={p.t_end:g}: d omega={2 * math.pi / p.t_end:.4f} "
        f"omega_pe; {n_steps} steps",
        f"analysis window k <= {p.k_window:g}, omega <= {w_max:.2f}: lineouts sampled at its Nyquist rate "
        "(time-averaged, anti-aliased)",
        f"particles: {nx * nt * nt * p.ppc * 2:,} ({p.ppc} per cell x 2 species)",
    ]
    return cfg, "\n".join(lines), info
