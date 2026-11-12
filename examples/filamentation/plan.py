"""Planner for the current-filamentation (Weibel) instability of counter-streaming beams.

Beams drift along +-x; the box is thin along x (its smallest k_x is far above the two-stream range,
so the electrostatic modes are excluded) and resolves k in the (y, z) plane. Linear theory:
pypic.theory.streaming.filamentation_growth (cold relativistic fluid, k perpendicular to the drift).
Saturation by magnetic trapping (omega_B ~ Gamma, Davidson et al. 1972): B_sat ~ gamma0 m Gamma^2 / (q beta0 k).
A magnetic field along the drift (b0, in units of the field of marginal stability) lowers the growth
rate and removes the long wavelengths; from b0 = 1 on every k is stable.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

from pypic.config import Beam, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides
from pypic.theory.streaming import filamentation_species

MIXTURES = ("e-p", "e-i")


@dataclass
class FilamentationParams:
    mixture: str = "e-p"  # e-p: electrons and positrons both counter-stream; e-i: electrons only, ions at rest
    u0: float = 0.5  # beam proper velocity (gamma beta)
    theta: float = 1e-4  # beam temperature T / m c^2 (cold theory valid for k v_th << Gamma)
    mass_ratio: float = 1836.15267
    L_perp_wavelengths: float = 4.0  # L_y = L_z = this many 2 pi / k_ref, k_ref = 1 c^-1 omega_p,tot
    cells_per_skin_depth: float = 10.0
    nx: int = 4  # thin along the drift
    ppc: int = 16
    growth_times: float = 25.0
    n_modes: int = 16  # k_y, k_z mode numbers analysed
    b0: float = 0.0  # magnetic field along the drift over the field of marginal stability (Omega^2 = A)
    gate_m_min: int = 2  # the gate compares the growth rates of the modes gate_m_min..10
    seed: int = 31
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> FilamentationParams:
    return FilamentationParams(**d)


def species_for(p: FilamentationParams):
    th = p.theta
    if p.mixture == "e-p":
        e = [Beam(0.5, (p.u0, 0, 0), th), Beam(0.5, (-p.u0, 0, 0), th)]
        return [SpeciesConfig("electron", -1.0, 1.0, e),
                SpeciesConfig("positron", 1.0, 1.0, [Beam(0.5, (p.u0, 0, 0), th), Beam(0.5, (-p.u0, 0, 0), th)])]
    if p.mixture == "e-i":
        return [SpeciesConfig("electron", -1.0, 1.0, [Beam(0.5, (p.u0, 0, 0), th), Beam(0.5, (-p.u0, 0, 0), th)]),
                SpeciesConfig("ion", 1.0, p.mass_ratio, [Beam(1.0, (0, 0, 0), th / p.mass_ratio)])]
    raise ValueError(f"mixture must be one of {MIXTURES}")


def plan_filamentation(p: FilamentationParams, n_gpus: int = 1):
    species = species_for(p)
    wp_tot = math.sqrt(sum(b.density * s.charge ** 2 / s.mass for s in species for b in s.beams))
    k_ref = wp_tot  # c = 1
    Lp = p.L_perp_wavelengths * 2 * math.pi / k_ref
    dx = 1.0 / (p.cells_per_skin_depth * wp_tot)
    nperp = int(math.ceil(Lp / dx / 8)) * 8
    dx = Lp / nperp
    dt = 0.95 / math.sqrt(3.0) * dx
    ks = 2 * math.pi * np.arange(1, p.n_modes + 1) / Lp
    g0 = math.sqrt(1 + p.u0 ** 2)
    # field of marginal stability: Omega^2 = A = sum wp^2 beta^2 / gamma over the drifting beams (|q| = m = 1)
    w_crit = math.sqrt(sum(b.density * s.charge ** 2 / s.mass for s in species for b in s.beams if b.drift[0])
                       * (p.u0 / g0) ** 2 / g0)
    omega_c = p.b0 * w_crit
    B0 = omega_c * g0
    gam = filamentation_species(ks, species, omega_c)
    gam_free = filamentation_species(ks, species)
    g_max = float(gam.max())
    # stable at every k, not only at the analysed ones: Omega^2 >= A (filamentation_growth)
    stable = bool(omega_c >= w_crit)
    if not stable and g_max == 0.0:
        b = sum(bm.density * s.charge ** 2 / s.mass / (1 + bm.drift[0] ** 2) ** 1.5 for s in species for bm in s.beams)
        m_first = int(math.floor(math.sqrt(b * omega_c ** 2 / (w_crit ** 2 - omega_c ** 2)) / ks[0])) + 1
        raise ValueError(f"b0={p.b0:g} is below the field of marginal stability, but the unstable modes start at "
                         f"m={m_first}, above the n_modes={p.n_modes} analysed ones: raise n_modes or lower b0")
    # a stable case has no growth time of its own: it runs as long as the case without the field
    g_run = g_max if not stable else float(gam_free.max())
    t_end = p.growth_times / g_run
    n_steps = int(math.ceil(t_end / dt))
    n_cells = (p.nx, nperp, nperp)
    modes = [[0, m, 0] for m in range(1, p.n_modes + 1)] + [[0, 0, m] for m in range(1, p.n_modes + 1)]
    diags = [
        {"name": "modes", "kind": "modes", "fields": ["Bz", "By"], "modes": modes},
        {"name": "slice_yz", "kind": "slice", "fields": ["Bz", "Jx", "n:electron"], "axis": "x", "coarsen": [1, 2, 2],
         "every": {"count": 120}},
        {"name": "kspec", "kind": "kspectrum", "fields": ["Bx", "By", "Bz"], "every": {"count": 40}},
        {"name": "phase", "kind": "histogram", "species": ["electron"], "quantities": ["y", "ux"],
         "bins": [256, 128], "ranges": [None, [-3 * p.u0, 3 * p.u0]], "every": {"count": 100}},
        {"name": "spectrum", "kind": "histogram", "quantities": ["ekin"], "bins": [120], "log": [True],
         "ranges": [[1e-6, 10.0]], "every": {"count": 50}},
    ]
    try:
        dims = auto_dims(n_gpus, n_cells)
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=n_cells, dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc, seed=p.seed,
                    dims=dims, diag_interval=max(1, n_steps // 500), gauss_interval=max(1, n_steps // 50),
                    sort_interval=50, B0=(B0, 0.0, 0.0),
                    output=OutputConfig(diagnostics=diags, timescales={"plasma_period": 2 * math.pi / wp_tot,
                                                                       "growth_time": 1.0 / g_run, "t_end": t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    beta0 = p.u0 / g0
    i_max = int(np.argmax(gam))
    b_sat = g0 * g_max ** 2 / (beta0 * ks[i_max])
    info = {"mixture": p.mixture, "u0": p.u0, "theta": p.theta, "wp_tot": wp_tot, "L_perp": Lp,
            "k_modes": ks.tolist(), "gamma_theory": gam.tolist(), "gamma_max": g_max,
            "gamma_limit": float(filamentation_species([1e6], species, omega_c)[0]), "B_sat_trapping": b_sat,
            "B0": B0, "omega_c": omega_c, "omega_c_critical": w_crit, "stable": stable,
            "gate_modes": [p.gate_m_min, 10], "params": asdict(p)}
    lines = ["=== pypic filamentation planner ===",
             f"{p.mixture}: beams u0=+-{p.u0:g} (beta={beta0:.3f}, gamma={g0:.3f}), theta={p.theta:g}, "
             f"omega_p,tot={wp_tot:.3f}",
             f"grid {n_cells}, dx={dx:.4f} ({p.cells_per_skin_depth:g}/skin depth), L_perp={Lp:.2f}; dt={dt:.4f}; "
             f"t_end={t_end:.1f} ({p.growth_times:g}/Gamma_max), {n_steps} steps",
             "cold theory Gamma(k_perp): " + ", ".join(f"m={m + 1}: {g:.3f}" for m, g in enumerate(gam[:8]))
             + f" ... -> {info['gamma_limit']:.3f} (k -> inf)",
             *([f"B0={B0:.4f} along the drift: Omega={omega_c:.4f} = {p.b0:g} of the marginal {w_crit:.4f}"
                + (f"; stable at every k, run length from the case without the field ({g_run:.3f})"
                   if stable else "")] if B0 else []),
             f"magnetic-trapping saturation estimate for the fastest mode: B ~ {b_sat:.3f} m c omega_pe / e"]
    return cfg, "\n".join(lines), info
