"""Planner for the isothermal expansion of a plasma slab into vacuum (Mora 2003; TNSA ion acceleration).

A slab of hot electrons (T_e) and cold ions fills 0 <= x < d (reflecting wall at x = 0 = symmetry
plane); vacuum fills d <= x < L (reflecting far wall). The expansion depends on
tau = omega_pi t / sqrt(2 e) only (lengths in lambda_D0), so a reduced mass ratio reaches the same
tau in fewer steps. The slab is thick compared with the rarefaction (c_s t) so the electrons stay
close to isothermal (Mora's assumption); the thin_foil case shows the adiabatic cooling instead.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import asdict, dataclass, field

from pypic.config import Beam, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory as T  # noqa: E402


@dataclass
class ExpansionParams:
    Te: float = 0.01  # T_e / m_e c^2 (lambda_D0 = 0.1 c/omega_pe)
    Ti_over_Te: float = 0.0
    mass_ratio: float = 100.0
    Z: float = 1.0
    tau_end: float = 10.0
    slab_lambda_D: float = 200.0  # slab thickness d in lambda_D0
    vacuum_front_margin: float = 1.6  # vacuum length = margin x Mora front position at tau_end
    dx_over_lambda_D: float = 1.0
    n_transverse: int = 4
    ppc: int = 256
    seed: int = 99
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> ExpansionParams:
    return ExpansionParams(**d)


def plan_expansion(p: ExpansionParams, n_gpus: int = 1):
    s = T.scales(p.Te, p.mass_ratio, p.Z)
    lam = s["lambda_D0"]
    t_end = p.tau_end * math.sqrt(2 * T.EULER) / s["omega_pi"]
    x_front_end = float(T.front_position(t_end, p.Te, p.mass_ratio, p.Z))
    d = p.slab_lambda_D * lam
    L_vac = max(p.vacuum_front_margin * x_front_end, 20 * lam)
    dx = p.dx_over_lambda_D * lam
    n_ranks = max(n_gpus, 1)
    nx = int(math.ceil((d + L_vac) / dx / (8 * n_ranks))) * 8 * n_ranks
    L = nx * dx
    dt = 0.95 / math.sqrt(3.0) * dx
    n_steps = int(math.ceil(t_end / dt))
    slab = {"axis": "x", "shape": "step", "x0": 0.0, "x1": d}
    ni = 1.0 / p.Z
    species = [SpeciesConfig("electron", -1.0, 1.0, [Beam(1.0, (0, 0, 0), p.Te, profile=slab)]),
               SpeciesConfig("ion", p.Z, p.mass_ratio,
                             [Beam(ni, (0, 0, 0), p.Te * p.Ti_over_Te / p.mass_ratio, profile=slab)])]
    nt = p.n_transverse
    A_perp = (nt * dx) ** 2
    e_max_end = float(T.max_energy(t_end, p.Te, p.mass_ratio, p.Z))
    diags = [
        {"name": "lineout", "kind": "lineout", "fields": ["Ex", "n:ion", "n:electron"], "axis": "x",
         "every": {"count": 200}},
        {"name": "ion_phase", "kind": "histogram", "species": ["ion"], "quantities": ["x", "ux"], "bins": [512, 256],
         "ranges": [[0.0, L], [-0.1 * math.sqrt(2 * 3 * e_max_end / p.mass_ratio),
                                math.sqrt(2 * 3 * e_max_end / p.mass_ratio)]], "every": {"count": 100}},
        {"name": "ion_spectrum", "kind": "histogram", "species": ["ion"], "quantities": ["ekin"], "bins": [600],
         "ranges": [[1e-4 * p.Te, 4 * e_max_end]], "log": [True], "every": {"count": 100}},
        {"name": "electron_spectrum", "kind": "histogram", "species": ["electron"], "quantities": ["ekin"],
         "bins": [160], "ranges": [[1e-4 * p.Te, 100 * p.Te]], "log": [True], "every": {"count": 50}},
    ]
    try:
        dims = auto_dims(n_gpus, (nx, nt, nt))
    except ValueError:
        dims = None
    cfg = SimConfig(n_cells=(nx, nt, nt), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=species, ppc=p.ppc,
                    seed=p.seed, dims=dims,
                    boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    diag_interval=max(1, n_steps // 500), gauss_interval=max(1, n_steps // 50), sort_interval=50,
                    output=OutputConfig(diagnostics=diags, timescales={"plasma_period": 2 * math.pi,
                                                                       "ion_period": 2 * math.pi / s["omega_pi"],
                                                                       "t_end": t_end}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    info = {"Te": p.Te, "mass_ratio": p.mass_ratio, "Z": p.Z, "edge": d, "A_perp": A_perp, "t_end": t_end,
            "x_front_end": x_front_end, "E_max_end": e_max_end, **s, "params": asdict(p)}
    lines = ["=== pypic plasma-expansion planner (Mora 2003) ===",
             f"T_e={p.Te:g} m_ec^2 ({p.Te * 511:.3g} keV), m_i/m_e={p.mass_ratio:g}, Z={p.Z:g}: lambda_D0={lam:.3f}, "
             f"c_s={s['cs']:.4g} c, omega_pi={s['omega_pi']:.4g}",
             f"slab d={d:.2f} ({p.slab_lambda_D:g} lambda_D0), vacuum {L - d:.2f}; grid {cfg.n_cells}, dx={dx:.3f}, "
             f"dt={dt:.4f}",
             f"run to tau={p.tau_end:g}: t_end={t_end:.1f} ({n_steps} steps); Mora: front x_f={x_front_end:.2f}, "
             f"E_max={e_max_end:.4g} m_ec^2 ({e_max_end * 511:.3g} keV)"]
    return cfg, "\n".join(lines), info
