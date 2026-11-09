"""Automatic choice of grid, time step, box size and run length for two-stream experiments.

The planner solves the linear dispersion relation for the requested plasma, places an integer
number of fastest-growing wavelengths in the box, picks the cell size / time step from the
stability and accuracy constraints of the selected solver, and reports every constraint margin.
Every chosen value can be overridden by the user (``overrides``).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

from ..config import ImplicitConfig, SimConfig
from ..parallel.decomp import auto_dims
from .dispersion import fastest_mode, filamentation_estimate
from .species import Composition, TwoStreamBeams, build_species

C_SI = 299792458.0
QE = 1.602176634e-19
ME = 9.1093837015e-31
EPS0 = 8.8541878128e-12
MEC2_EV = 510998.95


@dataclass
class TwoStreamParams:
    # --- plasma ---
    mixture: str = "e-i"  # e-i | e-p | e-i-p
    mass_ratio: float = 1836.15267
    Z: int = 1
    positron_fraction: float = 0.5
    mode: str = "electron"  # electron | symmetric | buneman
    u_drift: float = 0.2  # beam proper velocity gamma*v/c
    beam_fraction: float = 0.5
    T_e: float = 1e-4  # m_e c^2 units (use T_e_eV to give eV)
    T_i: float | None = None
    T_p: float | None = None
    T_e_eV: float | None = None
    T_i_eV: float | None = None
    T_p_eV: float | None = None
    n0_m3: float | None = None  # optional physical density for unit conversion in the report
    B0: tuple = (0.0, 0.0, 0.0)
    # --- numerics ---
    solver: str = "explicit"  # explicit | implicit
    n_wavelengths: int = 4  # fastest-growing wavelengths along x
    cells_per_wavelength: int | None = None  # default: 32 explicit / 16 implicit
    transverse_ratio: float = 0.25  # L_y = L_z = ratio * L_x
    ppc: int = 16  # particle positions per cell (each beam gets one particle per position)
    cfl: float | None = None  # c dt / dx (explicit: limited by 1/sqrt(3))
    wpe_dt_max: float | None = None  # default 0.2 explicit / 0.5 implicit
    growth_times: float = 30.0  # run length in units of 1/gamma_max
    t_end: float | None = None
    max_cells: int = 512 ** 3
    boundaries: tuple = (("periodic", "periodic"),) * 3
    n_records: int = 400  # number of diagnostic records over the run
    seed: int = 12345
    warm: bool = True
    overrides: dict = field(default_factory=dict)  # any SimConfig field

    def temperatures(self):
        Te = self.T_e if self.T_e_eV is None else self.T_e_eV / MEC2_EV
        Ti = self.T_i if self.T_i_eV is None else self.T_i_eV / MEC2_EV
        Tp = self.T_p if self.T_p_eV is None else self.T_p_eV / MEC2_EV
        return Te, Ti, Tp


def _round_cells(n: float, mult: int = 8) -> int:
    return max(mult, int(math.ceil(n / mult)) * mult)


def plan_two_stream(p: TwoStreamParams, n_gpus: int = 1) -> tuple[SimConfig, str, dict]:
    Te, Ti, Tp = p.temperatures()
    comp = Composition(p.mixture, p.mass_ratio, p.Z, p.positron_fraction)
    beams = TwoStreamBeams(p.mode, p.u_drift, p.beam_fraction, Te, Ti, Tp)
    species = build_species(comp, beams)
    lin = fastest_mode(species, warm=p.warm)
    lam = lin.wavelength
    implicit = p.solver == "implicit"
    cpw = p.cells_per_wavelength or (16 if implicit else 32)

    # ---- grid ----
    Lx = p.n_wavelengths * lam
    dx = lam / cpw
    nx = p.n_wavelengths * cpw
    ny = nz = _round_cells(p.transverse_ratio * nx)
    nx = _round_cells(nx)
    dx = Lx / nx
    # keep cubic cells; shrink the transverse size if the cell budget is exceeded
    while nx * ny * nz > p.max_cells and ny > 8:
        ny = nz = _round_cells(ny // 2)
    n_cells = (nx, ny, nz)
    lam_De = math.sqrt(max(Te, 1e-300))  # lambda_De = v_te/omega_pe = sqrt(T_e/m_e)/omega_pe

    # ---- time step ----
    cfl_max = 1.0 / math.sqrt(3.0)
    if implicit:
        cfl = p.cfl if p.cfl is not None else 1.0
        wdt = p.wpe_dt_max if p.wpe_dt_max is not None else 0.5
    else:
        cfl = p.cfl if p.cfl is not None else 0.95 * cfl_max
        if cfl >= cfl_max:
            raise ValueError(f"explicit scheme requires c dt/dx < 1/sqrt(3) = {cfl_max:.4f}")
        wdt = p.wpe_dt_max if p.wpe_dt_max is not None else 0.2
    B = math.sqrt(sum(b * b for b in p.B0))
    wce = B  # electron cyclotron frequency (normalised, m_e = e = 1)
    # highest plasma frequency of any species (electrons dominate): omega_pe = 1
    limits = {"CFL c dt/dx": cfl * dx, "omega_pe dt": wdt / 1.0}
    if wce > 0:
        limits["omega_ce dt"] = 0.2 / wce
    limits["growth resolution (gamma dt <= 0.05)"] = 0.05 / lin.gamma_max
    limits["oscillation resolution (|omega_r| dt <= 0.2)"] = 0.2 / max(abs(lin.omega_r), 1e-30)
    dt = min(limits.values())
    t_end = p.t_end if p.t_end is not None else p.growth_times / lin.gamma_max
    n_steps = int(math.ceil(t_end / dt))
    diag_interval = max(1, n_steps // p.n_records)

    dims = None
    try:
        dims = auto_dims(n_gpus, n_cells)
    except ValueError:
        pass

    cfg = SimConfig(
        n_cells=n_cells, dx=(dx, dx, dx), dt=dt, n_steps=n_steps, solver=p.solver,
        boundaries=tuple(tuple(b) for b in p.boundaries), species=species, ppc=p.ppc, seed=p.seed,
        dims=dims, B0=tuple(p.B0), diag_interval=diag_interval, gauss_interval=diag_interval * 10,
        sort_interval=50, mode_diag=min(64, nx // 2), implicit=ImplicitConfig(),
        phase_space_interval=max(1, n_steps // 6))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)

    # ---- report ----
    nbeams = sum(len(s.beams) for s in species)
    ncell = cfg.n_cells[0] * cfg.n_cells[1] * cfg.n_cells[2]
    npart = ncell * cfg.ppc * nbeams
    bytes_p = npart * 8 * (7 + (6 if implicit else 0)) * 1.3
    bytes_f = ncell * 8 * (9 + (12 if not implicit else 30))
    gpus = max(n_gpus, 1)
    info = {
        "k_max": lin.k_max, "gamma_max": lin.gamma_max, "omega_r": lin.omega_r, "wavelength": lam,
        "lambda_De": lam_De, "dx_over_lambdaDe": cfg.dx[0] / lam_De if lam_De > 0 else float("inf"),
        "dt": cfg.dt, "t_end": cfg.dt * cfg.n_steps, "particles": npart,
        "memory_per_gpu_GB": (bytes_p + bytes_f) / gpus / 1e9, "dims": cfg.dims,
        "filamentation_gamma_est": filamentation_estimate(species),
    }
    lines = [
        "=== pypic two-stream planner ===",
        f"mixture={p.mixture} mode={p.mode} m_i/m_e={p.mass_ratio:g} Z={p.Z} u_drift={p.u_drift} "
        f"T_e={Te:.3e} m_ec^2 ({Te * MEC2_EV:.3g} eV) solver={p.solver}",
        "species / beams (n in n0, u in c, theta = T/mc^2):",
    ]
    for s in species:
        for b in s.beams:
            lines.append(f"  {s.name:9s} q={s.charge:+.0f} m={s.mass:<9.4g} n={b.density:.4f} "
                         f"u_x={b.drift[0]:+.4f} theta={b.theta:.3e}")
    lines += [
        f"linear theory ({'warm' if p.warm else 'cold'} fluid): k_max={lin.k_max:.4f} omega_pe/c, "
        f"gamma_max={lin.gamma_max:.4f} omega_pe, omega_r={lin.omega_r:+.4f}, lambda={lam:.4f} c/omega_pe",
        f"filamentation (transverse) growth estimate: {info['filamentation_gamma_est']:.4f} omega_pe",
        f"grid: {cfg.n_cells} cells, dx={cfg.dx[0]:.5f} c/omega_pe ({lam / cfg.dx[0]:.1f} cells/wavelength), "
        f"L=({cfg.n_cells[0] * cfg.dx[0]:.3f}, {cfg.n_cells[1] * cfg.dx[1]:.3f}, {cfg.n_cells[2] * cfg.dx[2]:.3f})",
        f"dx/lambda_De = {info['dx_over_lambdaDe']:.3g}",
        f"dt={cfg.dt:.5f} /omega_pe  (c dt/dx={cfg.dt / cfg.dx[0]:.3f}, omega_pe dt={cfg.dt:.3f}, "
        f"gamma_max dt={lin.gamma_max * cfg.dt:.4f})",
        "  time-step limits: " + ", ".join(f"{k}: {v:.4g}" for k, v in limits.items()),
        f"run: t_end={cfg.dt * cfg.n_steps:.2f} /omega_pe = {cfg.dt * cfg.n_steps * lin.gamma_max:.1f}/gamma_max, "
        f"{cfg.n_steps} steps, diag every {cfg.diag_interval}",
        f"particles: {npart:,} ({cfg.ppc} positions/cell x {nbeams} beams), "
        f"memory ~{info['memory_per_gpu_GB']:.2f} GB/GPU on {gpus} GPU(s), rank grid {cfg.dims}",
    ]
    if not implicit and lam_De > 0 and cfg.dx[0] > math.pi * lam_De:
        lines.append(f"WARNING: dx = {cfg.dx[0] / lam_De:.1f} lambda_De > pi lambda_De: the explicit "
                     "momentum-conserving scheme may show finite-grid heating; consider solver='implicit' "
                     "(energy conserving) or more cells per wavelength.")
    if p.mode == "buneman":
        lines.append("NOTE: buneman mode carries a net current; the k=0 mode oscillates at omega_pe.")
    return cfg, "\n".join(lines), info


def apply_overrides(cfg: SimConfig, overrides: dict) -> SimConfig:
    """Return a new SimConfig with ``overrides`` applied. Nested values (``implicit`` as a dict,
    ``species`` as a list of dicts, tuples given as lists) are converted through
    :meth:`SimConfig.from_dict`; nested dicts are merged into the current values."""
    d = asdict(cfg)
    for k, v in overrides.items():
        if k not in d:
            raise KeyError(f"unknown SimConfig field {k!r}")
        if isinstance(d[k], dict) and isinstance(v, dict):
            unknown = set(v) - set(d[k])
            if unknown:
                raise KeyError(f"unknown {k} field(s) {sorted(unknown)}")
            d[k] = {**d[k], **v}
        else:
            d[k] = v
    return SimConfig.from_dict(d)


def params_from_dict(d: dict) -> TwoStreamParams:
    d = dict(d)
    if "B0" in d:
        d["B0"] = tuple(d["B0"])
    if "boundaries" in d:
        d["boundaries"] = tuple(tuple(b) for b in d["boundaries"])
    return TwoStreamParams(**d)


def params_to_dict(p: TwoStreamParams) -> dict:
    return asdict(p)
