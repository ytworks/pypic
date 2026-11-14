"""Planner of the topic's runs: a perpendicular piston shock moving into a slab of Alfvénic modes.

    python -m pypic.launch -n 8 research/topics/shock-alfvenic-turbulence/run.py \
        --config research/topics/shock-alfvenic-turbulence/configs/<name>.yaml --out <dir>

Geometry (README, Approach): x is the shock normal, y lies along the front and is periodic, B0 = B0 y
(theta_Bn = 90 deg). Electrons and ions flow at -V0 along x (motional field E0 = V0 B0 z).

- ``wall: true``: a reflecting wall at x- stops the flow; the perpendicular shock moves into the
  plasma. The face x+ injects the flow (``flux`` face, PIC-D13) and carries the modes (PIC-D26).
  The modes are loaded at t = 0 under an envelope h(x) that vanishes at the wall and the face
  (PIC-D24).
- ``wall: false``: x is periodic, no shock: the same modes in the drifting plasma (null test N2 and
  the grid heating check of the pilot).

Modes (design/setup.md). k_n = 2 pi n / L_y with k_n d_i <= kdi_max. Each mode is one circularly
polarized cold two-fluid wave that travels along +y (examples/alfven_inflow/theory.py): its
polarization does not rotate (a linear one would, design/estimates.md Sect. 1). dB_x and dB_z of a mode
have the same amplitude, so the rms of |dB| is sqrt(2) db. With k_x = 0 the upstream state depends on
(y, t) only. ``state`` chooses the amplitudes and phases:

- ``constant_b``: dB_x + i s dB_z = sqrt(2) db B0 f(e^{i k_1 y}) (s = +1 for ``helicity: R``, -1 for
  ``L``), where f is a finite Blaschke product, f(z) = z prod_j (z - a_j) / (1 - conj(a_j) z) with
  |a_j| = ``zero_radius``: |f| = 1 on the circle and f holds positive wavenumbers only, so |B| is
  uniform and every mode has the same helicity and direction. ``zeros: 0`` is one circular wave, an
  exact solution of the cold two-fluid equations with uniform |B| and density. With zeros, the
  dispersion of the modes changes their relative phases and |B| loses its uniformity
  (``b2_spread`` of the report: std / mean of |B|^2 over y from the linear phases, at t = 0 and at
  t_end). The angles of the zeros come from ``realization``.
- ``random``: weights of <dB_x^2> from the slab spectrum of design/estimates.py
  (P ~ 1 / (1 + (k L_c)^(5/3))), random phases, and ``helicity: random`` draws R or L per mode
  (``realization``). |B| is not uniform (a comparison state).

The plasma: beta_e and beta_i (temperatures from B0 = 1 / wpe_wce), mass_ratio, and the flow speed V0
chosen so that the perpendicular MHD shock of the piston has the Alfvén Mach number ``mach_alfven``
in its own frame (theory.piston_shock). Lengths of the parameters in d_i, times in 1 / Omega_ci.
"""
from __future__ import annotations

import importlib.util
import math
import os
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.optimize import brentq

from pypic.config import Beam, CurrentFilter, Injection, OutputConfig, SimConfig, SpeciesConfig
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(ROOT, "..", "..", ".."))


def _theory():
    """examples/alfven_inflow/theory.py (main-owned, read only), under a name of its own."""
    spec = importlib.util.spec_from_file_location("sat_alfven_theory",
                                                  os.path.join(REPO, "examples", "alfven_inflow", "theory.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T = _theory()
HELICITIES = ("random", "R", "L")
STATES = ("constant_b", "random")
L_DAMPED_KDI = 0.3  # above it an L mode damps within 40 / Omega_ci at beta_i = 1 (design/setup.md)


@dataclass
class TopicParams:
    mass_ratio: float = 100.0
    wpe_wce: float = 10.0  # omega_pe / Omega_ce; B0 = 1 / wpe_wce along y
    beta_e: float = 1.0
    beta_i: float = 1.0
    mach_alfven: float = 4.5  # of the piston shock in its own frame (perpendicular MHD jump conditions)
    db: float = 0.6  # rms dB_x / B0
    ly: float = 50.0  # d_i
    state: str = "constant_b"
    zeros: int = 0  # constant_b: zeros of the Blaschke product (0: one circular wave)
    zero_radius: float = 0.3  # constant_b: |a_j|; the amplitudes fall as zero_radius^n
    tail: float = 1e-4  # constant_b: share of <dB^2> that the dropped modes may hold
    lc: float | None = None  # random: correlation length of the spectrum (d_i); default ly / 4
    kdi_max: float = 0.5  # largest k d_i of the modes
    helicity: str = "R"
    realization: int = 0
    wall: bool = True
    lx: float = 200.0  # d_i
    wall_ramp: tuple = (0.4, 1.2)  # envelope rises from 0 to 1 over this x range (d_i)
    face_ramp: float = 0.8  # envelope falls to 0 over this distance before the face (d_i)
    nz: int = 2
    cells_per_de: float = 5.0
    cfl: float = 0.5
    ppc: int = 5  # per species and cell: 2 nz ppc particles per column of the plane
    t_end: float = 40.0  # 1 / Omega_ci
    filter_passes: int = 4
    marder_zone: float = 0.125
    npml: int = 12
    seed: int = 12345
    frames: int = 40  # slices of the plane z = 0
    profiles: int = 200  # lineouts along x (averaged over y)
    spectra: int = 40  # histograms (x, ekin)
    track_stride: int = 0  # electrons tracked: every n-th (0: none)
    track_every: int = 20  # steps between track samples
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> TopicParams:
    d = dict(d)
    if "wall_ramp" in d:
        d["wall_ramp"] = tuple(d["wall_ramp"])
    return TopicParams(**d)


def spectrum(ly: float, lc: float, kdi_max: float, db: float):
    """Mode numbers n, k_n d_i and the dB_x amplitudes a_n (sum a_n^2 / 2 = db^2) of the slab."""
    n_max = int(math.floor(kdi_max * ly / (2.0 * math.pi) + 1e-9))
    if n_max < 1:
        raise ValueError(f"no mode with k d_i <= {kdi_max:g} in L_y = {ly:g} d_i")
    n = np.arange(1, n_max + 1)
    k = 2.0 * math.pi * n / ly
    p = 1.0 / (1.0 + (k * lc) ** (5.0 / 3.0))
    return n, k, np.sqrt(2.0 * p / p.sum()) * db


def constant_b(ly: float, kdi_max: float, db: float, angles, radius: float, tail: float, nfft: int = 8192):
    """Mode numbers, k_n d_i, dB_x amplitudes (sum a_n^2 / 2 = db^2) and phases of the Blaschke state.

    f(z) = z prod_j (z - a_j) / (1 - conj(a_j) z), a_j = radius e^{i angles[j]}; its Fourier
    coefficients c_n (n >= 1) are kept up to kdi_max and until the rest holds less than ``tail``."""
    z = np.exp(2j * math.pi * np.arange(nfft) / nfft)
    f = z.copy()
    for a in radius * np.exp(1j * np.asarray(angles, dtype=float)):
        f *= (z - a) / (1.0 - np.conj(a) * z)
    c = np.fft.fft(f)[1:nfft // 2] / nfft
    n_max = int(math.floor(kdi_max * ly / (2.0 * math.pi) + 1e-9))
    if n_max < 1:
        raise ValueError(f"no mode with k d_i <= {kdi_max:g} in L_y = {ly:g} d_i")
    energy = np.cumsum(np.abs(c) ** 2)
    n_keep = min(n_max, int(np.searchsorted(energy, 1.0 - tail) + 1))
    c = c[:n_keep]
    kept = float(energy[n_keep - 1])
    n = np.arange(1, n_keep + 1)
    return n, 2.0 * math.pi * n / ly, np.abs(c) / math.sqrt(kept) * math.sqrt(2.0) * db, np.angle(c), kept


def b2_spread(k, a, phases, signs, ny: int = 4096) -> float:
    """std / mean over y of |B|^2 / B0^2 = 1 + dB_x^2 + dB_z^2 (dB_z of a mode: signs * a sin)."""
    y = np.linspace(0.0, 2.0 * math.pi / k[0], ny, endpoint=False)
    th = k[:, None] * y[None, :] + phases[:, None]
    b2 = 1.0 + np.sum(a[:, None] * np.cos(th), axis=0) ** 2 + np.sum((signs * a)[:, None] * np.sin(th), axis=0) ** 2
    return float(b2.std() / b2.mean())


def excursion(k, a, phases, ny: int = 4096):
    """Peak-to-peak field-line excursion xi_x(y) = int dB_x / B0 dy (d_i) of one realization at t = 0."""
    y = np.linspace(0.0, 2.0 * math.pi / k[0], ny, endpoint=False)
    xi = np.sum((a / k)[:, None] * np.sin(k[:, None] * y[None, :] + phases[:, None]), axis=0)
    return float(xi.max() - xi.min())


def flow_speed(mach: float, v_a: float, beta: float, mass: float) -> float:
    """V0 (box frame) of the piston whose shock has the Alfvén Mach number ``mach`` in its frame."""
    return brentq(lambda v: T.piston_shock(v, v_a, beta, mass=mass)["mach_alfven"] - mach,
                  0.05 * v_a, mach * v_a, xtol=1e-14)


def preload(table: dict, ly: float, dy: float, numbers, h: dict | None = None) -> tuple:
    """theory.preload for several modes: mode i of the table has k = 2 pi numbers[i] / ly."""
    fields, pert = [], {name: [] for name in table["modes"][0]["u"]}
    for m, n in zip(table["modes"], numbers):
        f, p = T.preload({"face": table["face"], "modes": [m]}, ly / n, dy, h)
        for spec in f:  # a profile or potential shape along y: cos of mode 1 over ly / n = mode n over ly
            sh = spec["shape"]["y"] if spec["kind"] == "potential" else spec
            sh["mode"], sh["length"] = int(n), ly
        for name, lst in p.items():
            for spec in lst:
                spec["shape"]["y"]["mode"], spec["shape"]["y"]["length"] = int(n), ly
            pert[name] += lst
        fields += f
    return fields, pert


def plan_topic(p: TopicParams, n_gpus: int = 1):
    if p.helicity not in HELICITIES:
        raise ValueError(f"helicity must be one of {HELICITIES}")
    if p.state not in STATES:
        raise ValueError(f"state must be one of {STATES}")
    if p.state == "constant_b" and p.helicity == "random":
        raise ValueError("state constant_b needs one helicity (R or L): both together do not keep |B| uniform "
                         "in a dispersive plasma")
    b = 1.0 / p.wpe_wce
    species = T.two_fluid(p.mass_ratio)
    d_i = T.ion_skin_depth(species)
    v_a = T.alfven_speed(b, species)
    omega_ci = b / p.mass_ratio
    theta_e = 0.5 * p.beta_e * b * b
    theta_i = 0.5 * p.beta_i * b * b / p.mass_ratio  # T_i / (m_i c^2)
    beta = p.beta_e + p.beta_i
    mass = p.mass_ratio + 1.0
    v0 = flow_speed(p.mach_alfven, v_a, beta, mass)
    shock = T.piston_shock(v0, v_a, beta, mass=mass)
    dx = 1.0 / p.cells_per_de
    dt = p.cfl * dx
    lx, ly = p.lx * d_i, p.ly * d_i
    nx, ny = int(round(lx / dx)), int(round(ly / dx))
    if abs(nx * dx - lx) > 1e-9 or abs(ny * dx - ly) > 1e-9:
        raise ValueError("lx and ly must be whole numbers of cells")
    rng = np.random.default_rng(p.realization)
    lc, kept = None, 1.0
    if p.state == "constant_b":
        angles = rng.uniform(0.0, 2.0 * math.pi, size=p.zeros)
        numbers, kdi, amp, phases, kept = constant_b(p.ly, p.kdi_max, p.db, angles, p.zero_radius, p.tail)
        branches = [p.helicity] * numbers.size
    else:
        lc = p.lc if p.lc is not None else 0.25 * p.ly
        numbers, kdi, amp = spectrum(p.ly, lc, p.kdi_max, p.db)
        phases = rng.uniform(0.0, 2.0 * math.pi, size=numbers.size)
        branches = (rng.choice(["R", "L"], size=numbers.size) if p.helicity == "random"
                    else [p.helicity] * numbers.size)
    table = {"face": "x+", "modes": []}
    waves = []
    for kd, a, ph, br in zip(kdi, amp, phases, branches):
        w = T.parallel_wave(kd / d_i, b, species, str(br), math.sqrt(2.0) * a)  # |B_x| = |B_z| = a b
        waves.append({"k_di": float(kd), "branch": str(br), "omega_over_omega_ci": w["omega"] / omega_ci,
                      "dbx": float(a), "phase": float(ph)})
        table["modes"] += T.fluctuation(w, v0, phase=float(ph))["modes"]
    # the state at t_end from the linear phases (plasma frame): dispersion changes |B|^2 and the field lines
    signs = np.array([1.0 if br == "R" else -1.0 for br in branches])
    late = phases - np.array([w["omega_over_omega_ci"] for w in waves]) * p.t_end
    spread = (b2_spread(kdi, amp, phases, signs), b2_spread(kdi, amp, late, signs)) if p.db > 0.0 else (0.0, 0.0)
    xi_pp, xi_end = excursion(kdi, amp, phases), excursion(kdi, amp, late)
    notes = []
    damped = [w["k_di"] for w in waves if w["branch"] == "L" and w["k_di"] > L_DAMPED_KDI]
    if damped and p.db > 0.0:
        notes.append(f"{len(damped)} L mode(s) above k d_i = {L_DAMPED_KDI:g} (largest {max(damped):.3f}): "
                     "ion cyclotron damping within the run (design/setup.md)")
    h = None
    if p.wall:
        w0, w1 = (r * d_i for r in p.wall_ramp)
        h = T.envelope(lx, (w0, w1), p.face_ramp * d_i)
    fields, pert = preload(table, ly, dx, numbers, h) if p.db > 0.0 else ([], {"e": [], "i": []})
    u0 = v0 / math.sqrt(1.0 - v0 * v0)
    sp = [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, (-u0, 0.0, 0.0), theta_e, perturbation=pert["e"])]),
          SpeciesConfig("i", 1.0, p.mass_ratio, [Beam(1.0, (-u0, 0.0, 0.0), theta_i, perturbation=pert["i"])])]
    n_steps = int(math.ceil(p.t_end / omega_ci / dt))
    plane = ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "n:e", "n:i", "flux:e", "flux:i", "ekin:e", "ekin:i"]
    diags = [
        {"name": "plane", "kind": "slice", "axis": "z", "position": 0.0, "fields": plane,
         "every": {"count": p.frames}},
        {"name": "profile", "kind": "lineout", "axis": "x", "average": False,
         "fields": ["n:e", "n:i", "Bx", "By", "Bz", "flux:e", "flux:i", "ekin:e", "ekin:i"],
         "every": {"count": p.profiles}},
        {"name": "spectrum", "kind": "histogram", "quantities": ["x", "ekin"],
         "bins": [max(1, nx // 50), 160], "ranges": [[0.0, lx], [1e-5, 10.0]], "log": [False, True],
         "every": {"count": p.spectra}},
    ]
    if p.track_stride > 0:
        diags.append({"name": "tracks", "kind": "tracks", "species": ["e"], "stride": p.track_stride,
                      "every": p.track_every, "with_fields": True})
    if p.wall:
        bc = (("reflecting", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
        injection = Injection(faces={"x+": "flux"}, fluctuation=table if p.db > 0.0 else
                              T.fluctuation(T.parallel_wave(kdi[0] / d_i, b, species, "R", 0.0), v0))
        capacity = 1.2 * shock["compression"]
    else:
        bc = (("periodic", "periodic"),) * 3
        injection = Injection()  # no face injects
        capacity = 1.3
    try:
        dims = auto_dims(n_gpus, (nx, ny, p.nz))
    except ValueError:
        dims = None
    cfg = SimConfig(
        n_cells=(nx, ny, p.nz), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=sp, ppc=p.ppc, seed=p.seed, dims=dims,
        boundaries=bc, npml=p.npml, current_filter=CurrentFilter(p.filter_passes, True),
        marder_zone=p.marder_zone if p.wall else 0.0,
        B0=(0.0, b, 0.0), E0=(0.0, 0.0, v0 * b), initial_fields=fields, injection=injection,
        capacity_factor=capacity, diag_interval=max(1, n_steps // 1000), gauss_interval=max(1, n_steps // 100),
        sort_interval=50,
        output=OutputConfig(diagnostics=diags, timescales={"ion_gyration": 2.0 * math.pi / omega_ci}))
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)
    l_sh = 0.68 * shock["mach_alfven"]  # d_i, design/estimates.py
    # when the plasma that entered through the face (PIC-D26, not preloaded) reaches the shock, and when a
    # fast-mode signal from the face does (speeds of the uniform state, B0 only); 1 / Omega_ci
    u1 = v0 + shock["v_shock"]
    t_face = {"flow": p.lx * v_a / u1, "fast": p.lx * v_a / (u1 + u1 / shock["mach_fast"])}
    if p.wall and p.t_end > t_face["flow"]:
        notes.append(f"plasma that entered through the face reaches the shock at {t_face['flow']:.1f} / Omega_ci, "
                     f"before t_end = {p.t_end:g}")
    b_total = math.sqrt(1.0 + float(np.sum(amp ** 2)) if p.db > 0.0 else 1.0)  # rms |B| / B0
    info = {
        "b0": b, "v_a": v_a, "v0": v0, "u0": u0, "d_i": d_i, "omega_ci": omega_ci, "theta_e": theta_e,
        "theta_i": theta_i, "beta": beta, "shock": shock, "lc": lc, "modes": waves, "xi_pp": xi_pp,
        "Xi": xi_pp / l_sh, "xi_pp_end": xi_end, "l_sh": l_sh, "b2_spread": {"t0": spread[0], "t_end": spread[1]},
        "kept": kept, "b_total": b_total, "t_face": t_face, "notes": notes,
        "debye_cells": math.sqrt(theta_e) / dx, "wall": p.wall,
        "lx": lx, "ly": ly, "envelope": h, "omega_ci_t_end": omega_ci * n_steps * dt,
        "particles": 2 * nx * ny * p.nz * p.ppc, "dims": cfg.dims, "params": asdict(p)}
    lines = [
        "=== shock-alfvenic-turbulence planner ===",
        f"plasma: m_i/m_e={p.mass_ratio:g}, omega_pe/Omega_ce={p.wpe_wce:g}, beta_e={p.beta_e:g}, beta_i={p.beta_i:g}, "
        f"v_A={v_a:.5f} c, d_i={d_i:g} c/omega_pe, Debye length {info['debye_cells']:.2f} dx",
        f"flow V0={v0:.5f} c = {v0 / v_a:.3f} v_A; piston shock: r={shock['compression']:.3f}, "
        f"V_sh={shock['v_shock'] / v_a:.3f} v_A, M_A={shock['mach_alfven']:.2f} (fast {shock['mach_fast']:.2f})"
        if p.wall else f"flow V0={v0:.5f} c = {v0 / v_a:.3f} v_A; periodic box, no shock",
        f"state {p.state}" + (f" (zeros {p.zeros}, radius {p.zero_radius:g}, kept {kept:.5f})"
                              if p.state == "constant_b" else f" (L_c={lc:g} d_i)")
        + f": {numbers.size} mode(s) (k d_i {kdi[0]:.3f}..{kdi[-1]:.3f}), rms dB_x/B0={p.db:g}, "
        f"helicity {p.helicity}, realization {p.realization}",
        f"  field lines: xi_pp={xi_pp:.2f} d_i, Xi={info['Xi']:.2f} (at t_end {xi_end:.2f} d_i); |B|^2 std/mean "
        f"{spread[0]:.4f} at t=0, {spread[1]:.4f} at t_end (linear phases); rms |B|/B0={b_total:.3f}, "
        f"M_A on it {shock['mach_alfven'] / b_total:.2f}",
        f"grid {cfg.n_cells} ({p.lx:g} x {p.ly:g} d_i), dx={dx:g}, dt={dt:g}, {n_steps} steps to "
        f"{info['omega_ci_t_end']:.1f} / Omega_ci; {info['particles']:,} particles at t = 0",
    ]
    if p.wall:
        lines.append(f"face: its plasma reaches the shock at {t_face['flow']:.1f} / Omega_ci, a fast-mode signal "
                     f"from it at {t_face['fast']:.1f} / Omega_ci")
    lines += [f"NOTE: {n}" for n in notes]
    return cfg, "\n".join(lines), info
