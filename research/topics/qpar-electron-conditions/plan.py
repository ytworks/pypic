"""Planner of the topic's runs: a piston shock whose magnetic field makes the angle theta_Bn with the normal.

    python -m pypic.launch -n 8 research/topics/qpar-electron-conditions/run.py \
        --config research/topics/qpar-electron-conditions/configs/<name>.yaml --out <dir>

Needs the engine of PR #52 (PIC-D12 / PIC-D18 amended: `external_fields` together with `B0` / `E0`; PIC-D13 /
PIC-D9 amended: a listed test species enters through an injecting face). V1 also needs the engine of PR #53
(PIC-D13 amended: the ids of a step cover the planes of the injecting faces); without it the planner refuses V1.

Geometry (README, Approach): quasi-1-D, x is the shock normal, y and z are periodic and a few cells wide. The
field lies in the x-y plane, B = b (cos theta_Bn, sin theta_Bn, 0). Electrons and ions flow at -v0 along x, so
the box is the downstream frame and v0 is the piston speed v_pt of `gupta2024electron`. The motional field is
E0 = (0, 0, v0 b sin theta_Bn).

- ``wall: true``: a reflecting wall at x- stops the flow and the shock moves into the plasma. The parts of the
  field have one place each (docs/plans/2026-11-14-oblique-piston-test-injection.md):
  - the normal component is an external field (the wall zeroes a normal Yee B);
  - the tangential component and E0 are Yee fields (a wall does not short an external E);
  - the face x+ injects the flow (PIC-D13) and carries the uniform B0 / E0 in through a PIC-D26 table of zero
    amplitude. The absorbing layer there is ``npml`` = 192 cells deep. The tangential field of the whole
    upstream (B_y, and E_z = v0 B_y with it) rises at a rate set by that depth: in flow-through pilots of V1's
    plasma (no wall, 40 d_i, to 8 / Omega_ci) B_y / B0y rose by 0.31 per 1 / Omega_ci with 12 cells, 0.019 with
    24 (linear to the end), 0.0010 with 48; with 96 and 192 the fitted rate was within +-0.0003, of either sign
    with the window, and the step time was the same at every depth (notebook, 2026-11-15). Not measured: runs
    longer than 8 / Omega_ci, other plasmas.
  - ``face`` is the injection rule of that face. ``paired`` (default): the inflow face of a flow, the
    electrons are the neutraliser. ``flux``: every species is launched on its own. With ``flux`` the tangential
    field at the face decayed from 0.4-0.6 / Omega_ci on in the step-4 pilots (flow at 3 electron thermal
    speeds, 30 and 8 degrees); with ``paired`` the field kept its shape (B_z within 0.01 B0y to
    8 / Omega_ci) and B_y changed only by the rise above (notebook, 2026-11-15). Seeds enter as at a flux face
    with either rule.
- ``through: true`` (with ``wall: true``): x- is an open face instead of the wall. The flow leaves through it
  and no shock forms: the control of the upstream state that the face x+ sets.
- ``wall: false``: every axis periodic, the whole field in the Yee arrays, no shock: the uniform flow of the
  null test N1 (numerical heating).

Units: lengths of the parameters in d_i = sqrt(mass_ratio) c / omega_pe, times in 1 / Omega_ci, speeds in c.
v_A = b / sqrt(mass_ratio) (ion mass only), beta_s = 2 T_s / b^2 per species, so the beta of
`gupta2024electron` (electrons and ions together, T_e = T_i) is beta_e + beta_i.

``mach_alfven`` is v_sh / v_A with v_sh the shock speed in the upstream frame. ``mach_convention`` says how
v_sh follows from v0:
- ``r4``: v_sh = 4 v0 / 3, the compression 4 that `gupta2024electron` assumes for its M_A and M_s;
- ``mhd``: v_sh = v0 r / (r - 1) with r the compression of the fast shock of ideal MHD at theta_Bn
  (gamma = 5/3, beta = beta_e + beta_i).
The report gives the Mach numbers of both conventions.

``seeds``: passive test electrons (PIC-D9), e.g. ``[{"name": "seed", "theta": 0.02}]``: charge -1, mass 1,
the drift of the flow, temperature ``theta`` (T / m_e c^2). They are loaded with the plasma and enter through
the face with it. Each seed species holds as many particles as the electrons.
"""
from __future__ import annotations

import math
import os
from dataclasses import asdict, dataclass, field

from pypic.config import (
    ID_LIMIT,
    Beam,
    Checkpoint,
    CurrentFilter,
    ExternalFields,
    Injection,
    OutputConfig,
    SimConfig,
    SpeciesConfig,
)
from pypic.injection import ids_per_step
from pypic.parallel.decomp import auto_dims
from pypic.setup.planner import apply_overrides

ROOT = os.path.dirname(os.path.abspath(__file__))
GAMMA = 5.0 / 3.0
CONVENTIONS = ("r4", "mhd")
FACES = ("paired", "flux")
BYTES_PER_SLOT = {"compact": 56}  # PIC-D19; the other precisions are not estimated here


@dataclass
class TopicParams:
    mass_ratio: float = 100.0
    v0: float = 0.1  # flow speed in the box, c (the piston speed v_pt of gupta2024electron)
    mach_alfven: float = 20.0  # v_sh / v_A, v_sh in the upstream frame
    mach_convention: str = "r4"
    theta_bn: float = 30.0  # degrees
    beta_e: float = 0.5
    beta_i: float = 0.5
    seeds: list = field(default_factory=list)  # [{"name": str, "theta": float}], passive test electrons
    wall: bool = True
    through: bool = False  # with a wall configuration: x- open instead of the wall (no shock; control run)
    face: str = "paired"  # injection rule of the x+ face with a wall (PIC-D13): "paired" or "flux"
    lx: float = 100.0  # d_i
    t_end: float = 10.0  # 1 / Omega_ci
    ny: int = 2
    nz: int = 2
    cells_per_de: float = 10.0
    cfl: float = 0.45  # c dt / dx
    ppc: int = 125  # per species and cell: 2 ny nz ppc plasma particles per cell column
    filter_passes: int = 32
    filter_compensate: bool = False
    marder_zone: float = 0.125
    npml: int = 192  # cells of the absorbing layer of every open face (see the docstring)
    precision: str = "compact"
    capacity: float | None = None  # capacity_factor; default 1.2 times the expected compression
    seed: int = 12345
    profiles: int = 200  # lineouts along x
    spectra: int = 50  # histograms (x, u) per species
    spectrum_cells: int = 500  # cells of x per bin of the histograms
    velocity_slab: list = field(default_factory=list)  # [x0, x1] in d_i: electrons (u_x, u_y, u_z) there
    velocity_outputs: int = 40
    track_stride: int = 0  # electrons tracked: every n-th (0: none)
    track_every: int = 20  # steps between track samples
    checkpoint_every: int = 0  # steps (0: none)
    # Particle dumps (x, u, weight, id) for N12 and N6 (a) (README; notebook 2026-11-15, "V1's diagnostics").
    # Lengths in d_i, times in 1 / Omega_ci, momenta in m_i v0 (box frame). An empty list: no such dump.
    dump_p: list = field(default_factory=lambda: [0.1, 1.5])  # thresholds p_d,e and p_d,i of the dumps
    dump_origin: list = field(default_factory=list)  # [x0, x1, t0]: e above p_d,e, spectra outputs of t >= t0
    dump_wall: list = field(default_factory=list)  # [depth, period]: e above p_d,e within depth of the wall
    dump_edge: list = field(default_factory=list)  # [depth, period]: e, i above p_d from depth inside the layer to x+
    dump_sample: list = field(default_factory=list)  # [stride, x0, x1, t, ...]: every stride-th e id, no threshold
    overrides: dict = field(default_factory=dict)


def params_from_dict(d: dict) -> TopicParams:
    d = dict(d)
    d["seeds"] = [dict(s) for s in d.get("seeds", [])]
    for k in ("velocity_slab", "dump_origin", "dump_wall", "dump_edge", "dump_sample"):
        d[k] = list(d.get(k, []))
    if "dump_p" in d:
        d["dump_p"] = list(d["dump_p"])
    return TopicParams(**d)


def _jump(r: float, m: float, th: float, beta: float, g: float = GAMMA) -> float:
    """Zero at the compression r of an ideal MHD shock: Alfvén Mach number m of the normal flow, angle th
    between B and the normal, plasma beta (all species)."""
    c2, s2, m2 = math.cos(th) ** 2, math.sin(th) ** 2, m * m
    return ((m2 - r * c2) ** 2 * (r * g * beta / 2.0 + 0.5 * m2 * (r * (g - 1.0) - (g + 1.0)))
            + 0.5 * m2 * s2 * r * ((g + r * (2.0 - g)) * m2 - r * c2 * ((g + 1.0) - r * (g - 1.0))))


def compression(m: float, th: float, beta: float, g: float = GAMMA, n: int = 20000) -> float | None:
    """Compression of the fast shock: the largest root in (1, (g + 1) / (g - 1)] with m^2 > r cos^2 th.
    None when there is none (the flow is slower than the fast speed)."""
    rmax = (g + 1.0) / (g - 1.0)
    best, prev, step = None, _jump(1.0 + 1e-9, m, th, beta, g), (rmax - 1.0) / n
    for i in range(1, n + 1):
        r = 1.0 + step * i
        cur = _jump(r, m, th, beta, g)
        if prev * cur <= 0.0 and m * m > r * math.cos(th) ** 2:
            a, b = r - step, r
            for _ in range(80):
                mid = 0.5 * (a + b)
                if _jump(a, m, th, beta, g) * _jump(mid, m, th, beta, g) <= 0.0:
                    b = mid
                else:
                    a = mid
            best = 0.5 * (a + b)
        prev = cur
    return best


def shock_state(p: TopicParams) -> dict:
    """b, v_A and the two readings of the shock (conventions ``r4`` and ``mhd``) for the parameters."""
    th, beta = math.radians(p.theta_bn), p.beta_e + p.beta_i
    if p.mach_convention == "r4":
        r_set = 4.0
    else:
        # r depends on M_A only: v_sh = v0 r / (r - 1) fixes v_A afterwards
        r_set = compression(p.mach_alfven, th, beta)
        if r_set is None:
            raise ValueError(f"no fast shock at M_A = {p.mach_alfven:g}, theta_Bn = {p.theta_bn:g} deg, "
                             f"beta = {beta:g}")
    v_sh = p.v0 * r_set / (r_set - 1.0)
    v_a = v_sh / p.mach_alfven
    # the other reading of the same plasma: the MHD compression at the flow's own Mach number, found by iteration
    r_mhd = r_set
    for _ in range(200):
        r_new = compression(p.v0 * r_mhd / (r_mhd - 1.0) / v_a, th, beta)
        if r_new is None or abs(r_new - r_mhd) < 1e-12:
            break
        r_mhd = r_new
    return {"v_a": v_a, "b": v_a * math.sqrt(p.mass_ratio), "r_set": r_set, "v_sh": v_sh,
            "r_mhd": r_mhd, "v_sh_mhd": p.v0 * r_mhd / (r_mhd - 1.0), "mach_r4": 4.0 * p.v0 / 3.0 / v_a,
            "mach_mhd": p.v0 * r_mhd / (r_mhd - 1.0) / v_a}


def plan_topic(p: TopicParams, n_gpus: int = 1):
    if p.mach_convention not in CONVENTIONS:
        raise ValueError(f"mach_convention must be one of {CONVENTIONS}")
    if not 0.0 <= p.theta_bn < 90.0:
        raise ValueError("theta_bn must be in [0, 90) degrees")
    if p.face not in FACES:
        raise ValueError(f"face must be one of {FACES}")
    sh = shock_state(p)
    b, v_a = sh["b"], sh["v_a"]
    th = math.radians(p.theta_bn)
    bx, by = b * math.cos(th), b * math.sin(th)
    d_i = math.sqrt(p.mass_ratio)
    omega_ci = b / p.mass_ratio
    theta_e = 0.5 * p.beta_e * b * b
    theta_i = 0.5 * p.beta_i * b * b / p.mass_ratio  # T_i / (m_i c^2)
    dx = 1.0 / p.cells_per_de
    dt = p.cfl * dx
    lx = p.lx * d_i
    nx = int(round(lx / dx))
    if abs(nx * dx - lx) > 1e-9:
        raise ValueError("lx must be a whole number of cells")
    n_steps = int(math.ceil(p.t_end / omega_ci / dt))
    u0 = p.v0 / math.sqrt(1.0 - p.v0 * p.v0)
    drift = (-u0, 0.0, 0.0)
    sp = [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, drift, theta_e)]),
          SpeciesConfig("i", 1.0, p.mass_ratio, [Beam(1.0, drift, theta_i)])]
    for s in p.seeds:
        if set(s) != {"name", "theta"}:
            raise ValueError("a seed is {'name': ..., 'theta': ...}")
        sp.append(SpeciesConfig(str(s["name"]), -1.0, 1.0, [Beam(1.0, drift, float(s["theta"]))], test=True))
    names = [s.name for s in sp]
    electrons = [n for n in names if n != "i"]

    u_e = 4.0 * math.sqrt(theta_e)  # thermal scale of the electron histogram
    diags = [
        {"name": "profile", "kind": "lineout", "axis": "x", "average": False,
         "fields": ["n:e", "n:i", "Bx", "By", "Bz", "Ex", "Ey", "Ez", "flux:e", "flux:i", "ekin:e", "ekin:i"],
         "every": {"count": p.profiles}},
    ]
    for n in names:
        lo, hi = (1e-2 * u_e, 1e3) if n != "i" else (1e-2 * math.sqrt(theta_i) + 1e-3 * p.v0, 1e2)
        diags.append({"name": f"spectrum_{n}", "kind": "histogram", "species": [n], "quantities": ["x", "u"],
                      "bins": [max(1, nx // p.spectrum_cells), 240], "ranges": [[0.0, lx], [lo, hi]],
                      "log": [False, True], "every": {"count": p.spectra}})
    if p.velocity_slab:
        # 64 bins of 0.2 sqrt(theta_e) per axis around the loaded drift (notebook, 2026-11-15, third review)
        x0, x1 = (float(v) * d_i for v in p.velocity_slab)
        w = 6.4 * math.sqrt(theta_e)
        diags.append({"name": "velocity_e", "kind": "histogram", "species": ["e"], "quantities": ["ux", "uy", "uz"],
                      "bins": [64, 64, 64], "ranges": [[-u0 - w, -u0 + w], [-w, w], [-w, w]],
                      "region": [[x0, x1], None, None], "every": {"count": p.velocity_outputs}})
    if p.track_stride > 0:
        diags.append({"name": "tracks", "kind": "tracks", "species": electrons, "stride": p.track_stride,
                      "every": p.track_every, "with_fields": True})
    # particle dumps of N12 and N6 (a); the thresholds are box-frame momenta, as gamma of each species
    p_e, p_i = (float(v) for v in p.dump_p)
    above = {"e": {"gamma_min": math.hypot(1.0, p_e * p.mass_ratio * p.v0)},
             "i": {"gamma_min": math.hypot(1.0, p_i * p.v0)}}
    per_gyration = 1.0 / (omega_ci * dt)  # steps per 1 / Omega_ci
    layer = lx - p.npml * dx  # inner edge of the absorbing layer of x+
    dumps = {"p_d": {"e": p_e, "i": p_i}}

    def dump(name, species, x0, x1, every, select=None, **more):
        diags.append({"name": name, "kind": "particles", "species": [species], "region": [[x0, x1], None, None],
                      "select": select or {}, "every": every, **more})

    if p.dump_origin:
        x0, x1, t0 = (float(v) for v in p.dump_origin)
        every = max(1, n_steps // p.spectra)  # the steps of the histograms
        dump("origin_e", "e", x0 * d_i, x1 * d_i, every, above["e"], t_start=t0 / omega_ci)
        dumps["origin"] = {"x_di": [x0, x1], "every": every, "t0": t0}
    if p.dump_wall:
        depth, period = (float(v) for v in p.dump_wall)
        every = max(1, int(round(period * per_gyration)))
        dump("wall_e", "e", 0.0, depth * d_i, every, above["e"])
        dumps["wall"] = {"depth_di": depth, "every": every}
    if p.dump_edge:
        depth, period = (float(v) for v in p.dump_edge)
        every = max(1, int(round(period * per_gyration)))
        for n in ("e", "i"):
            dump(f"edge_{n}", n, layer - depth * d_i, lx, every, above[n])
        dumps["edge"] = {"x_di": [layer / d_i - depth, p.lx], "layer_di": layer / d_i, "every": every}
    if p.dump_sample:
        stride, x0, x1 = int(p.dump_sample[0]), float(p.dump_sample[1]), float(p.dump_sample[2])
        if math.gcd(stride, p.ppc) != 1:
            raise ValueError("the stride of dump_sample must share no factor with ppc (ids are cell * ppc + slot)")
        steps = [int(round(float(t) * per_gyration)) for t in p.dump_sample[3:]]
        if not steps or min(steps) < 1 or max(steps) > n_steps:
            raise ValueError("dump_sample needs times within the run")
        for k, s in enumerate(steps):  # one output each: the only multiple of s inside the window is s
            dump(f"sample_e{k}", "e", x0 * d_i, x1 * d_i, s, stride=stride, t_start=(s - 0.5) * dt,
                 t_end=(s + 0.5) * dt)
        dumps["sample"] = {"stride": stride, "x_di": [x0, x1], "steps": steps}

    r_cap = max(sh["r_set"], sh["r_mhd"])
    zero = {"cos": [0.0, 0.0, 0.0], "sin": [0.0, 0.0, 0.0]}
    if p.wall:
        bc = (("open" if p.through else "reflecting", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
        table = {"face": "x+", "modes": [{"k": [0.0, 0.0, 0.0], "omega": 0.0, "phase": 0.0, "B": zero, "E": zero,
                                          "u": {n: zero for n in names}}]}
        fields = dict(B0=(0.0, by, 0.0), E0=(0.0, 0.0, p.v0 * by), external_fields=ExternalFields(B=(bx, 0.0, 0.0)),
                      injection=Injection(faces={"x+": p.face}, species=tuple(names), fluctuation=table,
                                          neutralizer="e" if p.face == "paired" else None),
                      marder_zone=p.marder_zone)
        capacity = p.capacity if p.capacity is not None else (1.3 if p.through else 1.2 * r_cap)
    else:
        bc = (("periodic", "periodic"),) * 3
        fields = dict(B0=(bx, by, 0.0), E0=(0.0, 0.0, p.v0 * by))
        capacity = p.capacity if p.capacity is not None else 1.3
    try:
        dims = auto_dims(n_gpus, (nx, p.ny, p.nz))
    except ValueError:
        dims = None
    diag = max(1, n_steps // 1000)  # the engine evaluates Gauss and div B only at steps that are multiples of both
    cfg = SimConfig(
        n_cells=(nx, p.ny, p.nz), dx=(dx,) * 3, dt=dt, n_steps=n_steps, species=sp, ppc=p.ppc, seed=p.seed,
        dims=dims, boundaries=bc, npml=p.npml, precision=p.precision,
        current_filter=CurrentFilter(p.filter_passes, p.filter_compensate), capacity_factor=capacity,
        diag_interval=diag, gauss_interval=10 * diag, sort_interval=50,
        checkpoint=Checkpoint(every=p.checkpoint_every),
        output=OutputConfig(diagnostics=diags, timescales={"ion_gyration": 2.0 * math.pi / omega_ci}), **fields)
    if p.overrides:
        cfg = apply_overrides(cfg, p.overrides)

    t_run = n_steps * dt
    v_box = p.v0 / (sh["r_mhd"] - 1.0)  # shock speed in the box, MHD reading
    x_shock = v_box * t_run / d_i
    columns = nx * p.ny * p.nz
    particles = len(sp) * columns * p.ppc
    # particle ids (PIC-D8, PIC-D13): the loaded ones, then one block per step for the entering particles; the
    # engine stops a run at the step where an id would reach 2^53, so such a run is refused here
    ids = sum(len(s.beams) for s in cfg.species) * math.prod(cfg.n_cells) * cfg.ppc
    if p.wall:
        per_step = ids_per_step(cfg)
        ids += cfg.n_steps * per_step
        if ids >= ID_LIMIT:
            stop = (ID_LIMIT - (ids - cfg.n_steps * per_step)) // per_step * cfg.dt * omega_ci
            raise ValueError(f"the run needs {ids:.3g} particle ids, {ids / ID_LIMIT:.2f} times the limit 2^53: the "
                             f"engine would stop it at {stop:.1f} / Omega_ci of {omega_ci * cfg.n_steps * cfg.dt:.1f}")
    notes = []
    shock = p.wall and not p.through
    if shock and x_shock > 0.8 * p.lx:
        notes.append(f"the shock reaches {x_shock:.0f} d_i of the {p.lx:g} d_i box (MHD reading)")
    debye = math.sqrt(theta_e) / dx
    if debye < 1.0:
        notes.append(f"Debye length {debye:.2f} dx: the numerical heating of this grid is measured by N1")
    slot = BYTES_PER_SLOT.get(cfg.precision)
    info = {
        "b": b, "bx": bx, "by": by, "v_a": v_a, "v0": p.v0, "u0": u0, "d_i": d_i, "omega_ci": omega_ci,
        "theta_e": theta_e, "theta_i": theta_i, "beta": p.beta_e + p.beta_i, "shock": sh,
        "mach_sonic_i": sh["v_sh"] / math.sqrt(theta_i), "mach_sonic_e": sh["v_sh"] / math.sqrt(theta_e),
        "debye_cells": debye, "wall": shock, "through": p.wall and p.through,
        "face": p.face if p.wall else None, "lx": lx,
        "omega_ci_t_end": omega_ci * t_run,
        "x_shock_end": x_shock, "particles": particles, "per_column": 2 * p.ny * p.nz * p.ppc,
        "capacity_bytes": None if slot is None else slot * particles * cfg.capacity_factor,
        "id_fraction": ids / ID_LIMIT, "dumps": dumps,
        "dims": cfg.dims, "notes": notes, "params": asdict(p)}
    lines = [
        "=== qpar-electron-conditions planner ===",
        f"plasma: m_i/m_e={p.mass_ratio:g}, b={b:.5f} (omega_pe/Omega_ce={1.0 / b:.2f}), theta_Bn={p.theta_bn:g} deg, "
        f"beta_e={p.beta_e:g}, beta_i={p.beta_i:g}, v_A={v_a:.6f} c, Debye length {debye:.2f} dx",
        f"flow v0={p.v0:g} c = {p.v0 / v_a:.2f} v_A; M_A={p.mach_alfven:g} set by convention {p.mach_convention}: "
        f"r4 reading M_A={sh['mach_r4']:.2f}; MHD reading r={sh['r_mhd']:.3f}, M_A={sh['mach_mhd']:.2f}, "
        f"shock at {v_box:.5f} c in the box",
        f"Mach numbers on v_sh={sh['v_sh']:.5f} c: sonic (ions, sqrt(T/m_i)) {info['mach_sonic_i']:.1f}, "
        f"electron sonic {info['mach_sonic_e']:.2f}",
        (f"wall at x-, {p.face} face at x+ (npml {cfg.npml}); shock at {x_shock:.0f} d_i at t_end (MHD reading)"
         if shock else f"open x-, {p.face} face at x+ (npml {cfg.npml}); flow-through control, no shock"
         if p.wall else "periodic box, no shock (N1)"),
        f"grid {cfg.n_cells} ({p.lx:g} d_i), dx={dx:g}, dt={dt:g}, {n_steps} steps to "
        f"{info['omega_ci_t_end']:.1f} / Omega_ci; filter {cfg.current_filter.passes} passes"
        f"{' + compensator' if cfg.current_filter.compensate else ''}; precision {cfg.precision}",
        f"particles at t = 0: {particles:,} ({len(sp)} species, {info['per_column']} plasma particles per column); "
        f"capacity factor {cfg.capacity_factor:.2f}"
        + ("" if slot is None else f", {info['capacity_bytes'] / 2 ** 30:.1f} GiB of particle storage"),
    ]
    if len(dumps) > 1:
        parts = [f"{k} {v.get('x_di', [0.0, v.get('depth_di')])} d_i "
                 + (f"at steps {v['steps']}" if "steps" in v else f"every {v['every']} steps")
                 for k, v in dumps.items() if k != "p_d"]
        lines.append(f"particle dumps above p_d = {p_e:g} (e), {p_i:g} (i) m_i v0 (sample: every "
                     f"{dumps.get('sample', {}).get('stride', '-')}-th id, no threshold): " + "; ".join(parts)
                     + ". The size of these in the output plan below ignores the region.")
    lines += [f"NOTE: {n}" for n in notes]
    return cfg, "\n".join(lines), info
