"""Simulation configuration (normalised units: c = eps0 = 1, lengths in c/omega_pe,
time in 1/omega_pe, charge in e, mass in m_e, density in n0 = reference electron density)."""
from __future__ import annotations

from dataclasses import MISSING, asdict, dataclass, field, fields

PERTURBATION_KEYS = {"axis", "mode", "phase", "density", "displacement", "velocity", "velocity_component"}


@dataclass
class Beam:
    """One drifting Maxwell-Juttner component of a species.

    ``perturbation``: list of sinusoidal modes applied after loading (PIC-D9), each a dict with
    ``axis`` (x|y|z, periodic), ``mode`` (k = 2 pi mode / L), ``phase`` and one or more of
      density: alpha       n = density * (1 + alpha cos(k x + phase)) exactly (particles moved)
      displacement: xi     x = x0 + xi sin(k x0 + phase) (Lagrangian displacement)
      velocity: u1         u_c += u1 sin(k x0 + phase), c = velocity_component (default: axis)
    Moves are deposited as currents, so the initial E satisfies Gauss's law exactly."""

    density: float  # in units of the reference density n0
    drift: tuple = (0.0, 0.0, 0.0)  # proper velocity u = gamma v / c
    theta: float = 0.0  # temperature T / (m c^2)
    perturbation: list = field(default_factory=list)
    # density profile along one axis: n(x) = density * shape(x) >= 0 (pypic.profiles), e.g.
    # {"axis": "x", "shape": "step", "x0": 0.0, "x1": 50.0}; beams of all species share positions,
    # so equal profiles keep the plasma exactly neutral; particles where shape < 1e-10 are dropped
    profile: dict | None = None

    def __post_init__(self):
        self.drift = tuple(self.drift)
        if self.profile is not None:
            from .profiles import validate

            self.profile = validate(self.profile, extra=("axis",))
            if self.profile.get("axis", "x") not in ("x", "y", "z"):
                raise ValueError("profile axis must be x, y or z")
        self.perturbation = [dict(p) for p in self.perturbation]
        for p in self.perturbation:
            unknown = set(p) - PERTURBATION_KEYS
            if unknown:
                raise KeyError(f"unknown perturbation key(s) {sorted(unknown)}")
            if p.get("axis", "x") not in ("x", "y", "z") or p.get("velocity_component", "x") not in ("x", "y", "z"):
                raise ValueError("perturbation axis / velocity_component must be x, y or z")
            if int(p.get("mode", 1)) < 1:
                raise ValueError("perturbation mode must be >= 1")
            if "density" in p and "displacement" in p:
                raise ValueError("give either density or displacement in one perturbation")
            if "density" in p and not abs(float(p["density"])) < 1.0:
                raise ValueError("density perturbation amplitude must satisfy |alpha| < 1")
            if not ({"density", "displacement", "velocity"} & set(p)):
                raise ValueError("perturbation needs density, displacement or velocity")


@dataclass
class SpeciesConfig:
    name: str
    charge: float  # in units of e
    mass: float  # in units of m_e
    beams: list = field(default_factory=list)  # list[Beam]
    # test species (PIC-D9): feel the fields, deposit no charge or current, excluded from the
    # energy / momentum budget and from Gauss's law; may be given particle by particle
    test: bool = False
    particles: list = field(default_factory=list)  # test species: [[x, y, z, ux, uy, uz(, w)], ...]

    def __post_init__(self):
        if self.particles and not self.test:
            raise ValueError(f"species {self.name!r}: explicit particle lists are for test species only")
        self.particles = [list(map(float, p)) for p in self.particles]
        for p in self.particles:
            if len(p) not in (6, 7) or (len(p) == 7 and p[6] <= 0):
                raise ValueError(f"species {self.name!r}: particles are [x, y, z, ux, uy, uz(, w > 0)]")


@dataclass
class CurrentFilter:
    """Digital filter of the deposited current (PIC-D10): ``passes`` passes of the binomial
    (1, 2, 1)/4 stencil along every axis, optionally followed by the compensator
    (-N/4, 1 + N/2, -N/4) that restores the long-wavelength response to fourth order. Suppresses
    the numerical Cherenkov instability of relativistic flows. Gauss's law stays exact (the charge
    density in the check is filtered identically) and the implicit energy stays exact (particles
    gather the filtered field). Periodic and reflecting faces only."""

    passes: int = 0
    compensate: bool = False


@dataclass
class ImplicitConfig:
    picard_tol: float = 1e-12  # relative change of E^{n+1} between Picard iterations
    picard_maxit: int = 100
    particle_tol: float = 1e-13  # per-particle fixed-point tolerance (in cells)
    particle_maxit: int = 100
    cg_tol: float = 1e-14
    cg_maxit: int = 1000


# ----------------------------------------------------------------------------- output
DIAG_KINDS = ("fields", "slice", "lineout", "modes", "kspectrum", "histogram", "particles", "tracks")
_COMMON = {"name", "kind", "every", "t_start", "t_end", "dtype", "keep_bits"}
# keys that are meaningful for each kind (validated in Diagnostic.from_dict / __post_init__)
DIAG_KEYS = {
    "fields": _COMMON | {"fields", "coarsen", "region"},
    "slice": _COMMON | {"fields", "axis", "position", "coarsen"},
    "lineout": _COMMON | {"fields", "axis", "average"},
    "modes": _COMMON | {"fields", "modes", "average"},
    "kspectrum": _COMMON | {"fields", "coarsen"},
    "histogram": _COMMON | {"species", "quantities", "bins", "ranges", "log", "region"},
    "particles": _COMMON | {"species", "stride", "select", "region"},
    "tracks": _COMMON | {"species", "stride", "with_fields"},
}


@dataclass
class Diagnostic:
    """One output diagnostic. Cadence ``every``: an int (steps), ``{"t": dt}``, ``{"count": n}``,
    ``{"per": <timescale>, "samples": n}``, ``{"nyquist": <timescale>, "oversample": s}`` or
    ``"auto"`` (resolved from ``OutputConfig.timescales`` by :mod:`pypic.output.schedule`).
    Field names: ``E``, ``B``, ``J`` (vectors) or components ``Ex``..``Jz``, ``rho``, and species
    moments ``n:<species>``, ``flux:<species>`` (vector n<v>), ``ekin:<species>``.
    See docs/design/output-analysis.md for the meaning of every key."""

    name: str
    kind: str
    every: object = "auto"
    t_start: float = 0.0
    t_end: float | None = None
    dtype: str = "float32"  # float32 | float64
    keep_bits: int | None = None  # lossy: keep this many mantissa bits (relative error <= 2^-(n+1))
    fields: tuple = ()
    species: tuple = ()  # empty: all species
    axis: str = "x"  # lineout: profile axis; slice: plane normal
    position: float | None = None  # slice plane coordinate in c/omega_pe (None: box centre)
    average: bool | None = None  # lineout/modes: mean over the output interval (None: auto)
    coarsen: tuple = (1, 1, 1)  # block-average factors
    region: tuple | None = None  # ((x0, x1), (y0, y1), (z0, z1)) in c/omega_pe; None = whole box
    modes: tuple = ()  # ((mx, my, mz), ...) with k = 2 pi m / L
    quantities: tuple = ()  # histogram axes: x y z ux uy uz u gamma ekin
    bins: tuple = ()
    ranges: tuple = ()  # per axis (lo, hi) or None (automatic)
    log: tuple = ()  # per axis: logarithmic bins
    stride: int = 1  # particles / tracks: keep particles with id % stride == 0
    select: dict = field(default_factory=dict)  # particles: {"gamma_min": g}
    with_fields: bool = False  # tracks: also record E and B at the particle

    def __post_init__(self):
        if self.kind not in DIAG_KINDS:
            raise ValueError(f"diagnostic {self.name!r}: kind must be one of {DIAG_KINDS}")
        for f in fields(Diagnostic):
            if f.name in DIAG_KEYS[self.kind]:
                continue
            dv = f.default if f.default_factory is MISSING else f.default_factory()
            if getattr(self, f.name) != dv:
                raise ValueError(f"diagnostic {self.name!r}: key {f.name!r} does not apply to kind {self.kind!r}")
        for k in ("fields", "species", "coarsen", "modes", "quantities", "bins", "ranges", "log"):
            v = getattr(self, k)
            if isinstance(v, str):
                v = (v,)
            setattr(self, k, tuple(tuple(x) if isinstance(x, list) else x for x in v))
        if self.region is not None:
            self.region = tuple(tuple(r) if r is not None else None for r in self.region)
        if self.dtype not in ("float32", "float64"):
            raise ValueError(f"diagnostic {self.name!r}: dtype must be float32 or float64")
        if self.keep_bits is not None and not 1 <= int(self.keep_bits) <= (23 if self.dtype == "float32" else 52):
            raise ValueError(f"diagnostic {self.name!r}: keep_bits out of range")
        if self.axis not in ("x", "y", "z"):
            raise ValueError(f"diagnostic {self.name!r}: axis must be x, y or z")
        if len(self.coarsen) != 3 or any(int(c) < 1 for c in self.coarsen):
            raise ValueError(f"diagnostic {self.name!r}: coarsen needs three factors >= 1")
        if self.stride < 1:
            raise ValueError(f"diagnostic {self.name!r}: stride must be >= 1")

    @classmethod
    def from_dict(cls, d) -> "Diagnostic":
        if isinstance(d, Diagnostic):
            return d
        d = dict(d)
        known = {f.name for f in fields(cls)}
        unknown = set(d) - known
        if unknown:
            raise KeyError(f"diagnostic {d.get('name')!r}: unknown key(s) {sorted(unknown)}")
        return cls(**d)


@dataclass
class OutputConfig:
    """File output (docs/design/output-analysis.md). ``dir=None`` writes nothing: the run keeps
    its in-memory history only (tests, quick experiments)."""

    dir: str | None = None
    codec: str = "zstd"  # zstd (default, needs hdf5plugin) | deflate (built into every HDF5) | none
    level: int = 3
    budget_GB: float | None = None  # storage budget: "auto" cadences are widened to fit, fixed ones raise
    n0_m3: float | None = None  # reference electron density for SI units in the files (None: normalised)
    timescales: dict = field(default_factory=dict)  # name -> time in 1/omega_pe (set by the planners)
    scalars: bool = True  # write the energy / conservation history (every diag_interval steps)
    diagnostics: list = field(default_factory=list)  # list[Diagnostic]

    def __post_init__(self):
        if self.codec not in ("zstd", "deflate", "none"):
            raise ValueError("output.codec must be zstd, deflate or none")
        self.diagnostics = [Diagnostic.from_dict(x) for x in self.diagnostics]
        names = [x.name for x in self.diagnostics]
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate diagnostic names in {names}")

    @classmethod
    def from_dict(cls, d) -> "OutputConfig":
        if isinstance(d, OutputConfig):
            return d
        d = dict(d or {})
        known = {f.name for f in fields(cls)}
        unknown = set(d) - known
        if unknown:
            raise KeyError(f"unknown output field(s) {sorted(unknown)}")
        return cls(**d)


@dataclass
class SimConfig:
    n_cells: tuple = (64, 16, 16)
    dx: tuple = (0.1, 0.1, 0.1)
    dt: float = 0.05
    n_steps: int = 1000
    solver: str = "explicit"  # "explicit" | "implicit"
    boundaries: tuple = (("periodic", "periodic"),) * 3  # per axis (lo, hi): periodic | open | reflecting
    npml: int = 12
    species: list = field(default_factory=list)  # list[SpeciesConfig]
    ppc: int = 16  # particle positions per cell (every beam gets one particle per position)
    seed: int = 12345
    dims: tuple | None = None  # rank grid (px,py,pz); None = automatic
    B0: tuple = (0.0, 0.0, 0.0)  # uniform external/initial magnetic field
    E0: tuple = (0.0, 0.0, 0.0)  # uniform initial electric field (no conducting face tangential to it)
    # initial fields added to B0 / E0 (PIC-D9), each divergence-free by construction:
    #  {"kind": "plane_wave", "axis": "x", "mode": 1, "polarization": "y", "amplitude": a, "phase": 0,
    #   "direction": 1}   exact discrete eigenmode of the solver (Yee or Crank-Nicolson)
    #  {"kind": "profile", "component": "Bx", "axis": "y", "shape": "tanh", ...}  (pypic.profiles; a
    #   component may only vary along another axis)
    initial_fields: list = field(default_factory=list)
    # subtract the box-averaged current every step (= a uniform external neutralising current):
    # removes the k = 0 plasma oscillation driven by a net beam current (Buneman); periodic boxes only
    remove_mean_current: bool = False
    diag_interval: int = 10  # energy / momentum history (and the scalars output) every N steps
    gauss_interval: int = 100  # 0 disables the Gauss/divB check
    sort_interval: int = 50
    implicit: ImplicitConfig = field(default_factory=ImplicitConfig)
    current_filter: CurrentFilter = field(default_factory=CurrentFilter)
    output: OutputConfig = field(default_factory=OutputConfig)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SimConfig":
        d = dict(d)
        sp = []
        for s in d.pop("species", []):
            # accept SpeciesConfig instances (Python API overrides) as well as plain dicts (YAML)
            s = asdict(s) if isinstance(s, SpeciesConfig) else dict(s)
            beams = [b if isinstance(b, Beam) else Beam(**b) for b in s.pop("beams", [])]
            sp.append(SpeciesConfig(beams=beams, **s))
        imp = d.pop("implicit", {})
        imp = imp if isinstance(imp, ImplicitConfig) else ImplicitConfig(**imp)
        out = OutputConfig.from_dict(d.pop("output", {}))
        flt = d.pop("current_filter", {})
        flt = flt if isinstance(flt, CurrentFilter) else CurrentFilter(**flt)
        for k in ("n_cells", "dx", "B0", "E0", "dims"):
            if d.get(k) is not None:
                d[k] = tuple(d[k])
        if "boundaries" in d:
            d["boundaries"] = tuple(tuple(b) for b in d["boundaries"])
        return cls(species=sp, implicit=imp, output=out, current_filter=flt, **d)
