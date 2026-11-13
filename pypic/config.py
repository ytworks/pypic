"""Simulation configuration (normalised units: c = eps0 = 1, lengths in c/omega_pe,
time in 1/omega_pe, charge in e, mass in m_e, density in n0 = reference electron density)."""
from __future__ import annotations

from dataclasses import MISSING, asdict, dataclass, field, fields

PERTURBATION_KEYS = {"axis", "mode", "phase", "density", "displacement", "velocity", "velocity_component", "shape",
                     "displacement_component"}
SHAPED_KEYS = {"shape", "velocity", "velocity_component", "displacement", "displacement_component"}


def validate_shape_map(shape, what: str) -> dict:
    """A map {axis: 1-D shape (pypic.profiles)} with at least one of the axes x, y, z (PIC-D24)."""
    from .profiles import validate

    if not isinstance(shape, dict) or not shape:
        raise ValueError(f"{what}: 'shape' maps one or more of the axes x, y, z to a 1-D shape")
    out = {}
    for a, spec in shape.items():
        if a not in ("x", "y", "z"):
            raise ValueError(f"{what}: shape axis must be x, y or z (got {a!r})")
        if not isinstance(spec, dict):
            raise ValueError(f"{what}: the shape along {a} must be a dict such as {{'shape': 'cos', ...}}")
        out[a] = validate(spec)
    return out


def validate_shaped_perturbation(p: dict) -> None:
    """Shaped perturbation of a beam (PIC-D24): amplitude x product of 1-D shapes."""
    unknown = set(p) - SHAPED_KEYS
    if unknown:
        raise KeyError(f"shaped perturbation: key(s) {sorted(unknown)} do not go with 'shape' "
                       "(axis, mode and phase belong to the sinusoidal form; density has no shaped form)")
    p["shape"] = validate_shape_map(p["shape"], "shaped perturbation")
    if "velocity" not in p and "displacement" not in p:
        raise ValueError("shaped perturbation needs velocity or displacement")
    for key, comp in (("velocity", "velocity_component"), ("displacement", "displacement_component")):
        if key in p:
            if p.get(comp) not in ("x", "y", "z"):
                raise ValueError(f"shaped perturbation: {comp} (x, y or z) is required with {key}")
            p[key] = float(p[key])
        elif comp in p:
            raise KeyError(f"shaped perturbation: {comp} without {key}")


POTENTIAL_KEYS = {"kind", "component", "shape", "scale"}


def validate_potential(spec: dict) -> dict:
    """Initial field of kind 'potential' (PIC-D24): B += discrete curl A, A_c = scale x product of 1-D
    shapes evaluated on the c-edges (the places of E_c)."""
    spec = dict(spec)
    unknown = set(spec) - POTENTIAL_KEYS
    if unknown:
        raise KeyError(f"initial field potential: unknown key(s) {sorted(unknown)}")
    if spec.get("component") not in ("Ax", "Ay", "Az"):
        raise ValueError(f"initial field potential: component must be Ax, Ay or Az (got {spec.get('component')!r})")
    spec["shape"] = validate_shape_map(spec.get("shape"), "initial field potential")
    spec["scale"] = float(spec.get("scale", 1.0))
    return spec


@dataclass
class Beam:
    """One drifting Maxwell-Juttner component of a species.

    ``perturbation``: list of sinusoidal modes applied after loading (PIC-D9), each a dict with
    ``axis`` (x|y|z, periodic), ``mode`` (k = 2 pi mode / L), ``phase`` and one or more of
      density: alpha       n = density * (1 + alpha cos(k x + phase)) exactly (particles moved)
      displacement: xi     x = x0 + xi sin(k x0 + phase) (Lagrangian displacement)
      velocity: u1         u_c += u1 sin(k x0 + phase), c = velocity_component (default: axis)
    or a shaped perturbation (PIC-D24): ``shape`` maps axes to 1-D shapes (pypic.profiles) whose
    product f(x0) multiplies the amplitudes; the axes may be periodic or not:
      velocity: u1         u_c += u1 f(x0), c = velocity_component (required)
      displacement: xi     x_c = x0_c + xi f(x0), c = displacement_component (required); refused for
                           particles closer than one cell plus the largest move to a non-periodic face
                           of axis c
    Moves are deposited as currents, so the initial E satisfies Gauss's law exactly; displacing one
    species by xi = -E / (q n) loads a target field E (PIC-D24)."""

    density: float  # in units of the reference density n0
    drift: tuple = (0.0, 0.0, 0.0)  # proper velocity u = gamma v / c
    theta: float = 0.0  # temperature T / (m c^2)
    perturbation: list = field(default_factory=list)
    # density profile along one axis: n(x) = density * shape(x) >= 0 (pypic.profiles), e.g.
    # {"axis": "x", "shape": "step", "x0": 0.0, "x1": 50.0}; beams of all species share positions,
    # so equal profiles keep the plasma exactly neutral; particles where shape < 1e-10 are dropped.
    # A list of such profiles with different axes among x, y, z is their product (PIC-D21): a
    # bunch. It is loaded only into the cells in which every shape reaches 1e-10, and the beams of
    # the plasma species that load the same cells must be neutral together (as for a radial profile)
    profile: dict | list | None = None
    # radial beam (PIC-D17): ``origin`` is a point in c/omega_pe. A profile with "axis": "r" is a
    # shape of the distance from it, and ``drift_radial`` is a proper velocity along the unit vector
    # from it (0 at the origin itself). A radial profile is loaded only into the cells its support
    # reaches; the beams of the plasma species that load the same cells must be neutral together
    origin: tuple | None = None
    drift_radial: float = 0.0

    def __post_init__(self):
        self.drift = tuple(self.drift)
        self.drift_radial = float(self.drift_radial)
        if isinstance(self.profile, (list, tuple)):
            from .profiles import validate

            self.profile = [validate(p, extra=("axis",)) for p in self.profile]
            axes = [p.get("axis", "x") for p in self.profile]
            if not axes or any(a not in ("x", "y", "z") for a in axes) or len(set(axes)) != len(axes):
                raise ValueError("a list of profiles needs one profile per axis among x, y, z")
        elif self.profile is not None:
            from .profiles import validate

            self.profile = validate(self.profile, extra=("axis",))
            if self.profile.get("axis", "x") not in ("x", "y", "z", "r"):
                raise ValueError("profile axis must be x, y, z or r")
        radial_profile = isinstance(self.profile, dict) and self.profile.get("axis", "x") == "r"
        if self.origin is None:
            if radial_profile or self.drift_radial:
                raise ValueError("a radial profile or drift_radial needs an origin")
        else:
            self.origin = tuple(float(v) for v in self.origin)
            if len(self.origin) != 3:
                raise ValueError("beam origin needs three components")
            if not (radial_profile or self.drift_radial):
                raise ValueError("a beam origin is only used by a radial profile or drift_radial")
            if self.drift_radial and any(self.drift):
                raise ValueError("give either drift or drift_radial in one beam")
            if self.perturbation:
                raise ValueError("a radial beam must not have a perturbation")
        self.perturbation = [dict(p) for p in self.perturbation]
        for p in self.perturbation:
            unknown = set(p) - PERTURBATION_KEYS
            if unknown:
                raise KeyError(f"unknown perturbation key(s) {sorted(unknown)}")
            if "shape" in p:
                validate_shaped_perturbation(p)
                continue
            if "displacement_component" in p:
                raise KeyError("displacement_component belongs to a shaped perturbation (with 'shape')")
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
    gather the filtered field). An open face is a wall of the filter at the face of the box and
    the Gauss check leaves out the nodes the charge left there reaches (PIC-D15, explicit solver);
    an open flow with a filter needs ``marder_zone`` to stay stationary."""

    passes: int = 0
    compensate: bool = False


@dataclass
class Dipole:
    """Magnetic dipole of the external field (PIC-D12): B = (3 (m.n) n - m) / r^3 for r >> core,
    regular inside ``core`` (B = 2 m / core^3 at the centre). ``center`` in c/omega_pe."""

    moment: tuple
    center: tuple
    core: float

    def __post_init__(self):
        self.moment = tuple(float(v) for v in self.moment)
        self.center = tuple(float(v) for v in self.center)
        self.core = float(self.core)
        if len(self.moment) != 3 or len(self.center) != 3:
            raise ValueError("dipole moment and center need three components")
        if not self.core > 0.0:
            raise ValueError("dipole core must be > 0")


@dataclass
class Laser:
    """Laser antenna (PIC-D20): a magnetisation on the cell layer that holds ``plane_x`` whose
    current, a discrete curl, emits the same pulse towards +x and towards -x. The pulse towards +x
    is a Gaussian beam with the waist ``waist`` at ``focus_x`` (``waist: None``: a plane wave), with
    E along ``polarization`` and the peak amplitude E = a0 * omega at the focus. Explicit solver only.

    ``envelope``: ``{"shape": "gaussian", "duration": tau, "t_peak": t0}`` is exp(-(t - t0)^2 /
    tau^2) for the field at the plane, from t = 0 on (choose t0 of several tau); ``{"shape":
    "flat", "ramp": T}`` rises as sin^2 over T and stays at 1. ``center``: (y, z) of the axis,
    default the middle of the box. ``phase``: the carrier is cos(omega (t - x / v) + phase) on the
    axis at the focus, t counted like the envelope at the plane."""

    plane_x: float
    a0: float
    omega: float
    envelope: dict
    polarization: str = "y"
    waist: float | None = None
    focus_x: float | None = None
    center: tuple | None = None
    phase: float = 0.0

    def __post_init__(self):
        self.plane_x, self.a0, self.omega, self.phase = (
            float(self.plane_x), float(self.a0), float(self.omega), float(self.phase))
        if self.polarization not in ("y", "z"):
            raise ValueError("laser polarization must be 'y' or 'z'")
        if not self.omega > 0.0:
            raise ValueError("laser omega must be > 0")
        if self.waist is not None:
            self.waist = float(self.waist)
            if not self.waist > 0.0:
                raise ValueError("laser waist must be > 0 (None: a plane wave)")
        self.focus_x = self.plane_x if self.focus_x is None else float(self.focus_x)
        if self.center is not None:
            self.center = tuple(float(v) for v in self.center)
            if len(self.center) != 2:
                raise ValueError("laser center is (y, z)")
        env = dict(self.envelope)
        shape = env.pop("shape", None)
        need = {"gaussian": ("duration", "t_peak"), "flat": ("ramp",)}
        if shape not in need:
            raise ValueError("laser envelope shape must be 'gaussian' or 'flat'")
        if set(env) != set(need[shape]):
            raise ValueError(f"laser envelope '{shape}' takes exactly {need[shape]}")
        env = {k: float(v) for k, v in env.items()}
        if not env[need[shape][0]] > 0.0:
            raise ValueError(f"laser envelope {need[shape][0]} must be > 0")
        self.envelope = {"shape": shape, **env}


@dataclass
class Window:
    """Moving window (PIC-D22): from the time ``start`` on the box follows a driver along +x at the
    speed of light, one cell at a time; positions stay those of the laboratory. x- must be open
    and x+ reflecting. The plasma is loaded up to ``margin`` (a length, c/omega_pe) before the x+
    wall, and every cell that enters the box behind this vacuum margin is loaded as at the start
    of the run. Explicit solver only."""

    start: float = 0.0
    margin: float = 1.0

    def __post_init__(self):
        self.start, self.margin = float(self.start), float(self.margin)
        if not self.start >= 0.0:
            raise ValueError("window.start must be >= 0")
        if not self.margin > 0.0:
            raise ValueError("window.margin must be > 0")


@dataclass
class ExternalFields:
    """Static fields the particles feel in addition to the Yee fields (PIC-D12): a uniform E, a
    uniform B and dipoles. The Yee fields are then the difference to them, so an external field
    may be normal or tangential to any face and no wall shorts it (E0 at a reflecting face: PIC-D18).
    Explicit solver only."""

    E: tuple = (0.0, 0.0, 0.0)
    B: tuple = (0.0, 0.0, 0.0)
    dipoles: list = field(default_factory=list)  # list[Dipole]

    def __post_init__(self):
        self.E = tuple(float(v) for v in self.E)
        self.B = tuple(float(v) for v in self.B)
        if len(self.E) != 3 or len(self.B) != 3:
            raise ValueError("external_fields.E and .B need three components")
        self.dipoles = [d if isinstance(d, Dipole) else Dipole(**d) for d in self.dipoles]

    @property
    def active(self) -> bool:
        return any(self.E) or any(self.B) or bool(self.dipoles)


@dataclass
class Body:
    """Absorbing body inside the box (PIC-D14): a particle that ends a step inside it is removed and
    its charge stays where the deposit left it, so the body charges up like an insulator; no
    particle is loaded inside. ``center`` and ``radius`` in c/omega_pe. ``relax`` (in omega_pe,
    0 = off) is the rate at which the electric field of the grid decays on the edges inside the
    body: the charge of the body then moves to its surface and the inside becomes an
    equipotential, with the total charge unchanged; not for use with a dipole (PIC-D14, measured).
    Explicit solver only."""

    center: tuple
    radius: float
    shape: str = "sphere"
    relax: float = 0.0

    def __post_init__(self):
        self.center = tuple(float(v) for v in self.center)
        self.radius = float(self.radius)
        if self.shape != "sphere":
            raise ValueError("body.shape must be sphere")
        if len(self.center) != 3:
            raise ValueError("body.center needs three components")
        if not self.radius > 0.0:
            raise ValueError("body.radius must be > 0")
        self.relax = float(self.relax)
        if not self.relax >= 0.0:
            raise ValueError("body.relax must be >= 0")


FACES = ("x-", "x+", "y-", "y+", "z-", "z+")  # face index = 2 * axis + side
ID_LIMIT = 2 ** 53  # particle ids (PIC-D8): float64 represents every integer below this exactly
# Particle slots of one species on one GPU: the kernels index particles with int (thread index, n);
# row offsets row * cap are 64-bit.
MAX_CAPACITY = 2 ** 31 - 1
# Device bytes by precision, for the memory estimates of the planners. SLOT_BYTES: one particle slot
# (fp64, mixed: 8 rows of float64; compact, PIC-D19: weight and id in float64, the cell in int32,
# three uint32 and three float64). SORT_BYTES: the temporaries of the sort per particle of the
# species being sorted: 40 was measured at the bow-shock gate; one sort of 3.4e7 and of 1.7e7
# particles on one GPU holds 40.2 (fp64, mixed) and 24.2 (compact: the key is the stored int32 cell,
# and the buffer of one row fits a block the argsort has freed).
SLOT_BYTES = {"fp64": 64, "mixed": 64, "compact": 56}
SORT_BYTES = {"fp64": 40, "mixed": 40, "compact": 24}


@dataclass
class Injection:
    """Plasma entering through open faces every step (PIC-D13). ``faces`` maps a face name
    (``x-`` .. ``z+``) to its mode: ``flux`` launches every species from the cell layer outside
    the face; ``paired`` is for an inflow face of a flow: the ``neutralizer`` species (the
    electrons) is returned at the face instead of leaving, and every other injected species
    enters together with a neutralizer particle at the same point, so the face collects no
    charge. The distribution is that of the species' own beams. Explicit solver only."""

    faces: dict = field(default_factory=dict)  # {"x-": "paired", "y-": "flux", ...}
    species: tuple | None = None  # injected species names; None: every plasma species with beams
    neutralizer: str | None = None  # needed by a paired face
    # a prescribed fluctuation carried in through one injecting face (PIC-D26; pypic.injection.check_fluctuation):
    # {"face": "x+", "modes": [{"k": [0, ky, kz], "omega": w, "phase": p, "B": {"cos": [..], "sin": [..]},
    #   "E": {...}, "u": {species: {"cos": [..], "sin": [..]}}}], "tolerance": 1e-3}
    fluctuation: dict | None = None

    def __post_init__(self):
        if not isinstance(self.faces, dict):
            raise ValueError("injection.faces maps a face name to its mode, e.g. {'x-': 'paired', 'y-': 'flux'}")
        self.faces = {str(k): str(v) for k, v in self.faces.items()}
        for k, v in self.faces.items():
            if k not in FACES:
                raise ValueError(f"injection.faces: unknown face {k!r} (one of {FACES})")
            if v not in ("flux", "paired"):
                raise ValueError(f"injection.faces[{k!r}]: mode must be flux or paired")
        if self.species is not None:
            self.species = (self.species,) if isinstance(self.species, str) else tuple(self.species)

    @property
    def active(self) -> bool:
        return bool(self.faces)


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
    # a directory that already holds run output is refused, so the files of two runs never mix;
    # overwrite=True removes that run's output files first (only those, see OutputManager)
    overwrite: bool = False

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
    # places of the ppc particles of a cell (PIC-D23): None = random (seed, global cell, slot);
    # (a, b, c) with a b c = ppc = a regular lattice of a x b x c places, the same in every cell
    lattice: tuple | None = None
    dims: tuple | None = None  # rank grid (px,py,pz); None = automatic
    B0: tuple = (0.0, 0.0, 0.0)  # uniform external/initial magnetic field
    E0: tuple = (0.0, 0.0, 0.0)  # uniform initial electric field; tangential to reflecting faces only (PIC-D18)
    # initial fields added to B0 / E0 (PIC-D9), each divergence-free by construction:
    #  {"kind": "plane_wave", "axis": "x", "mode": 1, "polarization": "y", "amplitude": a, "phase": 0,
    #   "direction": 1}   exact discrete eigenmode of the solver (Yee or Crank-Nicolson)
    #  {"kind": "profile", "component": "Bx", "axis": "y", "shape": "tanh", ...}  (pypic.profiles; a
    #   component may only vary along another axis)
    #  {"kind": "potential", "component": "Az", "shape": {"x": {...}, "y": {...}}, "scale": s}  B += the
    #   discrete curl of A (PIC-D24): div B = 0 bitwise; A tangential to a non-periodic face must be zero
    #   on it
    initial_fields: list = field(default_factory=list)
    # subtract the box-averaged current every step (= a uniform external neutralising current):
    # removes the k = 0 plasma oscillation driven by a net beam current (Buneman); periodic boxes only
    remove_mean_current: bool = False
    diag_interval: int = 10  # energy / momentum history (and the scalars output) every N steps
    gauss_interval: int = 100  # 0 disables the Gauss/divB check
    sort_interval: int = 50
    # "fp64" | "mixed" | "compact". mixed gathers E and B for the explicit push from float32 copies
    # (PIC-D11, Proposed); positions, momenta, the pusher and the current deposit stay FP64.
    # compact (PIC-D19, Proposed; explicit solver) also stores a particle in 56 bytes in place of
    # 64: cell and position within the cell in 2^-32 of a cell; momenta, pusher and current deposit
    # stay FP64
    precision: str = "fp64"
    implicit: ImplicitConfig = field(default_factory=ImplicitConfig)
    current_filter: CurrentFilter = field(default_factory=CurrentFilter)
    # static fields felt by the particles on top of the Yee fields (PIC-D12); not together with B0 / E0
    external_fields: ExternalFields = field(default_factory=ExternalFields)
    injection: Injection = field(default_factory=Injection)  # plasma entering through open faces (PIC-D13)
    body: Body | None = None  # absorbing sphere (PIC-D14)
    lasers: list = field(default_factory=list)  # list[Laser]: laser antennas (PIC-D20)
    window: Window | None = None  # moving window along +x (PIC-D22)
    # Marder number (0 < nu <= 1/6; 0 = off) of the zone at the open faces (PIC-D16): every step
    # div E - rho diffuses in the absorbing layer and the first nodes inside, where the removed and
    # the entering particles leave charge that no deposit balances. Explicit solver only.
    marder_zone: float = 0.0
    # particle array capacity = capacity_factor x the particles loaded on the rank: the room for the
    # density to pile up on one rank. Growing the arrays in a run holds the old and the new array
    # at once, so set this (and watch the ``fill`` history) rather than rely on it.
    capacity_factor: float = 1.3
    output: OutputConfig = field(default_factory=OutputConfig)

    def __post_init__(self):
        if not float(self.capacity_factor) >= 1.0:
            raise ValueError("capacity_factor must be >= 1")
        if not 0.0 <= float(self.marder_zone) <= 1.0 / 6.0:
            raise ValueError("marder_zone is a Marder number: 0 (off) or 0 < nu <= 1/6")
        fields = []
        for f in self.initial_fields:
            f = dict(f)
            if f.get("kind") == "potential":
                f = validate_potential(f)
            fields.append(f)
        self.initial_fields = fields
        if self.lattice is not None:
            lat = tuple(self.lattice)
            if len(lat) != 3 or any(int(v) != v or v < 1 for v in lat) or lat[0] * lat[1] * lat[2] != self.ppc:
                raise ValueError(f"lattice {lat} must be three positive integers with the product ppc = {self.ppc}")
            self.lattice = tuple(int(v) for v in lat)

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
        ext = d.pop("external_fields", {})
        ext = ext if isinstance(ext, ExternalFields) else ExternalFields(**ext)
        inj = d.pop("injection", {})
        inj = inj if isinstance(inj, Injection) else Injection(**inj)
        body = d.pop("body", None)
        body = body if body is None or isinstance(body, Body) else Body(**body)
        lasers = [v if isinstance(v, Laser) else Laser(**v) for v in d.pop("lasers", [])]
        window = d.pop("window", None)
        window = window if window is None or isinstance(window, Window) else Window(**window)
        for k in ("n_cells", "dx", "B0", "E0", "dims", "lattice"):
            if d.get(k) is not None:
                d[k] = tuple(d[k])
        if "boundaries" in d:
            d["boundaries"] = tuple(tuple(b) for b in d["boundaries"])
        return cls(species=sp, implicit=imp, output=out, current_filter=flt, external_fields=ext,
                   injection=inj, body=body, lasers=lasers, window=window, **d)
