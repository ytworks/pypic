"""Simulation configuration (normalised units: c = eps0 = 1, lengths in c/omega_pe,
time in 1/omega_pe, charge in e, mass in m_e, density in n0 = reference electron density)."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass
class Beam:
    """One drifting Maxwell-Juttner component of a species."""

    density: float  # in units of the reference density n0
    drift: tuple = (0.0, 0.0, 0.0)  # proper velocity u = gamma v / c
    theta: float = 0.0  # temperature T / (m c^2)


@dataclass
class SpeciesConfig:
    name: str
    charge: float  # in units of e
    mass: float  # in units of m_e
    beams: list = field(default_factory=list)  # list[Beam]


@dataclass
class ImplicitConfig:
    picard_tol: float = 1e-12  # relative change of E^{n+1} between Picard iterations
    picard_maxit: int = 100
    particle_tol: float = 1e-13  # per-particle fixed-point tolerance (in cells)
    particle_maxit: int = 100
    cg_tol: float = 1e-14
    cg_maxit: int = 1000


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
    diag_interval: int = 10
    gauss_interval: int = 100  # 0 disables the Gauss/divB check
    sort_interval: int = 50
    mode_diag: int = 32  # number of k_x Fourier modes of Ex recorded (0 disables)
    phase_space_interval: int = 0  # (x, u_x) histogram of each species every N steps (0 disables)
    phase_space_bins: tuple = (256, 128)
    phase_space_umax: float = 0.0  # |u_x| range of the histogram (0: automatic from the beams)
    output_dir: str = "output"
    implicit: ImplicitConfig = field(default_factory=ImplicitConfig)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SimConfig":
        d = dict(d)
        sp = []
        for s in d.pop("species", []):
            # accept SpeciesConfig instances (Python API overrides) as well as plain dicts (YAML)
            s = asdict(s) if isinstance(s, SpeciesConfig) else dict(s)
            beams = [b if isinstance(b, Beam) else Beam(**{**b, "drift": tuple(b.get("drift", (0, 0, 0)))})
                     for b in s.pop("beams", [])]
            sp.append(SpeciesConfig(beams=beams, **s))
        imp = d.pop("implicit", {})
        imp = imp if isinstance(imp, ImplicitConfig) else ImplicitConfig(**imp)
        for k in ("n_cells", "dx", "B0", "dims", "phase_space_bins"):
            if d.get(k) is not None:
                d[k] = tuple(d[k])
        if "boundaries" in d:
            d["boundaries"] = tuple(tuple(b) for b in d["boundaries"])
        return cls(species=sp, implicit=imp, **d)
