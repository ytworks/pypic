"""Initial electromagnetic fields (SimConfig.initial_fields, PIC-D9), divergence-free by construction.

plane_wave: {"kind": "plane_wave", "axis": "x", "mode": m, "polarization": "y", "amplitude": A,
             "phase": phi, "direction": +1 | -1}
    E_p = A cos(k x + phi) along `axis` (k = 2 pi m / L), B_q such that E x B points along
    direction * axis. It is the exact discrete eigenmode of the field solver:
      explicit (Yee leapfrog): sin(omega dt / 2) / dt = sin(k dx / 2) / dx; the stored B at integer
          time is the mean of B^(-1/2) and B^(+1/2): amplitude A cos(omega dt / 2)
      implicit (Crank-Nicolson): (2/dt) tan(omega dt / 2) = (2/dx) sin(k dx / 2); amplitude A
profile:    {"kind": "profile", "component": "Bx", "axis": "y", "shape": ..., ...} (pypic.profiles):
    the component varies only along an axis other than its own, so div E / div B stay zero.
"""
from __future__ import annotations

import math

import cupy as cp
import numpy as np

from .output.names import stagger_of
from .profiles import evaluate, validate

AXES = "xyz"


def omega_of(k: float, dx: float, dt: float, implicit: bool) -> float:
    s = dt / dx * math.sin(0.5 * k * dx)
    if implicit:
        return 2.0 / dt * math.atan(s)
    if s > 1.0:
        raise ValueError("plane wave beyond the Yee CFL limit")
    return 2.0 / dt * math.asin(s)


def _line(sim, comp: str, axis: int):
    """Global coordinate (c/omega_pe) of every local array index along `axis` for a component."""
    g = sim.grid
    return g.global_positions(axis, stagger_of(comp)[axis]) * g.dx[axis]


def _add(sim, comp: str, axis: int, values: np.ndarray):
    shape = [1, 1, 1]
    shape[axis] = -1
    getattr(sim.grid, comp)[...] += cp.asarray(values.reshape(shape))


def _check_tangential_E(cfg, c: int, what: str):
    for a in range(3):
        if a != c and cfg.boundaries[a][0] != "periodic":
            raise ValueError(f"{what}: E_{AXES[c]} is tangential to the conducting faces of axis {AXES[a]}")


def plane_wave(sim, spec: dict):
    cfg, g = sim.cfg, sim.grid
    a = AXES.index(spec.get("axis", "x"))
    p = AXES.index(spec.get("polarization", "y"))
    if p == a:
        raise ValueError("plane wave: polarization must be perpendicular to the axis")
    if cfg.boundaries[a][0] != "periodic":
        raise ValueError("plane wave: the propagation axis must be periodic")
    _check_tangential_E(cfg, p, "plane wave")
    q = 3 - a - p
    cyclic = (a, p, q) in ((0, 1, 2), (1, 2, 0), (2, 0, 1))
    d = int(spec.get("direction", 1))
    if d not in (1, -1):
        raise ValueError("plane wave: direction must be +1 or -1")
    A = float(spec["amplitude"])
    phi = float(spec.get("phase", 0.0))
    k = 2 * math.pi * int(spec.get("mode", 1)) / g.L[a]
    w = omega_of(k, g.dx[a], cfg.dt, sim.implicit)
    bfac = 1.0 if sim.implicit else math.cos(0.5 * w * cfg.dt)
    Ec, Bc = "E" + AXES[p], "B" + AXES[q]
    _add(sim, Ec, a, A * np.cos(k * _line(sim, Ec, a) + phi))
    _add(sim, Bc, a, (1 if cyclic else -1) * d * bfac * A * np.cos(k * _line(sim, Bc, a) + phi))
    return {"k": k, "omega": w, "phase_velocity": w / k}


def profile(sim, spec: dict):
    cfg = sim.cfg
    comp = spec["component"]
    if comp not in ("Ex", "Ey", "Ez", "Bx", "By", "Bz"):
        raise ValueError(f"profile component must be one of Ex..Bz (got {comp!r})")
    b = AXES.index(spec.get("axis", "x"))
    c = AXES.index(comp[1])
    if b == c:
        raise ValueError(f"profile: {comp} may not vary along its own axis (divergence)")
    if comp[0] == "E":
        _check_tangential_E(cfg, c, "field profile")
    shape = validate({k: v for k, v in spec.items() if k not in ("kind", "component", "axis")})
    _add(sim, comp, b, evaluate(shape, _line(sim, comp, b)))


def apply_initial_fields(sim, specs) -> list:
    info = []
    for spec in specs:
        kind = spec.get("kind")
        if kind == "plane_wave":
            info.append(plane_wave(sim, spec))
        elif kind == "profile":
            profile(sim, spec)
            info.append({})
        else:
            raise ValueError(f"initial field kind must be plane_wave or profile (got {kind!r})")
    return info
