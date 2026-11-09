"""Plasma composition (e-i, e-p, e-i-p) and two-stream beam set-up, always charge neutral.

Normalisation: n0 = total electron density (so omega_pe = 1), temperatures in m_e c^2,
drifts as proper velocity u = gamma v / c.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..config import SpeciesConfig
from ..particles import Beam

MIXTURES = ("e-i", "e-p", "e-i-p")
MODES = ("electron", "symmetric", "buneman")


def u_to_v(u: float) -> float:
    return u / math.sqrt(1.0 + u * u)


def v_to_u(v: float) -> float:
    return v / math.sqrt(1.0 - v * v)


@dataclass
class Composition:
    mixture: str = "e-i"
    mass_ratio: float = 1836.15267  # m_i / m_e
    Z: int = 1  # ion charge number
    positron_fraction: float = 0.5  # e-i-p: fraction of the positive charge carried by positrons

    def densities(self):
        """Number densities (n_e, n_i, n_p) with n_e = 1 and n_e = Z n_i + n_p (neutral)."""
        if self.mixture not in MIXTURES:
            raise ValueError(f"mixture must be one of {MIXTURES}")
        if self.mixture == "e-i":
            return 1.0, 1.0 / self.Z, 0.0
        if self.mixture == "e-p":
            return 1.0, 0.0, 1.0
        fp = float(self.positron_fraction)
        if not 0.0 < fp < 1.0:
            raise ValueError("positron_fraction must be in (0, 1) for e-i-p")
        return 1.0, (1.0 - fp) / self.Z, fp


@dataclass
class TwoStreamBeams:
    """How the streams are arranged.

    mode = "electron":  electrons split into a beam (+u_drift, density fraction beam_fraction) and
                        a return beam with zero net electron current; positive species at rest.
    mode = "symmetric": every species split into two neutral counter-streaming beams (+-u_drift),
                        e.g. colliding plasma shells / pair-plasma two-stream.
    mode = "buneman":   all electrons drift with +u_drift relative to positive species at rest
                        (carries a net current; the k=0 mode oscillates at omega_pe).
    """

    mode: str = "electron"
    u_drift: float = 0.2
    beam_fraction: float = 0.5
    T_e: float = 1e-4  # in m_e c^2
    T_i: float | None = None  # in m_e c^2 (default: T_e)
    T_p: float | None = None


def build_species(comp: Composition, beams: TwoStreamBeams) -> list[SpeciesConfig]:
    if beams.mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    n_e, n_i, n_p = comp.densities()
    m_i = float(comp.mass_ratio)
    Te = beams.T_e
    Ti = beams.T_e if beams.T_i is None else beams.T_i
    Tp = beams.T_e if beams.T_p is None else beams.T_p
    u = float(beams.u_drift)
    f = float(beams.beam_fraction)
    out = []

    def split(n, theta):
        return [Beam(0.5 * n, (u, 0.0, 0.0), theta), Beam(0.5 * n, (-u, 0.0, 0.0), theta)]

    # electrons
    if beams.mode == "electron":
        if not 0.0 < f < 1.0:
            raise ValueError("beam_fraction must be in (0, 1)")
        vb = u_to_v(u)
        vr = -f * vb / (1.0 - f)  # zero net electron current
        e_beams = [Beam(f * n_e, (u, 0.0, 0.0), Te), Beam((1.0 - f) * n_e, (v_to_u(vr), 0.0, 0.0), Te)]
    elif beams.mode == "symmetric":
        e_beams = split(n_e, Te)
    else:  # buneman
        e_beams = [Beam(n_e, (u, 0.0, 0.0), Te)]
    out.append(SpeciesConfig("electron", -1.0, 1.0, e_beams))

    def positive(name, q, m, n, T):
        if n <= 0:
            return
        theta = T / m
        if beams.mode == "symmetric":
            bl = split(n, theta)
        else:
            bl = [Beam(n, (0.0, 0.0, 0.0), theta)]
        out.append(SpeciesConfig(name, q, m, bl))

    positive("ion", float(comp.Z), m_i, n_i, Ti)
    positive("positron", 1.0, 1.0, n_p, Tp)
    net = sum(s.charge * b.density for s in out for b in s.beams)
    assert abs(net) < 1e-12, "composition is not neutral"
    return out


def species_summary(species: list[SpeciesConfig]) -> str:
    lines = []
    for s in species:
        for b in s.beams:
            lines.append(f"  {s.name:9s} q={s.charge:+.0f} m={s.m if hasattr(s, 'm') else s.mass:<10.4g}"
                         f" n={b.density:.4f} u_x={b.drift[0]:+.4f} theta={b.theta:.3e}")
    return "\n".join(lines)
