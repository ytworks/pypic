"""Laser antenna (PIC-D20): the current of a magnetisation on one layer of cells, ``J = curl M``.

For E along y the magnetisation is M_z on the faces (i_a + 1/2, j + 1/2, k) of the layer of cells
``i_a = floor(plane_x / dx)``; with ``S = -M_z`` its curl on the Yee edges is

    J_y(i_a) = S / dx,    J_y(i_a + 1) = -S / dx,    J_x(i_a + 1/2, j) = -(S(j + 1/2) - S(j - 1/2)) / dy

and the same with y -> z for E along z (S = M_y). The discrete divergence of this current is zero
at every node, so the antenna leaves no charge and Gauss's law is untouched. The update of E
becomes dE/dt = curl (B - M): the layer emits the same wave towards +x and towards -x (E odd, B
even about the layer). On the Yee grid a magnetisation ``M cos(omega t)`` gives the plane wave
``E = M tan(k dx / 2) sin(omega t - k |x - x_M|)`` on the + side, with ``x_M = (i_a + 1/2) dx`` and
the wavenumber k of the grid at omega; this fixes the amplitude of S exactly.

A Gaussian beam is emitted as its field on the plane, at the distance d = focus_x - x_M before
the focus: the amplitude w0 / w(d) exp(-r^2 / w(d)^2), the curved front as an advance in time
k_d r^2 / (2 omega R(d)) of the carrier and of the envelope, and the Gouy phase atan(d / z_R).
The beam is the paraxial beam of the grid: a wave with a small transverse wavenumber q has
k_x = k - q^2 / (2 k_d) on the Yee grid with ``k_d = sin(k dx) / dx``, so the Rayleigh length is
z_R = k_d w0^2 / 2 (k_d -> k for dx -> 0).

Everything is computed from global cell indices, so a rank's block equals the same part of the
whole box bit for bit; on a periodic axis the index is wrapped first, which makes M periodic.

Pure Python / array code: ``xp`` is NumPy (CPU tests) or CuPy (the engine).
"""
from __future__ import annotations

import math

import numpy as np

from .body import gauss_radius


def wavenumber(omega: float, dx: float, dt: float) -> float:
    """Wavenumber along x of the Yee grid at ``omega``: sin(omega dt / 2) / dt = sin(k dx / 2) / dx."""
    s = dx / dt * math.sin(0.5 * omega * dt)
    if not (omega * dt < math.pi and s < 1.0):
        raise ValueError(f"the grid carries no wave at the laser omega = {omega}: "
                         "(dx / dt) sin(omega dt / 2) must be < 1")
    return 2.0 / dx * math.asin(s)


def diffraction_wavenumber(k: float, dx: float) -> float:
    """k_d of the paraxial beam on the Yee grid: k_x = k - q^2 / (2 k_d) for a small transverse q."""
    return math.sin(k * dx) / dx


def plane_cell(laser, dx) -> int:
    """Index i_a of the layer of cells that holds the antenna; the nodes i_a and i_a + 1 carry J."""
    return math.floor(laser.plane_x / dx[0])


def first_cell(cfg, retreat: int = 0) -> int:
    """The lowest cell along x that may hold an antenna when x- is not periodic: off the wall and
    off the nodes next to an open face that the Gauss check leaves out."""
    return 1 + (retreat if cfg.boundaries[0][0] == "open" else 0)


def front_advance(laser, dx, dt) -> float:
    """c of the advance in time c r^2 of the curved front at the plane: k_d r^2 / (2 omega R(d)),
    d the distance from the plane to the focus (0 for a plane wave)."""
    if laser.waist is None:
        return 0.0
    d = laser.focus_x - (plane_cell(laser, dx) + 0.5) * dx[0]
    kd = diffraction_wavenumber(wavenumber(laser.omega, dx[0], dt), dx[0])
    zr = 0.5 * kd * laser.waist ** 2
    return 0.5 * kd / laser.omega * d / (d * d + zr * zr)


def check_config(cfg, retreat: int = 0, reach: int = 0) -> bool:
    """True if ``cfg`` has lasers; raises on what is not supported. ``retreat``: nodes next to an
    open face that are left out of the Gauss check; ``reach``: the reach of the current filter."""
    if not cfg.lasers:
        return False
    if cfg.solver != "explicit":
        raise ValueError("lasers are for the explicit solver only (PIC-D20)")
    n, dx = cfg.n_cells[0], cfg.dx
    for laser in cfg.lasers:
        wavenumber(laser.omega, dx[0], cfg.dt)
        ia = plane_cell(laser, dx)
        if not 0 <= ia < n:
            raise ValueError("laser plane_x is outside the box")
        # the two nodes of the antenna are nodes where Gauss's law is checked: not on a wall and
        # not among the nodes next to an open face (PIC-D15, PIC-D16)
        lo, hi = cfg.boundaries[0]
        if lo != "periodic":
            first = first_cell(cfg, retreat)
            last = n - 1 - (retreat if hi == "open" else 0)
            if ia < first or ia + 1 > last:
                raise ValueError(f"laser plane_x must lie in the cells {first} .. {last - 1} along x: off the "
                                 "walls and off the nodes next to an open face (PIC-D20)")
        body = cfg.body
        if body is not None:
            if abs((ia + 0.5) * dx[0] - body.center[0]) <= gauss_radius(body, dx, reach) + dx[0]:
                raise ValueError("the laser plane crosses the body (PIC-D20)")
    return True


class Antenna:
    """One laser on one block. ``index``: per axis the global cell index of every index of the
    block's arrays (ghosts included)."""

    def __init__(self, xp, laser, cfg, index):
        dx, n = cfg.dx, cfg.n_cells
        periodic = [b[0] == "periodic" for b in cfg.boundaries]

        def wrap(a, v):
            return np.mod(v, float(n[a])) if periodic[a] else v

        self.xp = xp
        self.pol = pol = 1 if laser.polarization == "y" else 2  # component of E = axis of the stagger of S
        self.dx = (dx[0], dx[pol])
        ia = float(plane_cell(laser, dx))
        ix = wrap(0, np.asarray(index[0], dtype=np.float64))
        self.rows = (np.nonzero(ix == ia)[0].tolist(), np.nonzero(ix == wrap(0, ia + 1.0))[0].tolist())
        self._ix, self._ia = ix, ia

        k = wavenumber(laser.omega, dx[0], cfg.dt)
        peak = laser.a0 * laser.omega / math.tan(0.5 * k * dx[0])
        center = laser.center if laser.center is not None else (0.5 * n[1] * dx[1], 0.5 * n[2] * dx[2])
        pos = []
        for a in (1, 2):
            v = np.asarray(index[a], dtype=np.float64)
            if a == pol:  # one more value below: J_x differences S along this axis
                v = np.concatenate([v[:1] - 1.0, v])
            pos.append((wrap(a, v) + (0.5 if a == pol else 0.0)) * dx[a] - center[a - 1])
        y, z = xp.asarray(pos[0]).reshape(-1, 1), xp.asarray(pos[1]).reshape(1, -1)
        r2 = y * y + z * z
        self.omega = laser.omega
        if laser.waist is None:
            self.amp = xp.full(r2.shape, peak)
            self.advance = xp.zeros(r2.shape)
            self.phase = laser.phase
        else:
            d = laser.focus_x - (ia + 0.5) * dx[0]
            zr = 0.5 * diffraction_wavenumber(k, dx[0]) * laser.waist ** 2
            w2 = laser.waist ** 2 * (1.0 + (d / zr) ** 2)
            self.amp = peak * math.sqrt(laser.waist ** 2 / w2) * xp.exp(-r2 / w2)
            self.advance = front_advance(laser, dx, cfg.dt) * r2  # k_d r^2 / (2 omega R)
            self.phase = laser.phase - math.atan(d / zr)
        self.envelope = dict(laser.envelope)

    def set_shift(self, shift: int, on: bool = True):
        """The rows of the plane after the moving window (PIC-D22) has gone ``shift`` cells: the
        plane stays at its lab cell. ``on`` False switches the antenna off."""
        lab = self._ix + float(shift)
        self.rows = ((np.nonzero(lab == self._ia)[0].tolist(), np.nonzero(lab == self._ia + 1.0)[0].tolist())
                     if on else ([], []))

    def _envelope(self, t):
        xp, env = self.xp, self.envelope
        if env["shape"] == "gaussian":
            u = (t - env["t_peak"]) / env["duration"]
            return xp.exp(-u * u)
        s = xp.sin(0.5 * math.pi * xp.clip(t / env["ramp"], 0.0, 1.0))
        return s * s

    def pieces(self, t: float):
        """The current at time ``t`` as [(component, index along x, 2-D array over y and z)]."""
        tr = t + self.advance
        S = self.amp * self._envelope(tr) * self.xp.sin(self.omega * tr + self.phase)
        if self.pol == 1:
            hi, lo = S[1:, :], S[:-1, :]
        else:
            hi, lo = S[:, 1:], S[:, :-1]
        out = []
        if self.rows[0]:
            jp, jx = hi / self.dx[0], (lo - hi) / self.dx[1]
            out += [(c, i, a) for i in self.rows[0] for c, a in ((self.pol, jp), (0, jx))]
        if self.rows[1]:
            jm = -hi / self.dx[0]
            out += [(self.pol, i, jm) for i in self.rows[1]]
        return out
