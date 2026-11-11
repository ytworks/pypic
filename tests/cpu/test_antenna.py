"""CPU tests of the laser antenna (PIC-D20): the current is a discrete curl computed from global
indices, its amplitude is exact on the Yee grid, and the combinations that are rejected."""
from __future__ import annotations

import math

import numpy as np
import pytest

from pypic import SimConfig
from pypic.antenna import Antenna, check_config, wavenumber
from pypic.config import Body, CurrentFilter, Laser

FLAT = {"shape": "flat", "ramp": 8.0}


def _current(ant, shape, t):
    J = np.zeros((3, *shape))
    for c, i, a in ant.pieces(t):
        J[c, i] += a
    return J


@pytest.mark.parametrize("pol", ["y", "z"])
def test_antenna_current_has_no_divergence_and_does_not_depend_on_the_block(pol):
    """The current of a focused pulse on a block with ghosts, y periodic and z not: its discrete
    divergence vanishes to round-off at every node, it is periodic along y bit for bit (so the
    divergence vanishes at the seam of the box too), and a block cut out of the grid holds the
    same bits as the same edges of the larger block."""
    n, dx = (40, 24, 28), (0.1, 0.25, 0.2)
    cfg = SimConfig(n_cells=n, dx=dx, dt=0.05,
                    boundaries=(("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting")),
                    lasers=[Laser(plane_x=1.234, a0=0.5, omega=6.0, polarization=pol, waist=1.2, focus_x=4.0,
                                  center=(2.6, 3.1), phase=0.3,
                                  envelope={"shape": "gaussian", "duration": 2.0, "t_peak": 5.0})])
    lo = (-2, -2, -2)
    shape = tuple(m + 4 for m in n)
    index = [np.arange(a, a + m, dtype=np.float64) for a, m in zip(lo, shape)]
    ant = Antenna(np, cfg.lasers[0], cfg, index)
    J = _current(ant, shape, 4.3)
    peak = np.max(np.abs(J[1 if pol == "y" else 2]))
    assert peak > 1.0 and np.max(np.abs(J[0])) > 0.01 * peak  # the pulse is on, and the closing current too
    div = ((J[0][1:, 1:, 1:] - J[0][:-1, 1:, 1:]) / dx[0] + (J[1][1:, 1:, 1:] - J[1][1:, :-1, 1:]) / dx[1]
           + (J[2][1:, 1:, 1:] - J[2][1:, 1:, :-1]) / dx[2])
    assert np.max(np.abs(div)) < 1e-12 * peak / dx[0]
    assert np.array_equal(J[:, :, :4], J[:, :, n[1]:]) and np.max(np.abs(J[:, :, :4])) > 1e-3 * peak
    off, m = (9, 0, 5), (7, 26, 11)  # the part holds both nodes of the antenna (cells 12 and 13 along x)
    part = Antenna(np, cfg.lasers[0], cfg, [v[o:o + k] for v, o, k in zip(index, off, m)])
    cut = (slice(None), *(slice(o, o + k) for o, k in zip(off, m)))
    assert np.array_equal(_current(part, m, 4.3), J[cut])


@pytest.mark.parametrize("pol", ["y", "z"])
def test_plane_wave_of_the_antenna_on_the_1d_yee_grid(pol):
    """A 1-D leapfrog with the antenna current, 5.8 cells per wavelength: after the ramp the field
    is E = +-a0 omega cos(omega t - k |x - x_M| + phase) on the two sides, with the wavenumber of
    the grid. The amplitude factor tan(k dx / 2) differs from sin(k dx / 2) by 14 % here."""
    nx, dx, dt, omega, a0, phase = 600, 0.1, 0.05, 2.0 * math.pi / 0.6, 0.3, 0.4
    cfg = SimConfig(n_cells=(nx, 1, 1), dx=(dx, 0.2, 0.3), dt=dt,
                    boundaries=(("reflecting", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic")),
                    lasers=[Laser(plane_x=(nx // 2 + 0.3) * dx, a0=a0, omega=omega, polarization=pol, phase=phase,
                                  envelope=FLAT)])
    ant = Antenna(np, cfg.lasers[0], cfg, [np.arange(nx + 1.0), np.arange(1.0), np.arange(1.0)])
    comp = 1 if pol == "y" else 2
    s = -1.0 if pol == "y" else 1.0  # dE/dt = s dB/dx - J, dB/dt = s dE/dx; B[i] at i + 1/2
    E, B = np.zeros(nx + 1), np.zeros(nx + 1)
    steps = int(0.45 * nx * dx / dt)  # the wave has not reached a wall
    for step in range(steps):
        J = _current(ant, (nx + 1, 1, 1), (step + 0.5) * dt)
        assert not J[0].any()
        B[:-1] += 0.5 * s * dt / dx * (E[1:] - E[:-1])
        E[1:-1] += dt * (s * (B[1:-1] - B[:-2]) / dx - J[comp, 1:-1, 0, 0])
        B[:-1] += 0.5 * s * dt / dx * (E[1:] - E[:-1])
    t = steps * dt
    x = (np.arange(nx + 1) - (nx // 2 + 0.5)) * dx
    behind_the_ramp = (np.abs(x) > dx) & (np.abs(x) < 0.8 * (t - FLAT["ramp"]))
    k = wavenumber(omega, dx, dt)
    ref = np.sign(x) * a0 * omega * np.cos(omega * t - k * np.abs(x) + phase)
    assert np.max(np.abs(E - ref)[behind_the_ramp]) < 2e-3 * a0 * omega


def test_lasers_reject_what_the_antenna_does_not_support():
    def cfg(plane_x=2.0, omega=6.0, **kw):
        laser = Laser(plane_x=plane_x, a0=0.1, omega=omega, envelope=FLAT)
        kw.setdefault("boundaries", (("open", "open"), ("periodic", "periodic"), ("periodic", "periodic")))
        return SimConfig(n_cells=(40, 16, 16), dx=(0.1, 0.1, 0.1), dt=0.05, lasers=[laser], **kw)

    assert check_config(cfg()) and not check_config(SimConfig())
    with pytest.raises(ValueError, match="explicit"):
        check_config(cfg(solver="implicit"))
    with pytest.raises(ValueError, match="no wave"):  # fewer than pi cells per wavelength
        check_config(cfg(omega=21.0))
    # the nodes of the antenna are checked nodes: with a filter of reach 3 the cells 4 .. 35
    flt = CurrentFilter(passes=2, compensate=True)
    assert check_config(cfg(plane_x=0.45, current_filter=flt), retreat=3, reach=3)
    for x in (0.35, 3.65):
        with pytest.raises(ValueError, match="open face"):
            check_config(cfg(plane_x=x, current_filter=flt), retreat=3, reach=3)
    with pytest.raises(ValueError, match="open face"):  # a node of the antenna on a wall
        check_config(cfg(plane_x=0.05, boundaries=(("reflecting", "reflecting"),) + (("periodic",) * 2,) * 2))
    with pytest.raises(ValueError, match="body"):
        check_config(cfg(body=Body(center=(2.5, 0.8, 0.8), radius=0.3)))
    with pytest.raises(ValueError, match="envelope"):
        Laser(plane_x=1.0, a0=0.1, omega=6.0, envelope={"shape": "flat", "duration": 1.0})
