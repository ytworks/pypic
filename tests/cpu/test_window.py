"""CPU tests of the moving window (PIC-D22): the schedule of the shifts and the configurations
that are rejected."""
from __future__ import annotations

import math

import pytest

from pypic import SimConfig
from pypic.config import Beam, Body, CurrentFilter, Laser, SpeciesConfig, Window
from pypic.window import antenna_off_shift, check_config, margin_cells, shift_after, steps_since_shift

BC = (("open", "reflecting"), ("periodic", "periodic"), ("periodic", "periodic"))


def cfg(start=2.0, margin=0.8, **kw):
    kw.setdefault("boundaries", BC)
    return SimConfig(n_cells=(40, 16, 16), dx=(0.1, 0.1, 0.1), dt=0.045, window=Window(start, margin), **kw)


def test_the_window_moves_one_cell_per_cell_crossing_time():
    """dt / dx = 0.45: the shift never grows by more than one cell per step, it follows the light
    (within a cell), and the steps between two shifts add up to the steps of the run."""
    c = cfg(start=2.0)
    s = [shift_after(c, j) for j in range(2000)]
    assert s[0] == 0 and set(b - a for a, b in zip(s, s[1:])) == {0, 1}
    assert all(0.0 <= (j * c.dt - c.window.start) / 0.1 - v < 1.0 for j, v in enumerate(s) if v > 0)
    assert all(v == 0 for j, v in enumerate(s) if j * c.dt < c.window.start)
    first = [j for j in range(1, 2000) if s[j] > s[j - 1]]
    since = [steps_since_shift(c, j) for j in first]
    # the first interval is counted from the time ``start`` (the first step at or after it)
    assert first[0] - since[0] == math.ceil(c.window.start / c.dt)
    assert since[1:] == [b - a for a, b in zip(first, first[1:])] and set(since) <= {2, 3}


def test_window_from_dict_and_margin_in_cells():
    c = SimConfig.from_dict({"n_cells": [40, 16, 16], "dx": [0.1, 0.1, 0.1], "dt": 0.045, "boundaries": BC,
                             "window": {"start": 3.0, "margin": 0.62}})
    assert isinstance(c.window, Window) and margin_cells(c) == 6 and check_config(c)
    assert not check_config(SimConfig())
    with pytest.raises(ValueError, match="start"):
        Window(start=-1.0)


@pytest.mark.parametrize("kw, match", [
    ({"solver": "implicit"}, "explicit"),
    ({"boundaries": (("open", "open"),) + BC[1:]}, "open, reflecting"),
    ({"boundaries": (("periodic", "periodic"),) * 3}, "open, reflecting"),
    ({"injection": {"faces": {"x-": "flux"}}}, "injection"),
    ({"body": Body(center=(2.0, 0.8, 0.8), radius=0.3)}, "body"),
    ({"external_fields": {"B": (0.0, 0.0, 0.1)}}, "external_fields"),
    ({"B0": (0.0, 0.0, 0.1)}, "B0"),
    ({"initial_fields": [{"kind": "plane_wave", "axis": "y", "mode": 1, "polarization": "z", "amplitude": 0.1}]},
     "initial_fields"),
    ({"species": [SpeciesConfig("e", -1.0, 1.0, [Beam(1.0, perturbation=[{"axis": "y", "mode": 1,
                                                                          "density": 0.1}])])]}, "perturbation"),
    ({"margin": 0.3}, "at least 4 cells"),
    ({"margin": 3.9}, "no cell for the plasma"),
])
def test_window_rejects_what_it_does_not_support(kw, match):
    kw = dict(kw)
    margin = kw.pop("margin", 0.8)
    c = SimConfig.from_dict({**cfg(margin=margin).to_dict(), **kw})
    with pytest.raises(ValueError, match=match):
        check_config(c)


def test_window_margin_grows_with_the_reach_of_the_filter():
    c = cfg(margin=0.6, current_filter=CurrentFilter(passes=2, compensate=True))
    with pytest.raises(ValueError, match="at least 7 cells"):
        check_config(c, retreat=3, reach=3)
    assert check_config(cfg(margin=0.7, current_filter=CurrentFilter(passes=2, compensate=True)), retreat=3, reach=3)


def test_window_keeps_the_antenna_until_the_pulse_has_ended():
    """The antenna in the cell 20 is switched off at the shift 20 (its cell would be the cell 0 of
    the window), which the window reaches at start + 2.0: a pulse that ends at t_peak + 5 duration
    = 6 needs start >= 4. A front that is curved away from the focus is late off the axis."""
    def with_laser(start, **kw):
        env = kw.pop("envelope", {"shape": "gaussian", "duration": 0.4, "t_peak": 4.0})
        return cfg(start=start, lasers=[Laser(plane_x=2.05, a0=0.1, omega=6.0, envelope=env, **kw)])

    c = with_laser(4.05)
    assert antenna_off_shift(c, c.lasers[0]) == 20 and check_config(c)
    with pytest.raises(ValueError, match="before the laser pulse has ended"):
        check_config(with_laser(3.95))
    # with a filter of reach 3 the cells 0 .. 3 are off: the antenna leaves 3 shifts earlier
    flt = CurrentFilter(passes=2, compensate=True)
    assert check_config(cfg(start=4.35, lasers=c.lasers, current_filter=flt), retreat=3, reach=3)
    with pytest.raises(ValueError, match="before the laser pulse has ended"):
        check_config(cfg(start=4.25, lasers=c.lasers, current_filter=flt), retreat=3, reach=3)
    with pytest.raises(ValueError, match="gaussian envelope"):
        check_config(with_laser(4.05, envelope={"shape": "flat", "ramp": 1.0}))
    assert check_config(with_laser(4.05, waist=0.5, focus_x=3.0))  # converging: the edge of the beam is early
    with pytest.raises(ValueError, match="curved front"):
        check_config(with_laser(4.05, waist=0.5, focus_x=1.0))
