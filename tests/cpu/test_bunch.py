"""CPU tests of the bunch driver (PIC-D21): a list of profiles, the interval a shape loads, and the
field of the Gaussian bunch (examples/wakefield/theory.py)."""
from __future__ import annotations

import importlib.util
import math
import os

import numpy as np
import pytest
from scipy.special import erf

from pypic import Beam
from pypic.profiles import support_interval

GAUSS = {"shape": "gauss", "center": 6.0, "width": 0.5}


def test_a_list_of_profiles_needs_one_profile_per_axis():
    b = Beam(1.0, drift=(3.0, 0.0, 0.0), profile=[{**GAUSS, "axis": a} for a in "xyz"])
    assert [p["axis"] for p in b.profile] == ["x", "y", "z"]
    Beam(1.0, profile=[{**GAUSS, "axis": "z"}])  # a slab loaded into the cells it reaches
    bad = [
        [],
        [{**GAUSS, "axis": "x"}, {**GAUSS, "axis": "x"}],
        [GAUSS, {**GAUSS, "axis": "x"}],  # the default axis is x
        [{**GAUSS, "axis": "r"}],
    ]
    for profile in bad:
        with pytest.raises(ValueError):
            Beam(1.0, profile=profile)


def test_support_interval_bounds_where_a_shape_is_cut():
    """Particles are loaded where the shape is >= 1e-10: both bounds lie outside that interval by
    less than two samples of the shape; a shape that is not cut reaches the ends, and a shape that
    is cut everywhere gives an empty interval."""
    length, sample = 20.0, 20.0 / 2 ** 16
    half = GAUSS["width"] * math.sqrt(math.log(1e10))
    lo, hi = support_interval(GAUSS, length)
    assert GAUSS["center"] - half - 2 * sample <= lo <= GAUSS["center"] - half
    assert GAUSS["center"] + half <= hi <= GAUSS["center"] + half + 2 * sample
    assert support_interval({"shape": "uniform"}, length) == (0.0, length)
    assert support_interval({**GAUSS, "center": 0.5}, length)[0] == 0.0
    assert support_interval({**GAUSS, "center": 40.0}, length) == (0.0, 0.0)


def _theory():
    spec = importlib.util.spec_from_file_location(
        "wakefield_theory", os.path.join(os.path.dirname(__file__), "..", "..", "examples", "wakefield", "theory.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_bunch_field_has_the_limits_of_a_point_a_sphere_and_a_line():
    """E_r in the mid-plane of the Gaussian bunch: a small bunch at rest is a point charge, a round
    bunch at rest is the Gaussian sphere, a bunch that is long in its rest frame is the line charge."""
    T = _theory()
    r = np.array([0.5, 1.0, 2.0, 3.0])
    assert np.allclose(T.bunch_field(r, 1e-4, 1e-4, 1.0) * 4.0 * math.pi * r ** 2, 1.0, rtol=1e-6)
    sphere = (erf(r / math.sqrt(2.0)) - math.sqrt(2.0 / math.pi) * r * np.exp(-0.5 * r ** 2)) / (4.0 * math.pi * r ** 2)
    assert np.allclose(T.bunch_field(r, 1.0, 1.0, 1.0), sphere, rtol=1e-8)
    line = (1.0 - np.exp(-0.5 * r ** 2)) / (2.0 * math.pi * r) / math.sqrt(2.0 * math.pi)
    assert np.allclose(T.bunch_field(r, 1.0, 1.0, 2000.0), line, rtol=1e-4)
    assert T.bunch_charge(0.05, 1.0, 2.0) == pytest.approx(0.05 * (2.0 * math.pi) ** 1.5 * 4.0)


def test_bunch_startup_goes_from_no_field_to_the_steady_field():
    """A bunch that starts on its partner: no field at t = 0, the steady field once the start has
    left (gamma = 3: after 60 sigma), and on the way B_theta and E_r differ."""
    T = _theory()
    r = np.array([0.5, 1.0, 2.0])
    e, b = T.bunch_startup(r, 0.0, 1.0, 1.0, 3.0, nk=400, nx=401)
    assert np.all(np.abs(e) < 1e-12) and np.all(np.abs(b) < 1e-12)
    e, b = T.bunch_startup(r, 60.0, 1.0, 1.0, 3.0)
    assert np.allclose(e, 1.0, atol=1e-3) and np.allclose(b, 1.0, atol=1e-3)
    e, b = T.bunch_startup(r, 2.0, 1.0, 1.0, 3.0, nk=400, nx=401)
    assert np.all(e > 0.0) and np.all(b > 0.0) and np.all(np.abs(e - b) > 0.01)
