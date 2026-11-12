"""CPU tests of the bunch driver (PIC-D21): a list of profiles and the interval a shape loads."""
from __future__ import annotations

import math

import pytest

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
