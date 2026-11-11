"""CPU tests of radial beams (PIC-D17): the combinations a beam rejects and the support of a shape."""
from __future__ import annotations

import math

import pytest

from pypic import Beam
from pypic.profiles import support

SHELL = {"axis": "r", "shape": "gauss", "center": 3.0, "width": 0.5}
ORIGIN = (4.0, 4.0, 4.0)


def test_radial_beam_combinations():
    Beam(1.0, profile=SHELL, origin=ORIGIN, drift_radial=0.2)
    Beam(1.0, origin=ORIGIN, drift_radial=-0.2)  # a uniform plasma that converges on the origin
    Beam(1.0, drift=(0.1, 0.0, 0.0), profile=SHELL, origin=ORIGIN)  # a shell in a uniform drift
    bad = [
        dict(profile=SHELL),  # a radial profile without an origin
        dict(drift_radial=0.2),  # a radial drift without an origin
        dict(origin=ORIGIN),  # an origin that nothing uses
        dict(origin=ORIGIN[:2], drift_radial=0.2),
        dict(origin=ORIGIN, drift_radial=0.2, drift=(0.1, 0.0, 0.0)),
        dict(origin=ORIGIN, drift_radial=0.2, perturbation=[{"axis": "x", "mode": 1, "velocity": 0.01}]),
        dict(profile={**SHELL, "axis": "w"}, origin=ORIGIN),
    ]
    for kw in bad:
        with pytest.raises(ValueError):
            Beam(1.0, **kw)


def test_support_bounds_the_radius_at_which_a_shape_is_cut():
    """Particles are loaded where the shape is >= 1e-10: the support is at least that radius and
    exceeds it by less than one sample of the shape; a shape that is not cut reaches rmax and a
    shape that is cut everywhere gives 0."""
    rmax = 20.0
    edge = SHELL["center"] + SHELL["width"] * math.sqrt(math.log(1e10))
    assert edge <= support(SHELL, rmax) <= edge + 2 * rmax / 2 ** 16
    assert support({"shape": "step", "x0": 1.0, "x1": 2.5}, rmax) >= 2.5
    assert support({"shape": "uniform"}, rmax) == rmax
    assert support(SHELL, 2.0) == 2.0  # the shell begins inside rmax and is not cut at it
    assert support({"shape": "gauss", "center": 40.0, "width": 0.5}, rmax) == 0.0
