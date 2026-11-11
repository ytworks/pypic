"""CPU tests of examples/wakefield/theory.py: the field of the Gaussian bunch and its start
(PIC-D21), the linear wakes of a pulse and of a bunch, the group velocity on the Yee grid."""
from __future__ import annotations

import importlib.util
import math
import os

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import erf, exp1, k0


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


def test_laser_wake_is_the_oscillator_driven_by_the_pulse():
    """phi'' + phi = a_env^2 / 4 integrated through the pulse leaves the amplitude of `laser_wake`
    (E_x = -d phi / d zeta has the amplitude of phi at k_w = 1); the wake is largest for
    duration = 2, and E_r / E_x at r is 4 r / w^2."""
    T = _theory()
    a0, tau = 0.1, 2.7
    sol = solve_ivp(lambda t, y: [y[1], 0.25 * a0 ** 2 * math.exp(-2.0 * (t / tau) ** 2) - y[0]], (-30.0, 30.0),
                    [0.0, 0.0], rtol=1e-11, atol=1e-16)
    ex, er = T.laser_wake(a0, tau, 4.0)
    assert math.hypot(*sol.y[:, -1]) == pytest.approx(ex, rel=1e-7) and er == 0.0
    best = max(np.linspace(1.0, 3.0, 201), key=lambda d: T.laser_wake(a0, d, 4.0)[0])
    assert best == pytest.approx(2.0, abs=0.011)
    ex_r, er_r = T.laser_wake(a0, tau, 4.0, r=2.0)
    assert ex_r == pytest.approx(ex * math.exp(-0.5)) and er_r / ex_r == pytest.approx(0.5)
    assert T.beam_at(3.0, 0.1, 4.0, 3.0) == pytest.approx((0.1 / math.sqrt(2.0), 4.0 * math.sqrt(2.0)))


def test_yee_group_velocity_has_the_limits_of_the_plasma_and_of_the_grid():
    """A fine grid gives sqrt(1 - omega_p^2 / omega^2); without a plasma dt = dx gives 1; and it is
    the slope of the dispersion relation of the grid."""
    T = _theory()
    assert T.yee_group_velocity(10.0, 1e-4, 5e-5) == pytest.approx(math.sqrt(0.99), rel=1e-6)
    assert T.yee_group_velocity(10.0, 0.02, 0.02, omega_p=0.0) == pytest.approx(1.0, rel=1e-12)
    dx, dt = 0.02, 0.01

    def k(w):
        return 2.0 / dx * math.asin(math.sqrt((dx / dt * math.sin(0.5 * w * dt)) ** 2 - 0.25 * dx * dx))

    slope = 2e-4 / (k(10.0 + 1e-4) - k(10.0 - 1e-4))
    assert T.yee_group_velocity(10.0, dx, dt) == pytest.approx(slope, rel=1e-7)
    assert slope < math.sqrt(0.99)
    with pytest.raises(ValueError):
        T.yee_group_velocity(0.9, dx, dt)


def test_beam_wake_has_the_limits_of_a_sheet_and_of_a_thin_bunch():
    """On the axis R is closed; a wide bunch has the wake of a sheet (the field behind a thin sheet
    is its charge per area) with the radial shape of the bunch; a thin bunch has K_0."""
    T = _theory()
    for s in (0.3, 1.0, 2.5):
        assert T.beam_wake_profile(0.0, s)[0] == pytest.approx(0.5 * s * s * math.exp(0.5 * s * s) * exp1(0.5 * s * s),
                                                               rel=1e-8)
    r = np.array([0.0, 20.0, 40.0])
    assert np.allclose(T.beam_wake_profile(r, 40.0), np.exp(-0.5 * (r / 40.0) ** 2), rtol=3e-3)
    r = np.array([0.5, 1.0, 2.0])
    assert np.allclose(T.beam_wake_profile(r, 0.01, 1.3), (1.3 * 0.01) ** 2 * k0(1.3 * r), rtol=1e-3)
    gamma, sx = 1e4, 1e-3
    assert T.beam_wake(0.05, sx, 40.0, gamma)[0] == pytest.approx(0.05 * math.sqrt(2.0 * math.pi) * sx, rel=2e-3)
    kw = 1.0 / math.sqrt(1.0 - 1.0 / 9.0)
    ratio = T.beam_wake(0.05, 1.0, 1.0, 3.0, [0.0, 1.0]) / (0.05 * math.sqrt(2.0 * math.pi) * math.exp(-0.5 * kw * kw))
    assert np.allclose(ratio, T.beam_wake_profile([0.0, 1.0], 1.0, kw), rtol=1e-12)
