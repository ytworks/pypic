"""Linear theory (pypic.theory) against literature values and exact limits."""
from __future__ import annotations

import math

import numpy as np

from pypic.theory import streaming as S
from pypic.theory import waves as W


def test_landau_damping_roots_match_literature():
    # Maxwellian Langmuir roots, omega / omega_pe at k lambda_De = 0.3 / 0.4 / 0.5
    # (e.g. Canosa 1973; Birdsall & Langdon): 1.1598-0.0126i, 1.2850-0.0661i, 1.4157-0.1534i
    vt = 0.05
    w = W.langmuir(np.array([0.3, 0.4, 0.5]) / vt, vt)
    ref = np.array([1.1598 - 0.0126j, 1.2850 - 0.0661j, 1.4157 - 0.1534j])
    assert np.max(np.abs(w - ref)) < 2e-4


def test_cold_magnetised_cutoffs_and_bernstein_limit():
    wce = 0.5
    p = W.cold_perpendicular(np.array([1e-6]), 1.0, wce)
    # X-mode cutoffs: omega_R,L = (+-wce + sqrt(wce^2 + 4 wpe^2)) / 2
    assert abs(p["X_upper"][0] - (wce + math.sqrt(wce ** 2 + 4)) / 2) < 1e-9
    assert abs(p["X_lower"][0] - (-wce + math.sqrt(wce ** 2 + 4)) / 2) < 1e-9
    # long-wavelength Bernstein branch tends to the upper hybrid frequency
    b = W.bernstein(np.array([1e-3]), 0.05, 1.0, wce, n_bands=2)
    assert abs(b[1, 0] - math.sqrt(1 + wce ** 2)) < 1e-6
    # whistler: omega -> k^2 c^2 wce / wpe^2 for k << wpe/c (corrections O(omega/wce))
    w = W.cold_parallel(np.array([0.05]), 1.0, 2.0)["whistler"][0]
    assert abs(w / (0.05 ** 2 * 2.0) - 1) < 0.02


def test_streaming_limits():
    assert abs(S.buneman_cold(1836.15267) - 0.0561) < 1e-4
    u = 3.0  # filamentation growth tends to beta wp / sqrt(gamma) for k >> wp/c
    g = math.sqrt(1 + u * u)
    assert abs(S.filamentation_cold(1e4, u) - (u / g) / math.sqrt(g)) < 1e-6


def _example_theory(name):
    import importlib.util
    import os

    path = os.path.join(os.path.dirname(__file__), "..", "..", "examples", name, "theory.py")
    spec = importlib.util.spec_from_file_location(f"theory_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_shock_jump_conditions():
    T = _example_theory("relativistic_shock")
    # constant Gamma = 3/2 reproduces Spitkovsky (2008): n2/n1 = 3.13, beta_sh = 0.47 for gamma0 = 15
    g = T.jump_conditions_gamma(15.0, 1.5)
    assert abs(g["compression"] - 3.1333) < 1e-3 and abs(g["beta_sh"] - 0.4677) < 1e-3
    # Synge gas: ultra-relativistic limits n2/n1 -> 4 (3-D, Gamma 4/3) and 3 (2-D, Gamma 3/2)
    assert abs(T.jump_conditions(1e4, 3)["compression"] - 4.0) < 2e-3
    assert abs(T.jump_conditions(1e4, 2)["compression"] - 3.0) < 2e-3


def test_mora_relations():
    T = _example_theory("plasma_expansion")
    Te, mi = 0.01, 1836.0
    t = np.linspace(50.0, 2000.0, 400)
    xf = T.front_position(t, Te, mi)
    vf = T.front_velocity(t, Te, mi)
    rel = (np.gradient(xf, t) - vf) / vf
    assert np.max(np.abs(rel[1:-1])) < 2e-3  # x_f is the integral of v_f (interior differences)
    assert np.allclose(T.max_energy(t, Te, mi), 0.5 * mi * vf ** 2, rtol=1e-12)  # E_max = m_i v_f^2 / 2


def test_profiles():
    from pypic.profiles import evaluate, validate

    y = np.array([0.0, 25.0, 50.0, 75.0, 99.0])
    f = evaluate(validate({"shape": "double_tanh", "centers": [25.0, 75.0], "width": 1.0}), y)
    assert np.allclose(f, [-1.0, 0.0, 1.0, 0.0, -1.0], atol=1e-9)
    assert np.allclose(evaluate({"shape": "step", "x0": 1.0, "x1": 2.0}, np.array([0.5, 1.0, 1.99, 2.0])),
                       [0, 1, 1, 0])
