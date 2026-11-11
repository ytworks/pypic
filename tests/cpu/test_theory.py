"""Linear theory (pypic.theory) against literature values and exact limits."""
from __future__ import annotations

import math

import numpy as np

from pypic.setup import mhd
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


def test_bow_shock_scales():
    T = _example_theory("bow_shock")
    g = mhd.GAMMA
    # the compression comes from the mass, momentum and induction equations; the energy flux
    # rho v^3 / 2 + gamma / (gamma - 1) p v + B^2 v is the equation they leave over
    for ma, beta in ((1.5, 0.1), (4.0, 0.64), (5.0, 0.64), (20.0, 3.0)):
        d = T.downstream(ma, beta)
        x, v2 = d["compression"], d["v"]
        up = 0.5 * ma ** 3 + g / (g - 1) * 0.5 * beta * ma + ma
        down = 0.5 * x * v2 ** 3 + g / (g - 1) * d["p"] * v2 + d["B"] ** 2 * v2
        assert x > 1.0 and abs(down / up - 1.0) < 1e-12, (ma, beta)
    # no shock at the fast magnetosonic Mach number 1, (gamma + 1) / (gamma - 1) for a strong shock
    beta = 0.64
    ma1 = math.sqrt(1.0 + 0.5 * g * beta)
    assert abs(mhd.fast_mach(ma1, beta) - 1.0) < 1e-14 and abs(mhd.compression(ma1, beta) - 1.0) < 1e-12
    assert mhd.compression(0.5 * ma1, beta) == 1.0
    assert abs(mhd.compression(1e4, beta) - 4.0) < 1e-6
    # oblique fast shock: the same test of the energy flux, with the tangential flow and the Poynting
    # flux B_z (v_x B_z - v_z B_x); (1.5, 0.1) along the field is a switch-on shock
    for ma, beta in ((1.5, 0.1), (2.0, 0.05), (4.0, 0.64), (20.0, 3.0)):
        for cos2 in (1e-6, 0.25, 0.5, 0.9, 1.0):
            d = T.oblique_downstream(ma, beta, cos2)
            x, bz1 = d["compression"], math.sqrt(1.0 - cos2)
            up = (0.5 * ma ** 2 + g / (g - 1) * 0.5 * beta) * ma + bz1 * ma * bz1
            down = ((0.5 * x * (d["v_x"] ** 2 + d["v_z"] ** 2) + g / (g - 1) * d["p"]) * d["v_x"]
                    + d["B_z"] * (d["v_x"] * d["B_z"] - d["v_z"] * d["B_x"]))
            assert x > 1.0 and d["p"] > 0.0 and abs(down / up - 1.0) < 1e-9, (ma, beta, cos2)
        # the limits: the perpendicular shock, and along the field the gas-dynamic shock of the sonic Mach number
        assert T.oblique_compression(ma, beta, 0.0) == mhd.compression(ma, beta)
        assert abs(T.oblique_compression(ma, beta, 1e-9) / mhd.compression(ma, beta) - 1.0) < 1e-6
    ms2 = 16.0 / (0.5 * g * 0.64)
    assert abs(T.oblique_compression(4.0, 0.64, 1.0) - (g + 1) * ms2 / ((g - 1) * ms2 + 2.0)) < 1e-12
    assert T.oblique_compression(1.5, 0.1, 1.0) == 2.25 and T.oblique_downstream(1.5, 0.1, 1.0)["B_z"] > 0.0
    assert T.oblique_compression(0.5 * ma1, beta, 0.5) == 1.0  # slower than the fast speed: no shock
    # magnetopause: the doubled dipole field of the planner's moment balances the ram pressure at R_mp
    r, v, mi = 20.0, 0.25, 16.0
    m = T.dipole_moment(r, v, mi)
    assert abs(0.5 * (2.0 * m / r ** 3) ** 2 / (mi * v * v) - 1.0) < 1e-14
    assert abs(T.standoff_distance(m, v, mi) / r - 1.0) < 1e-14


def test_blast_wave_scales():
    T = _example_theory("blast_wave")
    # Eq. 4 of Le et al. 2021 at M_A = 3: 1.09 / (11 / 9)^(1 / D)
    assert abs(T.debris_extent(1.0, 3.0, dim=2) - 0.98594) < 1e-5
    assert abs(T.debris_extent(1.0, 3.0) - 1.019475) < 1e-6
    # the part of a Gaussian shell inside a radius, against a quadrature of r^2 exp(-((r - r_s) / w)^2)
    r_s, w = 10.56, 3.2
    r = np.linspace(0.0, r_s + 12.0 * w, 400001)
    f = r * r * np.exp(-((r - r_s) / w) ** 2)
    for radius in (4.0, 8.0, r_s, 16.0):
        assert abs(T.shell_fraction_inside(r_s, w, radius) - f[r <= radius].sum() / f.sum()) < 1e-4
    assert T.shell_fraction_inside(r_s, w, 0.0) == 0.0 and abs(T.shell_fraction_inside(r_s, w, 1e3) - 1.0) < 1e-12


def test_profiles():
    from pypic.profiles import evaluate, validate

    y = np.array([0.0, 25.0, 50.0, 75.0, 99.0])
    f = evaluate(validate({"shape": "double_tanh", "centers": [25.0, 75.0], "width": 1.0}), y)
    assert np.allclose(f, [-1.0, 0.0, 1.0, 0.0, -1.0], atol=1e-9)
    assert np.allclose(evaluate({"shape": "step", "x0": 1.0, "x1": 2.0}, np.array([0.5, 1.0, 1.99, 2.0])),
                       [0, 1, 1, 0])


def test_density_perturbation_inverse_converges_where_the_density_nearly_vanishes():
    """Review of PR #2: with alpha = 0.99, Newton's method on the cumulative density left 185 of
    these 8192 positions unconverged (and the virtual move then needed ~1e16 rounds)."""
    from pypic.profiles import invert_cos_density

    L, alpha = 8.0, 0.99
    k = 2 * math.pi / L
    x0 = np.random.default_rng(12345).random(8192) * L
    x, res = invert_cos_density(x0, alpha, k, 0.0)
    assert res < 1e-12 * L
    assert np.max(np.abs(x + alpha / k * np.sin(k * x) - x0)) < 1e-12 * L
    assert np.all(np.diff(x[np.argsort(x0)]) >= 0)  # the map keeps the particle order


def test_drift_frame_orbits_cover_the_requested_lab_time():
    """Review of PR #2: exb (|E| -> |B|) and e_gt_b (|E| -> |B|+) mapped a fixed frame-time grid to the
    lab and stopped short of t_end. Compared with a direct integration of the Lorentz force."""
    from scipy.integrate import solve_ivp

    T = _example_theory("particle_orbits")

    def lorentz(E, B, t):
        def f(_, s):
            u = s[3:]
            g = math.sqrt(1 + u @ u)
            v = u / g
            return np.concatenate([v, [v[1] * B, E - v[0] * B, 0.0]])  # q = m = 1
        return solve_ivp(f, (0, t[-1]), np.zeros(6), t_eval=t, rtol=1e-11, atol=1e-11).y

    t = np.linspace(0, 10 * math.pi, 501)
    for fn, E, B in ((T.exb, 0.99, 1.0), (T.e_gt_b, 1.01, 1.0)):
        x, u = fn(t, (0.0, 0.0, 0.0), 1.0, 1.0, E, B)
        ref = lorentz(E, B, t)
        scale = np.max(np.abs(ref[:2]))
        assert np.max(np.abs(x[:2] - ref[:2])) < 1e-6 * scale, fn.__name__
        assert np.max(np.abs(u[:2] - ref[3:5])) < 1e-6 * np.max(np.abs(ref[3:5])), fn.__name__
