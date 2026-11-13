"""examples/alfven_inflow: the cold two-fluid parallel waves, the jump conditions of the piston shock,
the loaded wave and the planner. Pure Python (NumPy, SciPy): no CuPy needed."""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "examples" / "alfven_inflow"


def _module(name):
    """examples/alfven_inflow/<name>.py; its ``import theory`` must not meet another example's module."""
    spec = importlib.util.spec_from_file_location(f"alfven_inflow_{name}", FOLDER / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    other = sys.modules.pop("theory", None)
    sys.modules[spec.name] = mod  # the dataclasses of the planner look their module up
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("theory", None)
        if other is not None:
            sys.modules["theory"] = other
        if str(FOLDER) in sys.path:
            sys.path.remove(str(FOLDER))
    return mod


T = _module("theory")
B = 0.1


@pytest.mark.parametrize("mass_ratio, kdi", [(16, 0.785), (16, 3.0), (1836, 0.1)])
def test_the_roots_are_the_r_and_l_waves_of_stix(mass_ratio, kdi):
    """The roots of the 3x3 cold response (sigma E = i (w - k^2 / w) E on a circular vector) satisfy
    the textbook relations k^2 / w^2 = 1 - sum w_ps^2 / (w (w +/- Omega_s)); R above k v_A, L below."""
    sp = T.two_fluid(mass_ratio)
    k = kdi / T.ion_skin_depth(sp)
    roots = T.parallel_roots(k, B, sp)
    for branch in ("R", "L"):
        w = roots[branch][0]
        assert abs(T.stix(w, k, B, sp, branch)) < 1e-8 * k * k / (w * w), branch
    assert roots["L"][0] < k * T.alfven_speed(B, sp) < roots["R"][0]


@pytest.mark.parametrize("mass_ratio", [16, 1836])
def test_low_k_limit_and_hall_splitting(mass_ratio):
    """For k d_i -> 0 both waves tend to w = k v_A / sqrt(1 + v_A^2), and the Hall term splits them to
    first order: (w_R - w_L) / (k v_A) = k d_i (m_i - m_e) / sqrt(m_i (m_i + m_e)) / (1 + v_A^2)^2."""
    sp = T.two_fluid(mass_ratio)
    va, di = T.alfven_speed(B, sp), T.ion_skin_depth(sp)
    kdi = 0.01
    k = kdi / di
    r = T.parallel_roots(k, B, sp)
    wr, wl = r["R"][0] / (k * va), r["L"][0] / (k * va)
    assert 0.5 * (wr + wl) * math.sqrt(1 + va * va) == pytest.approx(1.0, abs=1e-4)
    split = kdi * (mass_ratio - 1) / math.sqrt(mass_ratio * (mass_ratio + 1)) / (1 + va * va) ** 2
    assert wr - wl == pytest.approx(split, rel=1e-3)


def test_the_table_of_the_wave_satisfies_maxwell_with_its_current():
    """The box-frame table (E = E_rest - v_d x B, the cold velocities) passes the engine's Maxwell and
    Ampère checks at a tolerance of 5e-3 (the displacement current of the motional field, omega v0 dB,
    is a relative 2.4e-3 that the rest-frame current does not carry); |dB| = 0.2 B0, circular."""
    from pypic import Beam, SimConfig, SpeciesConfig
    from pypic.config import Injection
    from pypic.injection import check_fluctuation

    sp = T.two_fluid(16)
    k = 2 * math.pi / 32
    v0 = 0.05
    table = T.fluctuation(T.parallel_wave(k, B, sp, "R", 0.2), v0)
    table["tolerance"] = 5e-3
    u0 = v0 / math.sqrt(1 - v0 * v0)
    cfg = SimConfig(n_cells=(16, 8, 2), species=[SpeciesConfig(n, q, m, [Beam(d, (-u0, 0, 0), 1e-3)])
                                                 for n, q, m, d in sp],
                    boundaries=(("open", "open"), ("periodic", "periodic"), ("periodic", "periodic")),
                    B0=(0, B, 0), E0=(0, 0, v0 * B), injection=Injection(faces={"x+": "flux"}, fluctuation=table))
    out = check_fluctuation(cfg)
    assert len(out["modes"]) == 1
    m = table["modes"][0]
    # |dB| = 0.2 B0; circular: equal amplitudes of B_x and B_z
    bx = math.hypot(m["B"]["cos"][0], m["B"]["sin"][0])
    bz = math.hypot(m["B"]["cos"][2], m["B"]["sin"][2])
    assert math.hypot(bx, bz) == pytest.approx(0.2 * B, rel=1e-12) and bx == pytest.approx(bz, rel=1e-12)


def test_the_potential_of_the_envelope_gives_the_wave_on_the_yee_grid():
    """The vector potential of the envelope load (theory.preload) differenced as the engine does,
    (A_z[j + 1] - A_z[j]) / dy for B_x and -(A_x[j + 1] - A_x[j]) / dy for B_z, gives h(x) times the
    wave's B on the staggered positions; h vanishes on the wall and the face and holds no
    displacement within one cell of them."""
    from pypic.profiles import evaluate

    sp = T.two_fluid(16)
    ly, dy, lx = 32.0, 0.25, 96.0
    table = T.fluctuation(T.parallel_wave(2 * math.pi / ly, B, sp, "R", 0.2), 0.07)
    h = T.envelope(lx, (4.0, 12.0), 8.0)
    fields, pert = T.preload(table, ly, dy, h)
    m = table["modes"][0]
    y = np.arange(128) * dy
    x = np.linspace(0.0, lx, 385)
    hx = evaluate(h, x)
    assert hx[0] == 0.0 and hx[-1] == 0.0 and np.all(hx[x > lx - 1.0] == 0.0) and np.all(hx[x < 4.0] == 0.0)
    assert np.max(np.abs(hx[(x > 12.0) & (x < lx - 9.0)] - 1.0)) == 0.0
    for f in fields:
        assert f["kind"] == "potential"
        c = f["component"]
        a = f["scale"] * evaluate(f["shape"]["y"], y)
        db = (np.roll(a, -1) - a) / dy * (1.0 if c == "Az" else -1.0)
        comp = 0 if c == "Az" else 2
        th = m["k"][1] * (y + 0.5 * dy)  # B_x and B_z sit half a cell up along y
        want = m["B"]["cos"][comp] * np.cos(th) + m["B"]["sin"][comp] * np.sin(th)
        assert np.max(np.abs(db - want)) < 1e-12 * B, c
        assert f["shape"]["x"] == h
    assert {p["displacement_component"] for p in pert["e"] if "displacement" in p} == {"x", "z"}


@pytest.mark.parametrize("v0", [0.03, 0.0728, 0.5])
def test_piston_shock_conserves_mass_momentum_and_energy(v0):
    """The piston shock in its own frame: n U, rho U^2 + p + B^2 / 2 (used for T2) and the energy flux
    rho U^3 / 2 + gamma p U / (gamma - 1) + B^2 U (with the Poynting flux of E = U B) are the same on
    both sides; V_sh (r - 1) = v0."""
    va, beta, mass = 0.0243, 0.4, 17.0
    s = T.piston_shock(v0, va, beta, mass=mass)
    r, u1 = s["compression"], v0 + s["v_shock"]
    u2 = u1 / r
    assert s["v_shock"] * (r - 1) == pytest.approx(v0, rel=1e-12)
    rho1 = mass
    b1 = va * math.sqrt(rho1)
    p1, p2 = 0.5 * beta * b1 * b1, r * s["temperature"]
    g = 5.0 / 3.0

    def flux(rho, u, p, b):
        return rho * u, rho * u * u + p + 0.5 * b * b, 0.5 * rho * u ** 3 + g / (g - 1) * p * u + b * b * u

    for a, b in zip(flux(rho1, u1, p1, b1), flux(rho1 * r, u2, p2, r * b1)):
        assert a == pytest.approx(b, rel=1e-9)


def test_piston_compression_tends_to_four():
    va, beta = 0.0243, 0.4
    assert T.piston_shock(100 * va, va, beta)["compression"] == pytest.approx(4.0, abs=0.01)
    assert T.piston_shock(2 * va, va, beta)["compression"] < T.piston_shock(4 * va, va, beta)["compression"]


def test_planner_configurations_are_accepted_and_the_wall_needs_the_envelope():
    import yaml

    from pypic import injection

    P = _module("plan")
    for name, wall in (("gate", False), ("piston", True)):
        with open(FOLDER / "configs" / f"{name}.yaml") as f:
            d = yaml.safe_load(f)
        d.pop("checkpoint", None)
        cfg, report, info = P.plan_alfven_inflow(P.params_from_dict(d), n_gpus=2)
        assert injection.check_config(cfg) and injection.check_fluctuation(cfg)
        assert info["wall"] == wall and (cfg.boundaries[0][0] == "reflecting") == wall
        if wall:
            assert info["shock"]["mach_alfven"] > 3.0
            null, _, _ = P.plan_alfven_inflow(P.params_from_dict({**d, "wave": False}), n_gpus=2)
            assert not null.initial_fields and not any(null.injection.fluctuation["modes"][0]["B"]["cos"])
            with pytest.raises(ValueError, match="envelope"):
                P.plan_alfven_inflow(P.params_from_dict({**d, "load": "uniform"}))
        else:
            with pytest.raises(ValueError, match="piston"):
                P.plan_alfven_inflow(P.params_from_dict({**d, "load": "envelope"}))
