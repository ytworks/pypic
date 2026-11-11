"""CPU tests of the set-up layer: mixtures, linear theory and the parameter planner."""
from __future__ import annotations

import math

import pytest

from pypic import Beam, SpeciesConfig


# ------------------------------------------------------------- setup / theory
def test_compositions_neutral():
    from pypic.setup.species import MIXTURES, MODES, Composition, TwoStreamBeams, build_species

    for mix in MIXTURES:
        for mode in MODES:
            sp = build_species(Composition(mix, 400.0, Z=2, positron_fraction=0.3), TwoStreamBeams(mode=mode))
            net = sum(s.charge * b.density for s in sp for b in s.beams)
            assert abs(net) < 1e-12
            ne = sum(b.density for s in sp if s.name == "electron" for b in s.beams)
            assert abs(ne - 1.0) < 1e-12


def test_symmetric_two_stream_growth_rate():
    from pypic.setup.species import Composition, TwoStreamBeams, build_species
    from pypic.theory.streaming import fastest_mode

    u = 0.05
    sp = build_species(Composition("e-i", 1e8), TwoStreamBeams("electron", u_drift=u, T_e=0.0))
    r = fastest_mode(sp, warm=False)
    g = math.sqrt(1 + u * u)
    assert abs(r.gamma_max - 1 / (2 * math.sqrt(2)) / g ** 1.5) < 1e-3
    assert abs(r.k_max * u / g - math.sqrt(3 / 8) / g ** 1.5) < 2e-3


def test_nested_overrides_are_converted():
    from pypic.config import ImplicitConfig as IC
    from pypic.setup.planner import TwoStreamParams, plan_two_stream

    cfg, _, _ = plan_two_stream(TwoStreamParams(solver="implicit", overrides={
        "implicit": {"picard_tol": 1e-10}, "dims": [1, 1, 1],
        "species": [{"name": "e", "charge": -1, "mass": 1,
                     "beams": [{"density": 1.0, "drift": [0.1, 0, 0], "theta": 1e-3}]}]}))
    assert isinstance(cfg.implicit, IC) and cfg.implicit.picard_tol == 1e-10
    assert cfg.implicit.picard_maxit == IC().picard_maxit
    assert isinstance(cfg.species[0], SpeciesConfig) and isinstance(cfg.species[0].beams[0], Beam)
    # typed objects from the Python API are accepted as well
    cfg, _, _ = plan_two_stream(TwoStreamParams(overrides={
        "implicit": IC(picard_tol=1e-9), "dims": (1, 1, 1),
        "species": [SpeciesConfig("e", -1, 1, [Beam(1.0, (0.1, 0, 0), 1e-3)])]}))
    assert cfg.implicit.picard_tol == 1e-9
    assert isinstance(cfg.species[0], SpeciesConfig) and cfg.species[0].beams[0].drift == (0.1, 0, 0)
    with pytest.raises(KeyError):
        plan_two_stream(TwoStreamParams(overrides={"implicit": {"no_such_field": 1}}))


def test_runner_output_overrides_are_typed(capsys):
    """Review of PR #2: output.* from the command line (or the YAML) was assigned unconverted, so
    diagnostics given as dicts failed in the scheduler."""
    import os

    from pypic.setup.planner import params_from_dict, plan_two_stream
    from pypic.setup.runner import run_example

    cfg = os.path.join(os.path.dirname(__file__), "..", "..", "examples", "two_stream", "configs", "ei.yaml")
    diags = "output.diagnostics=[{name: fd, kind: fields, fields: [E], every: 10}]"
    rc = run_example("two_stream", plan_two_stream, params_from_dict,
                     argv=["--config", cfg, "--plan-only", "--set", diags])
    assert rc == 0
    row = [ln.split() for ln in capsys.readouterr().out.splitlines() if ln.strip().startswith("fd ")]
    assert row and row[0][1:3] == ["fields", "10"]
    with pytest.raises(KeyError):
        run_example("two_stream", plan_two_stream, params_from_dict,
                    argv=["--config", cfg, "--plan-only", "--set", "output.no_such_field=1"])


def _example_planner(name):
    """examples/<name>/plan.py as a module of its own name; its ``import theory`` must not meet the
    theory module of another example."""
    import importlib.util
    import os
    import sys

    folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "examples", name))
    spec = importlib.util.spec_from_file_location(f"plan_{name}", os.path.join(folder, "plan.py"))
    mod = importlib.util.module_from_spec(spec)
    other = sys.modules.pop("theory", None)
    sys.modules[spec.name] = mod  # the dataclasses of the planner look their module up
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop("theory", None)
        if other is not None:
            sys.modules["theory"] = other
        if folder in sys.path:
            sys.path.remove(folder)
    return mod


def test_bow_shock_planner_gives_a_configuration_the_engine_accepts():
    import os

    import yaml

    from pypic import body, external, injection

    P = _example_planner("bow_shock")
    path = os.path.join(os.path.dirname(P.__file__), "configs", "gate.yaml")
    with open(path) as f:
        params = yaml.safe_load(f)
    # an oblique IMF: the wind must not feel a force, E + v x B = 0, whatever the direction of B
    for imf in ((0.0, 0.0, 1.0), (0.3, -1.0, 0.5)):
        cfg, report, info = P.plan_bow_shock(P.params_from_dict({**params, "imf": imf}), n_gpus=8)
        E, B = cfg.external_fields.E, cfg.external_fields.B
        for s in cfg.species:
            u = s.beams[0].drift
            v = [c / math.sqrt(1.0 + sum(c * c for c in u)) for c in u]
            force = [E[a] + v[(a + 1) % 3] * B[(a + 2) % 3] - v[(a + 2) % 3] * B[(a + 1) % 3] for a in range(3)]
            assert max(abs(c) for c in force) < 1e-15, (imf, s.name)
        # the checks that Simulation makes at construction (no GPU needed for them)
        reach = cfg.current_filter.passes + int(cfg.current_filter.compensate)
        assert external.check_config(cfg) and injection.check_config(cfg) and body.check_config(cfg, reach)
        assert math.prod(cfg.dims) == 8 and "bow-shock planner" in report
    assert abs(info["r_mp"] - 20.0) < 1e-12  # 5 d_i with m_i / m_e = 16
    # a body so small that the dipole field at its surface cannot be resolved in time is refused
    with pytest.raises(ValueError, match="Omega_ce dt"):
        P.plan_bow_shock(P.params_from_dict({**params, "body_radius": 0.4}))
    # output blocks so large that none lies inside the upstream region: its gate could not be measured
    with pytest.raises(ValueError, match="upstream region"):
        P.plan_bow_shock(P.params_from_dict({**params, "coarse": 0.5}))


def test_blast_wave_planner_gives_a_configuration_the_engine_accepts():
    import os

    import numpy as np
    import yaml

    from pypic import body, external, injection

    P = _example_planner("blast_wave")
    with open(os.path.join(os.path.dirname(P.__file__), "configs", "gate.yaml")) as f:
        params = yaml.safe_load(f)

    def plan(n_gpus=8, **over):
        return P.plan_blast_wave(P.params_from_dict({**params, **over}), n_gpus=n_gpus)

    cfg, report, info = plan()
    # the checks that Simulation makes at construction (no GPU needed for them)
    reach = cfg.current_filter.passes + int(cfg.current_filter.compensate)
    assert external.check_config(cfg) and injection.check_config(cfg) and body.check_config(cfg, reach)
    assert math.prod(cfg.dims) == 8 and "blast-wave planner" in report
    # the debris: as many ions as the ambient plasma holds inside R_m, electrons with the same shape
    # (neutral), every beam with the radial speed V_d = M_A v_A around the centre of the body
    sp = {s.name: s for s in cfg.species}
    r = np.linspace(0.0, 4.0 * info["r_m"], 400001)
    for name in ("de", "di"):
        number = 0.0
        for b in sp[name].beams:
            assert b.origin == tuple(info["center"]) and b.profile["axis"] == "r"
            assert b.drift_radial / math.sqrt(1.0 + b.drift_radial ** 2) == pytest.approx(info["v_d"], rel=1e-12)
            shape = np.exp(-((r - b.profile["center"]) / b.profile["width"]) ** 2)
            number += b.density * float(np.sum(4.0 * math.pi * r * r * shape) * (r[1] - r[0]))
        assert number == pytest.approx(4.0 / 3.0 * math.pi * info["r_m"] ** 3, rel=1e-6), name
    assert sorted(cfg.injection.species) == ["e", "i"]  # the debris does not enter through the faces
    # without the shell and the body: the null tests
    cfg0, _, _ = plan(debris=False, body=False)
    assert [s.name for s in cfg0.species] == ["e", "i"] and cfg0.body is None
    for match, n_gpus, over in (
            ("debris_beams", 8, {"debris_beams": 0}),
            ("omega_pe dt", 8, {"shell_width": 0.05}),  # a thin shell is dense
            ("overlaps the body", 8, {"body_radius": 1.2}),
            ("particle slots", 1, {"ppc": 32}),  # 2.7e9 slots of one species on one GPU
            ("half_width", 8, {"half_width": 1.2}),
            ("body radius", 8, {"dipole_radius": 0.1}),
            ("Omega_ce dt", 8, {"body_radius": 0.3, "dipole_radius": 0.5})):
        with pytest.raises(ValueError, match=match):
            plan(n_gpus, **over)
