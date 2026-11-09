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
