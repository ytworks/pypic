"""Configuration of frozen-in initial states (PIC-D24): shaped perturbations and the vector
potential. Pure Python: no CuPy needed."""
import pytest

from pypic.config import Beam, SimConfig

COS = {"shape": "cos", "mode": 1, "length": 8.0}


def test_shaped_perturbation_is_accepted_and_validated():
    b = Beam(1.0, perturbation=[{"shape": {"y": COS, "x": {"shape": "ramp", "x0": 1.0, "x1": 2.0}},
                                  "velocity": 0.1, "velocity_component": "x"},
                                 {"shape": {"y": COS}, "displacement": 0.01, "displacement_component": "z"}])
    assert b.perturbation[0]["shape"]["x"]["shape"] == "ramp"


@pytest.mark.parametrize("pert, error", [
    ({"shape": {"y": COS}, "velocity": 0.1}, "velocity_component"),
    ({"shape": {"y": COS}, "displacement": 0.1}, "displacement_component"),
    ({"shape": {"y": COS}, "velocity": 0.1, "velocity_component": "x", "axis": "y"}, "do not go with 'shape'"),
    ({"shape": {"y": COS}, "density": 0.1}, "do not go with 'shape'"),
    ({"shape": {"w": COS}, "velocity": 0.1, "velocity_component": "x"}, "shape axis"),
    ({"shape": {}, "velocity": 0.1, "velocity_component": "x"}, "one or more"),
    ({"shape": {"y": COS}, "velocity_component": "x"}, "velocity or displacement"),
    ({"axis": "x", "displacement": 0.1, "displacement_component": "y"}, "belongs to a shaped"),
])
def test_shaped_perturbation_refusals(pert, error):
    with pytest.raises((KeyError, ValueError), match=error):
        Beam(1.0, perturbation=[pert])


def test_potential_is_validated():
    cfg = SimConfig(initial_fields=[{"kind": "potential", "component": "Az", "shape": {"y": COS}, "scale": 2}])
    assert cfg.initial_fields[0]["scale"] == 2.0
    with pytest.raises(ValueError, match="Ax, Ay or Az"):
        SimConfig(initial_fields=[{"kind": "potential", "component": "Bz", "shape": {"y": COS}}])
    with pytest.raises(KeyError, match="unknown"):
        SimConfig(initial_fields=[{"kind": "potential", "component": "Az", "shape": {"y": COS}, "axis": "y"}])
    with pytest.raises(ValueError, match="one or more"):
        SimConfig(initial_fields=[{"kind": "potential", "component": "Az"}])
