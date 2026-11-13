"""Inflow fluctuation tables (PIC-D26): the Maxwell checks of injection.check_fluctuation. Pure
Python: no CuPy needed."""
import copy

import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import Injection, Window
from pypic.injection import check_fluctuation

K, W, B0 = 0.5, 0.2, 1.0


def _mode():
    """A transverse mode along y in a plasma at rest: B = b z cos, E from omega B = k x E, the
    current of Ampere carried by the ions (q n = 1)."""
    b = 0.1
    E = [-W * b / K, 0.0, 0.0]  # k x E = (0, K, 0) x (Ex, 0, 0) = (0, 0, -K Ex) = omega B
    Bc = [0.0, 0.0, b]
    # J_sin = -k x B_cos - omega E_cos = -(K b, 0, 0) - omega E
    Js = [-K * b - W * E[0], 0.0, 0.0]
    return {"k": [0.0, K, 0.0], "omega": W, "phase": 0.0, "B": {"cos": Bc}, "E": {"cos": E},
            "u": {"i": {"sin": Js}}}


def _cfg(mode, **kw):
    bc = (("open", "open"), ("periodic", "periodic"), ("periodic", "periodic"))
    sp = [SpeciesConfig("e", -1, 1, [Beam(1.0, (-0.05, 0, 0), 1e-3)]),
          SpeciesConfig("i", 1, 16, [Beam(1.0, (-0.05, 0, 0), 1e-4)])]
    return SimConfig(n_cells=(16, 8, 2), boundaries=bc, species=sp, B0=(0, B0, 0),
                     injection=Injection(faces={"x+": "flux"}, fluctuation={"face": "x+", "modes": [mode]}), **kw)


def test_a_consistent_mode_is_accepted():
    out = check_fluctuation(_cfg(_mode()))
    assert out["face"] == 1 and len(out["rows"][1]) == 1
    assert out["du_n"][1] == pytest.approx(abs(_mode()["u"]["i"]["sin"][0]))


@pytest.mark.parametrize("change, error", [
    (lambda m: m.update(k=[0.1, K, 0.0]), "normal of the face"),
    (lambda m: m["B"]["cos"].__setitem__(1, 0.05), "k . B"),
    (lambda m: m.update(omega=2 * W), "omega B - k x E"),
    (lambda m: m["u"]["i"]["sin"].__setitem__(0, 0.0), "Ampere"),
    (lambda m: m["E"].update(cos=[m["E"]["cos"][0], 0.5, 0.0]), "charge density"),
    (lambda m: m["u"].update(p={"cos": [0, 0, 0]}), "injected species"),
])
def test_inconsistent_modes_are_refused(change, error):
    m = copy.deepcopy(_mode())
    change(m)
    with pytest.raises(ValueError, match=error):
        check_fluctuation(_cfg(m))


def test_unsupported_combinations_are_refused():
    with pytest.raises(ValueError, match="moving window"):
        check_fluctuation(_cfg(_mode(), window=Window()))
    cfg = _cfg(_mode())
    cfg.injection.fluctuation["face"] = "y+"
    with pytest.raises(ValueError, match="injecting faces"):
        check_fluctuation(cfg)
