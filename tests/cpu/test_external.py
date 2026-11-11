"""CPU tests of the external static fields (PIC-D12): the dipole on the Yee faces and the
combinations that are rejected."""
from __future__ import annotations

import numpy as np
import pytest

from pypic import SimConfig
from pypic.config import ExternalFields
from pypic.external import check_config, dipole_B

DX = (0.5, 0.4, 0.3)
MOMENT, CENTER, CORE = (0.3, -0.2, 1.0), (6.1, 5.3, 4.7), 1.0


def _index(lo, n):
    return [np.arange(a, a + m, dtype=np.float64) for a, m in zip(lo, n)]


def test_dipole_has_no_discrete_divergence_and_does_not_depend_on_the_block():
    """div B on the cell centres vanishes to round-off (the field is a discrete curl), and a block
    cut out of the grid holds the same bits as the same faces of a larger block: the field of a
    rank does not depend on the decomposition."""
    lo, n = (-2, -2, -2), (28, 30, 34)
    B = dipole_B(np, _index(lo, n), DX, MOMENT, CENTER, CORE)
    div = (np.diff(B[0], axis=0)[:, :-1, :-1] / DX[0] + np.diff(B[1], axis=1)[:-1, :, :-1] / DX[1]
           + np.diff(B[2], axis=2)[:-1, :-1, :] / DX[2])
    assert np.max(np.abs(B[2])) > 1.0  # the block holds the core
    assert np.max(np.abs(div)) < 1e-12
    off, m = (9, 3, 14), (7, 20, 11)
    part = dipole_B(np, _index([a + o for a, o in zip(lo, off)], m), DX, MOMENT, CENTER, CORE)
    cut = tuple(slice(o, o + k) for o, k in zip(off, m))
    for c in range(3):
        assert np.array_equal(part[c], B[c][cut])


def test_dipole_matches_the_regularised_dipole_formula():
    """B on the Yee faces against B = [3 r (m.r) - m r^2 + 2 m a^2] / (r^2 + a^2)^(5/2) at the face
    centres: within 1 % of |B| beyond 6 core radii, and 2 m / a^3 at the centre of the dipole."""
    lo, n = (-12, -12, -12), (52, 56, 64)
    index = _index(lo, n)
    B = dipole_B(np, index, DX, MOMENT, CENTER, CORE)
    m = np.asarray(MOMENT)
    worst = 0.0
    for c in range(3):  # face of B_c: node along c, cell centre along the other two axes
        r = np.stack(np.meshgrid(*[(index[a] + (0.0 if a == c else 0.5)) * DX[a] - CENTER[a] for a in range(3)],
                                 indexing="ij"))
        r2 = np.sum(r * r, axis=0)
        mr = np.tensordot(m, r, axes=1)
        ref = (3.0 * r * mr - m[:, None, None, None] * (r2 - 2.0 * CORE**2)) / (r2 + CORE**2) ** 2.5
        far = r2 > (6.0 * CORE) ** 2
        worst = max(worst, np.max(np.abs(B[c] - ref[c])[far] / np.sqrt(np.sum(ref * ref, axis=0))[far]))
    assert worst < 0.01
    # a dipole centred on a Bz face of a fine grid: the face value is the central field
    fine = dipole_B(np, _index((-2, -2, -2), (5, 5, 5)), (0.01,) * 3, (0.0, 0.0, 1.0), (0.005, 0.005, 0.0), 1.0)
    assert abs(fine[2][2, 2, 2] - 2.0) < 1e-3


def test_external_fields_reject_the_implicit_solver_and_B0():
    ext = ExternalFields(B=(0.0, 0.0, 0.1))
    assert check_config(SimConfig(external_fields=ext)) and not check_config(SimConfig(B0=(0.0, 0.0, 0.1)))
    with pytest.raises(ValueError, match="explicit"):
        check_config(SimConfig(solver="implicit", external_fields=ext))
    with pytest.raises(ValueError, match="B0"):
        check_config(SimConfig(B0=(0.0, 0.0, 0.1), external_fields=ext))
