"""CPU tests of the injection at open faces (PIC-D13): the combinations that are rejected and the
cell layers the ranks fill."""
from __future__ import annotations

import itertools

import numpy as np
import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import FACES, Injection
from pypic.injection import check_config, face_slabs
from pypic.parallel.decomp import Decomposition

OPEN = (("open", "open"),) * 3


def config(faces, boundaries=OPEN, electrons=None, ions=None, solver="explicit", **kw):
    electrons = electrons or [Beam(1.0, (0.2, 0.0, 0.0), 0.01)]
    ions = ions or [Beam(1.0, (0.2, 0.0, 0.0), 0.001)]
    return SimConfig(n_cells=(8, 8, 8), boundaries=boundaries, solver=solver, dt=0.05,
                     species=[SpeciesConfig("e", -1, 1, electrons), SpeciesConfig("i", 1, 16, ions),
                              SpeciesConfig("t", -1, 1, test=True, particles=[[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]])],
                     injection=Injection(faces=faces, **kw))


def test_supported_and_rejected_combinations():
    assert not check_config(config({}))
    assert check_config(config({"x-": "flux", "y+": "flux"}))
    assert check_config(config({"x-": "paired", "y-": "flux"}, neutralizer="e"))
    side = (("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting"))
    out = [Beam(1.0, (-0.2, 0.0, 0.0), 0.01)]
    bad = [
        config({"x-": "flux"}, solver="implicit"),
        config({"y-": "flux"}, boundaries=side),  # not an open face
        config({"x-": "flux"}, boundaries=side),  # meets a reflecting face
        config({"x-": "flux"}, species=("e", "nobody")),
        config({"x-": "flux"}, species=("e", "t")),  # a test species
        config({"x-": "flux"}, ions=[Beam(1.0, profile={"axis": "x", "shape": "step", "x0": 0.0, "x1": 1.0})]),
        config({"x-": "flux"}, neutralizer="e"),  # no paired face uses it
        config({"x-": "paired"}),  # no neutralizer
        config({"x-": "paired"}, species=("e",), neutralizer="e"),  # nothing to pair with
        config({"x+": "paired"}, neutralizer="e"),  # the flow leaves through x+
        config({"x-": "paired"}, neutralizer="e", electrons=out, ions=out),
    ]
    for cfg in bad:
        with pytest.raises(ValueError):
            check_config(cfg)
    with pytest.raises(ValueError):
        Injection(faces={"w-": "flux"})
    with pytest.raises(ValueError):
        Injection(faces={"x-": "inflow"})


@pytest.mark.parametrize("dims", [(1, 1, 1), (2, 2, 1), (1, 3, 4)])
def test_every_cell_of_the_layer_is_filled_once(dims):
    """Over all ranks, a cell outside the box is filled exactly once if every face it lies outside
    of injects, and never otherwise; an edge or corner cell belongs to the first of its faces."""
    n = (4, 6, 8)
    faces = {"x-": "paired", "y-": "flux", "y+": "flux", "z-": "flux"}  # x+ and z+ only remove
    count = np.zeros(tuple(m + 2 for m in n), dtype=int)
    owner = np.full(count.shape, -1)
    for rank in range(dims[0] * dims[1] * dims[2]):
        axes = Decomposition.build(n, dims, rank, OPEN, 4).axes
        for face, g0, nl in face_slabs(axes, n, faces):
            block = tuple(slice(g + 1, g + 1 + m) for g, m in zip(g0, nl))
            count[block] += 1
            owner[block] = face
    for cell in itertools.product(*(range(-1, m + 1) for m in n)):
        outside = [2 * a + (c == n[a]) for a, c in enumerate(cell) if c in (-1, n[a])]
        filled = bool(outside) and all(FACES[f] in faces for f in outside)
        at = tuple(c + 1 for c in cell)
        assert count[at] == int(filled), cell
        assert owner[at] == (outside[0] if filled else -1), cell
