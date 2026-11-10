"""CPU tests of the absorbing body (PIC-D14): the combinations that are rejected and the nodes
that the Gauss check leaves out."""
from __future__ import annotations

import numpy as np
import pytest

from pypic import SimConfig
from pypic.body import check_config, gauss_mask, gauss_radius, relax_edges
from pypic.config import Body, Dipole, ExternalFields

OPEN_X = (("open", "open"), ("periodic", "periodic"), ("reflecting", "reflecting"))


def config(body, boundaries=OPEN_X, dipoles=(), **kw):
    # a box of 24 x 12 x 12 at dx = 0.5
    return SimConfig(n_cells=(48, 24, 24), dx=(0.5, 0.5, 0.5), dt=0.2, boundaries=boundaries, body=body,
                     external_fields=ExternalFields(dipoles=list(dipoles)), **kw)


def test_supported_and_rejected_combinations():
    sphere = Body(center=(12.0, 6.0, 6.0), radius=2.0)
    dipole = Dipole(moment=(0.0, 0.0, 1.0), center=(12.5, 6.0, 6.0), core=1.0)
    assert not check_config(config(None))
    assert check_config(config(sphere)) and check_config(config(sphere), reach=2)
    assert check_config(config(sphere, dipoles=[dipole]))
    # the nodes that hold the charge: 2 + sqrt(3)/2 from the centre, 2 + 3 sqrt(3)/2 with a reach of 2
    assert check_config(config(Body(center=(3.4, 6.0, 6.0), radius=2.0)))
    assert check_config(config(Body(center=(12.0, 2.9, 6.0), radius=2.0)))
    bad = [
        (config(sphere, solver="implicit"), 0),
        (config(Body(center=(3.3, 6.0, 6.0), radius=2.0)), 0),  # a cell from the open face x-
        (config(Body(center=(20.7, 6.0, 6.0), radius=2.0)), 0),  # and from x+
        (config(Body(center=(6.0, 6.0, 6.0), radius=2.0)), 2),  # clear without a filter, not with one
        (config(Body(center=(12.0, 2.8, 6.0), radius=2.0)), 0),  # through the periodic face y-
        (config(Body(center=(12.0, 6.0, 9.2), radius=2.0)), 0),  # through the wall z+
        (config(None, dipoles=[dipole]), 0),  # a dipole without a body
        (config(Body(center=(12.0, 6.0, 6.0), radius=1.4), dipoles=[dipole]), 0),  # the core sticks out
    ]
    for cfg, reach in bad:
        with pytest.raises(ValueError):
            check_config(cfg, reach)
    for kw in ({"radius": 0.0}, {"radius": 1.0, "shape": "cube"}, {"radius": 1.0, "center": (1.0, 2.0)},
               {"radius": 1.0, "relax": -0.1}):
        with pytest.raises(ValueError):
            Body(**{"center": (1.0, 1.0, 1.0), **kw})


def test_body_survives_the_round_trip_of_the_configuration():
    """run.json and the planner overrides rebuild the configuration from its dict."""
    cfg = config(Body(center=(12.0, 6.0, 6.0), radius=2.0))
    assert SimConfig.from_dict(cfg.to_dict()).body == cfg.body
    assert SimConfig.from_dict(config(None).to_dict()).body is None


@pytest.mark.parametrize("reach", [0, 2])
def test_gauss_mask_covers_every_node_an_absorbed_particle_charges(reach):
    """A particle inside the sphere leaves its charge on the nodes of its cell, and the filter
    carries it ``reach`` nodes further along every axis: none of those nodes is checked, while the
    mask stays a ball (it leaves out fewer nodes than the cube around it holds)."""
    dx, body = (0.5, 0.4, 0.3), Body(center=(6.1, 5.3, 4.7), radius=1.7)
    n = (30, 32, 36)
    keep = gauss_mask(np, [np.arange(m) * h for m, h in zip(n, dx)], body, dx, reach)
    assert keep.shape == n
    rng = np.random.default_rng(7)
    v = rng.normal(size=(4000, 3))
    x = np.asarray(body.center) + body.radius * v / np.linalg.norm(v, axis=1)[:, None] * rng.random((4000, 1)) ** 0.2
    cell = np.floor(x / np.asarray(dx)).astype(int)
    for o in np.ndindex(*(2 + 2 * reach,) * 3):
        node = cell + np.asarray(o) - reach
        assert not keep[node[:, 0], node[:, 1], node[:, 2]].any()
    rg = gauss_radius(body, dx, reach)
    assert np.count_nonzero(~keep) < 0.6 * np.prod([2.0 * rg / h for h in dx])  # a ball is 0.52 of its cube


def test_relaxed_edges_move_charge_only_between_nodes_inside_the_body():
    """``body.relax`` scales E on the edges it selects. Any field on those edges has zero divergence
    at every node outside the body and zero total divergence, so the relaxation cannot reach a node
    of the Gauss check and cannot change the charge of the body; and it selects every edge between
    two nodes inside."""
    dx, body = (0.5, 0.4, 0.3), Body(center=(6.1, 5.3, 4.7), radius=1.7, relax=0.05)
    n = (30, 32, 36)
    pos = [np.arange(m) * h for m, h in zip(n, dx)]
    edges = relax_edges(np, pos, body, dx)
    r2 = sum(((p - c) ** 2).reshape([-1 if b == a else 1 for b in range(3)])
             for a, (p, c) in enumerate(zip(pos, body.center)))
    node = r2 < body.radius**2
    rng = np.random.default_rng(3)
    div = np.zeros(n)
    for a in range(3):
        assert edges[a].sum() == (node & np.roll(node, -1, axis=a)).sum() > 100
        e = rng.normal(size=n) * edges[a]
        div += (e - np.roll(e, 1, axis=a)) / dx[a]
    assert np.abs(div[node]).max() > 1.0 and not div[~node].any()
    assert abs(div.sum()) < 1e-10
