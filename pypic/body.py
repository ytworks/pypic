"""Absorbing body inside the box (PIC-D14): the rules, the kernel constants and the Gauss mask.

A particle that ends a step inside the body is removed after its whole path was deposited, like a
particle that leaves through an open face. Its charge stays at the nodes its shape covered: the
body is an insulator that charges up to its floating potential, and Gauss's law does not hold at
those nodes for the charge density of the remaining particles. The current filter spreads that
charge by its reach, so the nodes left out of the check lie within ``gauss_radius`` of the centre.
With ``body.relax`` the electric field decays on the edges inside the body: its charge moves to the
surface and the total stays.

Pure Python / array code: the kernels are ``push_deposit``, ``mark_body`` and ``flag_dead`` in
``kernels/particles.cu``; ``xp`` is NumPy (CPU tests) or CuPy (the engine).
"""
from __future__ import annotations

import math


def gauss_radius(body, dx, reach: int) -> float:
    """Nodes closer than this to the centre can hold charge of absorbed particles. The shape of a
    particle at ``x`` covers the nodes with ``|node - x| < dx`` on every axis, and the current
    filter (``reach`` = passes + compensate) carries the charge that many nodes further."""
    return body.radius + (1 + reach) * math.sqrt(sum(h * h for h in dx))


def check_config(cfg, reach: int = 0) -> bool:
    """True if ``cfg`` has a body; raises on the combinations that are not supported.

    ``reach``: the reach of the current filter in nodes (0 without a filter). A dipole of
    ``external_fields`` needs a body that covers its core, with or without other use of the body:
    inside the core the field is as strong as 2 m / core^3 and is not that of a dipole."""
    body = cfg.body
    for d in cfg.external_fields.dipoles:
        if body is None or math.dist(d.center, body.center) + d.core > body.radius:
            raise ValueError("a dipole of external_fields needs a body that covers its core: "
                             "|dipole.center - body.center| + dipole.core <= body.radius (PIC-D14)")
    if body is None:
        return False
    if cfg.solver != "explicit":
        raise ValueError("body is for the explicit solver only (PIC-D14)")
    rg = gauss_radius(body, cfg.dx, reach)
    for a in range(3):
        length = cfg.n_cells[a] * cfg.dx[a]
        # at an open face the nodes of the Marder zone (PIC-D16) and those the filter fills from
        # the face (PIC-D15) must stay clear of the charge of the body
        clear = [(reach + 1) * cfg.dx[a] if b == "open" else 0.0 for b in cfg.boundaries[a]]
        if body.center[a] - rg < clear[0] or body.center[a] + rg > length - clear[1]:
            raise ValueError(
                f"body: axis {'xyz'[a]}: the nodes that hold its charge (radius {rg:.4g} around the centre) "
                f"must lie inside the box, {reach + 1} cell(s) from an open face")
    return True


def defines(body) -> dict:
    """Kernel constants of the body (``particle_defines``): none without a body, so the kernels of
    a run without a body are compiled from the same source as before."""
    if body is None:
        return {}
    return {"BODY_SPHERE": 1, "BODY_CX": body.center[0], "BODY_CY": body.center[1], "BODY_CZ": body.center[2],
            "BODY_R2": body.radius * body.radius}


def gauss_mask(xp, positions, body, dx, reach: int):
    """True at the nodes where Gauss's law is checked: those at least ``gauss_radius`` from the
    centre. ``positions``: per axis the node coordinates of a block (1-D arrays); the result has
    their shape product."""
    r2 = 0.0
    for a in range(3):
        d = (xp.asarray(positions[a], dtype=float) - body.center[a]).reshape([-1 if b == a else 1 for b in range(3)])
        r2 = r2 + d * d
    return r2 >= gauss_radius(body, dx, reach) ** 2


def relax_edges(xp, positions, body, dx):
    """Per field component, True at the edges on which ``body.relax`` acts: those whose two nodes
    are inside the body. Such an edge changes div E only at nodes inside the body, which the Gauss
    check leaves out, and takes from one of them what it gives to the other: the charge of the body
    is moved, not changed. ``positions`` as in ``gauss_mask``; component ``a`` of a block lies
    between the node and its neighbour ``dx[a]`` up along ``a``."""
    def inside(shift):
        r2 = 0.0
        for a in range(3):
            d = xp.asarray(positions[a], dtype=float) + shift[a] - body.center[a]
            r2 = r2 + (d * d).reshape([-1 if b == a else 1 for b in range(3)])
        return r2 < body.radius ** 2

    node = inside((0.0, 0.0, 0.0))
    return [node & inside(tuple(dx[a] if b == a else 0.0 for b in range(3))) for a in range(3)]
