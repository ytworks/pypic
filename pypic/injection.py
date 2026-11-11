"""Plasma entering through open faces (PIC-D13): the rules, the cell layers and the inflow tables.

Every step, after the particle push, the layer of cells just outside an injecting face is filled
with the species' own distribution; a particle of the layer that ends its step inside the box
enters, the others are dropped. One cell is enough: a particle moves less than a cell per step.
The random numbers are keyed by (seed, step, global cell, slot), so the entering particles do not
depend on the decomposition.

A ``flux`` face launches every species on its own. A ``paired`` face is the inflow face of a flow:
the neutraliser species is returned at the face with the inflow distribution instead of leaving,
and every other species (a carrier) enters together with a neutraliser particle born at the same
point, so the face collects no charge. A carrier that leaves through a paired face is paid for by
the next carrier entering in the same face column, which enters without a partner. The paired faces
of a run lie on one axis: a neutraliser returned near an edge where two of them met could cross
the other one again and again within its step.

Pure Python: the kernels are ``inject_face``, ``inject_partners`` and ``push_deposit`` in
``kernels/particles.cu``.
"""
from __future__ import annotations

import math

from .config import FACES


def _inward_velocity(beam, face):
    """Drift velocity of a beam along the inward normal of a face."""
    u = beam.drift
    v = u[face >> 1] / math.sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
    return -v if face & 1 else v


def inflow_table(beams, face):
    """Rows of the kernel's flux sampler for the beams of one species at one face:
    ``[cumulative density * vcap, theta, udx, udy, udz, vcap]``. ``vcap`` bounds the inward normal
    velocity of the beam (drift + 8 thermal velocities, at most c): the sampler picks a beam with
    probability proportional to ``density * vcap`` and accepts a velocity with probability
    ``v_n / vcap``, which gives the flux-weighted distribution ``v_n f(u)`` of the species."""
    rows, total = [], 0.0
    for b in beams:
        vcap = min(1.0, max(_inward_velocity(b, face), 0.0) + 8.0 * math.sqrt(b.theta))
        total += b.density * vcap
        rows.append([total, b.theta, *b.drift, vcap])
    return rows


def injected_species(cfg):
    """Names of the species that enter through the injecting faces."""
    if cfg.injection.species is not None:
        return list(cfg.injection.species)
    return [s.name for s in cfg.species if s.beams and not s.test]


def ids_per_step(cfg) -> int:
    """Particle ids (PIC-D8) that one step reserves for the entering particles: per stream (a beam
    of a species) one sub-block of the extended grid (cells -1..N) and one for the partners born
    with its particles at a paired face. A run of n steps needs loaded + n * ids_per_step < ID_LIMIT."""
    n_streams = sum(len(s.beams) for s in cfg.species)
    return 2 * n_streams * math.prod(n + 2 for n in cfg.n_cells) * cfg.ppc


def paired_faces(cfg) -> int:
    """Bit ``2 * axis + side`` per paired face (the kernel constant ``PAIRED_FACES``)."""
    return sum(1 << FACES.index(k) for k, v in cfg.injection.faces.items() if v == "paired")


def check_config(cfg) -> bool:
    """True if ``cfg`` injects plasma; raises on the combinations that are not supported."""
    inj = cfg.injection
    if not inj.active:
        return False
    if cfg.solver != "explicit":
        raise ValueError("injection is for the explicit solver only (PIC-D13)")
    for name in inj.faces:
        f = FACES.index(name)
        if cfg.boundaries[f >> 1][f & 1] != "open":
            raise ValueError(f"injection face {name!r} must be an open boundary")
        for a in range(3):
            if a != f >> 1 and "reflecting" in cfg.boundaries[a]:
                raise ValueError(f"injection face {name!r} must not meet a reflecting face (axis {'xyz'[a]})")
    by_name = {s.name: s for s in cfg.species}
    names = injected_species(cfg)
    if not names or len(set(names)) != len(names):
        raise ValueError("injection needs at least one species, each named once")
    for n in names:
        if n not in by_name:
            raise ValueError(f"injection: unknown species {n!r}")
        s = by_name[n]
        if s.test or not s.beams:
            raise ValueError(f"injection: species {n!r} must be a plasma species with beams")
        if any(b.profile is not None or b.perturbation for b in s.beams):
            raise ValueError(f"injection: the beams of species {n!r} must not have a profile or a perturbation")
    paired = [FACES.index(k) for k, v in inj.faces.items() if v == "paired"]
    if not paired:
        if inj.neutralizer is not None:
            raise ValueError("injection.neutralizer is only used by a paired face")
        return True
    if len({f >> 1 for f in paired}) > 1:
        raise ValueError("paired injection faces must lie on one axis: a neutralizer returned at one of two "
                         "paired faces that meet at an edge can leave through the other within its step")
    if inj.neutralizer not in names:
        raise ValueError("a paired injection face needs injection.neutralizer, one of the injected species")
    qn = by_name[inj.neutralizer].charge
    carriers = [n for n in names if n != inj.neutralizer]
    if not carriers or any(by_name[n].charge * qn >= 0.0 for n in carriers):
        raise ValueError("a paired injection face needs at least one other injected species, each with a charge "
                         "of the opposite sign to the neutralizer")
    for f in paired:
        if not inflow_table(by_name[inj.neutralizer].beams, f)[-1][0] > 0.0:
            raise ValueError(f"paired face {FACES[f]!r}: the neutralizer has no inward flux")
        for n in names:
            if sum(b.density * _inward_velocity(b, f) for b in by_name[n].beams) < 0.0:
                raise ValueError(f"paired face {FACES[f]!r} is an inflow face: species {n!r} drifts out of it")
    return True


def face_slabs(axes, n_cells, faces):
    """Blocks of cells outside the box that one rank fills: ``[(face, (gx0, gy0, gz0), (nx, ny, nz))]``
    in global cell indices (-1 or N along the normal of the face).

    ``axes``: the rank's ``AxisLayout`` per axis; ``faces``: the injecting face names. A rank holds
    the layer of a face over its own physical cells. A cell outside two or three faces (an edge or
    a corner of the box) is filled only if all of those faces inject, by the rank at that edge, and
    belongs to the first of its faces in the order x, y, z: over all ranks every cell of the layer
    is filled exactly once."""
    out = []
    for a in range(3):
        for side in (0, 1):
            if FACES[2 * a + side] not in faces or (axes[a].hi if side else axes[a].lo) is not None:
                continue
            g0, nl = [0, 0, 0], [1, 1, 1]
            g0[a] = n_cells[a] if side else -1
            for b in range(3):
                if b == a:
                    continue
                g0[b], nl[b] = axes[b].phys_lo, axes[b].n_phys
                if b > a and FACES[2 * b] in faces and axes[b].lo is None:
                    g0[b] -= 1
                    nl[b] += 1
                if b > a and FACES[2 * b + 1] in faces and axes[b].hi is None:
                    nl[b] += 1
            out.append((2 * a + side, tuple(g0), tuple(nl)))
    return out
