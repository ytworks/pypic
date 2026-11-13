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

A face can carry a prescribed fluctuation of the inflow (PIC-D26, ``check_fluctuation``): a sum of
modes uniform along the normal of the face, each a field and a velocity per species, c cos(th) +
s sin(th) with th = k . r - omega t + phase. The entering particles take its velocity at their
birth point; its fields enter through a total-field / scattered-field surface on the face.

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


def inflow_table(beams, face, extra=0.0):
    """Rows of the kernel's flux sampler for the beams of one species at one face:
    ``[cumulative density * vcap, theta, udx, udy, udz, vcap]``. ``vcap`` bounds the inward normal
    velocity of the beam (drift + 8 thermal velocities + ``extra``, the largest normal velocity of an
    inflow fluctuation (PIC-D26), at most c): the sampler picks a beam with probability proportional
    to ``density * vcap`` and accepts a velocity with probability ``v_n / vcap``, which gives the
    flux-weighted distribution ``v_n f(u)`` of the species."""
    rows, total = [], 0.0
    for b in beams:
        vcap = min(1.0, max(_inward_velocity(b, face), 0.0) + 8.0 * math.sqrt(b.theta) + float(extra))
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
        # shaped perturbations (PIC-D24) of the loaded plasma go with an inflow fluctuation (PIC-D26), which
        # carries the same pattern in; any other perturbation would differ from what enters
        shaped_ok = bool(inj.fluctuation)
        if any(b.profile is not None or b.origin is not None
               or (b.perturbation and not (shaped_ok and all("shape" in p for p in b.perturbation)))
               for b in s.beams):
            raise ValueError(f"injection: the beams of species {n!r} must not have a profile, a perturbation "
                             "(a shaped one only with injection.fluctuation) or a radial drift")
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


MODE_KEYS = {"k", "omega", "phase", "B", "E", "u"}
FLUCT_KEYS = {"face", "modes", "tolerance"}


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    return math.sqrt(_dot(a, a))


def _pair(d, what):
    d = d or {}
    if not isinstance(d, dict) or set(d) - {"cos", "sin"}:
        raise ValueError(f"inflow fluctuation: {what} is {{'cos': [3 numbers], 'sin': [3 numbers]}}")
    out = []
    for key in ("cos", "sin"):
        v = [float(x) for x in d.get(key, (0.0, 0.0, 0.0))]
        if len(v) != 3:
            raise ValueError(f"inflow fluctuation: {what}.{key} needs three components")
        out.append(v)
    return out


def check_fluctuation(cfg):
    """The inflow fluctuation of PIC-D26, checked; None without one. Returns ``{"face": f, "modes":
    [...], "rows": {species index: [[kx, ky, kz, omega, phase, ucx, ucy, ucz, usx, usy, usz], ...]},
    "du_n": {species index: bound of the normal velocity}}``.

    Every mode must satisfy Maxwell's equations with the current of the particles' velocities:
    k . B = 0 and omega B = k x E, and the current curl B - dE/dt equal to sum_s q_s n_s du_s / gamma_s
    (both for the cos and the sin part), within ``tolerance`` (default 1e-2) of the largest term. The
    entering particles carry no density fluctuation, so the charge density k . E that the mode needs
    must be below ``tolerance`` times sum_s |q_s n_s|: a wave in a drifting plasma has the charge
    -v_drift . J_rest (a relativistic effect of order v_drift dv / c^2), which is left out. The modes
    are uniform along the normal of the face (k_n = 0), so the plasma carries the pattern in
    unchanged."""
    inj = cfg.injection
    fl = inj.fluctuation
    if not fl:
        return None
    if not isinstance(fl, dict) or set(fl) - FLUCT_KEYS:
        raise ValueError(f"injection.fluctuation: keys are {sorted(FLUCT_KEYS)}")
    face = fl.get("face")
    if face not in inj.faces:
        raise ValueError("injection.fluctuation: 'face' must be one of the injecting faces")
    f = FACES.index(face)
    axis = f >> 1
    if cfg.window is not None:
        raise ValueError("injection.fluctuation is not supported with a moving window")
    tol = float(fl.get("tolerance", 1e-2))
    names = injected_species(cfg)
    index = {s.name: i for i, s in enumerate(cfg.species)}
    modes = fl.get("modes") or []
    if not modes:
        raise ValueError("injection.fluctuation needs at least one mode")
    rows = {index[n]: [] for n in names}
    qn = sum(abs(cfg.species[index[n]].charge) * sum(b.density for b in cfg.species[index[n]].beams) for n in names)
    du_n = {index[n]: 0.0 for n in names}
    out_modes = []
    for m in modes:
        if not isinstance(m, dict) or set(m) - MODE_KEYS:
            raise ValueError(f"inflow fluctuation mode: keys are {sorted(MODE_KEYS)}")
        k = [float(x) for x in m.get("k", ())]
        if len(k) != 3:
            raise ValueError("inflow fluctuation mode: k needs three components")
        if k[axis] != 0.0:
            raise ValueError("inflow fluctuation mode: k must have no component along the normal of the face")
        w, ph = float(m.get("omega", 0.0)), float(m.get("phase", 0.0))
        Bc, Bs = _pair(m.get("B"), "B")
        Ec, Es = _pair(m.get("E"), "E")
        us = m.get("u") or {}
        if not isinstance(us, dict) or set(us) - set(names):
            raise ValueError(f"inflow fluctuation mode: u maps injected species ({names}) to their velocity")
        J = [[0.0] * 3, [0.0] * 3]
        for n in names:
            sc = cfg.species[index[n]]
            uc, usn = _pair(us.get(n), f"u[{n}]")
            dens = sum(b.density for b in sc.beams)
            gam = sum(b.density * math.sqrt(1.0 + _dot(b.drift, b.drift)) for b in sc.beams) / max(dens, 1e-300)
            for part, v in ((0, uc), (1, usn)):
                for a in range(3):
                    J[part][a] += sc.charge * dens * v[a] / gam
            rows[index[n]].append([*k, w, ph, *uc, *usn])
            du_n[index[n]] += math.hypot(uc[axis], usn[axis])
        kk = _norm(k)
        scale = max(kk * max(_norm(Bc), _norm(Bs)), abs(w) * max(_norm(Ec), _norm(Es)),
                    _norm(J[0]), _norm(J[1]), 1e-300)
        checks = [
            ("k . B", [_dot(k, Bc), _dot(k, Bs)], kk * max(_norm(Bc), _norm(Bs), 1e-300)),
            ("the charge density k . E", [_dot(k, Ec), _dot(k, Es)], qn),
            ("omega B - k x E", [_norm([w * Bc[a] - _cross(k, Ec)[a] for a in range(3)]),
                                 _norm([w * Bs[a] - _cross(k, Es)[a] for a in range(3)])], scale),
            # Ampere: J_cos = k x B_sin + omega E_sin, J_sin = -k x B_cos - omega E_cos
            ("Ampere (cos)", [_norm([J[0][a] - _cross(k, Bs)[a] - w * Es[a] for a in range(3)])], scale),
            ("Ampere (sin)", [_norm([J[1][a] + _cross(k, Bc)[a] + w * Ec[a] for a in range(3)])], scale),
        ]
        for what, vals, ref in checks:
            if max(abs(v) for v in vals) > tol * ref:
                raise ValueError(f"inflow fluctuation mode k={k}: {what} is {max(abs(v) for v in vals):.3g}, above "
                                 f"{tol:g} of {ref:.3g} (PIC-D26)")
        out_modes.append({"k": k, "omega": w, "phase": ph, "B": (Bc, Bs), "E": (Ec, Es)})
    return {"face": f, "modes": out_modes, "rows": rows, "du_n": du_n}
