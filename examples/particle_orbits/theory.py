"""Exact relativistic orbits of a test particle (charge q, mass m; c = 1) used by the
particle-orbit example. All functions return arrays over the requested lab times t.

gyration     uniform B: u_perp rotates at omega = q B / (m gamma), radius u_perp m / |q B|
exb          E perpendicular to B, |E| < |B|: guiding centre drifts at E x B / B^2 (exact, any gamma)
hyperbolic   uniform E from rest: u = q E t / m, x = (m / qE) (sqrt(1 + (qEt/m)^2) - 1)
e_gt_b       E perpendicular to B, |E| > |B|: in the frame moving at B/E (along E x B) only
             E' = E / gamma_f remains -> hyperbolic motion there, boosted back to the lab
plane_wave   linearly polarised vacuum plane wave E_y = E0 cos(psi), B_z = E0 cos(psi),
             psi = omega t - k x + phi0, particle at rest at t = 0: light-front invariant
             gamma - u_x = 1, u_y = (q E0 / m omega)(sin psi - sin psi0), u_x = u_y^2 / 2
"""
from __future__ import annotations

import numpy as np


def gyration(t, x0, u0, q, m, B):
    """B along z; u0 = (ux, uy, uz). Returns positions (3, n) and momenta (3, n)."""
    t = np.asarray(t, dtype=float)
    ux, uy, uz = u0
    g = np.sqrt(1 + ux * ux + uy * uy + uz * uz)
    w = q * B / (m * g)  # signed angular frequency: du/dt = (q/m) u x B / gamma
    c, s = np.cos(w * t), np.sin(w * t)
    uxt = ux * c + uy * s
    uyt = -ux * s + uy * c
    # x(t) = x0 + int u/g dt
    x = x0[0] + (ux * s - uy * (c - 1)) / (w * g)
    y = x0[1] + (ux * (c - 1) + uy * s) / (w * g)
    z = x0[2] + uz / g * t
    return np.stack([x, y, z]), np.stack([uxt, uyt, np.full_like(t, uz)])


def hyperbolic(t, x0, q, m, E, axis=0):
    """Uniform E along `axis`, particle at rest at t = 0."""
    t = np.asarray(t, dtype=float)
    a = q * E / m
    u = a * t
    x = np.tile(np.asarray(x0, dtype=float)[:, None], (1, len(t)))
    x[axis] += (np.sqrt(1 + u * u) - 1) / a
    uu = np.zeros_like(x)
    uu[axis] = u
    return x, uu


def _frame_time(lab_time, t):
    """Frame time tp of every lab time t, given the increasing map lab_time(tp) with
    lab_time(0) = 0: the bracket is doubled until it covers max(t), then bisected (exact to
    float64; a fixed frame-time grid can stop short of t when the frame moves near c)."""
    t = np.asarray(t, dtype=float)
    hi = max(1.0, float(t.max()))
    for _ in range(200):
        if lab_time(np.array([hi]))[0] >= t.max():
            break
        hi *= 2.0
    else:  # pragma: no cover - lab_time grows at least like tp / gamma_frame
        raise RuntimeError("lab time does not reach the requested time")
    lo, up = np.zeros_like(t), np.full_like(t, hi)
    for _ in range(100):
        mid = 0.5 * (lo + up)
        below = lab_time(mid) < t
        lo, up = np.where(below, mid, lo), np.where(below, up, mid)
    return 0.5 * (lo + up)


def exb_drift(E, B):
    """Guiding-centre drift speed along x for E along y, B along z (|E| < |B|)."""
    return E / B


def e_gt_b(t, x0, q, m, E, B):
    """E along y, B along z with |E| > |B|, particle at rest at t = 0 at x0 (exact)."""
    t = np.asarray(t, dtype=float)
    beta = B / E
    gf = 1.0 / np.sqrt(1 - beta * beta)
    Ep = E / gf  # frame moving with +beta along x: B' = 0, E' = E / gamma_f along y
    upx = -gf * beta  # the particle, at rest in the lab, moves at -beta there
    a = q * Ep / m
    g0 = np.sqrt(1 + upx * upx)

    # orbit in the frame, parametrised by its time tp (the particle starts at the origin, tp = 0)
    def xp_of(tp):
        return upx / a * np.arcsinh(a * tp / g0)

    # back to the lab: t = gf (tp + beta xp), x = gf (xp + beta tp)
    tp = _frame_time(lambda s: gf * (s + beta * xp_of(s)), t)
    upy = a * tp
    gp = np.sqrt(1 + upx * upx + upy * upy)
    yp = (gp - g0) / a
    xp = xp_of(tp)
    x = gf * (xp + beta * tp) + x0[0]
    y = yp + x0[1]
    ux = gf * (upx + beta * gp)
    return np.stack([x, y, np.full_like(t, x0[2])]), np.stack([ux, upy, np.zeros_like(t)])


def plane_wave(t, x0, q, m, E0, omega, phi0=0.0, n_psi=400001):
    """Particle initially at rest at x0 in E_y = B_z = E0 cos(omega (t - x) + phi0) (vacuum, c = 1).
    Returns positions (3, n), momenta (3, n). With a0 = |q| E0 / (m omega), the drift velocity of
    a particle starting at a node of A (phi0 = 0 mod pi) is a0^2 / (4 + a0^2)."""
    t = np.asarray(t, dtype=float)
    a = q * E0 / (m * omega)
    psi0 = phi0 - omega * x0[0]
    # psi grows as d psi / dt = omega (1 - v_x) = omega / gamma (gamma - u_x = 1)
    psi_max = psi0 + omega * t.max() + 10.0
    psi = np.linspace(psi0, psi_max, n_psi)
    uy = a * (np.sin(psi) - np.sin(psi0))
    ux = 0.5 * uy * uy
    g = 1.0 + ux
    dpsi = psi[1] - psi[0]

    def cumint(f):  # trapezoid cumulative integral over psi
        out = np.zeros_like(f)
        out[1:] = np.cumsum(0.5 * (f[1:] + f[:-1])) * dpsi
        return out

    tt = cumint(g) / omega
    xx = x0[0] + cumint(ux) / omega
    yy = x0[1] + cumint(uy) / omega
    x = np.interp(t, tt, xx)
    y = np.interp(t, tt, yy)
    return np.stack([x, y, np.full_like(t, x0[2])]), np.stack([np.interp(t, tt, ux), np.interp(t, tt, uy),
                                                               np.zeros_like(t)])


def exb(t, x0, q, m, E, B):
    """E along y, B along z, |E| < |B|, particle at rest at t = 0 (exact): in the frame drifting at
    beta = E/B along x there is only B' = B / gamma_d, where the particle gyrates with u' = -gamma_d beta."""
    t = np.asarray(t, dtype=float)
    beta = E / B
    gd = 1.0 / np.sqrt(1 - beta * beta)

    def frame(tp):
        return gyration(tp, (0.0, 0.0, 0.0), (-gd * beta, 0.0, 0.0), q, m, B / gd)

    # back to the lab: t = gd (tp + beta x'), x = gd (x' + beta tp)
    tp = _frame_time(lambda s: gd * (s + beta * frame(s)[0][0]), t)
    Xp, Up = frame(tp)
    gp = np.sqrt(1 + np.sum(Up ** 2, axis=0))
    ux = gd * (Up[0] + beta * gp)
    return (np.stack([gd * (Xp[0] + beta * tp) + x0[0], Xp[1] + x0[1], np.full_like(t, x0[2])]),
            np.stack([ux, Up[1], np.zeros_like(t)]))
