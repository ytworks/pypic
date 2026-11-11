"""External static fields (PIC-D12): fields the particles feel in addition to the Yee fields.

With external fields the Yee arrays hold the difference to them. Maxwell's equations are linear
and the external fields are static solutions (a uniform E, a uniform B, dipoles kept up by a
current in their core that is not part of the plasma current), so the field solver, the absorbing
layer and the wall conditions act on the difference fields unchanged, and Gauss's law and div B
hold for them as before.

Pure array code: ``xp`` is NumPy (CPU tests) or CuPy (the engine).
"""
from __future__ import annotations


def check_config(cfg) -> bool:
    """True if ``cfg`` sets an external field; raises on the combinations that are not supported."""
    if not cfg.external_fields.active:
        return False
    if cfg.solver != "explicit":
        raise ValueError("external_fields are for the explicit solver only (PIC-D12)")
    if any(cfg.B0) or any(cfg.E0):
        raise ValueError("external_fields cannot be combined with B0 / E0: B0 and E0 are initial values of the "
                         "Yee fields; give the uniform external fields in external_fields.B / .E")
    return True


def dipole_B(xp, index, dx, moment, center, core):
    """Dipole field on the Yee faces of a block: Bx at (i, j+1/2, k+1/2), By, Bz likewise.

    ``index``: per axis the global node indices of the block (1-D float arrays). B is the discrete
    curl of ``A = m x r / (r^2 + core^2)^(3/2)`` taken on the Yee edges, so its discrete divergence
    vanishes to round-off, and every value depends on global indices only (not on the block).
    For r >> core, B = (3 (m.n) n - m) / r^3; at the centre B = 2 m / core^3.
    """
    shape = tuple(len(i) for i in index)
    a2 = float(core) ** 2

    def edge_A(c):  # A_c on the edges along c, one more node than the block at the top of every axis
        r = []
        for a in range(3):
            node = xp.concatenate([index[a], index[a][-1:] + 1.0])
            x = (node + (0.5 if a == c else 0.0)) * dx[a] - center[a]
            r.append(x.reshape([-1 if b == a else 1 for b in range(3)]))
        s = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + a2
        c1, c2 = (c + 1) % 3, (c + 2) % 3
        return (moment[c1] * r[c2] - moment[c2] * r[c1]) / (s * xp.sqrt(s))

    def diff(F, a):  # forward difference along a, on the block
        lo = [slice(0, n) for n in shape]
        hi = list(lo)
        hi[a] = slice(1, shape[a] + 1)
        return (F[tuple(hi)] - F[tuple(lo)]) / dx[a]

    A = [edge_A(c) for c in range(3)]
    return tuple(diff(A[(c + 2) % 3], (c + 1) % 3) - diff(A[(c + 1) % 3], (c + 2) % 3) for c in range(3))


def dipoles_B(xp, index, dx, dipoles):
    """Sum of the dipole fields (in list order) on the Yee faces of a block; None without dipoles."""
    total = None
    for d in dipoles:
        B = dipole_B(xp, index, dx, d.moment, d.center, d.core)
        total = B if total is None else tuple(t + b for t, b in zip(total, B))
    return total
