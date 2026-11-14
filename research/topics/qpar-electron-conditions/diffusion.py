"""Diffusion of test electrons along x in a frozen 1-D field (README, measurement 3; null test N7).

Method of gupta2025speed (S), appendix B: the fields of one output along x are cut to a window, mirrored so that
they are continuous and periodic, and held fixed; test particles are pushed in them with the Boris scheme and their
displacement along x is followed.

Conventions (README, "Reproduction"; they are the ones of gupta2025speed):

- D = <(x(t) - x(0))^2> / t, without a factor 1/2 (its Eq. 9). A process with <dx^2> = 2 kappa t has D = 2 kappa.
- D_B = v r_L = p^2 / (gamma m |q| B), without a factor 1/3 (its Eq. 2). Isotropic scattering at the rate nu along
  a field has kappa = v^2 / (3 nu), so D / D_B = 2 Omega / (3 nu) with Omega = |q| B / (gamma m).
- A value of D is taken only from a diffusive interval: a stretch of a factor ``factor`` in time over which the
  fitted slope of ln <dx^2> against ln t is within 1 +- ``tol`` (N7: factor 3, 0.1). The value is <dx^2> / t at the
  end of the last such interval. The slope of a straight line through <dx^2>(t) over the same interval is given
  beside it: the two differ by about (correlation time) / t, which shows how far the limit of Eq. 9 is reached.

Units are the engine's: lengths c / omega_pe, times 1 / omega_pe, u = gamma v in c, du/dt = (q/m) (E + v x B).
The arrays may be NumPy or CuPy arrays; ``xp`` is the module that made them.
"""
from __future__ import annotations

import numpy as np

# offset along x of each component in cells (pypic.output.names.STAGGER, x entry)
STAGGER_X = {"Ex": 0.5, "Ey": 0.0, "Ez": 0.0, "Bx": 0.0, "By": 0.5, "Bz": 0.5}


class Window:
    """Periodic 1-D fields: ``f[name]`` holds n samples at (i + stagger) dx, the period is n dx.
    A component that is absent is zero."""

    def __init__(self, fields: dict, dx: float, stagger: dict | None = None):
        self.f = dict(fields)
        self.dx = float(dx)
        self.stagger = dict(STAGGER_X if stagger is None else stagger)
        sizes = {len(v) for v in self.f.values()}
        if len(sizes) != 1:
            raise ValueError(f"components of different lengths: {sorted(sizes)}")
        self.n = sizes.pop()

    @property
    def length(self) -> float:
        return self.n * self.dx

    def at(self, name: str, x, xp=np):
        """Component ``name`` at the positions x (any real number; linear interpolation)."""
        s = x / self.dx - self.stagger[name]
        i = xp.floor(s)
        w = s - i
        i = i.astype(xp.int64) % self.n
        f = self.f[name]
        return f[i] * (1.0 - w) + f[(i + 1) % self.n] * w


def uniform(b, dx: float = 1.0, n: int = 8, xp=np) -> Window:
    """A window with the constant field ``b`` = (Bx, By, Bz) and no electric field."""
    return Window({c: xp.full(n, float(v)) for c, v in zip(("Bx", "By", "Bz"), b)}, dx)


def mirrored(fields: dict, i0: int, i1: int, dx: float, xp=np) -> Window:
    """The cells [i0, i1) of the lineouts ``fields`` (name -> samples at (i + stagger) dx), followed by their mirror
    image about x = i1 dx. The result is continuous and has the period 2 (i1 - i0) dx; x = 0 of the window is
    x = i0 dx of the lineout. A component at the nodes needs the sample i1."""
    out = {}
    for name, f in fields.items():
        f = xp.asarray(f, dtype=xp.float64)
        if STAGGER_X[name]:
            part = f[i0:i1]
            out[name] = xp.concatenate([part, part[::-1]])
        else:
            if i1 >= len(f):
                raise ValueError(f"{name}: the mirror needs the sample {i1}, the lineout has {len(f)}")
            out[name] = xp.concatenate([f[i0:i1], f[i0 + 1:i1 + 1][::-1]])
    return Window(out, dx)


def isotropic(n: int, u: float, rng, xp=np):
    """n proper velocities of magnitude u with directions uniform on the sphere; array (3, n)."""
    mu = 2.0 * rng.random(n) - 1.0
    phi = 2.0 * np.pi * rng.random(n)
    s = xp.sqrt(1.0 - mu * mu)
    return u * xp.stack([mu, s * xp.cos(phi), s * xp.sin(phi)])


def bohm(u: float, b: float, qm: float = -1.0) -> float:
    """D_B = v r_L of a particle of proper velocity u = gamma v in a field of strength b."""
    return u * u / (np.sqrt(1.0 + u * u) * abs(qm) * b)


def sample_steps(n_steps: int, per_decade: int = 24, first: int = 8) -> np.ndarray:
    """Steps at which the displacement is recorded: about ``per_decade`` per decade from ``first`` to n_steps."""
    k = np.unique(np.round(np.logspace(np.log10(first), np.log10(n_steps),
                                       int(per_decade * np.log10(n_steps / first)) + 1)).astype(np.int64))
    return k[k <= n_steps]


def push(win: Window, x, u, dt: float, n_steps: int, samples=None, qm: float = -1.0, noise: float = 0.0,
         rng=None, xp=np):
    """Push particles in the frozen fields of ``win`` with the Boris scheme.

    x (n,) and u (3, n) = gamma v are advanced in place; x is not wrapped, so x - x(0) is the displacement.
    ``qm`` = q / m (-1 for electrons).
    ``samples``: increasing step numbers; returns (t, <dx>, <dx^2>) at those steps.
    ``noise``: standard deviation of a magnetic field drawn anew for every particle, component and step (a field
    without memory, for N7). A direction then diffuses on the sphere, and the correlation of the velocity falls at
    the rate nu = (qm noise / gamma)^2 dt, up to terms of relative order nu dt.
    """
    x0 = x.copy()
    samples = np.asarray([] if samples is None else samples, dtype=np.int64)
    mean, msd = np.zeros(len(samples)), np.zeros(len(samples))
    has_e = any(c in win.f for c in ("Ex", "Ey", "Ez"))
    h = 0.5 * qm * dt
    j = 0
    for step in range(1, n_steps + 1):
        b = [win.at(c, x, xp) if c in win.f else 0.0 for c in ("Bx", "By", "Bz")]
        if noise:
            b = [bc + noise * rng.standard_normal(x.shape) for bc in b]
        if has_e:
            for k, c in enumerate(("Ex", "Ey", "Ez")):
                if c in win.f:
                    u[k] += h * win.at(c, x, xp)
        g = xp.sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
        tx, ty, tz = (h / g * bc for bc in b)
        f = 2.0 / (1.0 + tx * tx + ty * ty + tz * tz)
        vx = u[0] + u[1] * tz - u[2] * ty
        vy = u[1] + u[2] * tx - u[0] * tz
        vz = u[2] + u[0] * ty - u[1] * tx
        u[0] += f * (vy * tz - vz * ty)
        u[1] += f * (vz * tx - vx * tz)
        u[2] += f * (vx * ty - vy * tx)
        if has_e:
            for k, c in enumerate(("Ex", "Ey", "Ez")):
                if c in win.f:
                    u[k] += h * win.at(c, x, xp)
            g = xp.sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
        x += u[0] / g * dt
        if j < len(samples) and step == samples[j]:
            d = x - x0
            mean[j], msd[j] = float(d.mean()), float((d * d).mean())
            j += 1
    return samples * dt, mean, msd


def diffusive_interval(t, msd, factor: float = 3.0, tol: float = 0.1):
    """The last stretch of at least ``factor`` in time whose fitted log-slope is within 1 +- tol.

    Returns None if there is none, else a dict: ``t`` (start, end), ``log_slope``, ``d`` = <dx^2> / t at the end,
    ``d_line`` = the slope of a straight line through <dx^2>(t) over the stretch, ``n`` = samples in it.
    """
    t, msd = np.asarray(t, dtype=float), np.asarray(msd, dtype=float)
    for k in range(len(t) - 1, 1, -1):
        inside = np.nonzero(t <= t[k] / factor)[0]
        if not len(inside):
            break
        j = inside[-1]
        if k - j < 3 or np.any(msd[j:k + 1] <= 0):
            continue
        slope = float(np.polyfit(np.log(t[j:k + 1]), np.log(msd[j:k + 1]), 1)[0])
        if abs(slope - 1.0) <= tol:
            return {"t": [float(t[j]), float(t[k])], "log_slope": slope, "d": float(msd[k] / t[k]),
                    "d_line": float(np.polyfit(t[j:k + 1], msd[j:k + 1], 1)[0]), "n": int(k - j + 1)}
    return None
