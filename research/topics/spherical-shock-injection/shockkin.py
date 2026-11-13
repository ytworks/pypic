"""Electron reflection at a curved, expanding shock front: kinematics and a count model.

Pure numpy, no I/O. Units: c = 1, m_e = 1, lengths in d_e, times in 1/omega_pe.

Geometry. The upstream plasma is at rest in the lab frame (E = 0) in a uniform field B0 along z,
so a field line is a line of constant cylindrical radius rho from the axis through the centre.
The front is axisymmetric about that axis and symmetric about the equator: R(theta, t), with
theta the polar angle from B0. Everything below is written for the upper hemisphere
(0 <= theta <= pi/2); the lower one is its mirror image.

The foot point of a field line on the front, z_f(rho, t), moves along the field line with
    V_HT = dz_f/dt at constant rho = R_t R / (R_theta sin(theta) + R cos(theta)),
which is the speed of the de Hoffmann-Teller (HT) frame, V_n / cos(theta_Bn). In that frame the
motional electric field vanishes and an electron keeps its energy (up to the potential phi) and
its magnetic moment, so it is reflected when its pitch angle lies outside the loss cone. This is
the known one-step reflection ("fast Fermi", shock drift acceleration in the HT picture); the
protocol calls it S0 and uses it as the baseline that hypotheses have to exceed.
"""
from __future__ import annotations

import numpy as np

GAMMA_AD = 5.0 / 3.0


# ------------------------------------------------------------------ plasma speeds
def sound_speed(v_a: float, beta: float) -> float:
    """Sound speed for a total plasma beta (the definition of pypic.setup.mhd.fast_speed)."""
    return v_a * np.sqrt(0.5 * GAMMA_AD * beta)


def sound_speed_of(temp_e, temp_i, m_i: float):
    """The same sound speed from the temperatures (m_e c^2) and the ion mass (m_e): beta = 2 n (T_e + T_i) / B^2."""
    return np.sqrt(GAMMA_AD * (np.asarray(temp_e, float) + np.asarray(temp_i, float)) / m_i)


def fast_speed(v_a, c_s, theta_bn):
    """Fast magnetosonic speed for propagation at the angle theta_bn (radians) to the field."""
    s = v_a ** 2 + c_s ** 2
    return np.sqrt(0.5 * (s + np.sqrt(np.maximum(s ** 2 - 4.0 * v_a ** 2 * c_s ** 2 * np.cos(theta_bn) ** 2, 0.0))))


# ------------------------------------------------------------------ front models
def sphere(r0: float, speed: float):
    """Front model of a sphere: R = r0 + speed * t. Returns front(theta, t) -> (R, R_t, R_theta)."""

    def front(theta, t):
        th, tt = np.broadcast_arrays(np.asarray(theta, float), np.asarray(t, float))
        return r0 + speed * tt, np.full(th.shape, float(speed)), np.zeros(th.shape)

    return front


def _basis(theta, tau, order_t):
    cols, dt_cols, dth_cols = [], [], []
    for k in range(3):  # 1, cos 2 theta, cos 4 theta: symmetric about the axis and the equator
        ang, dang = np.cos(2 * k * theta), -2 * k * np.sin(2 * k * theta)
        for j in range(order_t + 1):
            cols.append(ang * tau ** j)
            dt_cols.append(ang * (j * tau ** (j - 1) if j else 0.0 * tau))
            dth_cols.append(dang * tau ** j)
    return np.stack(cols, -1), np.stack(dt_cols, -1), np.stack(dth_cols, -1)


class FittedFront:
    """A fitted front, callable as front(theta, t) -> (R, R_t, R_theta). It can be saved as a dict.

    The fit is a polynomial in time: it is not to be used outside t_range."""

    def __init__(self, coef, t0: float, span: float, order_t: int):
        self.coef, self.t0, self.span, self.order_t = np.asarray(coef, float), float(t0), float(span), int(order_t)

    @property
    def t_range(self) -> tuple:
        return self.t0, self.t0 + self.span

    def __call__(self, theta, t):
        th, tt = np.broadcast_arrays(np.asarray(theta, float), np.asarray(t, float))
        b, bt, bth = _basis(th, (tt - self.t0) / self.span, self.order_t)
        return b @ self.coef, (bt @ self.coef) / self.span, bth @ self.coef

    def to_dict(self) -> dict:
        return {"coef": self.coef.tolist(), "t0": self.t0, "span": self.span, "order_t": self.order_t}

    @classmethod
    def from_dict(cls, d: dict) -> "FittedFront":
        return cls(d["coef"], d["t0"], d["span"], d["order_t"])


def fit_front(theta, t, radius, order_t: int = 3):
    """Least-squares fit of measured front radii by {1, cos 2 theta, cos 4 theta} x {1..tau^order_t}.

    theta, t, radius: 1-D arrays of the measured points (theta folded into [0, pi/2] by the caller
    or not: the basis is symmetric). Returns (front, rms residual, (t_min, t_max))."""
    theta, t, radius = (np.asarray(a, float) for a in (theta, t, radius))
    t0, t1 = float(t.min()), float(t.max())
    span = max(t1 - t0, 1e-30)
    a, _, _ = _basis(theta, (t - t0) / span, order_t)
    coef, *_ = np.linalg.lstsq(a, radius, rcond=None)
    rms = float(np.sqrt(np.mean((a @ coef - radius) ** 2)))
    return FittedFront(coef, t0, span, order_t), rms, (t0, t1)


def foot(front, theta, t) -> dict:
    """Foot-point quantities of the field line that meets the front at polar angle theta at time t.

    Keys: rho, z, v_ht (speed of the foot point along the field line), cos_tbn, tbn (radians),
    v_n (normal speed of the front), drho_dtheta. Where the front is locally parallel to the
    field or folds back (R_theta sin + R cos <= 0), v_ht is inf and cos_tbn is 0."""
    theta = np.asarray(theta, float)
    r, rt, rth = front(theta, t)
    den = rth * np.sin(theta) + r * np.cos(theta)
    ok = den > 1e-12 * np.maximum(r, 1e-300)
    cos_tbn = np.where(ok, den, 0.0) / np.sqrt(r ** 2 + rth ** 2)
    v_n = rt * r / np.sqrt(r ** 2 + rth ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        v_ht = np.where(ok, rt * r / np.where(ok, den, 1.0), np.inf)
    return {"rho": r * np.sin(theta), "z": r * np.cos(theta), "v_ht": v_ht, "cos_tbn": cos_tbn,
            "tbn": np.arccos(np.clip(cos_tbn, 0.0, 1.0)), "v_n": v_n, "drho_dtheta": np.where(ok, den, 0.0)}


# ------------------------------------------------------------------ the reflection map
def ht_reflect(u_par, u_perp, v_ht, b, phi=0.0):
    """One adiabatic reflection off a mirror that moves along the field at v_ht.

    u_par is the momentum along the motion of the foot point (away from the centre), u_perp the
    momentum across the field; b = B_max / B_0 (> 1); phi >= 0 is the rise of the potential in
    the HT frame in units of m_e c^2 / e (it pulls electrons through the front).

    Returns (encounter, reflected, u_par_r, gamma_r):
      encounter  the foot point catches up with the electron (v_par < v_ht < 1);
      reflected  encounter and outside the loss cone: (gamma' + phi)^2 - 1 - b u_perp^2 < 0;
      u_par_r, gamma_r  the state after a reflection (computed for every encounter).
    The map (u_par, u_perp) -> (u_par_r, u_perp) is its own inverse: applied to a reflected
    electron it returns the momentum before the reflection."""
    u_par, u_perp, v_ht, b, phi = np.broadcast_arrays(*(np.asarray(a, float) for a in (u_par, u_perp, v_ht, b, phi)))
    gam = np.sqrt(1.0 + u_par ** 2 + u_perp ** 2)
    sub = v_ht < 1.0
    v = np.where(sub, v_ht, 0.0)
    g = 1.0 / np.sqrt(1.0 - v ** 2)
    up = g * (u_par - v * gam)  # HT frame
    gp = g * (gam - v * u_par)
    encounter = sub & (up < 0.0)
    reflected = encounter & ((gp + phi) ** 2 - 1.0 - b * u_perp ** 2 < 0.0)
    return encounter, reflected, g * (-up + v * gp), g * (gp - v * up)


def before_reflection(u_par_r, u_perp, v_ht):
    """Momentum and energy before a reflection that left the electron with (u_par_r, u_perp).
    Returns (u_par, gamma, gamma_ht): gamma_ht is the energy in the HT frame (for the loss cone)."""
    u_par_r, u_perp, v_ht = np.broadcast_arrays(*(np.asarray(a, float) for a in (u_par_r, u_perp, v_ht)))
    gam = np.sqrt(1.0 + u_par_r ** 2 + u_perp ** 2)
    v = np.where(v_ht < 1.0, v_ht, np.nan)
    g = 1.0 / np.sqrt(1.0 - v ** 2)
    up = g * (u_par_r - v * gam)
    gp = g * (gam - v * u_par_r)
    return g * (-up + v * gp), g * (gp - v * up), gp


# ------------------------------------------------------------------ the upstream distribution
def mj_norm(temp: float) -> float:
    """Normalisation of exp(-(gamma - 1) / temp) over d^3u (Maxwell-Juttner, the engine's loader)."""
    u = np.linspace(0.0, 40.0 * np.sqrt(temp) + 20.0 * temp, 40001)
    f = 4.0 * np.pi * u ** 2 * np.exp(-(np.sqrt(1.0 + u ** 2) - 1.0) / temp)
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(u)))


def mj_tail_fraction(temp: float, gamma_min: float) -> float:
    """Fraction of a Maxwell-Juttner distribution at rest with gamma >= gamma_min."""
    u0 = np.sqrt(max(gamma_min ** 2 - 1.0, 0.0))
    u = np.linspace(u0, u0 + 40.0 * np.sqrt(temp) + 20.0 * temp, 40001)
    f = 4.0 * np.pi * u ** 2 * np.exp(-(np.sqrt(1.0 + u ** 2) - 1.0) / temp)
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(u))) / mj_norm(temp)


def reflected_flux(v_ht: float, temp: float, b: float, phi: float, gamma_edges, v_edges, n_par=320, n_perp=240,
                   u_max=12.0):
    """Reflected electrons per unit area across the field, per unit time and unit upstream density.

    Midpoint quadrature of f (v_ht - v_par) over the reflected part of momentum space for a
    Maxwell-Juttner distribution of temperature temp (m_e c^2). Returns (flux[gamma bin,
    v_r bin], encounter flux): the first is binned by the energy and the parallel speed after
    the reflection, the second counts every electron the foot point catches."""
    ng, nv = len(gamma_edges) - 1, len(v_edges) - 1
    if not v_ht < 1.0:
        return np.zeros((ng, nv)), 0.0
    umax = u_max * np.sqrt(temp)
    dpar, dperp = 2.0 * umax / n_par, umax / n_perp
    u_par = (-umax + dpar * (np.arange(n_par) + 0.5))[:, None]
    u_perp = (dperp * (np.arange(n_perp) + 0.5))[None, :]
    gam = np.sqrt(1.0 + u_par ** 2 + u_perp ** 2)
    weight = np.exp(-(gam - 1.0) / temp) / mj_norm(temp) * 2.0 * np.pi * u_perp * dpar * dperp
    weight = weight * np.maximum(v_ht - u_par / gam, 0.0)
    enc, refl, u_r, g_r = ht_reflect(u_par, u_perp, v_ht, b, phi)
    w = (weight * refl).ravel()
    ig = np.searchsorted(gamma_edges, g_r.ravel(), side="right") - 1
    iv = np.searchsorted(v_edges, (u_r / g_r).ravel(), side="right") - 1
    ok = (w > 0) & (ig >= 0) & (ig < ng) & (iv >= 0) & (iv < nv)
    flux = np.bincount(ig[ok] * nv + iv[ok], weights=w[ok], minlength=ng * nv).reshape(ng, nv)
    return flux, float(np.sum(weight * enc))


def reflection_counts(front, t_grid, theta_grid, t_dumps, temp, b, phi, density, z_half, gamma_edges,
                      band_edges_deg, valid=None, n_v=200):
    """Count model B2: reflected electrons that are still in the box at each dump time.

    front        front(theta, t) -> (R, R_t, R_theta)
    t_grid       uniform times of reflection; theta_grid: uniform foot angles in (0, pi/2)
    temp, b, phi scalars or arrays (len(t_grid), len(theta_grid)): upstream temperature,
                 compression B_max/B_0, HT potential
    density      upstream macro-electrons per unit volume (or 1 for physical numbers)
    z_half       distance from the centre to the open face along the field
    valid        optional boolean array (t, theta): reflection is counted only where it is True
    Returns dict(counts[dump, band, gamma bin], reflected[band] (all reflections up to the last
    dump, in or out of the box), encountered[band]); both hemispheres are included."""
    t_grid, theta_grid = np.asarray(t_grid, float), np.asarray(theta_grid, float)
    t_dumps = np.asarray(t_dumps, float)
    gamma_edges = np.asarray(gamma_edges, float)
    band = np.radians(np.asarray(band_edges_deg, float))
    v_edges = np.linspace(0.0, 1.0, n_v + 1)
    dt, dth = t_grid[1] - t_grid[0], theta_grid[1] - theta_grid[0]
    shape = (len(t_grid), len(theta_grid))
    temp, b, phi = (np.broadcast_to(np.asarray(a, float), shape) for a in (temp, b, phi))
    nb, ng = len(band) - 1, len(gamma_edges) - 1
    counts = np.zeros((len(t_dumps), nb, ng))
    reflected, encountered = np.zeros(nb), np.zeros(nb)
    cache: dict = {}
    for it, t in enumerate(t_grid):
        if t >= t_dumps[-1]:
            break
        f = foot(front, theta_grid, t)
        for ith in range(len(theta_grid)):
            v_ht = float(f["v_ht"][ith])
            if not (0.0 < v_ht < 1.0) or (valid is not None and not valid[it, ith]):
                continue
            ib = int(np.searchsorted(band, f["tbn"][ith], side="right") - 1)
            if not 0 <= ib < nb:
                continue
            # 2 hemispheres x ring area across the field x time step
            area = 2.0 * 2.0 * np.pi * f["rho"][ith] * f["drho_dtheta"][ith] * dth * dt * density
            key = (round(v_ht, 4), round(float(temp[it, ith]), 6), round(float(b[it, ith]), 3),
                   round(float(phi[it, ith]), 4))
            if key not in cache:
                cache[key] = reflected_flux(key[0], key[1], key[2], key[3], gamma_edges, v_edges)
            flux, enc = cache[key]
            encountered[ib] += enc * area
            reflected[ib] += flux.sum() * area
            cum = np.concatenate([np.zeros((ng, 1)), np.cumsum(flux, axis=1)], axis=1)  # below v_edges[k]
            for idump, td in enumerate(t_dumps):
                if td <= t:
                    continue
                v_max = (z_half - f["z"][ith]) / (td - t)  # slower electrons are still in the box
                x = np.clip(v_max, 0.0, 1.0) * n_v
                k = min(int(x), n_v - 1)
                counts[idump, ib] += (cum[:, k] + (x - k) * flux[:, k]) * area
    return {"counts": counts, "reflected": reflected, "encountered": encountered}


# ------------------------------------------------------------------ per-particle test (B1)
class FootTable:
    """Foot points of the field lines on a (time, rho) grid, for tracing electrons back."""

    def __init__(self, front, t_grid, rho_max: float, n_rho: int = 400, n_theta: int = 721):
        self.t = np.asarray(t_grid, float)
        self.rho = np.linspace(0.0, rho_max, n_rho)
        th = np.linspace(0.0, 0.5 * np.pi, n_theta)
        shape = (len(self.t), n_rho)
        self.z, self.v_ht, self.tbn = (np.full(shape, np.nan) for _ in range(3))
        for i, t in enumerate(self.t):
            f = foot(front, th, t)
            rho_f = np.maximum.accumulate(f["rho"])  # monotonic envelope (a fold is not a foot)
            keep = np.concatenate([[True], np.diff(rho_f) > 0])
            inside = self.rho <= rho_f[-1]
            for arr, key in ((self.z, "z"), (self.v_ht, "v_ht"), (self.tbn, "tbn")):
                arr[i, inside] = np.interp(self.rho[inside], rho_f[keep], np.minimum(f[key][keep], 1e30))

    def at(self, rho, t):
        """(z, v_ht, tbn) of the foot of the field lines rho at the times t (bilinear; NaN outside)."""
        rho, t = np.asarray(rho, float), np.asarray(t, float)
        fi = np.clip((rho - self.rho[0]) / (self.rho[1] - self.rho[0]), 0, len(self.rho) - 1.000001)
        ft = np.clip((t - self.t[0]) / (self.t[1] - self.t[0]), 0, len(self.t) - 1.000001)
        i, j = fi.astype(int), ft.astype(int)
        a, c = fi - i, ft - j
        out = []
        for arr in (self.z, self.v_ht, self.tbn):
            out.append((1 - c) * ((1 - a) * arr[j, i] + a * arr[j, i + 1])
                       + c * ((1 - a) * arr[j + 1, i] + a * arr[j + 1, i + 1]))
        bad = (rho > self.rho[-1]) | (t < self.t[0]) | (t > self.t[-1])
        return tuple(np.where(bad, np.nan, o) for o in out)

    def trace_back(self, rho, z, v_par, t_now, chunk: int = 20000):
        """Time at which an upstream electron left the front, moving freely along its field line.

        rho: guiding-centre radius; z >= 0: distance from the equatorial plane; v_par: speed
        along the field away from that plane; t_now: time of the observation.
        Returns (t_r, v_ht, tbn) at the last crossing of the foot point, NaN where the straight
        orbit never met the front (field line not reached, or the electron was always ahead)."""
        rho, z, v_par = (np.asarray(a, float) for a in (rho, z, v_par))
        use = self.t <= t_now
        tt = self.t[use]
        t_r = np.full(rho.shape, np.nan)
        if len(tt) < 2:
            return t_r, t_r.copy(), t_r.copy()
        fi = np.clip((rho - self.rho[0]) / (self.rho[1] - self.rho[0]), 0, len(self.rho) - 1.000001)
        for s in range(0, rho.size, chunk):
            sl = slice(s, s + chunk)
            i = fi[sl].astype(int)
            a = fi[sl] - i
            zf = (1 - a) * self.z[use][:, i] + a * self.z[use][:, i + 1]  # (nt, n)
            g = zf - (z[sl][None, :] - v_par[sl][None, :] * (t_now - tt[:, None]))
            g = np.where(np.isnan(g), -np.inf, g)
            cross = (g[:-1] >= 0) & (g[1:] < 0)
            k = cross.shape[0] - 1 - np.argmax(cross[::-1], axis=0)  # the last crossing
            has = cross.any(axis=0) & (rho[sl] <= self.rho[-1])
            col = np.arange(g.shape[1])
            g0, g1 = g[k, col], g[k + 1, col]
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(np.isfinite(g1), g0 / (g0 - g1), 1.0)
            t_r[sl] = np.where(has, tt[k] + frac * (tt[k + 1] - tt[k]), np.nan)
        _, v_ht, tbn = self.at(rho, np.where(np.isnan(t_r), self.t[0], t_r))
        return t_r, np.where(np.isnan(t_r), np.nan, v_ht), np.where(np.isnan(t_r), np.nan, tbn)


def guiding_centre(x, y, ux, uy, b0: float, charge: float = -1.0):
    """Guiding centre across a uniform field b0 along z: r + (u x B) / (q B^2)."""
    return x + uy / (charge * b0), y - ux / (charge * b0)


def canonical_radius(x, y, ux, uy, b0: float, charge: float = -1.0):
    """Radius rho_P of the flux surface that the canonical angular momentum names, for the uniform
    field b0 along z: rho_P^2 = 2 P_phi / (q b0) = rho^2 + 2 (x u_y - y u_x) / (q b0), x and y from
    the axis. It equals rho_gc^2 - rho_L^2 and is conserved where the fields are axisymmetric; it is
    a radius only where the field is the uniform one (upstream). The sign of the square is kept for
    orbits that enclose the axis."""
    sq = x ** 2 + y ** 2 + 2.0 * (x * uy - y * ux) / (charge * b0)
    return np.sign(sq) * np.sqrt(np.abs(sq))


def birth_cell(ids, n_cells, ppc: int):
    """Cell (ix, iy, iz) in which an ambient particle was loaded: id = ppc * cell + slot with
    cell = (ix * ny + iy) * nz + iz. Holds for ids below ppc * (number of cells) only."""
    lin = np.asarray(ids, dtype=np.int64) // int(ppc)
    ny, nz = int(n_cells[1]), int(n_cells[2])
    return lin // (ny * nz), (lin // nz) % ny, lin % nz


def birth_radius(ids, n_cells, dx, ppc: int, centre):
    """Cylindrical radius (about the axis along z through `centre`) of the centre of the birth cell."""
    ix, iy, _ = birth_cell(ids, n_cells, ppc)
    return np.hypot((ix + 0.5) * dx[0] - centre[0], (iy + 0.5) * dx[1] - centre[1])


def thermal_bound(temp: float, n_total: float, n_max: float = 10.0) -> float:
    """Kinetic energy above which a Maxwell-Juttner distribution of temperature temp holds fewer
    than n_max of n_total electrons (the bound E_th of the per-electron test B1)."""
    lo, hi = 0.0, 200.0 * temp + 5.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mj_tail_fraction(temp, 1.0 + mid) * n_total > n_max:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def s0_test(u_par_r, u_perp, v_ht, b_max, e_thermal: float):
    """Is a reflected electron explained by one reflection of a thermal electron (S0)?

    It is when (1) the foot point was subluminal, (2) the kinetic energy before the reflection is
    at most e_thermal (thermal_bound) and (3) the electron was outside the loss cone for the
    largest compression b_max (phi = 0, the most permissive case).
    Returns (explained, kinetic energy before, mirrored); the energy is NaN for v_ht >= 1."""
    u_par0, gam0, gam_ht = before_reflection(u_par_r, u_perp, v_ht)
    with np.errstate(invalid="ignore"):
        mirror = gam_ht ** 2 - 1.0 - np.asarray(b_max, float) * np.asarray(u_perp, float) ** 2 < 0.0
        explained = mirror & (gam0 - 1.0 <= e_thermal) & (np.asarray(v_ht, float) < 1.0)
    return explained, gam0 - 1.0, mirror


# ------------------------------------------------------------------ measures on tracked orbits
def pitch_scatter_rate(u, b, t, lags=(4, 16, 64)):
    """Pitch-angle scattering rate from tracked orbits, nu = <(d mu)^2 / (1 - mu^2)> / lag.

    u, b: arrays (3, n_t, n) of momentum and field at the particle; t: (n_t,). mu is the pitch
    cosine. Returns {lag in samples: (nu, number of pairs)}; a diffusive process gives the same
    nu at every lag, a coherent one gives nu rising with the lag."""
    bm = np.sqrt(np.sum(b ** 2, axis=0))
    um = np.sqrt(np.sum(u ** 2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.sum(u * b, axis=0) / (um * bm)
    out = {}
    for lag in lags:
        if lag >= mu.shape[0]:
            continue
        d = mu[lag:] - mu[:-lag]
        mid = 0.5 * (mu[lag:] + mu[:-lag])
        with np.errstate(divide="ignore", invalid="ignore"):
            q = d ** 2 / (1.0 - mid ** 2)
        ok = np.isfinite(q) & (np.abs(mid) < 0.95)
        out[int(lag)] = (float(np.mean(q[ok]) / (t[lag] - t[0])) if ok.any() else float("nan"), int(ok.sum()))
    return out


def centre_momentum(u, t, dt_step: float):
    """Momentum at the time of the recorded fields. The engine records u half a step (dt_step / 2)
    before the fields, so u is moved forward by linear interpolation between two samples.
    u: (3, n_t, n); returns (3, n_t - 1, n) for the first n_t - 1 samples."""
    frac = 0.5 * dt_step / np.diff(t)
    return u[:, :-1] + (u[:, 1:] - u[:, :-1]) * frac[None, :, None]


def work_closure(u, e, t, dt_step: float, charge: float = -1.0, dgamma_min: float = 0.1, b=None):
    """Compare the work of the recorded field with the recorded change of energy, per orbit.

    u, e (and b): (3, n_t, n); u is centred on the field time first (centre_momentum).
    Returns dict(work, dgamma, n_used, median_rel_error) over the orbits with |dgamma| >=
    dgamma_min and no missing sample; with b also work_par (the part done by E along B)."""
    uc = centre_momentum(u, t, dt_step)
    gam = np.sqrt(1.0 + np.sum(uc ** 2, axis=0))
    v = uc / gam
    ee = e[:, :-1]
    dt = np.diff(t[:-1])[:, None]

    def integral(p):
        return np.sum(0.5 * (p[1:] + p[:-1]) * dt, axis=0)

    work = integral(charge * np.sum(ee * v, axis=0))
    dgam = gam[-1] - gam[0]
    ok = np.isfinite(work) & np.isfinite(dgam) & (np.abs(dgam) >= dgamma_min)
    rel = np.abs(work[ok] - dgam[ok]) / np.abs(dgam[ok])
    out = {"work": work, "dgamma": dgam, "n_used": int(ok.sum()),
           "median_rel_error": float(np.median(rel)) if ok.any() else float("nan")}
    if b is not None:
        bb = b[:, :-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            bhat = bb / np.sqrt(np.sum(bb ** 2, axis=0))
        out["work_par"] = integral(charge * np.sum(ee * bhat, axis=0) * np.sum(v * bhat, axis=0))
    return out


def magnetic_moment(u, b):
    """u_perp^2 / (2 |B|) along tracked orbits; u, b: (3, n_t, n)."""
    bm = np.sqrt(np.sum(b ** 2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        u_par = np.sum(u * b, axis=0) / bm
        return (np.sum(u ** 2, axis=0) - u_par ** 2) / (2.0 * bm)


def count_crossings(s, hysteresis: float):
    """Number of times a signed distance s (n_t, n) to the front changes side, with hysteresis.
    Returns the number of side changes per orbit (NaN samples are skipped)."""
    state = np.zeros(s.shape[1], dtype=np.int8)
    n = np.zeros(s.shape[1], dtype=np.int32)
    for row in s:
        new = np.where(row > hysteresis, 1, np.where(row < -hysteresis, -1, 0)).astype(np.int8)
        change = (new != 0) & (state != 0) & (new != state)
        n += change
        state = np.where(new != 0, new, state)
    return n


# ------------------------------------------------------------------ the measured front
class PolarMap:
    """Polar bins about a centre on a plane through it, for fields given on that plane.

    across, along: node coordinates of the first and the second array axis; the angle alpha is
    measured from +along towards +across (for the meridional plane (x, z): the polar angle from
    B0, continued to 360 degrees). The four quadrants are four samples of one polar angle
    theta in [0, 90] degrees: fold[j] is the theta bin of the angle bin j."""

    def __init__(self, across, along, centre, dr: float = 2.0, r_max: float = 190.0, dangle_deg: float = 5.0):
        da = np.asarray(across, float)[:, None] - centre[0]
        db = np.asarray(along, float)[None, :] - centre[1]
        rad = np.hypot(da, db)
        alpha = np.mod(np.arctan2(da, db), 2.0 * np.pi)
        self.n_r, self.n_a = int(round(r_max / dr)), int(round(360.0 / dangle_deg))
        if self.n_a % 4:
            raise ValueError("the angle bins must divide a quadrant")
        self.n_q = self.n_a // 4
        step = 2.0 * np.pi / self.n_a
        ir = np.floor(rad / dr).astype(int)
        ia = np.minimum(np.floor(alpha / step).astype(int), self.n_a - 1)
        self.index = np.where(ir < self.n_r, ir * self.n_a + ia, -1).ravel()
        self.inside = self.index >= 0
        self.r = (np.arange(self.n_r) + 0.5) * dr
        self.dr = dr
        self.alpha = (np.arange(self.n_a) + 0.5) * step
        self.theta = (np.arange(self.n_q) + 0.5) * step
        j = np.arange(self.n_a)
        self.quadrant = j // self.n_q
        self.fold = np.where(self.quadrant % 2 == 0, j % self.n_q, self.n_q - 1 - j % self.n_q)
        self.count = np.bincount(self.index[self.inside], minlength=self.n_r * self.n_a).reshape(self.n_r, self.n_a)

    def _sum(self, field):
        f = np.asarray(field, float).ravel()[self.inside]
        return np.bincount(self.index[self.inside], weights=f, minlength=self.n_r * self.n_a).reshape(
            self.n_r, self.n_a)

    def mean(self, field):
        """Bin averages (n_r, n_alpha); NaN in an empty bin."""
        with np.errstate(invalid="ignore", divide="ignore"):
            return self._sum(field) / self.count

    def quadrants(self, binned):
        """(n_r, n_alpha) -> (4, n_r, n_theta): the four quadrants, each on the theta bins."""
        out = np.empty((4, self.n_r, self.n_q))
        for q in range(4):
            cols = np.flatnonzero(self.quadrant == q)
            out[q][:, self.fold[cols]] = binned[:, cols]
        return out

    def fold_mean(self, field):
        """Bin averages over the four quadrants together, (n_r, n_theta)."""
        tot = self.quadrants(self._sum(field)).sum(axis=0)
        cnt = self.quadrants(self.count.astype(float)).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return tot / cnt

    def max(self, field):
        """Largest single-cell value per bin, (n_r, n_alpha); -inf in an empty bin."""
        out = np.full(self.n_r * self.n_a, -np.inf)
        np.maximum.at(out, self.index[self.inside], np.asarray(field, float).ravel()[self.inside])
        return out.reshape(self.n_r, self.n_a)


def outermost_crossing(profile, r, threshold: float, r_min: float = 16.0, smooth: int = 2):
    """Outermost radius at which a radial profile falls below a threshold.

    profile: (n_r, ...) on the uniform radii r; it is averaged over `smooth` radial bins first.
    The crossing has to have the averaged profile at or above the threshold in two neighbouring
    bins on its inner side, at radii >= r_min; its position is interpolated linearly. NaN where
    there is none, or where the profile is still above the threshold at the outer edge."""
    p = np.asarray(profile, float)
    p = np.where(np.isfinite(p), p, -np.inf)
    sm = sum(p[i:p.shape[0] - smooth + 1 + i] for i in range(smooth)) / smooth
    rs = sum(r[i:len(r) - smooth + 1 + i] for i in range(smooth)) / smooth
    above = sm >= threshold
    use = (rs >= r_min).reshape((-1,) + (1,) * (p.ndim - 1))
    hit = above[1:-1] & above[:-2] & ~above[2:] & use[1:-1]  # index k = 1 .. n - 2: above at k - 1, k; below at k + 1
    n = hit.shape[0]
    k = n - 1 - np.argmax(hit[::-1], axis=0) + 1
    found = hit.any(axis=0)
    idx = np.expand_dims(k, 0)
    s0 = np.take_along_axis(sm, idx, axis=0)[0]
    s1 = np.take_along_axis(sm, idx + 1, axis=0)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(np.isfinite(s1), (s0 - threshold) / (s0 - s1), 1.0)
    pos = rs[k] + np.clip(frac, 0.0, 1.0) * (rs[1] - rs[0])
    # a profile that is still above the threshold at the outer edge has no outer crossing
    return np.where(found & ~above[-1] & ~above[-2], pos, np.nan)


def front_kinematics(theta, t, radius, window: int = 5) -> dict:
    """Speed and geometry of a measured front, without the fit.

    radius: (n_t, n_theta) on the uniform bin centres theta in (0, pi/2), NaN where no front.
    R_t is the slope of a straight line through `window` outputs (NaN unless all are present);
    R_theta is a central difference, with the symmetry about the axis and the equator at the
    two ends and a one-sided difference beside a missing bin. sigma_v is the standard error of
    that slope (from the residuals of the line), carried to the normal speed.
    Returns dict(r_t, r_theta, v_n, sigma_v, cos_tbn, tbn, v_ht), each (n_t, n_theta)."""
    theta, t, radius = np.asarray(theta, float), np.asarray(t, float), np.asarray(radius, float)
    n_t = len(t)
    h = window // 2
    r_t = np.full(radius.shape, np.nan)
    r_t_err = np.full(radius.shape, np.nan)
    for i in range(h, n_t - h):
        tt = t[i - h:i + h + 1] - t[i]
        rr = radius[i - h:i + h + 1]
        dt_c, dr_c = (tt - tt.mean())[:, None], rr - rr.mean(axis=0)
        stt = np.sum(dt_c ** 2)
        r_t[i] = np.sum(dt_c * dr_c, axis=0) / stt
        r_t_err[i] = np.sqrt(np.sum((dr_c - r_t[i] * dt_c) ** 2, axis=0) / (window - 2) / stt)
    dth = theta[1] - theta[0]
    ext = np.concatenate([radius[:, :1], radius, radius[:, -1:]], axis=1)  # mirror images at both ends
    left, right = ext[:, :-2], ext[:, 2:]
    both = np.isfinite(left) & np.isfinite(right)
    r_theta = np.where(both, (right - left) / (2 * dth),
                       np.where(np.isfinite(right), (right - radius) / dth, (radius - left) / dth))
    den = r_theta * np.sin(theta) + radius * np.cos(theta)
    norm = np.sqrt(radius ** 2 + r_theta ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_tbn = np.abs(den) / norm
        v_ht = np.where(den > 0, r_t * radius / den, np.inf)
        v_ht = np.where(np.isfinite(r_t) & np.isfinite(den), v_ht, np.nan)
    return {"r_t": r_t, "r_theta": r_theta, "v_n": r_t * radius / norm, "sigma_v": r_t_err * radius / norm,
            "cos_tbn": cos_tbn, "tbn": np.arccos(np.clip(cos_tbn, 0.0, 1.0)), "v_ht": v_ht}


def vht_agree(ref, values, tol: float):
    """Do foot-point speeds agree? values: (k, ...) against the reference ref (...).

    True where the largest difference among the values and the reference is at most tol x |ref|,
    and where all of them are superluminal (>= 1: no reflection is possible by any of them, and
    the size of an infinite speed has no meaning). A missing value (NaN) never agrees."""
    ref = np.asarray(ref, float)
    allv = np.concatenate([ref[None], np.asarray(values, float)], axis=0)
    with np.errstate(invalid="ignore"):
        sup = (allv >= 1.0).all(axis=0)
        spread = allv.max(axis=0) - allv.min(axis=0)
        return sup | (np.isfinite(allv).all(axis=0) & (spread <= tol * np.abs(ref)))


def outputs_in(duration: float, dt_out: float) -> int:
    """Number of outputs in a row that stands for a duration: duration / dt_out rounded up, where a
    tenth of an output interval is not counted (2 / Omega_ci: 6 in L, 11 in M, 15 in H)."""
    return int(np.ceil(duration / dt_out - 0.1))


def _in_a_row(mask, n: int):
    """True where a True of mask (n_t, ...) belongs to at least n True values in a row in time."""
    m = np.asarray(mask, bool)
    run = np.zeros(m.shape, int)
    acc = np.zeros(m.shape[1:], int)
    for i in range(m.shape[0]):
        acc = np.where(m[i], acc + 1, 0)
        run[i] = acc
    out = np.zeros(m.shape, bool)
    acc = np.zeros(m.shape[1:], int)
    for i in range(m.shape[0] - 1, -1, -1):  # carry the length of each run back to its start
        acc = np.where(m[i], np.maximum(acc, run[i]), 0)
        out[i] = acc >= n
    return out


def valid_domain(radius, r_debris, r_bmax, v_n, sigma_v, tbn, t, v_a, c_s, omega_ci, r_outer: float,
                 lead: float = 0.7, mach: float = 1.5, mach_low: float = 1.1, age: float = 2.0,
                 omega_row: float | None = None) -> dict:
    """The valid shock domain of the protocol (Sect. 4) on the (time, theta) bins.

    c_s: the upstream sound speed, a number or one value per output (n_t,). Returns the masks
    a (front exists), b (the debris edge lies behind the position of B_max and the front leads it
    by lead x V_n / Omega_ci; true where no debris edge is found), c (V_n >= mach x v_f at the
    local theta_Bn and V_n - 2 sigma_V > mach_low x v_f), d (a-c in a row over age / Omega_ci),
    e (t >= age / Omega_ci and R <= r_outer - V_n / Omega_ci), valid (d and e), and n_row."""
    t = np.asarray(t, float)
    c_s = np.broadcast_to(np.asarray(c_s, float), t.shape)[:, None]
    # v_a and omega_ci may be arrays on the (time, theta) bins (a dipole: the local field at the front); the
    # persistence then counts outputs with omega_row, the ion gyro-frequency at the reference radius of the design
    n_row = outputs_in(age / (omega_row if omega_row is not None else omega_ci), float(np.median(np.diff(t))))
    a = np.isfinite(radius)
    with np.errstate(invalid="ignore"):
        v_f = fast_speed(v_a, c_s, tbn)
        b = a & (~np.isfinite(r_debris) | ((r_debris < r_bmax) & (radius - r_debris >= lead * v_n / omega_ci)))
        c = a & (v_n >= mach * v_f) & (v_n - 2.0 * sigma_v > mach_low * v_f)
        e = a & (t[:, None] >= age / omega_ci) & (radius <= r_outer - v_n / omega_ci)
    d = _in_a_row(a & b & c, n_row)
    return {"a": a, "b": b, "c": c, "d": d, "e": e, "valid": d & e, "n_row": n_row, "v_f": v_f}


def fit_temperature(hist, edges, t_guess: float, lo: float = 0.5, hi: float = 5.0):
    """Temperature of the bulk from a histogram of the kinetic energy (m = 1, energies in m c^2).

    hist: (..., n_bins) counts in the bins `edges`. For a Maxwell-Juttner distribution
    dN/dE is proportional to gamma u exp(-E / T), so ln(hist / (dE gamma u)) is a straight line
    in E with the slope -1/T. Weighted least squares (weight = count) over the bins with centres
    in [lo, hi] x t_guess. Returns (T, number of electrons used); NaN with fewer than 3 bins."""
    edges = np.asarray(edges, float)
    cen = np.sqrt(edges[1:] * edges[:-1])
    gam = 1.0 + cen
    sel = (cen >= lo * t_guess) & (cen <= hi * t_guess)
    h = np.asarray(hist, float)[..., sel]
    x = cen[sel]
    w = np.where(h > 0, h, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.where(h > 0, np.log(h / (np.diff(edges)[sel] * gam[sel] * np.sqrt(gam[sel] ** 2 - 1.0))), 0.0)
    sw = w.sum(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        xm = (w * x).sum(axis=-1) / sw
        ym = (w * y).sum(axis=-1) / sw
        sxx = (w * (x - xm[..., None]) ** 2).sum(axis=-1)
        sxy = (w * (x - xm[..., None]) * (y - ym[..., None])).sum(axis=-1)
        temp = -sxx / sxy
    ok = ((h > 0).sum(axis=-1) >= 3) & (sxy < 0)
    return np.where(ok, temp, np.nan), sw


def dipole_field(x, y, z, moment, b0=0.0):
    """Field of a point dipole moment (along z, at the origin) plus a uniform b0 along z, at (x, y, z)."""
    r2 = x * x + y * y + z * z
    r5 = np.maximum(r2, 1e-30) ** 2.5
    bx = 3.0 * moment * x * z / r5
    by = 3.0 * moment * y * z / r5
    bz = moment * (3.0 * z * z - r2) / r5 + b0
    return bx, by, bz


def dipole_tbn(theta, radius, r_theta):
    """theta_Bn between the outward normal of a surface r(theta) and a dipole field along z (both hemispheres
    folded): the field direction in the (r, theta) basis is (2 cos theta, sin theta). Returns (cos_tbn, tbn)."""
    th = np.broadcast_to(np.asarray(theta, float), np.shape(radius))
    with np.errstate(invalid="ignore", divide="ignore"):
        k = np.asarray(r_theta, float) / np.asarray(radius, float)
        num = np.abs(2.0 * np.cos(th) - k * np.sin(th))
        cos_tbn = num / (np.sqrt(4.0 * np.cos(th) ** 2 + np.sin(th) ** 2) * np.sqrt(1.0 + k * k))
    return cos_tbn, np.arccos(np.clip(cos_tbn, 0.0, 1.0))
