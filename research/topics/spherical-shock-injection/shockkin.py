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


def fit_front(theta, t, radius, order_t: int = 3):
    """Least-squares fit of measured front radii by {1, cos 2 theta, cos 4 theta} x {1..tau^order_t}.

    theta, t, radius: 1-D arrays of the valid points (theta folded into [0, pi/2] by the caller
    or not: the basis is symmetric). Returns (front, rms residual, (t_min, t_max))."""
    theta, t, radius = (np.asarray(a, float) for a in (theta, t, radius))
    t0, t1 = float(t.min()), float(t.max())
    span = max(t1 - t0, 1e-30)
    a, _, _ = _basis(theta, (t - t0) / span, order_t)
    coef, *_ = np.linalg.lstsq(a, radius, rcond=None)
    rms = float(np.sqrt(np.mean((a @ coef - radius) ** 2)))

    def front(th, tt):
        th, tt = np.broadcast_arrays(np.asarray(th, float), np.asarray(tt, float))
        b, bt, bth = _basis(th, (tt - t0) / span, order_t)
        return b @ coef, (bt @ coef) / span, bth @ coef

    return front, rms, (t0, t1)


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


def s0_test(u_par_r, u_perp, v_ht, b_max: float, u_thermal: float, k_thermal: float = 5.0):
    """Is a reflected electron explained by one reflection of a thermal electron (S0)?

    It is when the momentum before the reflection lies within k_thermal thermal momenta and the
    electron was outside the loss cone for the largest compression b_max (phi = 0, the most
    permissive case). Returns (explained, |u| before the reflection)."""
    u_par0, gam0, gam_ht = before_reflection(u_par_r, u_perp, v_ht)
    u0 = np.sqrt(np.maximum(gam0 ** 2 - 1.0, 0.0))
    with np.errstate(invalid="ignore"):
        mirror = gam_ht ** 2 - 1.0 - b_max * np.asarray(u_perp, float) ** 2 < 0.0
        explained = mirror & (u0 <= k_thermal * u_thermal) & (np.asarray(v_ht, float) < 1.0)
    return explained, u0


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


def work_closure(u, e, t, charge: float = -1.0, dgamma_min: float = 0.1):
    """Compare the work of the recorded field with the recorded change of energy, per orbit.

    u, e: (3, n_t, n); u is recorded half a step before e, so the velocity is averaged over
    neighbouring samples before it is multiplied with e. Returns dict(work, dgamma, n_used,
    median_rel_error) over the orbits with |dgamma| >= dgamma_min and no missing sample."""
    gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    v = u / gam
    p = charge * np.sum(e * v, axis=0)  # d gamma / dt
    dt = np.diff(t)[:, None]
    work = np.sum(0.5 * (p[1:] + p[:-1]) * dt, axis=0)
    dgam = gam[-1] - gam[0]
    ok = np.isfinite(work) & np.isfinite(dgam) & (np.abs(dgam) >= dgamma_min)
    rel = np.abs(work[ok] - dgam[ok]) / np.abs(dgam[ok])
    return {"work": work, "dgamma": dgam, "n_used": int(ok.sum()),
            "median_rel_error": float(np.median(rel)) if ok.any() else float("nan")}


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
