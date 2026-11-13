"""Census of the energetic electrons of a run (exploratory-protocol.md Sect. 3.2, 3.3, 5 N3b, 6, 7 K1).

Every electron of the selected dump `hot_e` gets an origin, a position class against the fitted
front, a band of the foot-point theta_Bn, its flux-surface displacement and, when it is upstream,
the result of the test B1' (one reflection of a thermal electron, S0), of the leakage test S1 and of
the trapping test S2. Lengths in d_e, times in 1/omega_pe, c = 1, m_e = 1; z is the direction of B0.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import shockkin as sk

E_EDGES = np.array([1.25, 1.5, 2.0, 3.0, 5.0, np.inf])  # gamma; the first bin only serves histories and N3b
E_NAMES = ("1.25-1.5", "1.5-2", "2-3", "3-5", ">5")
ORIGINS = ("ambient, born inside the shell", "ambient, swept by the front", "ambient, not reached",
           "injected at a face", "other", "debris")
CLASSES = ("upstream", "upstream, line not reached", "layer", "downstream", "no front")
RESULTS = ("not upstream", "S0 at exit", "S0 in range", "S1 leak", "S2 trapped", "undetermined", "residual",
           "residual, no encounter")
FLAGS = ("no front in the bin", "front detected", "valid domain", "usable for B1'")
UND = ("crossing may precede the fit", "no mid-ramp front at the foot", "|v_par - V_HT| < 0.1",
       "flight time above 0.1 / nu", "3-10 gyro-periods in the layer", "no mid-ramp front at the dump")
NO_FOOT = 9  # band index of an electron without a foot point
SHELL_RADIUS = 27.5  # outer edge of the initial debris shell [d_e]
N_RANGE = 16  # values of V_HT tried between the exit value and the value one dwell time earlier


def _band(tbn):
    """10-degree band of theta_Bn (radians); NO_FOOT where it is not defined."""
    deg = np.degrees(np.asarray(tbn, float))
    return np.where(np.isfinite(deg), np.clip(np.nan_to_num(deg) // 10, 0, 8), NO_FOOT).astype(np.int8)


def e_bin(gamma):
    """Index of the energy bin; -1 below the first edge."""
    return (np.searchsorted(E_EDGES, np.asarray(gamma, float), side="right") - 1).astype(np.int8)


class FrontModel:
    """The measured front of a run in the form the census needs."""

    def __init__(self, front, front_b, t_range, t_range_b, t, theta, detected, valid, usable, b_mean, b_cell,
                 has_b, t_up, temp_up, rho_debris, d_i, b0, centre, n_cells, dx, ppc, rho_max=200.0):
        self.front, self.front_b = front, front_b
        self.t_range, self.t_range_b = t_range, t_range_b
        self.t, self.theta = np.asarray(t, float), np.asarray(theta, float)
        self.dt, self.dtheta = self.t[1] - self.t[0], self.theta[1] - self.theta[0]
        self.dt_dump = self.dt  # spacing of the particle dumps; set by cmd_census when it is a multiple of dt
        self.flag_map = (np.asarray(detected, int) + np.asarray(valid, int) + np.asarray(usable, int)).astype(np.int8)
        self.b_mean, self.b_cell, self.has_b = (np.asarray(a) for a in (b_mean, b_cell, has_b))
        det = np.asarray(detected, bool)
        self.b_max = float(np.nanmax(np.where(det, self.b_mean, np.nan)))
        self.b_cell_max = float(np.nanmax(np.where(det, self.b_cell, np.nan)))
        self.t_up, self.temp_up = np.asarray(t_up, float), np.asarray(temp_up, float)
        self.rho_debris_t = np.asarray(rho_debris, float)
        self.d_i, self.b0, self.centre = float(d_i), float(b0), np.asarray(centre, float)
        self.n_cells, self.dx, self.ppc = tuple(int(n) for n in n_cells), tuple(float(d) for d in dx), int(ppc)
        self.n_total = float(self.ppc * np.prod(self.n_cells))
        self.tau_dwell = 3.0 * 2.0 * np.pi / self.b0
        step = self.dt / 8.0  # the dumps fall on the table times

        def table(fr, rng):
            if fr is None:
                return None
            n = int(np.floor((rng[1] - rng[0]) / step + 1e-9)) + 1
            return sk.FootTable(fr, rng[0] + step * np.arange(n), rho_max)
        self.tab, self.tab_b = table(front, t_range), table(front_b, t_range_b)

    @classmethod
    def from_run(cls, run, run_dir):
        z = np.load(Path(run_dir) / "analysis" / "front.npz")
        j = json.loads((Path(run_dir) / "analysis" / "front.json").read_text())
        plan = run.meta["planner"]
        fit = sk.FittedFront.from_dict(j["fit"]) if j.get("fit") else None
        fit_b = sk.FittedFront.from_dict(j["fit_B"]) if j.get("fit_B") else None
        with np.errstate(invalid="ignore"):
            deb = z["r_debris_thr"][:, 0, :] * np.sin(z["theta"])[None, :]
            rho_deb = np.where(np.isfinite(deb).any(axis=1), np.nanmax(np.nan_to_num(deb, nan=0.0), axis=1), 0.0)
            has_b = np.isfinite(z["r_b"]) & (z["b_mean"] >= 1.5)
        up = j["upstream"]
        return cls(fit, fit_b, fit.t_range if fit else None, fit_b.t_range if fit_b else None, z["t"], z["theta"],
                   np.isfinite(z["r"]), z["valid"], z["usable"], z["b_mean"], z["b_cell"], has_b, up["t"], up["T_e"],
                   rho_deb, plan["d_i"], plan["b0"], plan["center"], run.config["n_cells"], run.config["dx"],
                   run.config["ppc"])

    def _bin(self, theta, t):
        ok = np.isfinite(theta) & np.isfinite(t)
        it = np.clip(np.rint((np.where(ok, t, self.t[0]) - self.t[0]) / self.dt), 0, len(self.t) - 1).astype(int)
        ith = np.clip(np.where(ok, theta, 0.0) / self.dtheta, 0, len(self.theta) - 1).astype(int)
        return ok, it, ith

    def flag(self, theta, t):
        """0 no front in the (theta, t) bin, 1 detected, 2 in the valid domain, 3 usable for B1'."""
        theta, t = np.broadcast_arrays(np.asarray(theta, float), np.asarray(t, float))
        ok, it, ith = self._bin(theta, t)
        return np.where(ok, self.flag_map[it, ith], 0).astype(np.int8)

    def local(self, theta, t):
        """(mid-ramp front exists, bin-average compression) in the nearest (theta, t) bin."""
        theta, t = np.broadcast_arrays(np.asarray(theta, float), np.asarray(t, float))
        ok, it, ith = self._bin(theta, t)
        return ok & self.has_b[it, ith], np.where(ok, self.b_mean[it, ith], np.nan)

    def temp_upstream(self, t):
        return float(np.interp(t, self.t_up, self.temp_up))

    def rho_debris(self, t):
        return float(np.interp(t, self.t, self.rho_debris_t))


def origin_of(model, ids, t):
    """Origin code and the birth position (rho, z from the centre; NaN unless ambient) of the ids."""
    n = np.prod(model.n_cells)
    ids = np.asarray(ids, dtype=np.uint64)
    amb = ids < np.uint64(model.ppc * n)
    ix, iy, iz = sk.birth_cell(np.where(amb, ids, 0), model.n_cells, model.ppc)
    xb, yb, zb = ((i + 0.5) * d - c for i, d, c in zip((ix, iy, iz), model.dx, model.centre))
    rho_b, r_b = np.hypot(xb, yb), np.sqrt(xb ** 2 + yb ** 2 + zb ** 2)
    if model.front is not None and model.t_range[0] <= t <= model.t_range[1]:
        swept = r_b <= model.front(np.arctan2(rho_b, np.abs(zb)), t)[0]
    else:
        swept = np.zeros(ids.shape, bool)
    org = np.where(r_b <= SHELL_RADIUS, 0, np.where(swept, 1, 2))
    org = np.where(amb, org, np.where(ids >= np.uint64(20 * n), 3, 4)).astype(np.int8)
    return org, np.where(amb, rho_b, np.nan), np.where(amb, zb, np.nan)


def classify(model, part, t, e_th):
    """All per-electron quantities of one dump that need no other dump. Returns a dict of arrays."""
    x, y, z = (np.asarray(part[c], float) - model.centre[i] for i, c in enumerate("xyz"))
    ux, uy, uz = (np.asarray(part[c], float) for c in ("ux", "uy", "uz"))
    n = len(x)
    gam = np.sqrt(1.0 + ux ** 2 + uy ** 2 + uz ** 2)
    xg, yg = sk.guiding_centre(x, y, ux, uy, model.b0)
    rho, z_abs = np.hypot(xg, yg), np.abs(z)
    theta = np.arctan2(rho, z_abs)
    org, rho_b, _ = origin_of(model, part["id"], t)
    rho_p = sk.canonical_radius(x, y, ux, uy, model.b0)
    rec = {"id": np.asarray(part["id"], np.uint64), "origin": org, "gamma": gam, "rho_gc": rho, "z_abs": z_abs,
           "hemi": (z >= 0).astype(np.int8), "quad": ((xg >= 0) + 2 * (yg >= 0)).astype(np.int8),
           "mu": uz / np.sqrt(ux ** 2 + uy ** 2 + uz ** 2), "rho_p": rho_p, "drho": rho_p - rho_b}
    nan = np.full(n, np.nan)
    for k in ("t_r", "v_exit", "v_hi", "e_before", "e_before_min", "gam_ht", "tbn_r"):
        rec[k] = nan.copy()
    rec["explained"], rec["und"] = np.zeros(n, np.uint8), np.zeros(n, np.uint8)
    rec["mirror_local"], rec["luminal"] = np.zeros(n, bool), np.zeros(n, bool)
    if model.front is None or not (model.t_range[0] <= t <= model.t_range[1]):
        rec.update(cls=np.full(n, 4, np.int8), delta=nan, band=np.full(n, NO_FOOT, np.int8), flag=np.zeros(n, np.int8))
        return rec
    delta = np.hypot(rho, z_abs) - model.front(theta, t)[0]  # radial distance ahead of the front
    z_foot, v_foot, tbn_foot = model.tab.at(rho, t)
    cls = np.where(delta > model.d_i, np.where(np.isfinite(z_foot), 0, 1), np.where(delta >= -3 * model.d_i, 2, 3))
    here = sk.foot(model.front, theta, t)  # the front along the electron's own ray (layer, downstream)
    up = cls <= 1
    tbn = np.where(up, tbn_foot, here["tbn"])
    v_at = np.where(up, v_foot, here["v_ht"])
    flag = model.flag(theta, t)
    tb = model.tab_b
    if tb is None or not (tb.t[0] < t <= tb.t[-1] + 1e-9):
        rec["und"][up] |= 1 << 5
    elif up.any():
        s = np.flatnonzero(up)
        u_par, u_perp = np.sign(z[s]) * uz[s], np.hypot(ux[s], uy[s])
        u_par = np.where(z[s] == 0, np.abs(uz[s]), u_par)
        v_par = u_par / gam[s]
        t_r, v_exit, tbn_r = tb.trace_back(rho[s], z_abs[s], v_par, t)
        met = np.isfinite(t_r)
        tr0 = np.where(met, t_r, t)
        samples = [np.where(met, v_exit, np.nan)]
        for k in range(1, 5):
            v_k = tb.at(rho[s], tr0 - model.tau_dwell * k / 4)[1]
            samples.append(np.where(np.isfinite(v_k), v_k, 1.0))  # the line was touched later: up to c
        samples = np.minimum(np.array(samples), 0.999)
        lo, hi = samples.min(axis=0), samples.max(axis=0)
        v_rng = lo[:, None] + (hi - lo)[:, None] * np.linspace(0.0, 1.0, N_RANGE)[None, :]
        theta_f = np.arctan2(rho[s], tb.at(rho[s], tr0)[0])  # polar angle of the foot at the exit (NaN: none)
        exists, b_loc = model.local(theta_f, tr0)
        b_loc = np.where(np.isfinite(b_loc), np.maximum(b_loc, 1.0), 1.0)
        bits = np.zeros(len(s), np.uint8)
        e_min = np.full(len(s), np.nan)
        for j, b in enumerate((model.b_max, b_loc, model.b_cell_max)):
            b_arr = np.broadcast_to(np.asarray(b, float), (len(s),))
            ex0, e0, mir0 = sk.s0_test(u_par, u_perp, np.where(met, v_exit, np.nan), b_arr, e_th)
            ex, e_all, mir = sk.s0_test(u_par[:, None], u_perp[:, None], v_rng, b_arr[:, None], e_th)
            bits |= ((ex0 & met).astype(np.uint8) << (2 * j)) | ((ex.any(axis=1) & met).astype(np.uint8) << (2 * j + 1))
            if j == 0:
                with np.errstate(invalid="ignore"):
                    e_min = np.where(mir.any(axis=1) & met, np.nanmin(np.where(mir, e_all, np.inf), axis=1), np.nan)
                rec["e_before"][s] = np.where(met, e0, np.nan)
            if j == 1:
                rec["mirror_local"][s] = mir0 & met
        rec["explained"][s] = bits
        rec["e_before_min"][s] = e_min
        rec["gam_ht"][s] = sk.before_reflection(u_par, u_perp, np.where(met, v_exit, np.nan))[2]
        rec["t_r"][s], rec["v_exit"][s], rec["tbn_r"][s] = t_r, v_exit, tbn_r
        rec["v_hi"][s] = np.where(met, hi, np.nan)
        z_b_now = tb.at(rho[s], t)[0]
        early = ~met & (np.isfinite(tb.z[0][np.clip(np.searchsorted(tb.rho, rho[s]), 0, len(tb.rho) - 1)])
                        | (z_abs[s] < np.nan_to_num(z_b_now, nan=-np.inf)))
        und = early.astype(np.uint8) | ((met & ~exists).astype(np.uint8) << 1)
        with np.errstate(invalid="ignore"):
            und |= (met & (np.abs(v_par - v_exit) < 0.1)).astype(np.uint8) << 2
        rec["und"][s] = und
        tbn[s] = np.where(met, tbn_r, tbn[s])
        v_at[s] = np.where(met, v_exit, v_at[s])
        flag[s] = np.where(met, model.flag(theta_f, tr0), flag[s])
    with np.errstate(invalid="ignore"):
        rec["luminal"] = np.nan_to_num(v_at, nan=0.0) >= 1.0
    rec.update(cls=cls.astype(np.int8), delta=delta, band=_band(tbn), flag=flag)
    return rec


class History:
    """What the earlier dumps say about each electron, by id."""

    FIELDS = {"first_dump": np.int16, "first_cls": np.int8, "first_flag": np.int8, "ever_down": bool,
              "run": np.int16, "run_max": np.int16, "last": np.int16}

    def __init__(self):
        self.ids = np.empty(0, np.uint64)
        self.f = {k: np.empty(0, d) for k, d in self.FIELDS.items()}

    def update(self, k, ids, cls, flag):
        """Take dump k in. Returns per electron: first_dump, first_cls, first_flag, ever_down (in an
        earlier dump) and run_max (largest number of successive dumps in the layer class so far)."""
        n = len(ids)
        if len(np.unique(ids)) != n:
            raise ValueError("ids of a dump are not unique")
        out = {"first_dump": np.full(n, k, np.int16), "first_cls": cls.astype(np.int8), "first_flag":
               flag.astype(np.int8), "ever_down": np.zeros(n, bool), "run_max": (cls == 2).astype(np.int16)}
        found = np.zeros(n, bool)
        if len(self.ids):
            idx = np.minimum(np.searchsorted(self.ids, ids), len(self.ids) - 1)
            found = self.ids[idx] == ids
            j = idx[found]
            for key in ("first_dump", "first_cls", "first_flag", "ever_down"):
                out[key][found] = self.f[key][j]
            run = np.where(cls[found] == 2, np.where(self.f["last"][j] == k - 1, self.f["run"][j], 0) + 1, 0)
            self.f["run"][j] = run
            self.f["run_max"][j] = np.maximum(self.f["run_max"][j], run)
            self.f["ever_down"][j] |= cls[found] == 3
            self.f["last"][j] = k
            out["run_max"][found] = self.f["run_max"][j]
        new = ~found
        add = {"first_dump": out["first_dump"][new], "first_cls": out["first_cls"][new], "first_flag":
               out["first_flag"][new], "ever_down": cls[new] == 3, "run": (cls[new] == 2).astype(np.int16),
               "run_max": (cls[new] == 2).astype(np.int16), "last": np.full(int(new.sum()), k, np.int16)}
        all_ids = np.concatenate([self.ids, ids[new]])
        order = np.argsort(all_ids, kind="stable")
        self.ids = all_ids[order]
        self.f = {key: np.concatenate([self.f[key], add[key].astype(d)])[order] for key, d in self.FIELDS.items()}
        return out


N3B_KEYS = ("n", "n_mu", "sum_q", "sum_dgamma", "sum_dgamma2", "sum_drho", "sum_drho2")


def _pairs(prev, rec, d_ahead):
    """N3b: changes between two successive dumps of the electrons that stayed upstream, far enough ahead
    of the front (d_ahead) that they cannot have met it in between. Sums per energy bin."""
    acc = np.zeros((len(E_EDGES) - 1, len(N3B_KEYS)))
    if prev is None or not len(prev["id"]) or not len(rec["id"]):
        return acc
    _, i0, i1 = np.intersect1d(prev["id"], rec["id"], assume_unique=True, return_indices=True)
    with np.errstate(invalid="ignore"):
        ok = (prev["cls"][i0] <= 1) & (rec["cls"][i1] <= 1) & (prev["delta"][i0] > d_ahead)
    i0, i1 = i0[ok], i1[ok]
    mid = 0.5 * (prev["mu"][i0] + rec["mu"][i1])
    with np.errstate(divide="ignore", invalid="ignore"):
        q = (rec["mu"][i1] - prev["mu"][i0]) ** 2 / (1.0 - mid ** 2)
    use = np.isfinite(q) & (np.abs(mid) < 0.95)
    dg, dr = rec["gamma"][i1] - prev["gamma"][i0], rec["rho_p"][i1] - prev["rho_p"][i0]
    eb = e_bin(0.5 * (prev["gamma"][i0] + rec["gamma"][i1]))
    for b in range(acc.shape[0]):
        m = eb == b
        acc[b] = (m.sum(), (m & use).sum(), q[m & use].sum(), dg[m].sum(), (dg[m] ** 2).sum(), dr[m].sum(),
                  (dr[m] ** 2).sum())
    return acc


def run_census(model, dumps, temp_down_of=None, trace=True):
    """Classify the electrons of successive dumps. dumps: iterable of (t, part) in time order; part
    holds x, y, z, ux, uy, uz, id. Returns (records: dict of arrays over all electrons, n3b sums
    (dump, energy bin, N3B_KEYS))."""
    hist, prev, recs, acc, t_prev = History(), None, [], [], None
    for k, (t, part) in enumerate(dumps):
        e_th = sk.thermal_bound(model.temp_upstream(t), model.n_total)
        rec = classify(model, part, t, e_th) if trace else classify_only(model, part, t)
        rec.update(hist.update(k, rec["id"], rec["cls"], rec["flag"]))
        n = len(rec["id"])
        rec["dump"], rec["t"], rec["e_th"] = np.full(n, k, np.int16), np.full(n, float(t)), np.full(n, e_th)
        t_dn = float(temp_down_of(t)) if temp_down_of is not None else np.nan
        bound = sk.thermal_bound(t_dn, model.n_total) if np.isfinite(t_dn) else np.nan
        with np.errstate(invalid="ignore"):
            rec["s1"] = ((rec["cls"] <= 1) & ((rec["explained"] & 3) == 0) & (rec["v_exit"] < 1.0)
                         & ~rec["mirror_local"] & rec["ever_down"] & (rec["gam_ht"] - 1.0 <= bound))
        d_ahead = model.d_i + 1.25 * (t - t_prev) if t_prev is not None else np.inf
        acc.append(_pairs(prev, rec, d_ahead))
        prev, t_prev = rec, t
        recs.append(rec)
    keys = recs[0].keys() if recs else ()
    return {key: np.concatenate([r[key] for r in recs]) for key in keys}, np.array(acc)


def classify_only(model, part, t):
    """classify without the back-tracing (for the debris electrons, which are counted by weight)."""
    tab_b, model.tab_b = model.tab_b, None
    try:
        return classify(model, part, t, np.nan)
    finally:
        model.tab_b = tab_b


def dump_spacing(times, slice_t, rtol=1e-6):
    """The spacing of the particle dumps. They must be evenly spaced and fall on the times of the slices (so the
    spacing is a whole multiple of that of the slices); otherwise SystemExit."""
    times, slice_t = np.asarray(times, float), np.asarray(slice_t, float)
    dt = slice_t[1] - slice_t[0]
    if len(times) < 2:
        raise SystemExit("fewer than two particle dumps")
    step = times[1] - times[0]
    k = np.rint((times - slice_t[0]) / dt)
    on_slices = np.allclose(times, slice_t[0] + k * dt, rtol=0, atol=rtol * dt) and k.max() < len(slice_t)
    even = np.allclose(np.diff(times), step, rtol=0, atol=rtol * dt) and step > 0
    if not (on_slices and even):
        raise SystemExit("the particle dumps are not evenly spaced on the times of the slices")
    return float(step)


def scatter_rates(acc, dt_dump, nu_fallback, n_min=100):
    """Pool the N3b sums over the dumps. Returns (rows per energy bin, nu used per energy bin)."""
    tot = acc.sum(axis=0) if len(acc) else np.zeros((len(E_EDGES) - 1, len(N3B_KEYS)))
    rows, nu_use = [], []
    for b, name in enumerate(E_NAMES):
        n, n_mu, s_q, s_g, s_g2, s_r, s_r2 = tot[b]
        nu = s_q / n_mu / dt_dump if n_mu > 0 else np.nan
        measured = n_mu >= n_min and np.isfinite(nu) and nu > 0
        nu_use.append(nu if measured else nu_fallback)
        rows.append({"gamma": name, "pairs": int(n), "pairs_mu": int(n_mu), "nu": float(nu), "nu_used": nu_use[-1],
                     "nu_from": "this run" if measured else "null run", "flight_limit": 0.1 / nu_use[-1],
                     "dgamma_mean": float(s_g / n) if n else None,
                     "energy_diffusion": float(s_g2 / n / (2 * dt_dump)) if n else None,
                     "drho_mean": float(s_r / n) if n else None, "drho_rms": float(np.sqrt(s_r2 / n)) if n else None})
    return rows, np.array(nu_use)


def assign(rec, model, nu_use, p999_of=None):
    """The result of every electron (RESULTS) and the K1 flags; adds result, flight, k1, beyond_debris."""
    up = rec["cls"] <= 1
    period = 2.0 * np.pi / model.b0
    n_s2, n_und = int(np.ceil(10 * period / model.dt_dump)), int(np.ceil(3 * period / model.dt_dump))
    flight = rec["t"] - rec["t_r"]
    limit = 0.1 / np.asarray(nu_use)[np.clip(e_bin(rec["gamma"]), 0, len(nu_use) - 1)]
    und = rec["und"].copy()
    with np.errstate(invalid="ignore"):
        und |= (up & (flight > limit)).astype(np.uint8) << 3
    hard = (und & 0b101111) != 0
    s0x, s0r = (rec["explained"] & 1) != 0, (rec["explained"] & 2) != 0
    wait = up & ~hard & ~s0x & ~s0r & ~rec["s1"] & (rec["run_max"] >= n_und) & (rec["run_max"] < n_s2)
    und |= wait.astype(np.uint8) << 4
    rec["result"] = np.select([~up, hard, s0x, s0r, rec["s1"], rec["run_max"] >= n_s2, wait, ~np.isfinite(rec["t_r"])],
                              [0, 5, 1, 2, 3, 4, 5, 7], default=6).astype(np.int8)
    rec["und"], rec["flight"] = und, flight
    times = np.unique(rec["t"])
    thr = np.array([max(2.0, float(p999_of(t)) if p999_of is not None else 0.0) for t in times])
    deb = np.array([model.rho_debris(t) for t in times])
    it = np.searchsorted(times, rec["t"])
    with np.errstate(invalid="ignore"):
        rec["k1"] = (rec["origin"] <= 2) & up & (rec["gamma"] >= 1.5) & (np.abs(rec["drho"]) > thr[it])
        rec["beyond_debris"] = rec["rho_p"] > deb[it]
    rec["k1_threshold"] = thr[it]
    return rec


def tables(rec, n_dumps):
    """Counts: all[dump, origin, class, energy] and up[dump, flag, band, energy, result, hemisphere, quadrant]."""
    eb = e_bin(rec["gamma"])
    ok = eb >= 0
    all_t = np.zeros((n_dumps, len(ORIGINS), len(CLASSES), len(E_NAMES)), np.int64)
    np.add.at(all_t, (rec["dump"][ok], rec["origin"][ok], rec["cls"][ok], eb[ok]), 1)
    up = ok & (rec["cls"] <= 1)
    up_t = np.zeros((n_dumps, len(FLAGS), NO_FOOT + 1, len(E_NAMES), len(RESULTS), 2, 4), np.int64)
    np.add.at(up_t, (rec["dump"][up], rec["flag"][up], rec["band"][up], eb[up], rec["result"][up], rec["hemi"][up],
                     rec["quad"][up]), 1)
    return all_t, up_t


# ------------------------------------------------------------------ downstream reference (Sect. 6)
def fit_t(hist, edges, t_guess, lo=0.5, hi=5.0):
    """sk.fit_temperature with one guess per histogram: hist (..., n_bins), t_guess (...)."""
    edges = np.asarray(edges, float)
    cen = np.sqrt(edges[1:] * edges[:-1])
    gam = 1.0 + cen
    h = np.asarray(hist, float)
    tg = np.asarray(t_guess, float)[..., None]
    with np.errstate(invalid="ignore"):
        w = np.where((h > 0) & (cen >= lo * tg) & (cen <= hi * tg), h, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.where(w > 0, np.log(h / (np.diff(edges) * gam * np.sqrt(gam ** 2 - 1.0))), 0.0)
        sw = w.sum(axis=-1)
        xm, ym = (w * cen).sum(axis=-1) / sw, (w * y).sum(axis=-1) / sw
        sxx = (w * (cen - xm[..., None]) ** 2).sum(axis=-1)
        sxy = (w * (cen - xm[..., None]) * (y - ym[..., None])).sum(axis=-1)
        temp = -sxx / sxy
        ok = ((w > 0).sum(axis=-1) >= 3) & (sxy < 0)
    return np.where(ok, temp, np.nan), sw


def bulk_temperature(hist, edges, n_iter=3):
    """Temperature of the bulk of each histogram: three fits, starting from the mean energy / 1.5."""
    cen = np.sqrt(np.asarray(edges[1:], float) * np.asarray(edges[:-1], float))
    h = np.asarray(hist, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        temp = (h * cen).sum(axis=-1) / h.sum(axis=-1) / 1.5
    for _ in range(n_iter):
        new, _ = fit_t(h, edges, temp)
        temp = np.where(np.isfinite(new), new, temp)
    return temp


class TailTable:
    """Fraction of a Maxwell-Juttner distribution above a kinetic energy, against the temperature."""

    def __init__(self, e_cut, t_lo=0.002, t_hi=1.5, n=80):
        self.lt = np.log(np.geomspace(t_lo, t_hi, n))
        self.lf = np.log(np.maximum([sk.mj_tail_fraction(float(np.exp(v)), 1.0 + e_cut) for v in self.lt], 1e-300))

    def __call__(self, temp):
        return np.exp(np.interp(np.log(np.clip(temp, np.exp(self.lt[0]), np.exp(self.lt[-1]))), self.lt, self.lf))


def downstream_cells(hist, edges, temp, tail, e_cut_index):
    """Measured and expected (bulk fit at 0.9, 1.0, 1.1 T) weight above the energy edge e_cut_index."""
    n = np.asarray(hist, float).sum(axis=-1)
    measured = np.asarray(hist, float)[..., e_cut_index:].sum(axis=-1)
    return measured, [n * tail(f * temp) for f in (0.9, 1.0, 1.1)]


def downstream_reference(run, model, r_debris=None, e_cut=0.5):
    """Sect. 6: energetic electrons behind the front against the tail of the local bulk temperature,
    from the meridional slab spectra. Returns (rows per output, function T_down(t))."""
    s = run.series("spec_meridional")
    x_e, z_e, ek = (np.asarray(e, float) for e in s.edges)
    xc = 0.5 * (x_e[1:] + x_e[:-1]) - model.centre[0]
    zc = 0.5 * (z_e[1:] + z_e[:-1]) - model.centre[2]
    xx, zz = np.meshgrid(xc, zc, indexing="ij")
    r, theta = np.hypot(xx, zz), np.arctan2(np.abs(xx), np.abs(zz))
    k_cut = int(np.searchsorted(ek, e_cut - 1e-12))
    tail = TailTable(float(ek[k_cut]))
    band = np.clip(np.degrees(theta) // 10, 0, 8).astype(int)
    rows = []
    for i, t in enumerate(np.asarray(s.t, float)):
        if model.front is None or not (model.t_range[0] <= t <= model.t_range[1]):
            continue
        h = np.asarray(s.read("e", i), float)
        down = (r - model.front(theta, t)[0] < -3 * model.d_i) & (h.sum(axis=-1) > 0)
        if r_debris is not None:
            ok, it, ith = model._bin(theta, np.full(theta.shape, t))
            clear = down & (r > np.nan_to_num(r_debris[it, ith], nan=0.0))
        else:
            clear = down
        temp = bulk_temperature(h, ek)
        measured, expected = downstream_cells(h, ek, temp, tail, k_cut)
        row = {"t": float(t), "e_cut": float(ek[k_cut])}
        for name, sel in (("behind the front", down), ("behind the front, outside the debris", clear)):
            sel = sel & np.isfinite(temp)
            if not sel.any():
                row[name] = {"cells": 0}
                continue
            q = np.percentile(temp[sel], (10, 50, 90))
            ent = {"cells": int(sel.sum()), "T_p10": float(q[0]), "T_median": float(q[1]), "T_p90": float(q[2]),
                   "weight": float(h[sel].sum()), "measured": float(measured[sel].sum()),
                   "expected": [float(e[sel].sum()) for e in expected], "bands": []}
            for b in range(9):
                m = sel & (band == b)
                ent["bands"].append({"cells": int(m.sum()), "measured": float(measured[m].sum()),
                                     "expected": [float(e[m].sum()) for e in expected]})
            row[name] = ent
        rows.append(row)

    def temp_down(t):
        pts = [(w["t"], (w["behind the front, outside the debris"].get("T_median")
                         or w["behind the front"].get("T_median"))) for w in rows]
        pts = [(a, b) for a, b in pts if b is not None]
        return float(np.interp(t, [a for a, _ in pts], [b for _, b in pts])) if pts else np.nan
    return rows, temp_down


# ------------------------------------------------------------------ synthetic checks (Sect. 12)
def _synthetic_model():
    t, theta = 18.5 * np.arange(41), np.radians(2.5 + 5.0 * np.arange(18))
    yes = np.ones((41, 18), bool)
    front = sk.sphere(20.0, 0.2)
    return FrontModel(front, front, (0.0, 740.0), (0.0, 740.0), t, theta, yes, yes, yes, 3.0 * yes, 3.0 * yes, yes,
                      [0.0, 740.0], [0.01, 0.01], np.zeros(41), 4.0, 0.2857142857, (192.0, 192.0, 192.0),
                      (768, 768, 768), (0.5, 0.5, 0.5), 2)


def _place(model, rng, rho, z_abs, u_par, u_perp, rho_birth, iz):
    """Electrons with the guiding-centre radius rho, the height z_abs and the momentum (u_par away from the
    equatorial plane, u_perp), at random azimuth, gyro-phase and hemisphere; ambient ids of cells at rho_birth."""
    n = len(rho)
    az, ph, sgn = rng.uniform(0, 2 * np.pi, n), rng.uniform(0, 2 * np.pi, n), rng.choice([-1.0, 1.0], n)
    ux, uy = u_perp * np.cos(ph), u_perp * np.sin(ph)
    x = rho * np.cos(az) + uy / model.b0  # inverse of sk.guiding_centre for the charge -1
    y = rho * np.sin(az) - ux / model.b0
    c, d = model.centre, model.dx
    ix = np.floor((c[0] + rho_birth * np.cos(az)) / d[0]).astype(np.int64)
    iy = np.floor((c[1] + rho_birth * np.sin(az)) / d[1]).astype(np.int64)
    ids = model.ppc * ((ix * model.n_cells[1] + iy) * model.n_cells[2] + iz)
    return {"x": x + c[0], "y": y + c[1], "z": sgn * z_abs + c[2], "ux": ux, "uy": uy, "uz": sgn * u_par,
            "id": ids.astype(np.uint64)}


def _join(parts):
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def selfcheck(check):
    """Checks of the census on constructed electrons; `check(name, ok, detail)` records one result."""
    rng = np.random.default_rng(7)
    model = _synthetic_model()
    b0 = model.b0
    times = 18.5 * np.arange(14, 28)
    t1 = times[-1]

    def reflected_state(n, lo_deg, hi_deg):
        t_r, th = rng.uniform(200.0, 480.0, n), np.radians(rng.uniform(lo_deg, hi_deg, n))
        rad = 20.0 + 0.2 * t_r
        return t_r, th, rad * np.sin(th), rad * np.cos(th), 0.2 / np.cos(th)

    # (1) thermal electrons after one reflection; half of them born 6 d_e away from their flux surface
    n = 4000
    t_r, th, rho, z_f, v_ht = reflected_state(n, 50.0, 75.0)
    u_par0, u_perp = rng.normal(0.0, 0.1, n), rng.uniform(0.0, 0.5, n)  # not a thermal sample: every energy <= 0.15
    _, refl, u_par_r, gam_r = sk.ht_reflect(u_par0, u_perp, v_ht, 3.0)
    keep = refl & (np.sqrt(1.0 + u_par0 ** 2 + u_perp ** 2) - 1.0 <= 0.15)
    t_r, th, rho, z_f, v_ht, u_perp, u_par_r, gam_r = (a[keep][:400] for a in
                                                       (t_r, th, rho, z_f, v_ht, u_perp, u_par_r, gam_r))
    n_a = len(t_r)
    moved = np.arange(n_a) % 2 == 0
    rho_p = np.sqrt(rho ** 2 - (u_perp / b0) ** 2)
    a = _place(model, rng, rho, z_f + u_par_r / gam_r * (t1 - t_r), u_par_r, u_perp, rho_p + 6.0 * moved,
               250 + np.arange(n_a))
    # (2) electrons whose perpendicular energy alone is above the thermal bound
    n_b = 100
    t_b, th_b, rho_b, z_b, v_b = reflected_state(n_b, 50.0, 70.0)
    uperp_b, vp_b = rng.uniform(0.8, 1.2, n_b), v_b + rng.uniform(0.12, 0.2, n_b)
    b = _place(model, rng, rho_b, z_b + vp_b * (t1 - t_b), vp_b * np.sqrt((1 + uperp_b ** 2) / (1 - vp_b ** 2)),
               uperp_b, rho_b, 250 + n_a + np.arange(n_b))
    # (4) leak: inside the loss cone, cold in the frame of the foot point, downstream one dump earlier
    n_d = 50
    t_d = rng.uniform(300.0, 450.0, n_d)
    rad = 20.0 + 0.2 * t_d
    g = 1.0 / np.sqrt(1.0 - 0.4 ** 2)
    gp = np.sqrt(1.0 + 0.5 ** 2 + 0.05 ** 2)
    up_d, gam_d = g * (0.5 + 0.4 * gp), g * (gp + 0.4 * 0.5)
    iz_d = 250 + n_a + n_b + np.arange(n_d)
    d = _place(model, rng, rad * np.sin(np.pi / 3), rad * np.cos(np.pi / 3) + up_d / gam_d * (t1 - t_d),
               np.full(n_d, up_d), np.full(n_d, 0.05), rad * np.sin(np.pi / 3), iz_d)
    d_before = _place(model, rng, np.full(n_d, 20.0), np.full(n_d, 20.0), np.full(n_d, 0.3), np.full(n_d, 0.05),
                      rad * np.sin(np.pi / 3), iz_d)
    d_before["id"] = d["id"]
    # (5) trapped: 13 (and, for the undetermined ones, 5) successive dumps in the layer, then upstream
    n_e = 40
    iz_e = 250 + n_a + n_b + n_d + np.arange(n_e)
    t_e = np.full(n_e, 400.0)
    rad = 20.0 + 0.2 * t_e
    vp_e = 0.4 + 0.15
    e = _place(model, rng, rad * np.sin(np.pi / 3), rad * np.cos(np.pi / 3) + vp_e * (t1 - t_e),
               np.full(n_e, vp_e * np.sqrt(2.0 / (1 - vp_e ** 2))), np.ones(n_e), rad * np.sin(np.pi / 3), iz_e)
    # (6) far upstream in the last two dumps, with a known change of the pitch angle
    n_f = 300
    iz_f = np.arange(n_f) % 200
    mu0 = rng.uniform(-0.8, 0.8, n_f)
    f_parts = []
    for mu in (mu0, mu0 + 0.01):
        p = _place(model, rng, np.full(n_f, 150.0), np.full(n_f, 100.0), np.abs(mu), np.sqrt(1 - mu ** 2),
                   np.full(n_f, 150.0), iz_f)
        p["uz"] = mu.copy()
        p["id"] = np.uint64(20 * 768 ** 3) + np.arange(n_f, dtype=np.uint64)  # injected at a face
        f_parts.append(p)
    dumps = []
    for k, t in enumerate(times[:-1]):
        rad = 20.0 + 0.2 * t - 2.0
        lay = _place(model, rng, np.full(n_e, rad * np.sin(1.2)), np.full(n_e, rad * np.cos(1.2)),
                     np.full(n_e, 0.8), np.zeros(n_e), np.full(n_e, 100.0), iz_e)
        lay["id"] = e["id"]
        if k < 8:  # the second half is in the layer in the last five dumps only
            lay = {key: v[:n_e // 2] for key, v in lay.items()}
        parts = [lay] + ([d_before, f_parts[0]] if k == len(times) - 2 else [])
        dumps.append((t, _join(parts)))
    dumps.append((t1, _join([a, b, d, e, f_parts[1]])))
    rec, acc = run_census(model, dumps, temp_down_of=lambda t: 0.05)
    rows, nu_use = scatter_rates(acc, 18.5, 1e-5)
    rec = assign(rec, model, np.full(len(E_NAMES), 1e-5), p999_of=lambda t: 3.0)
    last = rec["dump"] == len(times) - 1
    res, cls = rec["result"][last], rec["cls"][last]
    sl_a, sl_b = slice(0, n_a), slice(n_a, n_a + n_b)
    sl_d, sl_e = slice(n_a + n_b, n_a + n_b + n_d), slice(n_a + n_b + n_d, n_a + n_b + n_d + n_e)
    up_a = cls[sl_a] <= 1
    v_par = (u_par_r / gam_r)
    near = np.abs(v_par - v_ht) < 0.1
    want = np.where(near, 5, 1)
    dt_r = np.abs(rec["t_r"][last][sl_a] - t_r)[up_a]
    dtbn = np.abs(np.degrees(rec["tbn_r"][last][sl_a]) - np.degrees(th))[up_a]
    check("census, one reflection", up_a.mean() > 0.8 and (res[sl_a][up_a] == want[up_a]).all()
          and dt_r.max() < 0.5 and dtbn.max() < 0.5,
          f"{int(up_a.sum())}/{n_a} upstream: {int((res[sl_a][up_a] == 1).sum())} S0 at exit, "
          f"{int((res[sl_a][up_a] == 5).sum())} undetermined (|v_par - V_HT| < 0.1); t_r within {dt_r.max():.2f}, "
          f"theta_Bn within {dtbn.max():.2f} deg")
    up_b = cls[sl_b] <= 1
    check("census, residual", up_b.mean() > 0.8 and (res[sl_b][up_b] == 6).all() and (res[sl_b][~up_b] == 0).all(),
          f"{int((res[sl_b] == 6).sum())}/{int(up_b.sum())} upstream electrons with a perpendicular energy above "
          f"E_th are residual; {int((~up_b).sum())} still within d_i of the front are not judged")
    hot = up_a & (rec["gamma"][last][sl_a] >= 1.5)
    k1 = rec["k1"][last][sl_a]
    check("census, K1", hot.sum() > 50 and (k1[hot] == moved[hot]).all() and not k1[~hot].any(),
          f"{int(k1[hot].sum())} of {int(hot.sum())} upstream electrons above gamma 1.5 flagged, "
          f"{int(moved[hot].sum())} were born 6 d_e away")
    check("census, leak and trapping", (res[sl_d] == 3).all() and (res[sl_e][:n_e // 2] == 4).all()
          and (res[sl_e][n_e // 2:] == 5).all() and (rec["und"][last][sl_e][n_e // 2:] == 16).all(),
          f"{int((res[sl_d] == 3).sum())}/{n_d} S1; {int((res[sl_e] == 4).sum())}/{n_e // 2} S2 after 13 dumps in "
          f"the layer, {int((res[sl_e] == 5).sum())}/{n_e // 2} undetermined after 5")
    with np.errstate(invalid="ignore"):
        nu_want = np.mean(0.01 ** 2 / (1 - (mu0 + 0.005) ** 2)) / 18.5
    pooled = acc.sum(axis=0)
    nu_got = pooled[:, 2].sum() / pooled[:, 1].sum() / 18.5
    sl_f = slice(n_a + n_b + n_d + n_e, None)
    check("census, scattering between dumps", pooled[:, 0].sum() == n_f and abs(nu_got / nu_want - 1) < 1e-6
          and (res[sl_f] == 7).all() and (rec["origin"][last][sl_f] == 3).all(),
          f"{int(pooled[:, 0].sum())} pairs, nu = {nu_got:.3e} (constructed {nu_want:.3e}); far-upstream electrons "
          f"are 'residual, no encounter'")
    # the per-cell temperature fit and the tail of a Maxwell-Juttner histogram
    edges = np.geomspace(1e-4, 100.0, 97)
    cen = np.sqrt(edges[1:] * edges[:-1])
    temps = np.array([0.01, 0.05, 0.1])
    h = np.diff(edges) * (1 + cen) * np.sqrt((1 + cen) ** 2 - 1) * np.exp(-cen / temps[:, None])
    h *= 1e6 / h.sum(axis=1, keepdims=True)
    got = bulk_temperature(h, edges)
    one, _ = fit_t(h[1], edges, np.float64(0.05))
    ref, _ = sk.fit_temperature(h[1], edges, 0.05)
    k_cut = int(np.searchsorted(edges, 0.5 - 1e-12))
    measured, expected = downstream_cells(h, edges, got, TailTable(float(edges[k_cut])), k_cut)
    ratio = measured[1:] / expected[1][1:]
    check("census, downstream reference", np.allclose(got, temps, rtol=0.01) and np.isclose(one, ref, rtol=1e-12)
          and np.all(np.abs(ratio - 1) < 0.1),
          f"bulk fit {got.round(4).tolist()} for {temps.tolist()}; tail above {edges[k_cut]:.2f} measured / expected "
          f"{ratio.round(3).tolist()} at T = 0.05, 0.1")

    # (8) the dump spacing: equal to the slices (L, S25) or a multiple of it (P5); refused off the slices
    sl = 9.0 * np.arange(41)

    def refused(times):
        try:
            dump_spacing(times, sl)
        except SystemExit:
            return True
        return False
    ok = (dump_spacing(sl, sl) == 9.0 and dump_spacing(18.0 * np.arange(21), sl) == 18.0
          and refused(4.5 + 18.0 * np.arange(20)) and refused(np.array([0.0, 18.0, 27.0])) and refused(sl[:1]))
    check("census, dump spacing", ok, "dumps on every slice: 9; on every second: 18; off the slices, uneven or "
          "single: refused")

    # (9) counts without back-tracing (large runs) agree with the census on the class and the energy bin
    last = rec["dump"] == len(dumps) - 1
    want = np.zeros((len(CLASSES), len(E_NAMES)), np.int64)
    eb = e_bin(rec["gamma"][last])
    np.add.at(want, (rec["cls"][last][eb >= 0], eb[eb >= 0]), 1)
    got_c, got_u, _ = count_dump(model, dict(dumps[-1][1], time=dumps[-1][0]), dumps[-1][0], chunk=997)
    check("census, counts without back-tracing", np.array_equal(got_c, want) and got_u.sum() == want[:2].sum(),
          f"{int(want.sum())} electrons in the last dump; class x energy table equal in chunks of 997")


# ------------------------------------------------------------------ the command
def _open(run_dir):
    import hdf5plugin  # noqa: F401  (the files are compressed with zstd)

    from pypic.analysis import open_run
    return open_run(str(run_dir))


def _nonzero(table, names_a, names_b):
    return {a: {b: int(table[i, j]) for j, b in enumerate(names_b) if table[i, j]}
            for i, a in enumerate(names_a) if table[i].any()}


def summarize(rec, all_t, up_t, model, times, null=None):
    """The tables of Sect. 6 and 7 (K1) as a dict; counts are macro-electrons with gamma >= 1.5 unless said."""
    n_dumps = len(times)
    hot = slice(1, None)  # energy bins from gamma = 1.5
    n1 = {round(r["t"], 3): r for r in null["N1"]} if null else {}
    per_dump, k1_rows, up_rows = [], [], []
    amb_up = (rec["origin"] <= 2) & (rec["cls"] <= 1) & (rec["gamma"] >= 1.5)
    for k, t in enumerate(times):
        floor = n1.get(round(float(t), 3))
        null_n = {"ambient": floor["ambient>=1.5"], "injected": floor["injected>=1.5"]} if floor else None
        per_dump.append({"t": float(t), "n_dumped": int(all_t[k].sum()), "n>=1.5": int(all_t[k][..., hot].sum()),
                         "by_origin_class": _nonzero(all_t[k][..., hot].sum(axis=-1), ORIGINS, CLASSES),
                         "null_run": null_n, "signal_threshold": 100 if null_n is None else
                         max(100, 10 * max(max(null_n.values()), 3))})
        m = amb_up & (rec["dump"] == k)
        row = {"t": float(t), "threshold": float(rec["k1_threshold"][m][0]) if m.any() else None}
        for name, sel in (("whole front", m), ("valid domain", m & (rec["flag"] >= 2))):
            k1 = sel & rec["k1"]
            row[name] = {"upstream ambient": int(sel.sum()), "k1": int(k1.sum()), "expected": 0.001 * int(sel.sum()),
                         "outward": int((k1 & (rec["drho"] > 0)).sum()), "inward": int((k1 & (rec["drho"] < 0)).sum()),
                         "beyond the debris": int((k1 & rec["beyond_debris"]).sum()),
                         "beyond the debris, first seen in the valid domain":
                             int((k1 & rec["beyond_debris"] & (rec["first_flag"] >= 2)).sum()),
                         "by result": {RESULTS[r]: int((k1 & (rec["result"] == r)).sum()) for r in range(len(RESULTS))
                                       if (k1 & (rec["result"] == r)).any()}}
        k1_rows.append(row)
        domains = (("detected front", slice(1, None)), ("valid domain", slice(2, None)), ("no front", slice(0, 1)))
        for name, fl in domains:
            tab = up_t[k, fl].sum(axis=(0, 4, 5))  # band, energy, result
            if not tab[:, hot].any():
                continue
            entry = {"t": float(t), "domain": name, "results": {RESULTS[r]: int(tab[:, hot, r].sum())
                                                                for r in range(1, len(RESULTS))}, "bands": []}
            for b in range(NO_FOOT + 1):
                for e in range(1, len(E_NAMES)):
                    if tab[b, e].any():
                        entry["bands"].append({"band": "none" if b == NO_FOOT else f"{10 * b}-{10 * b + 10}",
                                               "gamma": E_NAMES[e], "n": int(tab[b, e].sum()),
                                               "measurable": bool(tab[b, e].sum() >= 100),
                                               "counts": tab[b, e, 1:].tolist()})
            hq = up_t[k, fl][:, :, hot].sum(axis=(0, 1, 2))  # result, hemisphere, quadrant
            entry["residual by hemisphere"] = hq[6].sum(axis=1).tolist()
            entry["residual by quadrant"] = hq[6].sum(axis=0).tolist()
            entry["S0 by hemisphere"] = (hq[1] + hq[2]).sum(axis=1).tolist()
            entry["S0 by quadrant"] = (hq[1] + hq[2]).sum(axis=0).tolist()
            up_rows.append(entry)
    und = {}
    sel = (rec["result"] == 5) & (rec["gamma"] >= 1.5)
    for b, name in enumerate(UND):
        und[name] = int((sel & ((rec["und"] >> b) & 1 == 1)).sum())
    und["flight time only"] = int((sel & (rec["und"] == 8)).sum())
    return {"n_dumps": n_dumps, "result_names": list(RESULTS[1:]), "per_dump": per_dump, "upstream": up_rows,
            "K1": k1_rows, "undetermined_by_cause_electron_dumps": und}


def debris_table(run, model, w_ambient):
    """Debris electrons of the dump by class and energy bin, in units of one ambient macro-electron."""
    p = run.particles("hot_e")
    rows = []
    for it, t in zip(p.iterations, p.times):
        d = p.get(it, "de")
        if not len(d["id"]):
            continue
        rec = classify_only(model, d, float(t))
        eb = e_bin(rec["gamma"])
        w = np.asarray(d["w"], float) / w_ambient
        tab = np.zeros((len(CLASSES), len(E_NAMES)))
        np.add.at(tab, (rec["cls"][eb >= 0], eb[eb >= 0]), w[eb >= 0])
        rows.append({"t": float(t), "entries": int(len(w)), "weight>=1.25": float(tab.sum()),
                     "weight>=1.5": float(tab[:, 1:].sum()),
                     "by_class>=1.5": {c: float(tab[i, 1:].sum()) for i, c in enumerate(CLASSES) if tab[i, 1:].any()},
                     "upstream_by_gamma": {E_NAMES[e]: float(tab[:2, e].sum()) for e in range(len(E_NAMES))
                                           if tab[:2, e].any()}})
    return rows


def cmd_census(args) -> int:
    run_dir = Path(args.run)
    run = _open(run_dir)
    model = FrontModel.from_run(run, run_dir)
    if model.front is None:
        print("no fitted front in analysis/front.json: run the command front first; nothing to classify against")
        return 1
    null = json.loads(Path(args.null).read_text()) if args.null else None
    stem, tau_protocol = "census", model.tau_dwell
    if getattr(args, "tau", None):  # the replacement rule of protocol Sect. 3.2; the first result is kept
        model.tau_dwell, stem = float(args.tau), "census-tau"
    r_debris = np.load(run_dir / "analysis" / "front.npz")["r_debris_thr"][:, 0, :]
    down_rows, temp_down = downstream_reference(run, model, r_debris)
    p = run.particles("hot_e")
    times = np.asarray(p.times, float)
    model.dt_dump = dump_spacing(times, model.t)
    rec, acc = run_census(model, ((float(t), p.get(it, "e")) for it, t in zip(p.iterations, times)), temp_down)
    nu_null = float(null["tracks"]["nu_local"][-1][1]) if null else np.nan  # v > 0.3, lag 16
    nu_rows, nu_use = scatter_rates(acc, model.dt_dump, nu_null)
    p999 = None
    if null:
        rows = null["tracks"]["delta_rho"]
        tt, pp = [r["t"] for r in rows], [r["birth all"]["p99.9"] for r in rows]

        def p999(t):
            return np.interp(t, tt, pp)
    rec = assign(rec, model, nu_use, p999)
    all_t, up_t = tables(rec, len(times))
    out = summarize(rec, all_t, up_t, model, times, null)
    names = [s["name"] for s in run.config["species"]]
    w_ambient = float(np.prod(model.dx)) / model.ppc
    out.update({"run": str(run_dir), "commit": run.info.get("git_commit"),
                "null": str(args.null) if args.null else None,
                "b_max": model.b_max, "b_cell_max": model.b_cell_max, "tau_dwell": model.tau_dwell,
                "tau_dwell_protocol": tau_protocol,
                "dump_spacing": model.dt_dump,
                "layer_dumps_S2": int(np.ceil(10 * 2 * np.pi / model.b0 / model.dt_dump)),
                "layer_dumps_undetermined": int(np.ceil(3 * 2 * np.pi / model.b0 / model.dt_dump)),
                "E_th": [float(sk.thermal_bound(model.temp_upstream(t), model.n_total)) for t in times[::8]],
                "T_down": [{"t": float(t), "T": temp_down(float(t))} for t in times[::4]],
                "N3b": nu_rows, "downstream": down_rows,
                "debris": debris_table(run, model, w_ambient) if "de" in names else [],
                "notes": ["class: an electron within d_i of the front is in the layer even when its field line is "
                          "not reached (protocol Sect. 6 reads 'upstream' for it); the conservative reading",
                          "delta is the radial distance to the fitted density front, not the normal distance",
                          "downstream reference: slab |y| < 8 d_e, lab-frame energies (the bulk flow is included)",
                          "N3b: both dumps need gamma >= 1.25, so the energy change of the lowest bin is biased"]})
    dest = run_dir / "analysis"
    (dest / f"{stem}.json").write_text(json.dumps(out, indent=1))
    small = {k: (v.astype(np.float32) if v.dtype == np.float64 and k not in ("t",) else v) for k, v in rec.items()}
    np.savez_compressed(dest / f"{stem}.npz", all_table=all_t, up_table=up_t, n3b=acc, times=times, **small)
    report(out)
    print(f"wrote {dest / (stem + '.json')} and {stem}.npz ({len(rec['id'])} electron records)")
    return 0


def count_dump(model, part, t, chunk=2_000_000):
    """Counts of one dump without back-tracing: by_class[class, energy], upstream[flag, energy], and the
    radius (10/50/90 %) of the electrons with gamma >= 1.5 per class."""
    by_class = np.zeros((len(CLASSES), len(E_NAMES)), np.int64)
    upstream = np.zeros((4, len(E_NAMES)), np.int64)
    radii = [[] for _ in CLASSES]
    n = len(part["id"])
    for a in range(0, n, chunk):
        sub = {k: (np.asarray(v)[a:a + chunk] if np.ndim(v) == 1 and len(v) == n else v) for k, v in part.items()}
        rec = classify_only(model, sub, float(t))
        eb = e_bin(rec["gamma"])
        ok = eb >= 0
        np.add.at(by_class, (rec["cls"][ok], eb[ok]), 1)
        up = ok & (rec["cls"] <= 1)
        np.add.at(upstream, (rec["flag"][up], eb[up]), 1)
        r = np.hypot(rec["rho_gc"], rec["z_abs"])
        for c in range(len(CLASSES)):
            m = (rec["cls"] == c) & (rec["gamma"] >= 1.5)
            if m.any():
                radii[c].append(r[m])
    rad = {}
    for c, name in enumerate(CLASSES):
        if radii[c]:
            q = np.percentile(np.concatenate(radii[c]), (10, 50, 90))
            rad[name] = [float(v) for v in q]
    return by_class, upstream, rad


def cmd_counts(args) -> int:
    """Counts of the ambient electrons per dump by class and energy bin, without back-tracing (no result, no
    origin): for runs whose dumps are too large for the census (P5: millions of electrons per dump)."""
    run_dir = Path(args.run)
    run = _open(run_dir)
    model = FrontModel.from_run(run, run_dir)
    if model.front is None:
        print("no fitted front in analysis/front.json: run the command front first")
        return 1
    p = run.particles("hot_e")
    times = np.asarray(p.times, float)
    model.dt_dump = dump_spacing(times, model.t)
    rows = []
    for it, t in zip(p.iterations, times):
        part = p.get(it, "e")
        if not len(part["id"]):
            rows.append({"t": float(t), "n": 0})
            continue
        by_class, upstream, rad = count_dump(model, part, t)
        rows.append({"t": float(t), "n": int(len(part["id"])),
                     "by_class": {c: by_class[i].tolist() for i, c in enumerate(CLASSES) if by_class[i].any()},
                     "upstream_by_flag": {f: upstream[i].tolist() for i, f in enumerate(FLAGS) if upstream[i].any()},
                     "radius_gamma>=1.5": rad})
        print(f"t {t:6.1f}  n {len(part['id']):9d}  " + "; ".join(
            f"{c}: {int(by_class[i, 1:].sum())}" for i, c in enumerate(CLASSES) if by_class[i, 1:].any()), flush=True)
    out = {"run": str(run_dir), "commit": run.info.get("git_commit"), "energy_bins": list(E_NAMES),
           "classes": list(CLASSES), "flags": list(FLAGS), "dump_spacing": model.dt_dump, "rows": rows,
           "notes": ["no back-tracing: no result (S0, S1, ...), no K1, no origin",
                     "class: an electron within d_i of the front is in the layer",
                     "radius: of the guiding centre in the x-y plane and |z|, in d_e"]}
    dest = run_dir / "analysis" / "counts.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"wrote {dest}")
    return 0


def report(out):
    """A short text view of census.json."""
    print(f"census of {out['run']}: b_max {out['b_max']:.2f} (cell {out['b_cell_max']:.2f}), tau_dwell "
          f"{out['tau_dwell']:.0f}, E_th {out['E_th'][0]:.3f}-{out['E_th'][-1]:.3f}")
    print("\nelectrons with gamma >= 1.5 by origin and class (every fourth dump)")
    for row in out["per_dump"][::4] + out["per_dump"][-1:]:
        cells = "; ".join(f"{o}: " + ", ".join(f"{c} {n}" for c, n in d.items())
                          for o, d in row["by_origin_class"].items())
        print(f"  t {row['t']:6.1f}  n {row['n>=1.5']:8d}  null {row['null_run']}  {cells}")
    print("\nupstream electrons with gamma >= 1.5 by result " + str(out["result_names"]))
    for e in out["upstream"]:
        if e["domain"] != "no front" and round(e["t"] / out.get("dump_spacing", 18.5)) % 4 == 0:
            print(f"  t {e['t']:6.1f} {e['domain']:15s} {list(e['results'].values())}  residual by hemisphere "
                  f"{e['residual by hemisphere']} quadrant {e['residual by quadrant']}")
    last = [e for e in out["upstream"] if e["t"] == out["per_dump"][-1]["t"]]
    for e in last:
        print(f"\n  last dump, {e['domain']}: band, gamma, n, counts per result")
        for b in e["bands"]:
            mark = "" if b["measurable"] else "  (n.m.)"
            print(f"    {b['band']:>6s} {b['gamma']:>5s} {b['n']:7d} {b['counts']}{mark}")
    print("\nundetermined, by cause (electron-dumps):", out["undetermined_by_cause_electron_dumps"])
    print("\nN3b (pairs of successive dumps far upstream):")
    for r in out["N3b"]:
        print(f"  gamma {r['gamma']:>8s} pairs {r['pairs']:7d} nu {r['nu']:.3e} used {r['nu_used']:.3e} "
              f"({r['nu_from']}) flight limit {r['flight_limit']:.0f}")
    print("\nK1 (whole front / valid domain): t, threshold, upstream ambient, K1, expected, outward, inward, "
          "beyond debris")
    for r in out["K1"][::4] + out["K1"][-1:]:
        cells = " / ".join(f"{d['upstream ambient']} {d['k1']} {d['expected']:.1f} {d['outward']} {d['inward']} "
                           f"{d['beyond the debris']}" for d in (r["whole front"], r["valid domain"]))
        print(f"  t {r['t']:6.1f} thr {r['threshold']}  {cells}")
    print("\ndownstream reference: t, cells, T median (p10-p90), measured above e_cut, expected at 0.9/1.0/1.1 T")
    for r in out["downstream"]:
        for name in ("behind the front", "behind the front, outside the debris"):
            d = r[name]
            if d.get("cells"):
                print(f"  t {r['t']:6.1f} {name:38s} {d['cells']:6d} T {d['T_median']:.4f} ({d['T_p10']:.4f}-"
                      f"{d['T_p90']:.4f}) measured {d['measured']:.3e} expected "
                      + " ".join(f"{e:.3e}" for e in d["expected"]))
    if out["debris"]:
        print("\ndebris electrons (weight in ambient macro-electrons, gamma >= 1.5) by class")
        for r in out["debris"][::4] + out["debris"][-1:]:
            print(f"  t {r['t']:6.1f} entries {r['entries']:9d} weight {r['weight>=1.5']:.3e} {r['by_class>=1.5']}")
