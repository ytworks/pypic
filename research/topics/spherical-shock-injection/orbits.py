"""Measures on the tracked orbits (exploratory-protocol.md Sect. 7 and N4, f_FN of Sect. 8).

The tracked set is fixed at the start of a run and holds ambient electrons only. For every orbit:
the side of the fitted density front against time, the encounters with the layer, the dwell in
it; for the orbits that become energetic also the magnetic moment, the energy and the canonical
radius before and after the encounters, the work of the recorded field and its part along B.
The census (census.py) is joined by id: it gives the false-negative rate of the test B1'.
"""
from __future__ import annotations

import json
from pathlib import Path

import census
import numpy as np
import shockkin as sk

CATS = ("never left the upstream", "encounter, back upstream", "in the layer at the end", "transmitted",
        "returned from downstream", "starts behind the front", "no sample in the range of the fit")
GAMMA_HOT = 1.25  # the threshold of the particle dumps


def sides(model, pos, t):
    """Radial distance to the fitted density front (NaN outside the fit's time range or where the
    particle has left the box) and the polar angle, folded to one hemisphere. pos: (3, n_t, n)."""
    c = model.centre
    x, y, z = (np.asarray(pos[i], float) - c[i] for i in range(3))
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arctan2(np.hypot(x, y), np.abs(z))
    tt = np.broadcast_to(np.asarray(t, float)[:, None], r.shape)
    ok = np.isfinite(r) & (tt >= model.t_range[0]) & (tt <= model.t_range[1])
    delta = np.full(r.shape, np.nan)
    delta[ok] = r[ok] - model.front(theta[ok], tt[ok])[0]
    return delta, theta


def encounters(delta, d_i: float) -> dict:
    """Histories of the side of the front. Upstream is delta > d_i, downstream delta < -3 d_i, the layer
    between (Sect. 6). An encounter is an entry into delta <= d_i; the next one counts only after the
    orbit has been beyond 2 d_i again (hysteresis of 1 d_i). Returns per orbit: cat (CATS), n_enc,
    n_cross (side changes about the front with a hysteresis of 1 d_i), run_max (longest stay in the
    layer, samples), n_in (samples not upstream), i_enter, i_exit (first and last sample not upstream),
    i_min and d_min (deepest sample)."""
    n_t, n = delta.shape
    first, last = np.zeros(n, np.int8), np.zeros(n, np.int8)
    started, armed = np.zeros(n, bool), np.ones(n, bool)
    ever_down, ever_in = np.zeros(n, bool), np.zeros(n, bool)
    n_enc, run, run_max, n_in = (np.zeros(n, np.int32) for _ in range(4))
    i_enter, i_exit, i_min = (np.full(n, -1, np.int32) for _ in range(3))
    d_min = np.full(n, np.inf)
    for k in range(n_t):
        row = delta[k]
        ok = np.isfinite(row)
        with np.errstate(invalid="ignore"):
            side = np.where(row > d_i, 1, np.where(row < -3.0 * d_i, -1, 0)).astype(np.int8)
            far = ok & (row > 2.0 * d_i)
            deeper = ok & (row < d_min)
        first = np.where(ok & ~started, side, first)
        started |= ok
        last = np.where(ok, side, last)
        inside = ok & (side <= 0)
        hit = inside & armed
        n_enc += hit
        armed = (armed & ~hit) | far
        ever_down |= ok & (side < 0)
        i_enter = np.where(inside & ~ever_in, k, i_enter)
        ever_in |= inside
        i_exit = np.where(inside, k, i_exit)
        n_in += inside
        run = np.where(ok & (side == 0), run + 1, np.where(ok, 0, run))
        run_max = np.maximum(run_max, run)
        i_min = np.where(deeper, k, i_min)
        d_min = np.where(deeper, row, d_min)
    cat = np.select([~started, first != 1, ~ever_in, ever_down & (last == 1), ever_down, last == 1],
                    [6, 5, 0, 4, 3, 1], default=2).astype(np.int8)
    return {"cat": cat, "n_enc": n_enc, "n_cross": sk.count_crossings(delta, d_i), "run_max": run_max, "n_in": n_in,
            "i_enter": i_enter, "i_exit": i_exit, "i_min": i_min, "d_min": d_min}


def cumulative_work(u, e, b, t, dt_step: float, charge: float = -1.0):
    """Work of the recorded field along the orbits, cumulative in time, and its part along B. The momentum
    is centred on the field time first. Returns (gamma, work, work_par), each (n_t - 1, n); missing samples
    add nothing."""
    uc = sk.centre_momentum(u, t, dt_step)
    gam = np.sqrt(1.0 + np.sum(uc ** 2, axis=0))
    v = uc / gam
    ee, bb = e[:, :-1], b[:, :-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        bhat = bb / np.sqrt(np.sum(bb ** 2, axis=0))
    p = charge * np.sum(ee * v, axis=0)
    p_par = charge * np.sum(ee * bhat, axis=0) * np.sum(v * bhat, axis=0)
    dt = np.diff(np.asarray(t, float)[:-1])[:, None]
    out = []
    for q in (p, p_par):
        step = np.nan_to_num(0.5 * (q[1:] + q[:-1]) * dt)
        out.append(np.concatenate([np.zeros((1, q.shape[1])), np.cumsum(step, axis=0)]))
    return gam, out[0], out[1]


def _window(a, lo, hi, good):
    """Median of a[lo:hi] over the samples that are good; NaN with fewer than three."""
    lo, hi = max(lo, 0), min(hi, len(a))
    if hi <= lo:
        return np.nan
    v = a[lo:hi][good[lo:hi]]
    return float(np.median(v)) if len(v) >= 3 else np.nan


def measure(model, pos, u, b, e, t, dt_step: float) -> dict:
    """All measures of a set of orbits (arrays (3, n_t, n), float64). Quantities "in" are medians over one
    gyro-period upstream before the first entry into the layer, "out" the same after the last exit."""
    t = np.asarray(t, float)
    delta, theta = sides(model, pos, t)
    out = encounters(delta, model.d_i)
    n = delta.shape[1]
    dts = float(t[1] - t[0])
    n_w = max(4, int(round(2.0 * np.pi / model.b0 / dts)))
    mu = sk.magnetic_moment(u, b)
    gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    c = model.centre
    rho_p = sk.canonical_radius(pos[0] - c[0], pos[1] - c[1], u[0], u[1], model.b0)
    gam_c, work, work_par = cumulative_work(u, e, b, t, dt_step)
    up = np.isfinite(delta) & (delta > model.d_i)
    keys = ("mu_in", "mu_out", "gamma_in", "gamma_out", "rho_in", "rho_out", "work_enc", "work_par_enc",
            "dgamma_enc", "theta_min", "t_min", "t_enter", "t_exit")
    res = {k: np.full(n, np.nan) for k in keys}
    last = work.shape[0] - 1
    for j in np.flatnonzero(out["i_enter"] >= 0):
        i0, i1 = int(out["i_enter"][j]), int(out["i_exit"][j])
        res["t_enter"][j], res["t_exit"][j] = t[i0], t[i1]
        im = int(out["i_min"][j])
        res["theta_min"][j], res["t_min"][j] = theta[im, j], t[im]
        for name, a in (("mu", mu), ("gamma", gam), ("rho", rho_p)):
            res[f"{name}_in"][j] = _window(a[:, j], i0 - n_w, i0, up[:, j])
            res[f"{name}_out"][j] = _window(a[:, j], i1 + 1, i1 + 1 + n_w, up[:, j])
        a0, a1 = max(i0 - 1, 0), min(i1 + 1, last)
        res["work_enc"][j] = work[a1, j] - work[a0, j]
        res["work_par_enc"][j] = work_par[a1, j] - work_par[a0, j]
        res["dgamma_enc"][j] = gam_c[a1, j] - gam_c[a0, j]
    ok = np.isfinite(gam_c)
    any_ok = ok.any(axis=0)
    i_first = np.argmax(ok, axis=0)
    i_last = last - np.argmax(ok[::-1], axis=0)
    cols = np.arange(n)
    dgam = np.where(any_ok, gam_c[i_last, cols] - gam_c[i_first, cols], np.nan)
    w_all = np.where(any_ok, work[i_last, cols] - work[i_first, cols], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        res["closure"] = np.where(np.abs(dgam) >= 0.1, np.abs(w_all - dgam) / np.abs(dgam), np.nan)
        res["par_share"] = (work_par[i_last, cols] - work_par[i_first, cols]) / w_all
    res["dgamma_run"], res["gamma_max"] = dgam, np.nanmax(np.where(np.isfinite(gam), gam, -np.inf), axis=0)
    res["dwell"], res["time_in"] = out["run_max"] * dts, out["n_in"] * dts
    with np.errstate(invalid="ignore"):
        ok_f = (res["t_min"] >= model.t_range_b[0]) & (res["t_min"] <= model.t_range_b[1]) \
            if model.front_b is not None else np.zeros(n, bool)
    res["tbn_min"], res["vht_min"] = np.full(n, np.nan), np.full(n, np.nan)
    if ok_f.any():
        f = sk.foot(model.front_b, res["theta_min"][ok_f], res["t_min"][ok_f])
        res["tbn_min"][ok_f], res["vht_min"][ok_f] = f["tbn"], f["v_ht"]
    res["flag_min"] = np.zeros(n, np.int8)
    seen = np.isfinite(res["t_min"])
    if seen.any():
        res["flag_min"][seen] = model.flag(res["theta_min"][seen], res["t_min"][seen])
    res.update({k: out[k] for k in ("cat", "n_enc", "n_cross", "d_min")})
    return res


def binomial_upper(k: int, n: int, level: float = 0.95) -> float:
    """Upper bound of a binomial share at the given one-sided level (exact; 3/n for k = 0 at 95 %)."""
    if n <= 0:
        return 1.0
    from scipy.stats import beta
    return 1.0 if k >= n else float(beta.ppf(level, k + 1, n - k))


def _pct(a, qs=(10, 50, 90)):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return {"n": int(len(a)), **({f"p{q}": float(np.percentile(a, q)) for q in qs} if len(a) else {})}


def _by_cat(cat, mask=None):
    m = np.ones(len(cat), bool) if mask is None else mask
    return {CATS[c]: int((m & (cat == c)).sum()) for c in range(len(CATS)) if (m & (cat == c)).any()}


def join_census(ids, res, cen, t_dump) -> dict:
    """The census records of the tracked electrons, per orbit: the result of B1' at the first dump after the
    last exit in which the electron is upstream, and whether any dump calls it a residual."""
    order = np.argsort(ids)
    pos = np.searchsorted(ids[order], cen["id"])
    pos = np.clip(pos, 0, len(ids) - 1)
    hit = ids[order][pos] == cen["id"]
    j = order[pos[hit]]
    rec = {k: np.asarray(cen[k])[hit] for k in ("dump", "result", "cls", "gamma", "t_r", "tbn_r", "k1", "drho")}
    n = len(ids)
    first = {k: np.full(n, np.nan) for k in ("result", "gamma", "t_r", "tbn_r", "t")}
    resid, k1 = np.zeros(n, bool), np.zeros(n, bool)
    t_rec = np.asarray(t_dump, float)[rec["dump"]]
    t_exit = np.where(np.isfinite(res["t_exit"]), res["t_exit"], -np.inf)
    up_after = (rec["cls"] <= 1) & (t_rec > t_exit[j])
    for i in np.argsort(-t_rec):  # the earliest record is written last
        if up_after[i]:
            for k in ("result", "gamma", "t_r", "tbn_r"):
                first[k][j[i]] = rec[k][i]
            first["t"][j[i]] = t_rec[i]
    np.logical_or.at(resid, j, (rec["result"] >= 6) & (rec["gamma"] >= 1.5))
    np.logical_or.at(k1, j, rec["k1"])
    return {"first": first, "residual": resid, "k1": k1, "n_records": int(hit.sum())}


def summarize(res, joined, tau_protocol: float) -> dict:
    hot, hot15 = res["gamma_max"] >= GAMMA_HOT, res["gamma_max"] >= 1.5
    out = {"n_tracked": int(len(hot)), "all_by_category": _by_cat(res["cat"]),
           "n_gamma_max>=1.25": int(hot.sum()), "n_gamma_max>=1.5": int(hot15.sum()),
           "gamma_max>=1.25_by_category": _by_cat(res["cat"], hot),
           "gamma_max>=1.5_by_category": _by_cat(res["cat"], hot15)}
    done = np.abs(res["dgamma_run"]) >= 0.1
    with np.errstate(invalid="ignore"):
        good = res["closure"] < 0.1
    out["N4"] = {"orbits with |dgamma| >= 0.1": int(done.sum()), "within 10 %": int((done & good).sum()),
                 "closure": _pct(res["closure"][done], (50, 90, 99)),
                 "failing by category": _by_cat(res["cat"], done & ~good),
                 "parallel share of the work (orbits within 10 %)": _pct(res["par_share"][done & good])}
    # single reflections: one encounter, never downstream, back upstream, visible in the dumps afterwards
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = res["mu_out"] / res["mu_in"]
        refl = (res["cat"] == 1) & (res["n_enc"] == 1) & (res["gamma_out"] >= GAMMA_HOT)
        kept = refl & (np.abs(ratio - 1.0) <= 0.2)
        gain = res["gamma_out"] - res["gamma_in"]
    # every energy: one encounter that reached the front, never downstream, back upstream
    with np.errstate(invalid="ignore"):
        back = (res["cat"] == 1) & (res["n_enc"] == 1) & np.isfinite(ratio)
        deep = back & (res["d_min"] <= 0.0)
        kept_all = deep & (np.abs(ratio - 1.0) <= 0.2)
        valid_all = deep & (res["flag_min"] >= 2)
    out["reflections, all energies"] = {
        "one encounter, back upstream": int(back.sum()), "of these, reached the front (delta <= 0)": int(deep.sum()),
        "mu within 20 %": int(kept_all.sum()), "in the valid domain": int(valid_all.sum()),
        "in the valid domain, mu within 20 %": int((valid_all & kept_all).sum()),
        "mu_out / mu_in": _pct(ratio[deep]), "mu_out / mu_in, valid domain": _pct(ratio[valid_all]),
        "gain": _pct(gain[deep], (10, 50, 90, 99)), "gain, valid domain": _pct(gain[valid_all], (10, 50, 90, 99)),
        "gamma_out above": {str(g): int((deep & (res["gamma_out"] >= g)).sum()) for g in (1.05, 1.1, 1.25, 1.5)},
        "gamma_out above, valid domain": {str(g): int((valid_all & (res["gamma_out"] >= g)).sum())
                                          for g in (1.05, 1.1, 1.25, 1.5)},
        "dwell in the layer": _pct(res["dwell"][deep], (50, 90, 99)),
        "dwell in the layer, gamma_out >= 1.1": _pct(res["dwell"][deep & (res["gamma_out"] >= 1.1)], (50, 90, 99)),
        "foot theta_Bn at the deepest point (deg)": _pct(np.degrees(res["tbn_min"][deep])),
        "parallel share of the work in the encounter (|work| >= 0.02)":
            _pct((res["work_par_enc"] / res["work_enc"])[deep & (np.abs(res["work_enc"]) >= 0.02)]),
        "work - dgamma in the encounter": _pct((res["work_enc"] - res["dgamma_enc"])[deep], (10, 50, 90))}
    first = joined["first"]
    rows = {}
    for name, sel in (("gamma_out >= 1.25", kept), ("gamma_out >= 1.5", kept & (res["gamma_out"] >= 1.5))):
        r = first["result"][sel]
        expl, fail, und = int(np.isin(r, (1, 2)).sum()), int(np.isin(r, (6, 7)).sum()), int((r == 5).sum())
        judged = expl + fail
        rows[name] = {"reflections with mu within 20 %": int(sel.sum()), "explained": expl, "not explained": fail,
                      "undetermined": und, "other or not in a dump": int(sel.sum()) - judged - und,
                      "f_FN": fail / judged if judged else None, "f_FN upper 95 %": binomial_upper(fail, judged)}
    with np.errstate(invalid="ignore"):
        dt_r = first["t_r"] - res["t_min"]
        dtbn = np.degrees(first["tbn_r"] - res["tbn_min"])
        lum = refl & (res["vht_min"] >= 1.0)
    out["reflections"] = {
        "one encounter, back upstream, gamma_out >= 1.25": int(refl.sum()), "mu within 20 %": int(kept.sum()),
        "mu_out / mu_in": _pct(ratio[refl]), "gain": _pct(gain[refl]), "gamma_out": _pct(res["gamma_out"][refl]),
        "dwell in the layer": _pct(res["dwell"][refl]), "dwell protocol": tau_protocol,
        "time not upstream": _pct(res["time_in"][refl]),
        "foot theta_Bn at the deepest point (deg)": _pct(np.degrees(res["tbn_min"][refl])),
        "in the valid domain": int((refl & (res["flag_min"] >= 2)).sum()), "foot superluminal": int(lum.sum()),
        "canonical radius out - in": _pct((res["rho_out"] - res["rho_in"])[refl]),
        "B1'": rows, "census t_r - time of the deepest point": _pct(dt_r[kept]),
        "census theta_Bn - foot theta_Bn of the track (deg)": _pct(dtbn[kept])}
    with np.errstate(invalid="ignore"):
        multi = (res["cat"] == 1) & (res["n_enc"] >= 2) & (res["gamma_out"] >= 1.5)
        k3 = (res["cat"] == 1) & (res["gamma_out"] >= 1.5) & (res["vht_min"] >= 1.0) & ~(np.abs(ratio - 1.0) <= 0.2)
        leak = (res["cat"] == 4) & (res["gamma_out"] >= 1.5)
    out["K2 several encounters, gamma_out >= 1.5"] = {"n": int(multi.sum()), "n_enc": _pct(res["n_enc"][multi]),
                                                      "gain": _pct(gain[multi])}
    out["K3 superluminal foot, mu not conserved, gamma_out >= 1.5"] = int(k3.sum())
    out["returned from downstream, gamma_out >= 1.5"] = int(leak.sum())
    r15 = joined["residual"]
    out["census residuals on tracks (gamma >= 1.5)"] = {
        "n": int(r15.sum()), "by_category": _by_cat(res["cat"], r15), "n_enc": _pct(res["n_enc"][r15]),
        "mu_out / mu_in": _pct(ratio[r15]), "dwell": _pct(res["dwell"][r15]),
        "closure within 10 %": int((r15 & good).sum())}
    out["K1 on tracks"] = {"n": int(joined["k1"].sum()), "by_category": _by_cat(res["cat"], joined["k1"])}
    out["census records of tracked electrons"] = joined["n_records"]
    return out


def cmd_tracks(args) -> int:
    run_dir = Path(args.run)
    run = census._open(run_dir)
    model = census.FrontModel.from_run(run, run_dir)
    if model.front is None:
        print("no fitted front in analysis/front.json: run the command front first")
        return 1
    tr = run.tracks("tracks_e", "e")
    t, ids = np.asarray(tr.t, float), np.asarray(tr.ids)
    groups = (("x", "y", "z"), ("ux", "uy", "uz"), ("Bx", "By", "Bz"), ("Ex", "Ey", "Ez"))
    data = {k: tr.f[k][...] for g in groups for k in g}  # a row of the file is one chunk: read whole datasets
    parts = []
    for a in range(0, len(ids), args.chunk):
        sel = slice(a, min(len(ids), a + args.chunk))
        parts.append(measure(model, *(np.stack([data[k][:, sel] for k in g]).astype(float) for g in groups),
                             t, float(run.dt)))
    del data
    res = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    cen = np.load(run_dir / "analysis" / "census.npz")
    joined = join_census(ids, res, cen, cen["times"])
    out = summarize(res, joined, model.tau_dwell)
    out.update({"run": str(run_dir), "commit": run.info.get("git_commit"), "sample_step": float(t[1] - t[0]),
                "stride": int(ids[1] - ids[0]) if len(ids) > 1 else None})
    dest = run_dir / "analysis"
    (dest / "tracks.json").write_text(json.dumps(out, indent=1))
    np.savez_compressed(dest / "tracks.npz", id=ids, residual=joined["residual"],
                        k1=joined["k1"], **{f"first_{k}": v for k, v in joined["first"].items()},
                        **{k: (v.astype(np.float32) if v.dtype == np.float64 else v) for k, v in res.items()})
    print(json.dumps(out, indent=1))
    print(f"wrote {dest / 'tracks.json'} and tracks.npz")
    return 0


# ------------------------------------------------------------------ checks on synthetic orbits
def selfcheck(check) -> None:
    """Orbits with a prescribed distance to the front of the synthetic model of census.py."""
    model = census._synthetic_model()
    rng = np.random.default_rng(11)
    t = np.arange(0.0, 400.5, 0.5)
    n_each = 40
    d_i, b0 = model.d_i, model.b0

    def profile(kind):
        ramp = lambda a, b, t0, t1: a + (b - a) * np.clip((t - t0) / (t1 - t0), 0.0, 1.0)  # noqa: E731
        if kind == 0:  # stays upstream
            return 30.0 + 0.0 * t
        if kind == 1:  # one encounter: delta <= d_i from t = 125.2 to 194.9
            return np.where(t < 190, ramp(30.0, -1.0, 100, 130), ramp(-1.0, 40.0, 190, 230))
        if kind == 2:  # two encounters
            d = np.where(t < 120, ramp(30.0, 0.0, 60, 90), ramp(0.0, 12.0, 120, 150))
            return np.where(t < 250, d, np.where(t < 300, ramp(12.0, 0.0, 200, 230), ramp(0.0, 30.0, 260, 300)))
        if kind == 3:  # transmitted
            return ramp(30.0, -30.0, 100, 160)
        return np.where(t < 220, ramp(30.0, -20.0, 100, 150), ramp(-20.0, 30.0, 250, 300))  # returned

    kinds = np.repeat(np.arange(5), n_each)
    n = len(kinds)
    theta = np.radians(rng.uniform(35.0, 80.0, n)) * np.ones((len(t), 1))
    az = rng.uniform(0, 2 * np.pi, n)
    delta = np.stack([profile(k) for k in kinds], axis=1)
    r = model.front(theta, t[:, None] * np.ones((1, n)))[0] + delta
    c = model.centre
    pos = np.stack([c[0] + r * np.sin(theta) * np.cos(az), c[1] + r * np.sin(theta) * np.sin(az),
                    c[2] + r * np.cos(theta)])
    inside = delta <= d_i
    bz = b0 * np.where(inside, 2.0, 1.0)
    b = np.stack([0.0 * bz, 0.0 * bz, bz])
    after = np.cumsum(inside, axis=0) > 0
    u_perp = 0.1 * np.sqrt(bz / b0)  # the magnetic moment is conserved
    u_par = np.where(after & ~inside, 1.2, 0.05)
    phase = rng.uniform(0, 2 * np.pi, n) + b0 * t[:, None]
    u = np.stack([u_perp * np.cos(phase), u_perp * np.sin(phase), u_par])
    res = measure(model, pos, u, b, np.zeros_like(u), t, 0.25)
    want_cat = np.array([0, 1, 1, 3, 4])[kinds]
    want_enc = np.array([0, 1, 2, 1, 1])[kinds]
    one = kinds == 1
    with np.errstate(invalid="ignore"):
        ratio = res["mu_out"] / res["mu_in"]
    ok = (np.array_equal(res["cat"], want_cat) and np.array_equal(res["n_enc"], want_enc)
          and np.all(np.abs(res["dwell"][one] - 69.7) < 1.5) and np.all(np.abs(ratio[one] - 1.0) < 1e-6)
          and np.all(np.abs(res["gamma_out"][one] - np.sqrt(1.0 + 1.44 + 0.01)) < 1e-6)
          and np.all((res["t_min"][one] >= 129.0) & (res["t_min"][one] <= 191.0)))
    check("orbits, sides and encounters", ok,
          f"{int((res['cat'] == want_cat).sum())}/{n} categories, {int((res['n_enc'] == want_enc).sum())}/{n} "
          f"encounter counts; dwell {np.median(res['dwell'][one]):.1f} (69.7 between the two passages of d_i); "
          f"mu_out/mu_in - 1 at most {np.nanmax(np.abs(ratio[one] - 1.0)):.1e}")
    # the join and f_FN: the census calls 30 of the 40 reflections explained, 4 residual, 6 undetermined
    ids = np.arange(n, dtype=np.uint64) * 947
    j_one = np.flatnonzero(one)
    result = np.array([1] * 30 + [6] * 4 + [5] * 6)
    cen = {"id": np.concatenate([ids[j_one], ids[j_one]]), "dump": np.array([10] * 40 + [16] * 40),
           "result": np.concatenate([np.full(40, 6), result]), "cls": np.zeros(80, int), "gamma": np.full(80, 1.6),
           "t_r": np.full(80, 160.0), "tbn_r": np.full(80, 1.0), "k1": np.zeros(80, bool), "drho": np.zeros(80)}
    joined = join_census(ids, res, cen, 18.5 * np.arange(41))  # dump 10 (t = 185) is before the exit
    out = summarize(res, joined, 66.0)
    row = out["reflections"]["B1'"]["gamma_out >= 1.5"]
    ok = (row["explained"], row["not explained"], row["undetermined"]) == (30, 4, 6) \
        and abs(row["f_FN"] - 4 / 34) < 1e-12 and abs(binomial_upper(0, 30) - 0.095) < 0.002 \
        and out["K2 several encounters, gamma_out >= 1.5"]["n"] == n_each \
        and out["reflections, all energies"]["mu within 20 %"] == n_each
    check("orbits, join with the census and f_FN", ok,
          f"explained/not/undetermined {row['explained']}/{row['not explained']}/{row['undetermined']}, f_FN "
          f"{row['f_FN']:.3f} (4/34), upper bound for 0 of 30: {binomial_upper(0, 30):.3f} (0.095); "
          f"{out['K2 several encounters, gamma_out >= 1.5']['n']} orbits with two encounters")
