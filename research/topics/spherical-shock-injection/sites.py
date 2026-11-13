"""Sites of the driven shock (sites.md): boundaries along each polar-angle ray from the |B| profile, and the
change of the energy of the tracked electrons per site and per local theta_Bn.

Sites, from the centre outwards, per ray (polar angle bin of 5 deg) and output:
  0 cavity            r < r0
  1 wall foot         r0 <= r < r1
  2 steep wall        r1 <= r <= r_Bmax
  3 transition layer  r_Bmax < r <= r_f + d_i
  4 upstream          r > r_f + d_i
r_f is the front of `analyze.py front` (ambient ion density 1.5 n0). r_Bmax is the maximum of the |B| profile of
the ray within 40 d_e inside r_f + d_i. r1: from r_Bmax inwards, the innermost radius of the contiguous run with
|B| >= B0. r0: from r1 inwards, the innermost radius of the contiguous run with |B| > 1.1 B_cav, where B_cav is
the median |B| between r1 - 40 and r1 - 20 d_e.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

NAMES = ("cavity", "wall foot", "steep wall", "transition layer", "upstream")
DR = 2.0             # radial bin of the profiles [d_e]
CAV_FACTOR = 1.1     # the foot starts where |B| exceeds the cavity level by this factor
CAV_RANGE = (40.0, 20.0)  # the cavity level is the median |B| between r1 - 40 and r1 - 20 d_e
SEARCH = 40.0        # r_Bmax is searched within this distance inside r_f + d_i
TBN_EDGES = np.arange(0.0, 90.1, 15.0)


def boundaries(r, b, r_f, d_i):
    """r0, r1, r_Bmax of one ray: r ascending radii, b the |B| / B0 profile on them, r_f the front."""
    nan = (np.nan, np.nan, np.nan)
    if not np.isfinite(r_f):
        return nan
    win = (r >= r_f + d_i - SEARCH) & (r <= r_f + d_i) & np.isfinite(b)
    if not win.any():
        return nan
    k = np.flatnonzero(win)[np.argmax(b[win])]
    if b[k] < 1.0:
        return nan
    j = k
    while j > 0 and np.isfinite(b[j - 1]) and b[j - 1] >= 1.0:
        j -= 1
    r1 = r[j]
    cav_sel = (r >= r1 - CAV_RANGE[0]) & (r <= r1 - CAV_RANGE[1]) & np.isfinite(b)
    if not cav_sel.any():
        return (np.nan, r1, r[k])
    level = CAV_FACTOR * np.median(b[cav_sel])
    i = j
    while i > 0 and np.isfinite(b[i - 1]) and b[i - 1] > level:
        i -= 1
    return (r[i], r1, r[k])


def site_of(r, r0, r1, rb, rf, d_i):
    """Site index (NAMES) of radii r against the boundaries (arrays broadcast); -1 where a boundary is missing."""
    r, r0, r1, rb, rf = np.broadcast_arrays(*(np.asarray(a, float) for a in (r, r0, r1, rb, rf)))
    s = np.full(r.shape, -1, np.int8)
    ok = np.isfinite(r) & np.isfinite(r0) & np.isfinite(r1) & np.isfinite(rb) & np.isfinite(rf)
    s[ok & (r < r0)] = 0
    s[ok & (r >= r0) & (r < r1)] = 1
    s[ok & (r >= r1) & (r <= rb)] = 2
    s[ok & (r > rb) & (r <= rf + d_i)] = 3
    s[ok & (r > rf + d_i)] = 4
    return s


def ray_profiles(bmag, x, z, theta_edges, r_edges):
    """Median of a meridional map bmag[x, z] in (polar angle, radius) bins; x, z are cell centres from the
    centre. Both hemispheres and both signs of x are folded."""
    X, Z = np.meshgrid(np.abs(x), np.abs(z), indexing="ij")
    rr = np.hypot(X, Z).ravel()
    th = np.arctan2(X, Z).ravel()
    it = np.digitize(th, theta_edges) - 1
    ir = np.digitize(rr, r_edges) - 1
    ok = (it >= 0) & (it < len(theta_edges) - 1) & (ir >= 0) & (ir < len(r_edges) - 1)
    out = np.full((len(theta_edges) - 1, len(r_edges) - 1), np.nan)
    key = it[ok] * (len(r_edges) - 1) + ir[ok]
    val = bmag.ravel()[ok]
    order = np.argsort(key, kind="stable")
    key, val = key[order], val[order]
    starts = np.flatnonzero(np.r_[True, np.diff(key) != 0])
    ends = np.r_[starts[1:], len(key)]
    for a, e in zip(starts, ends):
        out.flat[key[a]] = np.median(val[a:e])
    return out


def cmd_sites(args) -> int:
    from census import _open
    run_dir = Path(args.run)
    run = _open(run_dir)
    pl = run.meta["planner"]
    b0, d_i, c = float(pl["b0"]), float(pl["d_i"]), np.asarray(pl["center"], float)
    fz = np.load(run_dir / "analysis" / "front.npz")
    theta = fz["theta"]
    dth = theta[1] - theta[0]
    theta_edges = np.r_[theta - 0.5 * dth, theta[-1] + 0.5 * dth]
    s = run.series("slice_y").f
    n = s["Bx"].shape[1]
    dx = float(run.config["dx"][0]) if hasattr(run.config["dx"], "__len__") else float(run.config["dx"])
    xc = (np.arange(n) + 0.5) * dx - c[0]
    zc = (np.arange(s["Bx"].shape[2]) + 0.5) * dx - c[2]
    r_edges = np.arange(0.0, float(np.hypot(xc.max(), zc.max())) + DR, DR)
    rc = 0.5 * (r_edges[1:] + r_edges[:-1])
    t = np.asarray(s["t"], float)
    if not np.allclose(t, fz["t"]):
        raise SystemExit("slice_y and front.npz are not at the same times")
    nt, nth = len(t), len(theta)
    prof = np.full((nt, nth, len(rc)), np.nan, np.float32)
    r0, r1, rb = (np.full((nt, nth), np.nan) for _ in range(3))
    for i in range(nt):
        bmag = np.sqrt(sum(np.asarray(s[k][i], float) ** 2 for k in ("Bx", "By", "Bz"))) / b0
        prof[i] = ray_profiles(bmag, xc, zc, theta_edges, r_edges)
        for j in range(nth):
            r0[i, j], r1[i, j], rb[i, j] = boundaries(rc, prof[i, j].astype(float), fz["r"][i, j], d_i)
    dest = run_dir / "analysis" / "sites.npz"
    np.savez_compressed(dest, t=t, theta=theta, r_centres=rc, profile=prof, r0=r0, r1=r1, r_bmax=rb, r_f=fz["r"],
                        d_i=d_i)
    eq = int(np.argmin(np.abs(np.degrees(theta) - 87.5)))
    print("equator (87.5 deg): t, r0, r1, r_Bmax, r_f [d_e]")
    for i in range(0, nt, 4):
        print(f"  {t[i]:6.1f}  {r0[i, eq]:6.1f} {r1[i, eq]:6.1f} {rb[i, eq]:6.1f} {fz['r'][i, eq]:6.1f}")
    print(f"wrote {dest}")
    return 0


def local_tbn(x, y, z, bx, by, bz, theta, r_f_row):
    """theta_Bn [deg] between the local B and the normal of the front surface r_f(theta) at the electron's
    polar angle (both hemispheres folded). r_f_row: r_f on the theta bins at that time."""
    rho = np.hypot(x, y)
    pol = np.arctan2(rho, np.abs(z))
    ok = np.isfinite(r_f_row)
    if ok.sum() < 2:
        return np.full(np.shape(x), np.nan)
    drdth = np.interp(pol, theta[ok], np.gradient(r_f_row[ok], theta[ok]))
    rf = np.interp(pol, theta[ok], r_f_row[ok])
    # normal in the (rho, |z|) half plane: r_hat - (1/r) dr/dtheta theta_hat
    s, c = np.sin(pol), np.cos(pol)
    k = drdth / np.maximum(rf, 1e-9)
    n_rho, n_z = s - k * c, c + k * s
    norm = np.hypot(n_rho, n_z)
    b_rho = (x * bx + y * by) / np.maximum(rho, 1e-9)
    b_z = bz * np.sign(z)
    bm = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
    cos = np.abs(b_rho * n_rho + b_z * n_z) / np.maximum(norm * bm, 1e-30)
    return np.degrees(np.arccos(np.clip(cos, 0.0, 1.0)))


def cmd_sitegain(args) -> int:
    """Net change of gamma of the tracked electrons per site and per local theta_Bn band, up to the time of the
    largest gamma of each electron; groups by gamma_max."""
    from census import _open
    run_dir = Path(args.run)
    run = _open(run_dir)
    pl = run.meta["planner"]
    c, d_i = np.asarray(pl["center"], float), float(pl["d_i"])
    st = np.load(run_dir / "analysis" / "sites.npz")
    theta, ts = st["theta"], st["t"]
    tr = run.tracks("tracks_e", "e")
    t = np.asarray(tr.t, float)
    n = len(tr.ids)
    # gamma_max of every tracked electron, in chunks
    gmax = np.full(n, np.nan)
    for a in range(0, n, args.chunk):
        sl = slice(a, min(n, a + args.chunk))
        u2 = sum(np.asarray(tr.f[k][:, sl], float) ** 2 for k in ("ux", "uy", "uz"))
        gmax[sl] = np.nanmax(np.sqrt(1.0 + u2), axis=0)
    rng = np.random.default_rng(1)
    groups = {"gamma_max >= 2": np.flatnonzero(gmax >= 2.0),
              "gamma_max 1.5-2": np.flatnonzero((gmax >= 1.5) & (gmax < 2.0)),
              "reference, gamma_max < 1.25 (20000 random)": rng.choice(np.flatnonzero(gmax < 1.25),
                                                                       min(20000, int((gmax < 1.25).sum())),
                                                                       replace=False)}
    i_out = np.array([int(np.argmin(np.abs(ts - q))) for q in t])
    out = {"run": str(run_dir), "sites": list(NAMES), "tbn_edges": TBN_EDGES.tolist(), "groups": {}}
    for name, cols in groups.items():
        cols = np.sort(cols)
        D = {k: np.asarray(tr.f[k][:, cols], float) for k in ("x", "y", "z", "ux", "uy", "uz", "Bx", "By", "Bz")}
        x, y, z = (D[k] - c[i] for i, k in enumerate("xyz"))
        g = np.sqrt(1.0 + D["ux"] ** 2 + D["uy"] ** 2 + D["uz"] ** 2)
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        pol = np.arctan2(np.hypot(x, y), np.abs(z))
        ith = np.clip(np.rint((np.nan_to_num(pol) - theta[0]) / (theta[1] - theta[0])).astype(int), 0, len(theta) - 1)
        rows = i_out[:, None]
        site = site_of(r, st["r0"][rows, ith], st["r1"][rows, ith], st["r_bmax"][rows, ith], st["r_f"][rows, ith], d_i)
        tbn = np.full(g.shape, np.nan)
        for k in range(len(t)):
            tbn[k] = local_tbn(x[k], y[k], z[k], D["Bx"][k], D["By"][k], D["Bz"][k], theta, st["r_f"][i_out[k]])
        tmax = np.argmax(np.where(np.isfinite(g), g, 0.0), axis=0)
        win = np.arange(len(t))[:, None] < tmax[None, :]
        dg = np.full(g.shape, np.nan)
        dg[:-1] = g[1:] - g[:-1]
        ok = win & np.isfinite(dg) & (site >= 0)
        net = np.zeros(len(NAMES))
        time = np.zeros(len(NAMES))
        grid = np.zeros((len(NAMES), len(TBN_EDGES) - 1))
        ib = np.clip(np.digitize(np.nan_to_num(tbn, nan=-1.0), TBN_EDGES) - 1, -1, len(TBN_EDGES) - 2)
        for s_ in range(len(NAMES)):
            m = ok & (site == s_)
            net[s_], time[s_] = dg[m].sum(), m.sum()
            for b_ in range(len(TBN_EDGES) - 1):
                grid[s_, b_] = dg[m & (ib == b_)].sum()
        unplaced = win & np.isfinite(dg) & (site < 0)
        out["groups"][name] = {"n": int(len(cols)), "net_dgamma": net.tolist(),
                               "time_share": (time / max(1.0, time.sum())).tolist(),
                               "net_dgamma_by_site_and_tbn": grid.tolist(),
                               "dgamma_unplaced": float(dg[unplaced].sum()),
                               "time_unplaced_share": float(unplaced.sum() / max(1, (win & np.isfinite(dg)).sum()))}
        print(f"\n== {name}: n = {len(cols)}; up to gamma_max; unplaced share of the time "
              f"{out['groups'][name]['time_unplaced_share']:.2f} (net {dg[unplaced].sum():+.1f})")
        for s_ in range(len(NAMES)):
            row = " ".join(f"{v:+7.1f}" for v in grid[s_])
            print(f"   {NAMES[s_]:17s} net {net[s_]:+8.1f}  time {time[s_] / max(1.0, time.sum()):.2f}  "
                  f"by local theta_Bn {TBN_EDGES[0]:.0f}-{TBN_EDGES[-1]:.0f} deg in 15: {row}")
    dest = run_dir / "analysis" / "sitegain.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"wrote {dest}")
    return 0


def selfcheck(check):
    """Boundaries and sites on a constructed profile; theta_Bn on a sphere with a uniform field."""
    r = np.arange(1.0, 160.0, 2.0)
    b = np.full_like(r, 0.5)
    b[(r >= 85) & (r < 99)] = np.linspace(0.6, 0.9, int(((r >= 85) & (r < 99)).sum()))
    up = (r >= 99) & (r <= 109)
    b[up] = np.linspace(1.2, 4.5, int(up.sum()))
    down = (r > 109) & (r <= 117)
    b[down] = np.linspace(3.5, 1.0, int(down.sum()))
    b[r > 117] = 1.0
    r0, r1, rb = boundaries(r, b, 114.0, 4.0)
    s = site_of(np.array([60.0, 90.0, 101.0, 112.0, 130.0]), r0, r1, rb, 114.0, 4.0)
    ok = (r0, r1, rb) == (85.0, 99.0, 109.0) and s.tolist() == [0, 1, 2, 3, 4] \
        and np.isnan(boundaries(r, b, np.nan, 4.0)[0])
    check("sites, boundaries of a constructed |B| profile", ok,
          f"r0, r1, r_Bmax = {r0}, {r1}, {rb} (85, 99, 109); sites of 60, 90, 101, 112, 130 d_e: {s.tolist()}")
    # uniform B along z on a sphere: theta_Bn equals the polar angle (folded to 0-90)
    theta = np.radians(2.5 + 5.0 * np.arange(18))
    pol = np.radians(np.array([10.0, 45.0, 80.0, 100.0]))
    x, z = 120 * np.sin(pol), 120 * np.cos(pol)
    tb = local_tbn(x, np.zeros(4), z, np.zeros(4), np.zeros(4), np.ones(4), theta, np.full(18, 120.0))
    check("sites, local theta_Bn on a sphere in a uniform field", np.allclose(tb, [10, 45, 80, 80], atol=1e-6),
          f"{tb.round(3).tolist()} for polar angles 10, 45, 80, 100 (10, 45, 80, 80)")
