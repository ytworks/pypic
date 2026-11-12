"""Count model B2 on the measured front of a run against its census (exploratory-protocol.md Sect. 3.4),
and the energies before the reflection of the electrons that B1' explains (Sect. 3.2).

Upper bound: the measured compression b(theta, t), phi = 0, the measured upstream temperature.
Lower bound: the initial temperature and b lower by its scatter; phi is not measured, so phi = 0 here too
and the lower bound is too high by the effect of the potential.
The energy bin gamma 1.25-1.5 lies below the threshold fixed before the runs; it is carried as an extra.
"""
from __future__ import annotations

import json
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import census
import numpy as np
import shockkin as sk

BAND_EDGES_DEG = list(range(0, 91, 10))
DOMAINS = ("detected front", "valid domain")
KINDS = ("explained by S0", "undetermined", "all upstream")
N_MIN = 100


def b_scatter(b, has, n_min: int = 6):
    """Scatter of the bin-average compression between successive outputs at one angle, per angle (n_theta,):
    1.48 x the median absolute deviation of the differences of neighbouring outputs about their median, over
    sqrt 2. A slow change of b drops out; the value is that of sigma for independent Gaussian noise. This
    stands for the standard error of b, which the front measurement does not store. Angles with fewer than
    n_min differences take the median of the others."""
    b = np.where(has, np.asarray(b, float), np.nan)
    d = b[1:] - b[:-1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # angles without a front: all NaN
        sig = 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d, axis=0)), axis=0) / np.sqrt(2.0)
    sig = np.where(np.isfinite(d).sum(axis=0) >= n_min, sig, np.nan)
    return np.where(np.isfinite(sig), sig, np.nanmedian(sig) if np.isfinite(sig).any() else 0.0)


def before_moments(v_ht: float, temp: float, b: float, gamma_edges, n_par=320, n_perp=240, u_max=12.0):
    """Per energy bin after one reflection (phi = 0): the reflected flux, the flux-weighted mean kinetic
    energy before the reflection, and the share of the flux that had more than 10 temp before it. The
    quadrature is that of shockkin.reflected_flux."""
    ng = len(gamma_edges) - 1
    if not v_ht < 1.0:
        return np.zeros(ng), np.full(ng, np.nan), np.full(ng, np.nan)
    umax = u_max * np.sqrt(temp)
    dpar, dperp = 2.0 * umax / n_par, umax / n_perp
    u_par = (-umax + dpar * (np.arange(n_par) + 0.5))[:, None]
    u_perp = (dperp * (np.arange(n_perp) + 0.5))[None, :]
    gam = np.sqrt(1.0 + u_par ** 2 + u_perp ** 2)
    weight = np.exp(-(gam - 1.0) / temp) * u_perp * np.maximum(v_ht - u_par / gam, 0.0)
    _, refl, _, g_r = sk.ht_reflect(u_par, u_perp, v_ht, b)
    w = (weight * refl).ravel()
    ig = np.searchsorted(gamma_edges, g_r.ravel(), side="right") - 1
    ok = (w > 0) & (ig >= 0) & (ig < ng)
    e0 = (gam - 1.0).ravel()
    flux = np.bincount(ig[ok], weights=w[ok], minlength=ng)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.bincount(ig[ok], weights=(w * e0)[ok], minlength=ng) / flux
        share = np.bincount(ig[ok], weights=(w * (e0 > 10.0 * temp))[ok], minlength=ng) / flux
    norm = 2.0 * np.pi * dpar * dperp / sk.mj_norm(temp)
    return flux * norm, mean, share


def _case(job):
    front, t_grid, theta, t_dumps, temp, b, density, z_half, mask = job
    res = sk.reflection_counts(front, t_grid, theta, t_dumps, temp, b, 0.0, density, z_half, census.E_EDGES,
                               BAND_EDGES_DEG, valid=mask)
    return res["counts"], res["reflected"]


def model_counts(model, t_dumps, sigma_b, t_step: float = 3.7, workers: int = 4):
    """Counts of B2 [bound (upper, lower), domain, dump, band, energy bin] and the reflections in all."""
    t0, t1 = model.t_range_b
    t_grid = np.arange(t0 + 0.5 * t_step, t1, t_step)
    theta = np.radians(np.arange(0.5, 90.0, 1.0))
    tt, th = np.meshgrid(t_grid, theta, indexing="ij")
    flag = model.flag(th, tt)
    has, b = model.local(th, tt)
    ith = np.clip(np.searchsorted(0.5 * (model.theta[1:] + model.theta[:-1]), theta), 0, len(model.theta) - 1)
    temp_up = np.array([model.temp_upstream(t) for t in t_grid])[:, None] * np.ones_like(tt)
    temp_0 = model.temp_upstream(0.0)
    b_low = np.maximum(b - sigma_b[ith][None, :], 1.0001)
    density = model.ppc / float(np.prod(model.dx))
    z_half = model.n_cells[2] * model.dx[2] - model.centre[2]
    masks = [has & (flag >= 1), has & (flag >= 2)]
    jobs = [(model.front_b, t_grid, theta, t_dumps, tp, np.where(has, bb, 1.0001), density, z_half, m)
            for tp, bb in ((temp_up, b), (temp_0 * np.ones_like(tt), b_low)) for m in masks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        res = list(pool.map(_case, jobs))
    counts = np.array([r[0] for r in res]).reshape(2, 2, len(t_dumps), len(BAND_EDGES_DEG) - 1, -1)
    reflected = np.array([r[1] for r in res]).reshape(2, 2, -1)
    return counts, reflected, {"t_grid": [float(t_grid[0]), float(t_grid[-1]), t_step], "T_initial": temp_0,
                               "density": density, "z_half": float(z_half)}


def measured_counts(cen, n_dumps: int):
    """Upstream electrons of the census [dump, domain, band, energy bin, kind (KINDS)], by the band of the
    reflection point; electrons without one are left out."""
    eb = np.searchsorted(census.E_EDGES, cen["gamma"], side="right") - 1
    ok = (cen["cls"] <= 1) & (cen["band"] < len(BAND_EDGES_DEG) - 1) & (eb >= 0)
    kinds = (np.isin(cen["result"], (1, 2)), cen["result"] == 5, np.ones(len(eb), bool))
    m = np.zeros((n_dumps, 2, len(BAND_EDGES_DEG) - 1, len(census.E_EDGES) - 1, 3))
    for d, f_min in enumerate((1, 2)):
        for k, kind in enumerate(kinds):
            s = ok & (cen["flag"] >= f_min) & kind
            np.add.at(m[:, d, :, :, k], (cen["dump"][s], cen["band"][s], eb[s]), 1)
    return m


def compare(model_c, meas, bins) -> dict:
    """The rule of Sect. 3.4 on the cells (dump, band, energy bin in `bins`) that hold at least N_MIN upstream
    electrons: is the count explained by S0 inside [lower / 2, 2 x upper]? Also with the undetermined added."""
    out = {}
    for d, name in enumerate(DOMAINS):
        up, lo = model_c[0, d][..., bins], model_c[1, d][..., bins]
        s0, und, tot = (meas[:, d][..., bins, k] for k in range(3))
        cells = tot >= N_MIN
        inside = (s0 >= 0.5 * lo) & (s0 <= 2.0 * up)
        inside_u = (s0 + und >= 0.5 * lo) & (s0 <= 2.0 * up)
        with np.errstate(divide="ignore", invalid="ignore"):
            r_up, r_lo = (s0 / up)[cells], (s0 / lo)[cells]
        out[name] = {"cells with >= 100 electrons": int(cells.sum()), "inside": int((cells & inside).sum()),
                     "inside with the undetermined added": int((cells & inside_u).sum()),
                     "below lower / 2": int((cells & (s0 < 0.5 * lo)).sum()),
                     "above 2 x upper": int((cells & (s0 > 2.0 * up)).sum()),
                     "S0 / upper": census_pct(r_up), "S0 / lower": census_pct(r_lo)}
    return out


def census_pct(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return {"n": int(len(a)), **({f"p{q}": float(np.percentile(a, q)) for q in (10, 50, 90)} if len(a) else {})}


def energies_before(cen, model, b_ref: dict) -> list:
    """Sect. 3.2: the energies before the reflection of the upstream electrons explained by S0, per domain
    and energy bin, against the flux-weighted Maxwell-Juttner mean for the same foot speeds and bins."""
    eb = np.searchsorted(census.E_EDGES, cen["gamma"], side="right") - 1
    s0 = (cen["cls"] <= 1) & np.isin(cen["result"], (1, 2)) & np.isfinite(cen["e_before"]) & (cen["v_exit"] < 1.0)
    temp = np.array([model.temp_upstream(t) for t in np.unique(cen["t_r"][s0])])
    temp_of = dict(zip(np.unique(cen["t_r"][s0]).tolist(), temp.tolist(), strict=True))
    cache: dict = {}
    rows = []
    for d, name in enumerate(DOMAINS):
        for ie in range(len(census.E_EDGES) - 1):
            sel = np.flatnonzero(s0 & (cen["flag"] >= d + 1) & (eb == ie))
            row = {"domain": name, "gamma": census.E_NAMES[ie], "n": int(len(sel)), "measurable": len(sel) >= N_MIN}
            if len(sel):
                tp = np.array([temp_of[float(t)] for t in cen["t_r"][sel]])
                e0 = cen["e_before"][sel].astype(float)
                row.update({"mean": float(e0.mean()), "share above 10 T": float((e0 > 10.0 * tp).mean()),
                            "mean T": float(tp.mean())})
                for label, b in b_ref.items():
                    mean, share = np.zeros(len(sel)), np.zeros(len(sel))
                    for j, i in enumerate(sel):
                        key = (round(float(cen["v_exit"][i]) / 0.005) * 0.005, round(tp[j] / 2.5e-4) * 2.5e-4, b)
                        if key not in cache:
                            cache[key] = before_moments(key[0], key[1], key[2], census.E_EDGES)
                        mean[j], share[j] = cache[key][1][ie], cache[key][2][ie]
                    m, s = float(np.nanmean(mean)), float(np.nanmean(share))
                    # a share that neither side resolves (fewer than 5 expected and fewer than 10 found) is not tested
                    n_hot = row["share above 10 T"] * len(sel)
                    tested = s * len(sel) >= 5.0 or n_hot >= 10.0
                    thermal = abs(row["mean"] / m - 1.0) <= 0.3 and (
                        not tested or 0.5 <= row["share above 10 T"] / max(s, 1e-30) <= 2.0)
                    row[label] = {"b": b, "model mean": m, "model share above 10 T": s, "share tested": bool(tested),
                                  "verdict": ("thermal seed" if thermal else "pre-heated seed population")
                                  if row["measurable"] else "not measurable"}
            rows.append(row)
    return rows


def cmd_baseline(args) -> int:
    run_dir = Path(args.run)
    run = census._open(run_dir)
    model = census.FrontModel.from_run(run, run_dir)
    if model.front_b is None:
        print("no fitted mid-ramp front in analysis/front.json: B2 cannot be evaluated")
        return 1
    z = np.load(run_dir / "analysis" / "front.npz")
    cen = dict(np.load(run_dir / "analysis" / "census.npz"))
    cj = json.loads((run_dir / "analysis" / "census.json").read_text())
    t_dumps = cen["times"]
    has = np.isfinite(z["r_b"]) & (z["b_mean"] >= 1.5)
    sigma_b = b_scatter(z["b_mean"], has)
    counts, reflected, info = model_counts(model, t_dumps, sigma_b, workers=args.workers)
    meas = measured_counts(cen, len(t_dumps))
    pre = list(range(1, len(census.E_EDGES) - 1))
    b_med = float(np.nanmedian(np.where(has, z["b_mean"], np.nan)))
    out = {"run": str(run_dir), "commit": run.info.get("git_commit"), "domains": DOMAINS, "kinds": KINDS,
           "gamma_bins": census.E_NAMES, "band_edges_deg": BAND_EDGES_DEG, "times": t_dumps.tolist(), **info,
           "sigma_b": {"definition": "1.48 MAD of the differences of the bin-average compression between "
                       "successive outputs, over sqrt 2, per angle",
                       "per_angle": sigma_b.tolist(), "median": float(np.median(sigma_b))},
           "phi": "not measured: phi = 0 in both bounds, so the lower bound is too high",
           "N3": [{"gamma": r["gamma"], "nu from": r["nu_from"], "nu x tau_dwell": r["nu_used"] * model.tau_dwell}
                  for r in cj["N3b"]],
           "reflected_in_all": {"upper": reflected[0].sum(axis=-1).tolist(),
                                "lower": reflected[1].sum(axis=-1).tolist()},
           "model_counts": counts.tolist(), "measured": meas.tolist(),
           "rule gamma >= 1.5": compare(counts, meas, pre), "rule gamma 1.25-1.5 (extra)": compare(counts, meas, [0]),
           "energies_before": energies_before(cen, model, {"b_max": float(model.b_max), "b_median": b_med})}
    (run_dir / "analysis" / "baseline.json").write_text(json.dumps(out, indent=1))
    report(out)
    print(f"wrote {run_dir / 'analysis' / 'baseline.json'}")
    return 0


def report(out) -> None:
    c, m, t = np.array(out["model_counts"]), np.array(out["measured"]), np.array(out["times"])
    print(f"B2 of {out['run']}: T initial {out['T_initial']:.4f}, sigma_b median {out['sigma_b']['median']:.3f}, "
          f"{out['density']:.0f} macro-electrons per d_e^3; {out['phi']}")
    print("N3, nu x tau_dwell:", "; ".join(f"gamma {r['gamma']} {r['nu x tau_dwell']:.3f} ({r['nu from']})"
                                          for r in out["N3"]))
    for d, name in enumerate(out["domains"]):
        print(f"\n{name}: electrons in the box, all bands; lower / upper bound | explained by S0 / undetermined / "
              "all upstream")
        for label, sl in (("gamma 1.25-1.5 (extra)", slice(0, 1)), ("gamma >= 1.5", slice(1, None))):
            print(f"  {label}")
            for i in range(4, len(t), 3):
                lo, up = c[1, d, i][:, sl].sum(), c[0, d, i][:, sl].sum()
                k = m[i, d][:, sl].sum(axis=(0, 1))
                print(f"    t {t[i]:6.1f}  {lo:9.1f} / {up:9.1f} | {k[0]:6.0f} / {k[1]:6.0f} / {k[2]:6.0f}")
        i = int(np.argmax(m[:, d, :, 0, 2].sum(axis=1)))
        print(f"  by band at t = {t[i]:.1f}, gamma 1.25-1.5: lower / upper | S0 / undetermined / all")
        for ib in range(c.shape[3]):
            if c[0, d, i, ib, 0] >= 0.5 or m[i, d, ib, 0, 2] > 0:
                k = m[i, d, ib, 0]
                print(f"    {BAND_EDGES_DEG[ib]:2d}-{BAND_EDGES_DEG[ib + 1]:2d}  {c[1, d, i, ib, 0]:8.1f} / "
                      f"{c[0, d, i, ib, 0]:8.1f} | {k[0]:5.0f} / {k[1]:5.0f} / {k[2]:5.0f}")
    for key in ("rule gamma >= 1.5", "rule gamma 1.25-1.5 (extra)"):
        print(f"\n{key}")
        for name, r in out[key].items():
            print(f"  {name}: {r}")
    print("\nenergies before the reflection (electrons explained by S0)")
    for r in out["energies_before"]:
        if r["n"]:
            ref = "; ".join(f"{k}: model {r[k]['model mean']:.4f}, {r[k]['model share above 10 T']:.3f} -> "
                            f"{r[k]['verdict']}" for k in ("b_max", "b_median"))
            print(f"  {r['domain']:15s} gamma {r['gamma']:8s} n {r['n']:6d} mean {r['mean']:.4f} "
                  f"share > 10 T {r['share above 10 T']:.3f} | {ref}")


# ------------------------------------------------------------------ checks on synthetic data
def selfcheck(check) -> None:
    rng = np.random.default_rng(5)
    # the scatter of b: a smooth decay with Gaussian noise of 0.05
    n_t, n_th = 41, 18
    smooth = 3.0 - 1.2 * np.linspace(0.0, 1.0, n_t)[:, None] + 0.2 * np.linspace(0.0, 1.0, n_th)[None, :]
    has = np.ones((n_t, n_th), bool)
    has[:, 0], has[:38, 1] = False, False
    sig = b_scatter(smooth + rng.normal(0.0, 0.05, (n_t, n_th)), has)
    ok = abs(np.median(sig[2:]) / 0.05 - 1.0) < 0.15 and np.all(sig[:2] == np.median(sig[2:]))
    check("baseline, scatter of the compression", ok,
          f"median {np.median(sig[2:]):.4f} for a noise of 0.05 (range {sig[2:].min():.3f}-{sig[2:].max():.3f}); "
          "angles with too few values take the median")
    # the energies before the reflection against a weighted Monte-Carlo sample
    v_ht, temp, b = 0.3, 0.01, 2.5
    edges = np.array([1.0, 1.1, 1.25, np.inf])
    flux, mean, share = before_moments(v_ht, temp, b, edges)
    s = np.sqrt(1.3 * temp)
    u = rng.normal(0.0, s, (3, 2_000_000))
    gam = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
    u_perp = np.hypot(u[0], u[1])
    w = np.exp(-(gam - 1.0) / temp + np.sum(u ** 2, axis=0) / (2.0 * s ** 2)) * np.maximum(v_ht - u[2] / gam, 0.0)
    _, refl, _, g_r = sk.ht_reflect(u[2], u_perp, v_ht, b)
    ig = np.searchsorted(edges, g_r, side="right") - 1
    mc_mean = np.array([np.sum((w * (gam - 1.0))[refl & (ig == k)]) / np.sum(w[refl & (ig == k)]) for k in range(3)])
    mc_share = np.array([np.sum(w[refl & (ig == k) & (gam - 1.0 > 10 * temp)]) / np.sum(w[refl & (ig == k)])
                         for k in range(3)])
    ref = sk.reflected_flux(v_ht, temp, b, 0.0, edges, np.linspace(0.0, 1.0, 201))[0].sum(axis=1)
    ok = np.allclose(flux, ref, rtol=1e-9) and np.all(np.abs(mean / mc_mean - 1.0) < 0.03) \
        and np.all(np.abs(share - mc_share) < 0.01)
    check("baseline, energies before the reflection", ok,
          f"mean {np.round(mean, 4).tolist()} (Monte Carlo {np.round(mc_mean, 4).tolist()}), share above 10 T "
          f"{np.round(share, 3).tolist()} ({np.round(mc_share, 3).tolist()}); flux equal to reflected_flux")
    # the rule of Sect. 3.4
    c = np.zeros((2, 2, 1, 1, 2))
    c[0], c[1] = 400.0, 150.0
    m = np.zeros((1, 2, 1, 2, 3))
    m[0, :, 0, 0] = (300.0, 20.0, 340.0)  # inside
    m[0, :, 0, 1] = (60.0, 5.0, 130.0)  # below lower / 2
    r = compare(c, m, [0, 1])["valid domain"]
    ok = (r["cells with >= 100 electrons"], r["inside"], r["below lower / 2"], r["above 2 x upper"]) == (2, 1, 1, 0)
    check("baseline, rule of the bounds", ok, f"{r['inside']} of {r['cells with >= 100 electrons']} cells inside, "
          f"{r['below lower / 2']} below half the lower bound")
