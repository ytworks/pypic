"""Count model B2 on the measured front of a run against its census (exploratory-protocol.md Sect. 3.4),
and the energies before the reflection of the electrons that B1' explains (Sect. 3.2).

Upper bound: the measured compression b(theta, t), phi = 0, the measured upstream temperature.
Lower bound: the initial temperature and b lower by its scatter, with the measured phi (analysis/phi.npz of the
command phi) where its error is below 50 %; elsewhere phi = 0 and the lower bound is too high there.
The share of all encounters that the model reflects is compared with that of the tracked orbits.
The energy bin gamma 1.25-1.5 lies below the threshold fixed before the runs; it is carried as an extra.
"""
from __future__ import annotations

import json
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import census
import numpy as np
import potential
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
    front, t_grid, theta, t_dumps, temp, b, phi, density, z_half, mask = job
    # a first bin from gamma = 1 holds the reflections below the bins of the census
    res = sk.reflection_counts(front, t_grid, theta, t_dumps, temp, b, phi, density, z_half,
                               [1.0, *census.E_EDGES], BAND_EDGES_DEG, valid=mask)
    return res["counts"][..., 1:], res["reflected"], res["encountered"]


def model_counts(model, t_dumps, sigma_b, phi_map=None, t_step: float = 3.7, workers: int = 4):
    """Counts of B2 [bound (upper, lower), domain, dump, band, energy bin]; the reflections of all energies and
    the encounters [bound, domain, band]. phi_map(theta, t) -> (phi, source): used in the lower bound."""
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
    phi, src = phi_map(th, tt) if phi_map is not None else (np.zeros_like(tt), np.zeros(tt.shape, int))
    jobs = [(model.front_b, t_grid, theta, t_dumps, tp, np.where(has, bb, 1.0001), ph, density, z_half, m)
            for tp, bb, ph in ((temp_up, b, 0.0), (temp_0 * np.ones_like(tt), b_low, phi)) for m in masks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        res = list(pool.map(_case, jobs))
    counts = np.array([r[0] for r in res]).reshape(2, 2, len(t_dumps), len(BAND_EDGES_DEG) - 1, -1)
    reflected, encountered = (np.array([r[k] for r in res]).reshape(2, 2, -1) for k in (1, 2))
    info = {"t_grid": [float(t_grid[0]), float(t_grid[-1]), t_step], "T_initial": temp_0, "density": density,
            "z_half": float(z_half),
            "phi_cells": [{"domain": name, "cells": int(m.sum()), "windowed": int((m & (src == 2)).sum()),
                           "mean of the bin": int((m & (src == 1)).sum()), "not measured": int((m & (src == 0)).sum()),
                           "median phi where measured": float(np.median(phi[m & (src > 0)]))
                           if (m & (src > 0)).any() else None} for name, m in zip(DOMAINS, masks)]}
    return counts, reflected, encountered, info


def track_shares(run_dir, model):
    """Tracked orbits that met the front, by the band of the foot of their field line when they entered the
    layer and by their end: [domain, band, (back upstream, transmitted, returned from downstream, in the layer at
    the end, back upstream after a sample behind the density front)]. None without analysis/tracks.npz."""
    path = Path(run_dir) / "analysis" / "tracks.npz"
    if not path.exists():
        return None
    z = np.load(path)
    cat, t_in, rho = z["cat"], z["t_enter"], z["rho_in"]
    deep = (cat == 1) & (z["d_min"] <= 0.0)
    ok = np.isin(cat, (1, 2, 3, 4)) & np.isfinite(t_in) & np.isfinite(rho)
    ok &= (t_in >= model.t_range_b[0]) & (t_in <= model.t_range_b[1])
    z_foot, _, tbn = model.tab_b.at(rho[ok], t_in[ok])
    flag = model.flag(np.arctan2(rho[ok], np.abs(z_foot)), t_in[ok])
    good = np.isfinite(tbn)
    band = np.minimum((np.degrees(np.where(good, tbn, 0.0)) // 10).astype(int), len(BAND_EDGES_DEG) - 2)
    out = np.zeros((2, len(BAND_EDGES_DEG) - 1, 5), int)
    kind = np.array([0, 0, 3, 1, 2])[cat[ok]]  # CATS of orbits.py: 1 back, 2 in the layer, 3 transmitted, 4 returned
    for d, fmin in enumerate((1, 2)):
        s = good & (flag >= fmin)
        np.add.at(out[d], (band[s], kind[s]), 1)
        np.add.at(out[d], (band[s & deep[ok]], 4), 1)
    return out


def share_rows(reflected, encountered, tracks) -> list:
    """Per domain and band: the share of the encounters that the model reflects (upper bound: phi = 0; lower
    bound: with phi) and the share of the tracked orbits that came back upstream, of those with a known end."""
    rows = []
    for d, dom in enumerate(DOMAINS):
        for ib in range(len(BAND_EDGES_DEG) - 1):
            row = {"domain": dom, "band": [BAND_EDGES_DEG[ib], BAND_EDGES_DEG[ib + 1]]}
            for k, bound in enumerate(("upper", "lower")):
                row[f"model {bound}"] = float(reflected[k, d, ib] / encountered[k, d, ib]) \
                    if encountered[k, d, ib] > 0 else None
            if tracks is not None:
                back, done = int(tracks[d, ib, 0]), int(tracks[d, ib, :3].sum())
                row.update({"tracks back": back, "tracks with a known end": done,
                            "tracks in the layer at the end": int(tracks[d, ib, 3]),
                            "tracks back from behind the front": int(tracks[d, ib, 4]),
                            "tracks share": back / done if done >= N_MIN else None,
                            "tracks share error": float(np.sqrt(max(back, 1)) / done) if done >= N_MIN else None})
            rows.append(row)
    return rows


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
    phi_map = potential.PhiMap.load(run_dir)
    counts, reflected, encountered, info = model_counts(model, t_dumps, sigma_b, phi_map, workers=args.workers)
    meas = measured_counts(cen, len(t_dumps))
    pre = list(range(1, len(census.E_EDGES) - 1))
    b_med = float(np.nanmedian(np.where(has, z["b_mean"], np.nan)))
    out = {"run": str(run_dir), "commit": run.info.get("git_commit"), "domains": DOMAINS, "kinds": KINDS,
           "gamma_bins": census.E_NAMES, "band_edges_deg": BAND_EDGES_DEG, "times": t_dumps.tolist(), **info,
           "sigma_b": {"definition": "1.48 MAD of the differences of the bin-average compression between "
                       "successive outputs, over sqrt 2, per angle",
                       "per_angle": sigma_b.tolist(), "median": float(np.median(sigma_b))},
           "phi": "not measured: phi = 0 in both bounds, so the lower bound is too high" if phi_map is None else
                  "measured phi in the lower bound where its error is below 50 %, 0 elsewhere",
           "N3": [{"gamma": r["gamma"], "nu from": r["nu_from"], "nu x tau_dwell": r["nu_used"] * model.tau_dwell}
                  for r in cj["N3b"]],
           "reflected, all energies": {"upper": reflected[0].tolist(), "lower": reflected[1].tolist()},
           "encountered": {"upper": encountered[0].tolist(), "lower": encountered[1].tolist()},
           "share reflected": share_rows(reflected, encountered, track_shares(run_dir, model)),
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
    for r in out.get("phi_cells", []):
        print(f"phi in the lower bound, {r['domain']}: {r['cells']} cells of the model; windowed value in "
              f"{r['windowed']}, mean of the bin in {r['mean of the bin']}, not measured in {r['not measured']}; "
              f"median where measured {r['median phi where measured']}")
    print("\nshare of the encounters that end in a reflection, all energies (model: feet slower than light only): "
          "model with phi = 0 / model of the lower bound | tracked orbits back upstream (of them with a sample "
          "behind the density front) / with a known end (in the layer at the end)")
    for r in out.get("share reflected", []):
        if r["model upper"] is None and not r.get("tracks with a known end"):
            continue
        fmt = "  ".join("   n/a" if r[k] is None else f"{r[k]:6.4f}" for k in ("model upper", "model lower"))
        line = f"  {r['domain']:14s} {r['band'][0]:2d}-{r['band'][1]:2d}  {fmt}"
        if "tracks back" in r:
            share = "n/a" if r["tracks share"] is None else f"{r['tracks share']:.4f} +- {r['tracks share error']:.4f}"
            line += (f" | {r['tracks back']:5d} ({r['tracks back from behind the front']:5d}) / "
                     f"{r['tracks with a known end']:6d} "
                     f"({r['tracks in the layer at the end']:6d})  {share}")
        print(line)
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
    # the share of the encounters that end in a reflection: the bin from gamma = 1 leaves the counts as they
    # were; a potential lowers the share; the shares of the tracked orbits by the band of their foot
    front, t_grid, theta = sk.sphere(20.0, 0.2), np.arange(0.5, 100, 1.0), np.radians(np.arange(0.5, 90, 1.0))
    t_dumps, mask = np.array([50.0, 100.0]), np.ones((100, 90), bool)
    c0, r0, e0 = _case((front, t_grid, theta, t_dumps, 0.01, 2.5, 0.0, 1.0, 60.0, mask))
    c1, r1, e1 = _case((front, t_grid, theta, t_dumps, 0.01, 2.5, np.full((100, 90), 0.02), 1.0, 60.0, mask))
    ref = sk.reflection_counts(front, t_grid, theta, t_dumps, 0.01, 2.5, 0.0, 1.0, 60.0, census.E_EDGES,
                               BAND_EDGES_DEG, valid=mask)
    one = sk.reflected_flux(0.2 / np.cos(np.radians(45.5)), 0.01, 2.5, 0.0, [1.0, np.inf], np.linspace(0, 1, 201))
    rows = share_rows(np.array([[r0, r0], [r1, r1]]), np.array([[e0, e0], [e1, e1]]), None)
    lo, hi = rows[4]["model lower"], rows[4]["model upper"]  # 40-50 degrees, detected front

    class _Model:
        t_range_b = (0.0, 100.0)

        class tab_b:  # noqa: N801
            @staticmethod
            def at(rho, t):
                return 10.0 * np.ones_like(rho), 0.3 * np.ones_like(rho), np.radians(rho)

        @staticmethod
        def flag(theta, t):
            return np.where(t > 50.0, 2, 1)

    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "analysis").mkdir()
        # 600 orbits with the foot at 45 degrees: 60 back, 440 transmitted, 100 in the layer; 300 of them late;
        # 5 that never met the front and 5 that entered before the range of the fit
        cat = np.array([1] * 60 + [3] * 440 + [2] * 100 + [0] * 5 + [1] * 5)
        t_in = np.concatenate([np.tile([20.0, 80.0], 300), np.full(5, np.nan), np.full(5, -5.0)])
        np.savez(Path(tmp) / "analysis" / "tracks.npz", cat=cat, t_enter=t_in, rho_in=np.full(610, 45.0),
                 d_min=np.where(np.arange(610) < 20, -1.0, 2.0))
        tr = track_shares(tmp, _Model)
    row = share_rows(np.array([[r0, r0], [r1, r1]]), np.array([[e0, e0], [e1, e1]]), tr)[4]
    ok = np.allclose(c0, ref["counts"], rtol=1e-12) and np.all(r0 >= c0[-1].sum(axis=-1)) and 0.0 < lo < 0.7 * hi \
        and abs(hi / (one[0].sum() / one[1]) - 1.0) < 0.1 and tr[0, 4].tolist() == [60, 440, 0, 100, 20] \
        and tr[1, 4, :4].sum() == 300 and tr[..., :4].sum() == 900 and row["tracks share"] == 60 / 500
    check("baseline, share of the encounters that end in a reflection", ok,
          f"40-50 degrees: {hi:.4f} without a potential ({one[0].sum() / one[1]:.4f} at 45.5 degrees), {lo:.4f} with "
          f"phi = 2 T; counts of the census bins unchanged; tracked orbits {tr[0, 4].tolist()} -> "
          f"{row['tracks share']:.3f}")
