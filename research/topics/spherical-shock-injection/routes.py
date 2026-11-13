"""M6 of global-plan.md: the route by which the most energetic electrons of a run reach their energy.

analyze.py routes <run>   writes analysis/routes.json (needs campaign_hist.npz, campaign_tracks.npz,
                          campaign_entries.npy and front.npz from analyze.py campaign)

Electrons: the top 1 % by the gyration-averaged drift-frame energy (level_drift_frame of campaign_tracks.npz) of
the tracked electrons born outside the initial debris shell that have a history. The gain G and the window run
from the start to the maximum of the lab gamma smoothed over SMOOTH_BLOCKS blocks (about one gyration), so that
a gyration phase in a moving downstream does not count as a gain. A *step* is a run of history blocks (every 4th
sample, work summed per block) at one site
and one 10-degree polar band; its gain is the summed work of all fields. Steps with a gain of 0.05 or more carry
the energy. Reflections are the entries of the entry table that return upstream with gamma_out > 1.02 gamma_in.
The signatures (first match in the order R5, R1, R3a, R3c, R3b, R2, R4; else "other"); G is the total gain,
gamma_max - gamma at the start:
  R5   one reflection gives 70 % of G or more (gamma_out - gamma_in), at V_HT >= 0.5 c
  R1   two or more reflections that together give 50 % of G or more, the later at a polar angle larger by 5
       degrees or more, with a larger gain
  R3a  two or more reflections that together give 50 % of G or more, within 5 degrees of polar angle
  R3c  two or more cycles upstream -> downstream (wall or cavity) -> upstream, and 50 % of G or more gained
       between the first and the last crossing
  R3b  50 % of G or more in transition-layer steps lasting 10 time units or more, with p_perp^2 / B changing by
       a factor 2 or more over them
  R2   50 % of G or more in wall-foot, steep-wall or cavity steps
The share of the gain done by E_phi, E_par, E_r, E_theta is reported per route (diagnostic, not a condition).
Revised 2026-11-13 on P5 (the pilot, not a campaign run) before any campaign run was classified: the first
version selected by the lab gamma_max and asked E_phi >= 50 % and a first step in the layer for R2/R3b; 40 % of
P5's top electrons fell in "other", most of them gaining on the steep wall. Then "unplaced" was added (50 % of G
or more where the sites are not defined: latitudes without a magnetic wall, mostly near the poles), and
the gain of a step was changed from the recorded work to the change of the smoothed gamma over the step: the
work recorded every 2 steps closes on only part of the energetic orbits (block work over the gamma change,
median 0.66 in P5), while the gamma changes add up to G exactly. The work by component stays a diagnostic,
reported for electrons whose work closes within 30 %. Frozen after this.
  R4   three or more gaining steps at two or more sites, losing steps that sum to 30 % of the gaining ones or
       more, and no step with 50 % of G
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROUTES = ("R5", "R1", "R3a", "R3c", "R3b", "R2", "R4", "unplaced", "other")
STEP_MIN = 0.05
SMOOTH_BLOCKS = 6


def steps_of(site, band, gamma, w_phi):
    """Steps of one electron: list of (site, band, start index, end index, gain, work by E_phi); the gain of a
    step is the change of gamma from its first sample to the first sample of the next step (to the last sample
    for the last step), so that the gains add up to gamma[-1] - gamma[0]."""
    out = []
    a = 0
    n = len(site)
    for b in range(1, n + 1):
        if b == n or site[b] != site[a] or band[b] != band[a]:
            e = min(b, n - 1)
            out.append((int(site[a]), int(band[a]), a, b, float(gamma[e] - gamma[a]), float(np.sum(w_phi[a:b]))))
            a = b
    return out


def classify(G, steps, refl, cycles, dt_block, mu_ratio_layer):
    """Route of one electron. refl: list of (polar_deg, gain, V_HT) in time order; cycles: (number of
    upstream -> downstream -> upstream cycles, gain between the first and last crossing)."""
    if G <= 0:
        return "other"
    gains = [r[1] for r in refl]
    if refl and max(gains) >= 0.7 * G:
        k = int(np.argmax(gains))
        if refl[k][2] >= 0.5:
            return "R5"
    if len(refl) >= 2 and sum(gains) >= 0.5 * G:
        later = [(i, j) for i in range(len(refl)) for j in range(i + 1, len(refl))
                 if refl[j][0] >= refl[i][0] + 5.0 and refl[j][1] > refl[i][1]]
        if later:
            return "R1"
        if max(r[0] for r in refl) - min(r[0] for r in refl) <= 5.0:
            return "R3a"
    if cycles[0] >= 2 and cycles[1] >= 0.5 * G:
        return "R3c"
    lay = [s for s in steps if s[0] == 3 and (s[3] - s[2]) * dt_block >= 10.0 and s[4] > 0]
    if sum(s[4] for s in lay) >= 0.5 * G and mu_ratio_layer >= 2.0:
        return "R3b"
    if sum(s[4] for s in steps if s[0] in (0, 1, 2)) >= 0.5 * G:
        return "R2"
    if sum(s[4] for s in steps if s[0] < 0) >= 0.5 * G:
        return "unplaced"
    gain_steps = [s for s in steps if s[4] >= STEP_MIN]
    loss = -sum(s[4] for s in steps if s[4] < 0)
    if (len(gain_steps) >= 3 and len({s[0] for s in gain_steps}) >= 2
            and loss >= 0.3 * sum(s[4] for s in gain_steps) and max(s[4] for s in gain_steps) < 0.5 * G):
        return "R4"
    return "other"


def cycles_of(site, gamma):
    """Upstream -> downstream (sites 0-2) -> upstream cycles; gain between the first and the last crossing."""
    idx = [i for i in range(len(site)) if site[i] in (0, 1, 2, 4)]
    seq = [(i, 4 if site[i] == 4 else 0) for i in idx]
    marks = [seq[0]] if seq else []
    for q in seq[1:]:
        if q[1] != marks[-1][1]:
            marks.append(q)
    n_cyc = sum(1 for k in range(len(marks) - 2) if marks[k][1] == 4 and marks[k + 1][1] == 0 and marks[k + 2][1] == 4)
    if len(marks) >= 3:
        return n_cyc, float(gamma[marks[-1][0]] - gamma[marks[1][0]])
    return n_cyc, 0.0


def cmd_routes(args) -> int:
    run_dir = Path(args.run)
    an = run_dir / "analysis"
    with np.load(an / "campaign_hist.npz") as f:
        h = {k: f[k] for k in f.files}      # decompress once; per-electron access of an npz decompresses again
    tk = np.load(an / "campaign_tracks.npz")
    ent = np.load(an / "campaign_entries.npy")
    fz = np.load(an / "front.npz")
    theta, rd0 = fz["theta"], fz["r_debris"][0]
    cols = h["col"]
    t = h["t"]
    dt_block = float(t[1] - t[0])
    inside = tk["birth_r"] < np.interp(np.radians(tk["birth_polar"]), theta, rd0)
    gmax = tk["gmax"]
    have = np.zeros(len(gmax), bool)
    have[cols] = True
    level = tk["level_drift_frame"]
    cand = np.flatnonzero(have & ~inside & np.isfinite(level))
    cut = np.percentile(level[cand], 99.0)
    top = cand[level[cand] >= cut]
    pos = {int(c): i for i, c in enumerate(cols)}
    refl_by = {}
    for row in ent:
        if row[7] == 1 and row[9] > 1.02 * row[8]:
            refl_by.setdefault(int(row[0]), []).append((row[1], row[2], row[9] - row[8], row[4]))
    counts = {r: 0 for r in ROUTES}
    per = []
    for c in top:
        j = pos[int(c)]
        site = h["site"][:, j]
        x, y, z = h["x"][:, j], h["y"][:, j], h["z"][:, j]
        pol = np.degrees(np.arctan2(np.hypot(x, y), np.abs(z)))
        band = np.clip(np.nan_to_num(pol, nan=-10.0) // 10, -1, 8).astype(int)
        u2 = h["ux"][:, j] ** 2 + h["uy"][:, j] ** 2 + h["uz"][:, j] ** 2
        g = np.sqrt(1.0 + u2)
        ok = np.isfinite(g)
        if ok.sum() < 3:
            continue
        kern = np.ones(SMOOTH_BLOCKS) / SMOOTH_BLOCKS
        gs = np.convolve(np.where(ok, g, 0.0), kern, mode="same") / np.maximum(
            np.convolve(ok.astype(float), kern, mode="same"), 1e-9)
        gs = np.where(ok, gs, np.nan)
        k0 = int(np.flatnonzero(ok)[0])
        km = int(np.nanargmax(np.where(ok, gs, -1.0)))
        G = float(gs[km] - gs[k0])
        sl = slice(k0, km + 1)
        steps = steps_of(site[sl], band[sl], gs[sl], h["W_phi"][sl, j])
        refl = [(r[1], r[2], r[3]) for r in sorted(refl_by.get(int(c), [])) if t[k0] <= r[0] <= t[km]]
        cyc = cycles_of(site[sl], g[sl])
        B = np.sqrt(h["Bx"][:, j] ** 2 + h["By"][:, j] ** 2 + h["Bz"][:, j] ** 2)
        upar = (h["ux"][:, j] * h["Bx"][:, j] + h["uy"][:, j] * h["By"][:, j] + h["uz"][:, j] * h["Bz"][:, j]) / B
        mu = (u2 - upar ** 2) / B
        lay = np.flatnonzero(site[sl] == 3) + k0
        mu_ratio = float(np.nanmax(mu[lay]) / max(np.nanmin(mu[lay]), 1e-12)) if len(lay) >= 2 else 1.0
        route = classify(G, steps, refl, cyc, dt_block, mu_ratio)
        counts[route] += 1
        comp = {k: float(np.nansum(h[k][sl, j])) for k in ("W_par", "W_phi", "W_r", "W_theta", "W_all")}
        comp["closes"] = bool(G > 0 and abs(comp["W_all"] - G) <= 0.3 * G)
        per.append({"electron": int(c), "gamma_smoothed_max": float(gs[km]), "level_drift_frame": float(level[c]),
                    "route": route, "G": G, "work": comp,
                    "n_reflections": len(refl), "cycles": cyc[0],
                    "gain_by_site": {str(s_): float(sum(q[4] for q in steps if q[0] == s_)) for s_ in range(-1, 5)},
                    "polar_at_max": float(pol[km])})
    comp_by_route = {}
    for r in ROUTES:
        E = [e for e in per if e["route"] == r and e["work"]["closes"]]
        if E:
            tot = sum(e["work"]["W_all"] for e in E)
            comp_by_route[r] = {k: sum(e["work"][k] for e in E) / tot if tot else None
                                for k in ("W_par", "W_phi", "W_r", "W_theta")}
    out = {"run": str(run_dir), "n_top": len(per), "level_cut": float(cut), "counts": counts,
           "shares": {r: counts[r] / max(1, len(per)) for r in ROUTES},
           "work_component_shares_by_route (electrons whose work closes within 30 %)": comp_by_route,
           "share_work_closes": float(np.mean([e["work"]["closes"] for e in per])) if per else None,
           "electrons": per, "rules": __doc__}
    (an / "routes.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in ("electrons", "rules")}, indent=1))
    return 0


def selfcheck(check):
    st = [(3, 5, 0, 10, 0.1, 0.02), (1, 8, 10, 30, 0.6, 0.1)]
    g_ = np.array([1.0, 1.1, 1.3, 1.2, 1.6])
    sg = steps_of(np.array([4, 3, 3, 2, 2]), np.zeros(5, int), g_, np.zeros(5))
    ok0 = abs(sum(q[4] for q in sg) - (g_[-1] - g_[0])) < 1e-12
    ok = ok0 and classify(0.7, st, [], (0, 0.0), 2.0, 1.0) == "R2"
    ok = ok and classify(0.5, [], [(20.0, 0.4, 0.6)], (0, 0.0), 2.0, 1.0) == "R5"
    ok = ok and classify(0.5, [], [(20.0, 0.1, 0.3), (30.0, 0.2, 0.4)], (0, 0.0), 2.0, 1.0) == "R1"
    ok = ok and classify(0.5, [], [(20.0, 0.15, 0.3), (22.0, 0.15, 0.3)], (0, 0.0), 2.0, 1.0) == "R3a"
    ok = ok and cycles_of(np.array([4, 4, 3, 2, 2, 3, 4, 4, 3, 1, 4]), np.linspace(1, 2, 11))[0] == 2
    check("routes, signatures of constructed electrons", ok, "R2, R5, R1, R3a and two upstream-downstream cycles")
