"""Campaign A (campaign-a.md): the analyses of one run that need the tracks or the particle dumps, in one pass
over the tracks, and the pruning of those large outputs once the results are written.

analyze.py campaign <run>   writes analysis/campaign.json, campaign_entries.npy, campaign_tracks.npz
analyze.py prune <run>      deletes tracks and hot_e once campaign.json and counts.json exist; analysis/pruned.json

Sites as in sites.md. An entry is a step from upstream (site 4) into the transition layer (site 3); it ends when
the electron leaves the layer: back upstream ("return", a reflection) or into the steep wall.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import sites

VHT_EDGES = np.array([0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 3.0])
X_EDGES = np.array([0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, np.inf])
POLAR_EDGES = np.arange(0.0, 90.1, 10.0)
SMOOTH = 12.5   # time units over which E and B are averaged for the drift frame (about one gyration)
AFTER = 8       # samples after leaving the layer at which gamma_out is read
HIST_EVERY = 4  # the saved histories keep every 4th sample, with the work summed over each block of 4
HIST_GAMMA = 1.25  # histories of every electron with gamma_max >= this, and of a random 2 %


def entries_of(site_row):
    """(k_in, k_out, outcome) for every entry of one electron: k_in the first sample in the layer, k_out the first
    sample after it, outcome 1 when that sample is upstream, 0 when it is in the steep wall. Entries that leave
    the layer into anything else or do not leave it before the end are dropped."""
    s = np.asarray(site_row)
    out = []
    for k in np.flatnonzero((s[:-1] == 4) & (s[1:] == 3)) + 1:
        e = k
        while e + 1 < len(s) and s[e + 1] == 3:
            e += 1
        if e + 1 >= len(s) or s[e + 1] not in (2, 4):
            continue
        out.append((int(k), int(e + 1), 1 if s[e + 1] == 4 else 0))
    return out


def adiabatic_reflection(u_par, u_perp, gamma, v_ht, b_c):
    """Adiabatic mirror in the de Hoffmann-Teller frame, no potential (prediction H_A of campaign-a.md).
    u_par: momentum along the field, oriented so that the frame moves at +v_ht; u_perp >= 0; b_c the compression
    B_max / B_1. Returns (reflected, gamma after reflection in the lab frame). No frame for v_ht >= 1."""
    u_par, u_perp, gamma, v_ht, b_c = np.broadcast_arrays(*(np.asarray(q, float) for q in (u_par, u_perp, gamma,
                                                                                       v_ht, b_c)))
    ok = v_ht < 1.0
    vh = np.where(ok, v_ht, 0.0)
    g_h = 1.0 / np.sqrt(1.0 - vh ** 2)
    up = g_h * (u_par - vh * gamma)          # momentum along B in the frame
    gp = g_h * (gamma - vh * u_par)          # energy in the frame
    sin2 = u_perp ** 2 / np.maximum(u_perp ** 2 + up ** 2, 1e-30)
    refl = ok & (up < 0.0) & (sin2 > 1.0 / np.maximum(b_c, 1.0 + 1e-9))
    g_out = g_h * (gp - vh * up)             # u_par -> -u_par in the frame, back to the lab
    return refl, np.where(refl, g_out, gamma)


def _smooth(a, w):
    k = np.ones(w) / w
    ok = np.isfinite(a)
    num = np.apply_along_axis(lambda q: np.convolve(q, k, mode="same"), 0, np.where(ok, a, 0.0))
    den = np.apply_along_axis(lambda q: np.convolve(q, k, mode="same"), 0, ok.astype(float))
    return np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)


def drift_frame_level(g, u, E, B, w):
    """Gyration-averaged energy in the local E x B frame (fields smoothed over w samples): gamma' level."""
    Es = np.stack([_smooth(q, w) for q in E])
    Bs = np.stack([_smooth(q, w) for q in B])
    vE = np.cross(Es, Bs, axis=0) / np.maximum(np.sum(Bs ** 2, axis=0), 1e-30)
    sp = np.linalg.norm(vE, axis=0)
    ok = np.isfinite(sp) & (sp < 0.95) & np.isfinite(g)
    v = u / g
    gp = np.where(ok, g * (1.0 - np.sum(vE * v, axis=0)) / np.sqrt(1.0 - np.where(ok, sp, 0.0) ** 2), np.nan)
    return _smooth(gp, w)


def _binned_mean(x, y, edges):
    out, n = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        n.append(int(m.sum()))
        out.append(float(np.mean(y[m])) if m.any() else None)
    return out, n


def cmd_campaign(args) -> int:
    from census import _open
    run_dir = Path(args.run)
    run = _open(run_dir)
    pl = run.meta["planner"]
    c, d_i = np.asarray(pl["center"], float), float(pl["d_i"])
    fz = np.load(run_dir / "analysis" / "front.npz")
    st = np.load(run_dir / "analysis" / "sites.npz")
    theta, ts = st["theta"], st["t"]
    tr = run.tracks("tracks_e", "e")
    t = np.asarray(tr.t, float)
    dt = t[1] - t[0]
    w = int(round(SMOOTH / dt)) | 1
    n = len(tr.ids)
    i_out = np.array([int(np.argmin(np.abs(ts - q))) for q in t])
    rows = i_out[:, None]
    i_front = np.array([int(np.argmin(np.abs(fz["t"] - q))) for q in t])
    rd0 = fz["r_debris"][0]
    rng = np.random.default_rng(int(args.seed))
    gmax = np.full(n, np.nan)
    kmax = np.zeros(n, np.int32)
    birth_r = np.full(n, np.nan)
    birth_pol = np.full(n, np.nan)
    level = np.full(n, np.nan)
    ent = []
    net = {k: np.zeros(len(sites.NAMES)) for k in ("random", "inside", "outside")}
    tim = {k: np.zeros(len(sites.NAMES)) for k in net}
    hist = {k: [] for k in ("col", "x", "y", "z", "ux", "uy", "uz", "Bx", "By", "Bz", "site", "W_par", "W_phi",
                            "W_r", "W_theta", "W_all")}
    for a in range(0, n, args.chunk):
        cols = np.arange(a, min(n, a + args.chunk))
        D = {k: np.asarray(tr.f[k][:, cols], float) for k in ("x", "y", "z", "ux", "uy", "uz", "Bx", "By", "Bz")}
        x, y, z = (D[k] - c[i] for i, k in enumerate("xyz"))
        u = np.stack([D["ux"], D["uy"], D["uz"]])
        g = np.sqrt(1.0 + np.sum(u ** 2, axis=0))
        B = np.stack([D["Bx"], D["By"], D["Bz"]])
        Bm = np.linalg.norm(B, axis=0)
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        pol = np.arctan2(np.hypot(x, y), np.abs(z))
        ith = np.clip(np.rint((np.nan_to_num(pol) - theta[0]) / (theta[1] - theta[0])).astype(int), 0, len(theta) - 1)
        s = sites.site_of(r, st["r0"][rows, ith], st["r1"][rows, ith], st["r_bmax"][rows, ith], st["r_f"][rows, ith],
                          d_i)
        gg = np.where(np.isfinite(g), g, 0.0)
        kmax[cols] = np.argmax(gg, axis=0)
        gmax[cols] = gg[kmax[cols], np.arange(len(cols))]
        birth_r[cols], birth_pol[cols] = r[0], np.degrees(pol[0])
        inside = r[0] < np.interp(pol[0], theta, rd0)
        # entries into the layer
        for j in range(len(cols)):
            for k, ko, outc in entries_of(s[:, j]):
                q = k - 1
                nrm = sites.local_normal(x[q:k, j], y[q:k, j], z[q:k, j], theta, st["r_f"][i_out[k]])[:, 0]
                bh = B[:, q, j] / Bm[q, j]
                cos_bn = float(np.dot(bh, nrm))
                if cos_bn < 0:
                    bh, cos_bn = -bh, -cos_bn
                tb = float(np.degrees(np.arccos(min(1.0, cos_bn))))
                vn = float(np.interp(pol[k, j], theta, np.nan_to_num(fz["v_n"][i_front[k]], nan=0.0)))
                vht = vn / max(cos_bn, 1e-3)
                gi = g[q, j]
                upar = float(np.dot(u[:, q, j], bh))
                uperp = float(np.sqrt(max(np.sum(u[:, q, j] ** 2) - upar ** 2, 0.0)))
                bc = float(np.interp(pol[k, j], theta, np.nan_to_num(fz["b_mean"][i_front[k]], nan=1.0)))
                ka = min(len(t) - 1, ko + AFTER) if outc == 1 else ko
                t_up, ahead = 0.0, np.nan
                if outc == 1:
                    e2 = ko
                    while e2 + 1 < len(t) and s[e2 + 1, j] == 4:
                        e2 += 1
                    t_up = (e2 - ko + 1) * dt
                    ahead = float(np.nanmax(r[ko:e2 + 1, j] - st["r_f"][i_out[ko:e2 + 1], ith[ko:e2 + 1, j]]))
                refl_a, g_a = adiabatic_reflection(upar, uperp, gi, vht, bc)
                ent.append((cols[j], t[k], np.degrees(pol[k, j]), tb, vht, np.sqrt(1.0 - 1.0 / gi ** 2),
                            abs(upar) / gi, outc, gi, g[ka, j], float(inside[j]), upar, uperp, bc, t_up, ahead, vn,
                            float(refl_a), float(g_a)))
        # net change of gamma per site: a random 10 % of the chunk, and inside / outside born electrons with
        # gamma_max >= 1.5 (whole run, no window)
        dg = np.full(g.shape, np.nan)
        dg[:-1] = g[1:] - g[:-1]
        pick = {"random": rng.random(len(cols)) < 0.1, "inside": inside & (gmax[cols] >= 1.5),
                "outside": ~inside & (gmax[cols] >= 1.5)}
        for key, sel in pick.items():
            for k_ in range(len(sites.NAMES)):
                m = (s[:, sel] == k_) & np.isfinite(dg[:, sel])
                net[key][k_] += float(dg[:, sel][m].sum())
                tim[key][k_] += float(m.sum())
        # drift-frame level and histories: electrons with gamma_max >= HIST_GAMMA and a random 2 % of the chunk
        sel = np.flatnonzero((gmax[cols] >= HIST_GAMMA) | (rng.random(len(cols)) < 0.02))
        if len(sel):
            E = np.stack([np.asarray(tr.f[k][:, cols[sel]], float) for k in ("Ex", "Ey", "Ez")])
            lv = drift_frame_level(g[:, sel], u[:, :, sel], E, B[:, :, sel], w)
            level[cols[sel]] = np.nanmax(np.where(np.isfinite(lv), lv, -1.0), axis=0)
            xs, ys, zs = x[:, sel], y[:, sel], z[:, sel]
            v = u[:, :, sel] / g[:, sel]
            Bs = B[:, :, sel]
            bh = Bs / np.maximum(np.linalg.norm(Bs, axis=0), 1e-30)
            rr = np.sqrt(xs ** 2 + ys ** 2 + zs ** 2)
            rho = np.hypot(xs, ys)
            ephi = np.stack([-ys, xs, np.zeros_like(xs)]) / np.maximum(rho, 1e-9)
            er = np.stack([xs, ys, zs]) / np.maximum(rr, 1e-9)
            eth = np.cross(ephi, er, axis=0)
            pw = {"W_par": -np.sum(E * bh, axis=0) * np.sum(v * bh, axis=0),
                  "W_phi": -np.sum(E * ephi, axis=0) * np.sum(v * ephi, axis=0),
                  "W_r": -np.sum(E * er, axis=0) * np.sum(v * er, axis=0),
                  "W_theta": -np.sum(E * eth, axis=0) * np.sum(v * eth, axis=0),
                  "W_all": -np.sum(E * v, axis=0)}
            nb = len(t) // HIST_EVERY
            for k_, val in pw.items():
                blk = np.nan_to_num(val[: nb * HIST_EVERY]).reshape(nb, HIST_EVERY, -1).sum(axis=1) * dt
                hist[k_].append(blk.astype(np.float32))
            idx = slice(0, nb * HIST_EVERY, HIST_EVERY)
            for k_, val in (("x", xs), ("y", ys), ("z", zs), ("ux", u[0][:, sel]), ("uy", u[1][:, sel]),
                            ("uz", u[2][:, sel]), ("Bx", Bs[0]), ("By", Bs[1]), ("Bz", Bs[2])):
                hist[k_].append(val[idx].astype(np.float32))
            hist["site"].append(s[idx][:, sel].astype(np.int8))
            hist["col"].append(cols[sel])
        print(f"campaign: {min(n, a + args.chunk)}/{n} tracked electrons; entries {len(ent)}", flush=True)
    ent = np.array(ent, float).reshape(-1, 19)
    dest = run_dir / "analysis"
    np.savez_compressed(dest / "campaign_hist.npz", t=t[0: (len(t) // HIST_EVERY) * HIST_EVERY: HIST_EVERY],
                        every=HIST_EVERY, gamma_cut=HIST_GAMMA,
                        **{k: (np.concatenate(v) if k == "col" else np.concatenate(v, axis=1))
                           for k, v in hist.items() if v})
    np.save(dest / "campaign_entries.npy", ent)
    np.savez_compressed(dest / "campaign_tracks.npz", ids=np.asarray(tr.ids), gmax=gmax, t_gmax=t[kmax],
                        birth_r=birth_r, birth_polar=birth_pol, level_drift_frame=level,
                        net=np.stack([net[k] for k in ("random", "inside", "outside")]),
                        tim=np.stack([tim[k] for k in ("random", "inside", "outside")]))
    return write_summary(run_dir)


def write_summary(run_dir) -> int:
    """campaign.json from campaign_entries.npy and campaign_tracks.npz (no access to the tracks)."""
    from census import _open
    run_dir = Path(run_dir)
    dest = run_dir / "analysis"
    pl = _open(run_dir).meta["planner"]
    d_i = float(pl["d_i"])
    fz = np.load(dest / "front.npz")
    theta, rd0 = fz["theta"], fz["r_debris"][0]
    ent = np.load(dest / "campaign_entries.npy")
    tk = np.load(dest / "campaign_tracks.npz")
    gmax, birth_r, birth_pol, level = tk["gmax"], tk["birth_r"], tk["birth_polar"], tk["level_drift_frame"]
    n = len(gmax)
    net = dict(zip(("random", "inside", "outside"), tk["net"])) if "net" in tk.files else {}
    tim = dict(zip(("random", "inside", "outside"), tk["tim"])) if "tim" in tk.files else {}
    out = {"run": str(run_dir), "n_tracked": int(n),
           "entry_columns": ["electron", "t", "polar_deg", "local_tbn_deg", "V_HT", "v", "v_par", "return",
                             "gamma_in", "gamma_out", "born_inside", "u_par (frame direction)", "u_perp",
                             "compression", "t_upstream_after_return", "max_ahead_of_front", "V_n",
                             "reflect_H_A", "gamma_out_H_A"],
           "smooth": SMOOTH, "after_samples": AFTER}
    if len(ent):
        ins = ent[:, 10]
        outside = ins == 0
        out["entries"] = int(len(ent))
        out["entries_born_outside"] = int(outside.sum())
        e_o = ent[outside]
        vht_o, v_o, ret_o, gin_o, gout_o, pol_o = e_o[:, 4], e_o[:, 5], e_o[:, 7], e_o[:, 8], e_o[:, 9], e_o[:, 2]
        out["V_HT_edges"] = VHT_EDGES.tolist()
        out["P_return_by_VHT"], out["n_by_VHT"] = _binned_mean(vht_o, ret_o, VHT_EDGES)
        out["X_edges"] = [float(q) for q in X_EDGES[:-1]] + [None]
        out["P_return_by_VHT_over_v"], out["n_by_VHT_over_v"] = _binned_mean(vht_o / v_o, ret_o, X_EDGES)
        out["polar_edges"] = POLAR_EDGES.tolist()
        out["P_return_by_polar"], out["n_by_polar"] = _binned_mean(pol_o, ret_o, POLAR_EDGES)
        b0 = float(pl["b0"])
        ref2 = (ret_o == 1) & (gout_o > 1.02 * gin_o)
        persist = (e_o[:, 14] >= 2 * np.pi * gin_o / b0) | (e_o[:, 15] >= d_i + 2 * np.sqrt(gin_o ** 2 - 1) / b0)
        pa = (e_o[:, 17] == 1) & (e_o[:, 18] > 1.02 * gin_o)
        out["P_ref_by_VHT"], _ = _binned_mean(vht_o, ref2.astype(float), VHT_EDGES)
        out["P_ref_persistent_by_VHT"], _ = _binned_mean(vht_o, (ref2 & persist).astype(float), VHT_EDGES)
        out["P_ref_H_A_by_VHT"], _ = _binned_mean(vht_o, pa.astype(float), VHT_EDGES)
        out["P_ref_by_VHT_over_v"], _ = _binned_mean(vht_o / v_o, ref2.astype(float), X_EDGES)
        out["P_ref_persistent_by_VHT_over_v"], _ = _binned_mean(vht_o / v_o, (ref2 & persist).astype(float),
                                                                X_EDGES)
        out["P_ref_by_polar"], _ = _binned_mean(pol_o, ref2.astype(float), POLAR_EDGES)
        out["P_ref_persistent_by_polar"], _ = _binned_mean(pol_o, (ref2 & persist).astype(float), POLAR_EDGES)
        out["P_ref_H_A_by_polar"], _ = _binned_mean(pol_o, pa.astype(float), POLAR_EDGES)
        iv = np.clip(np.digitize(vht_o, VHT_EDGES) - 1, -1, len(VHT_EDGES) - 2)
        ip = np.clip(np.digitize(pol_o, POLAR_EDGES) - 1, 0, len(POLAR_EDGES) - 2)
        tab = np.zeros((len(VHT_EDGES) - 1, len(POLAR_EDGES) - 1), int)
        np.add.at(tab, (iv[iv >= 0], ip[iv >= 0]), 1)
        out["n_by_VHT_and_polar"] = tab.tolist()
        rr = ref2 & persist
        if rr.any():
            out["persistent returns: gamma_out/gamma_in p50 p90 p99"] = \
                np.percentile(gout_o[rr] / gin_o[rr], (50, 90, 99)).tolist()
            okb = rr & (vht_o < 0.99)
            if okb.any():
                bnd = (1 + vht_o[okb] ** 2) / (1 - vht_o[okb] ** 2)
                out["persistent returns: gamma_out over the bound from rest p10 p50 p90"] = \
                    np.percentile(gout_o[okb] / bnd, (10, 50, 90)).tolist()
        hi = (ret_o == 1) & (gout_o >= 1.25)
        out["share_return_gamma_ge_1p25_by_polar"], _ = _binned_mean(pol_o, hi.astype(float), POLAR_EDGES)
        g_hi_ratio = (ret_o == 1) & (gout_o >= 1.25 * gin_o)
        out["share_return_gain_25pc_by_polar"], _ = _binned_mean(pol_o, g_hi_ratio.astype(float), POLAR_EDGES)
        ok = (ret_o == 1) & (vht_o < 0.99)
        if ok.any():
            bound = (1 + vht_o[ok] ** 2) / (1 - vht_o[ok] ** 2)
            out["return_gamma_out_over_bound_p10_50_90"] = np.percentile(gout_o[ok] / bound, (10, 50, 90)).tolist()
            out["return_gamma_out_p50_p90_p99_max"] = np.percentile(gout_o[ret_o == 1], (50, 90, 99)).tolist() + \
                [float(gout_o[ret_o == 1].max())]
    ins_all = birth_r < np.interp(np.radians(birth_pol), theta, rd0)
    for lab, m in (("born inside", ins_all), ("born outside", ~ins_all)):
        out[f"gamma_max {lab}: n, share >= 1.5, >= 2, >= 3"] = [int(m.sum())] + \
            [float(np.mean(gmax[m] >= q)) for q in (1.5, 2.0, 3.0)]
        lv = level[m & np.isfinite(level)]
        lvh = level[m & np.isfinite(level) & (gmax >= 1.5)]
        out[f"drift-frame level {lab}: p50 p90 p99 of gamma_max >= 1.5; n"] = \
            (np.percentile(lvh, (50, 90, 99)).tolist() if len(lvh) else []) + [int(len(lvh))]
        out[f"drift-frame level {lab}, random 2 %: p50 p90 p99"] = \
            np.percentile(lv, (50, 90, 99)).tolist() if len(lv) else []
    out["net_dgamma_by_site"] = {k: dict(zip(sites.NAMES, v.tolist())) for k, v in net.items()}
    out["time_share_by_site"] = {k: dict(zip(sites.NAMES, (v / max(1.0, v.sum())).tolist())) for k, v in tim.items()}
    (dest / "campaign.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if "edges" not in k}, indent=1))
    print(f"wrote {dest / 'campaign.json'}, campaign_entries.npy, campaign_tracks.npz")
    return 0


def cmd_prune(args) -> int:
    """Delete the tracks and the particle dumps of a run once the analyses that need them have written their
    results. The run is reproducible from its record (commit, config, seed, GPUs)."""
    run_dir = Path(args.run)
    an = run_dir / "analysis"
    need = ["null.json"] if args.null else ["campaign.json", "campaign_entries.npy", "campaign_tracks.npz",
                                             "campaign_hist.npz", "counts.json"]
    missing = [q for q in need if not (an / q).exists()]
    if missing:
        print(f"not pruned: missing {missing}")
        return 1
    gone = []
    for sub in ("tracks", "particles"):
        p = run_dir / sub
        if p.exists():
            size = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
            shutil.rmtree(p)
            gone.append({"path": str(p), "bytes": int(size)})
    (an / "pruned.json").write_text(json.dumps({"deleted": gone, "reason": "analysed; reproducible from the run "
                                                "record (commit, config, seed)"}, indent=1))
    print(f"pruned {sum(q['bytes'] for q in gone) / 1e9:.1f} GB: {[q['path'] for q in gone]}")
    return 0


def selfcheck(check):
    vh = 0.4
    gh = 1.0 / np.sqrt(1.0 - vh ** 2)
    refl, g_out = adiabatic_reflection([0.0, 0.0, 0.0], [1e-3, 0.5, 0.05], [1.0, np.sqrt(1.25), np.sqrt(1.0025)],
                                       [vh, vh, 1.2], [3.0, 3.0, 3.0])
    want = gh ** 2 * (1 + vh ** 2)
    ok = (not refl[0]) and refl[1] and (not refl[2]) and np.isclose(g_out[0], 1.0)
    # an electron at rest is in the loss cone (sin^2 = 0); one with u_perp = 0.5 at rest along B reflects to
    # gamma_HT^2 (gamma + v^2 gamma) since u_par' = -gamma_HT v gamma
    ok = ok and np.isclose(g_out[1], want * np.sqrt(1.25))
    check("campaign, adiabatic reflection in the de Hoffmann-Teller frame", ok,
          f"reflected {refl.tolist()} (False, True, False: rest in the loss cone, u_perp 0.5 out of it, no frame "
          f"above c); gamma_out {g_out[1]:.6f} against gamma_HT^2 (1 + V^2) gamma = {want * np.sqrt(1.25):.6f}")
    s = [4, 4, 3, 3, 4, 4, 3, 2, 2, 4, 3, 3, -1, 4, 3]
    got = entries_of(s)
    check("campaign, entries into the layer", got == [(2, 4, 1), (6, 7, 0)],
          f"{got} for sites {s} (expected a return, then an entry into the wall; the others dropped)")
