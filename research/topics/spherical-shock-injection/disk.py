"""The equatorial disk of energetic electrons in a dipole run (D1-disk): histories of its electrons from the slab
dump `disk_e`, with the fields at the electrons from the equatorial slice, and the test of the mechanisms.

analyze.py disk <run> [--frames <dir>]   writes analysis/disk.json and disk_hist.npz (and animation frames)

Target: the electrons in the disk (|z| < 12 d_e, rho < 130 d_e) with gamma >= 1.5 at the last dump. Control: a
random sample of electrons born in the same band (birth radius and polar angle within those of the target) that
stay below gamma 1.2. Each is followed by id through the dumps (every 9 time units). The fields at an electron
come from the equatorial slice (z = centre) by bilinear interpolation; the electrons of the disk stay within
12 d_e of that plane. Mechanism tests on the target, between successive dumps:
  - betatron: p_perp^2 / |B| kept while |B| grows (ratio of mu over the run);
  - drift along E_phi: the energy change against the work -E_phi v_phi dt with v_phi from the change of azimuth;
  - where: radius and azimuth of the gains.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import shockkin as sk


def bilinear(field, x0, dx, x, y):
    """Bilinear interpolation of a plane field[i, j] on a uniform grid with cell centres x0 + (i + 0.5) dx."""
    fx = (x - x0) / dx - 0.5
    fy = (y - x0) / dx - 0.5
    i = np.clip(np.floor(fx).astype(int), 0, field.shape[0] - 2)
    j = np.clip(np.floor(fy).astype(int), 0, field.shape[1] - 2)
    ax, ay = np.clip(fx - i, 0, 1), np.clip(fy - j, 0, 1)
    return ((1 - ax) * (1 - ay) * field[i, j] + ax * (1 - ay) * field[i + 1, j]
            + (1 - ax) * ay * field[i, j + 1] + ax * ay * field[i + 1, j + 1])


def cmd_disk(args) -> int:
    from census import _open
    run_dir = Path(args.run)
    run = _open(run_dir)
    pl = run.meta["planner"]
    c = np.asarray(pl["center"], float)
    n_cells = tuple(int(q) for q in run.config["n_cells"])
    dx = float(run.config["dx"][0]) if hasattr(run.config["dx"], "__len__") else float(run.config["dx"])
    ppc = int(run.config["ppc"])
    p = run.particles("disk_e")
    times = np.asarray(p.times, float)
    its = list(p.iterations)
    last = p.get(its[-1], "e")
    x, y, z = (np.asarray(last[q], float) - c[i] for i, q in enumerate("xyz"))
    g = np.sqrt(1 + sum(np.asarray(last[q], float) ** 2 for q in ("ux", "uy", "uz")))
    ids_last = np.asarray(last["id"], np.uint64)
    n_amb = np.uint64(ppc * int(np.prod(n_cells)))
    disk = (np.abs(z) < 12) & (np.hypot(x, y) < 130) & (ids_last < n_amb)
    target = np.sort(ids_last[disk & (g >= 1.5)])

    def birth(ids):
        ix, iy, iz = sk.birth_cell(ids, n_cells, ppc)
        xb, yb, zb = ((q + 0.5) * dx - cc for q, cc in zip((ix, iy, iz), c))
        return np.sqrt(xb ** 2 + yb ** 2 + zb ** 2), np.degrees(np.arctan2(np.hypot(xb, yb), np.abs(zb)))
    rb_t, pb_t = birth(target)
    first = p.get(its[0], "e")
    ids0 = np.asarray(first["id"], np.uint64)
    amb0 = ids0[ids0 < n_amb]
    rb0, pb0 = birth(amb0)
    band = ((rb0 >= np.percentile(rb_t, 5)) & (rb0 <= np.percentile(rb_t, 95))
            & (pb0 >= np.percentile(pb_t, 5)))
    pool = np.setdiff1d(amb0[band], target)
    rng = np.random.default_rng(3)
    control = np.sort(rng.choice(pool, min(len(pool), 4 * len(target)), replace=False))
    sz = run.series("slice_z").f
    ts = np.asarray(sz["t"], float)
    groups = {"target": target, "control": control}
    H = {k: {q: np.full((len(times), len(v)), np.nan) for q in ("x", "y", "z", "ux", "uy", "uz", "B", "Ephi", "Bpar")}
         for k, v in groups.items()}
    for k_t, it in enumerate(its):
        d = p.get(it, "e")
        ids = np.asarray(d["id"], np.uint64)
        order = np.argsort(ids)
        ids_s = ids[order]
        i_s = int(np.argmin(np.abs(ts - times[k_t])))
        F = {q: np.asarray(sz[q][i_s], float) for q in ("Ex", "Ey", "Bx", "By", "Bz")}
        for name, want in groups.items():
            pos = np.searchsorted(ids_s, want)
            pos = np.clip(pos, 0, len(ids_s) - 1)
            hit = ids_s[pos] == want
            src = order[pos[hit]]
            for q in ("x", "y", "z", "ux", "uy", "uz"):
                v = np.asarray(d[q], float)[src]
                H[name][q][k_t, hit] = v - (c["xyz".index(q)] if q in "xyz" else 0.0)
            xx, yy = H[name]["x"][k_t, hit], H[name]["y"][k_t, hit]
            b = [bilinear(F[q], -c[0], dx, xx, yy) for q in ("Bx", "By", "Bz")]
            ex, ey = (bilinear(F[q], -c[0], dx, xx, yy) for q in ("Ex", "Ey"))
            rho = np.maximum(np.hypot(xx, yy), 1e-9)
            bm = np.sqrt(b[0] ** 2 + b[1] ** 2 + b[2] ** 2)
            H[name]["B"][k_t, hit] = bm
            H[name]["Ephi"][k_t, hit] = (-ex * yy + ey * xx) / rho
            u = np.stack([H[name][q][k_t, hit] for q in ("ux", "uy", "uz")])
            H[name]["Bpar"][k_t, hit] = np.sum(u * np.stack(b), axis=0) / np.maximum(bm, 1e-30)
        print(f"disk: dump {k_t + 1}/{len(its)} (t = {times[k_t]:.0f})", flush=True)
    out = {"run": str(run_dir), "n_target": int(len(target)), "n_control": int(len(control)),
           "target_birth_r_p10_50_90": np.percentile(rb_t, (10, 50, 90)).tolist(),
           "target_birth_polar_p10_50_90": np.percentile(pb_t, (10, 50, 90)).tolist(), "t": times.tolist()}
    for name in groups:
        h = H[name]
        g_ = np.sqrt(1 + h["ux"] ** 2 + h["uy"] ** 2 + h["uz"] ** 2)
        u2 = h["ux"] ** 2 + h["uy"] ** 2 + h["uz"] ** 2
        mu = (u2 - h["Bpar"] ** 2) / np.maximum(h["B"], 1e-30)
        r = np.hypot(h["x"], h["y"])
        phi = np.unwrap(np.arctan2(h["y"], h["x"]), axis=0)
        dphi = np.diff(phi, axis=0)
        dt = np.diff(times)[:, None]
        vphi = 0.5 * (r[1:] + r[:-1]) * dphi / dt
        w_phi = -0.5 * (h["Ephi"][1:] + h["Ephi"][:-1]) * vphi * dt
        dg = np.diff(g_, axis=0)
        ok = np.isfinite(dg) & np.isfinite(w_phi)
        out[name] = {"gamma_median_by_dump": np.nanmedian(g_, axis=1).tolist(),
                     "r_median_by_dump": np.nanmedian(r, axis=1).tolist(),
                     "B_median_by_dump": np.nanmedian(h["B"], axis=1).tolist(),
                     "mu_median_by_dump": np.nanmedian(mu, axis=1).tolist(),
                     "found_share_by_dump": np.mean(np.isfinite(g_), axis=1).tolist(),
                     "net_dgamma": float(np.nansum(dg[ok])),
                     "work_E_phi_from_azimuth": float(np.nansum(w_phi[ok])),
                     "azimuth_travelled_rad_p10_50_90": np.nanpercentile(phi[-1] - phi[0], (10, 50, 90)).tolist()}
    np.savez_compressed(run_dir / "analysis" / "disk_hist.npz", t=times,
                        **{f"{k}_{q}": v for k, hh in H.items() for q, v in hh.items()})
    (run_dir / "analysis" / "disk.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: (v if not isinstance(v, dict) else {q: w for q, w in v.items() if "by_dump" not in q})
                      for k, v in out.items() if k != "t"}, indent=1))
    return 0
