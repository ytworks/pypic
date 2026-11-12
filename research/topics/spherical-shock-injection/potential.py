"""The rise of the potential along the field lines through the layer: phi of exploratory-protocol.md Sect. 3.1,
measured from the meridional slices.

In the frame of the foot point the electric field along B is the same as in the box to the order V_HT^2 (E.B is
an invariant), so phi = - integral of E_par dl along a field line at one time, from the upstream into the layer;
phi > 0 means that an electron gains energy on the way in. The planes hold all three components of E and B.
With the symmetry about the axis of B0 the field line is followed in the plane with the two in-plane components,
and the out-of-plane part of the path is in the weight |B| / B_pol: phi = - sum of (E.B) / B_pol ds_pol.

A path starts START d_i ahead of the fitted density front, where phi = 0, and runs towards the mid-plane to DEPTH
d_i behind the front. The values are taken where a path first passes each level of the distance to the front,
and averaged over the paths whose foot (the passage of the front) lies in a bin of 5 degrees. The error is the
scatter of the eight quadrants (two planes), and the same paths with the E of the null run give the noise.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import census
import numpy as np
from scipy.ndimage import map_coordinates

START, DEPTH = 3.0, 6.0  # d_i ahead of and behind the density front
LEVEL_STEP = 0.5  # d_e between the levels of the distance to the front
DS = 0.25  # step of a path in the plane
DTHETA0 = 0.25  # degrees between the starting points of a quadrant
N_BIN = 18  # bins of the foot angle, 5 degrees
CLIP = 0.1  # of B0: the smallest in-plane field used in the weight |B| / B_pol
PLANES = (("slice_y", "x", 0), ("slice_x", "y", 1))
WINDOW = 5  # outputs in the running mean that the count model uses
REL_ERR_MAX = 0.5  # protocol Sect. 3.4: phi is used if its error is below 50 %


def _sample(field, a, z):
    arr, (ax_a, ax_z) = field
    return map_coordinates(arr, [(a - ax_a[0]) / (ax_a[1] - ax_a[0]), (z - ax_z[0]) / (ax_z[1] - ax_z[0])],
                           order=1, mode="nearest")


def levels(d_i: float):
    return np.arange(START * d_i, -DEPTH * d_i - 1e-9, -LEVEL_STEP)


def trace_plane(fields, centre, radius, d_i: float, b0: float, across: str, e_null=None, r_outer: float = np.inf,
                n_steps: int = 2400) -> dict:
    """Paths along the field lines of one plane.

    fields   key (Ex ... Bz) -> (array [across, z], (axis across, axis z))
    centre   (across, z) of the centre; radius(theta) -> radius of the front at the polar angle from B0
    across   name of the in-plane component across B0 ("x" or "y")
    e_null   optional fields Ex, Ey, Ez of the null run: integrated along the same paths
    Returns per path: theta_foot, quadrant, and at each level phi, phi_null, b (|B| / B0); clipped: the share of
    the steps where the in-plane field was below CLIP B0."""
    lev = levels(d_i)
    k_zero = int(round(START * d_i / LEVEL_STEP))
    th0 = np.radians(np.arange(0.5 * DTHETA0, 90.0, DTHETA0))
    n0 = len(th0)
    sa, sz = np.repeat([1.0, -1.0, 1.0, -1.0], n0), np.repeat([1.0, 1.0, -1.0, -1.0], n0)
    th = np.tile(th0, 4)
    r0 = radius(th) + START * d_i
    a, z = centre[0] + sa * r0 * np.sin(th), centre[1] + sz * r0 * np.cos(th)
    n = len(a)
    step_dir = -sz * np.sign(b0)  # along B towards the mid-plane
    phi, phi_n = np.zeros(n), np.zeros(n)
    out = {k: np.full((n, len(lev)), np.nan) for k in ("phi", "phi_null", "b")}
    theta_foot = np.full(n, np.nan)
    k = np.zeros(n, int)
    delta_prev = np.full(n, START * d_i)
    alive = r0 <= r_outer
    n_clip = n_all = 0
    for _ in range(n_steps):
        idx = np.flatnonzero(alive)
        if not idx.size:
            break
        ai, zi, d = a[idx], z[idx], step_dir[idx]
        ba, bz = _sample(fields["B" + across], ai, zi), _sample(fields["Bz"], ai, zi)
        bp = np.maximum(np.hypot(ba, bz), 1e-30)
        am, zm = ai + 0.5 * DS * d * ba / bp, zi + 0.5 * DS * d * bz / bp
        bv = {c: _sample(fields["B" + c], am, zm) for c in "xyz"}
        bp = np.hypot(bv[across], bv["z"])
        weight = DS / np.maximum(bp, CLIP * abs(b0))
        n_clip, n_all = n_clip + int((bp < CLIP * abs(b0)).sum()), n_all + idx.size
        bp = np.maximum(bp, 1e-30)
        an, zn = ai + DS * d * bv[across] / bp, zi + DS * d * bv["z"] / bp
        before, before_n = phi[idx], phi_n[idx]
        phi[idx] = before - d * sum(_sample(fields["E" + c], am, zm) * bv[c] for c in "xyz") * weight
        if e_null is not None:
            phi_n[idx] = before_n - d * sum(_sample(e_null["E" + c], am, zm) * bv[c] for c in "xyz") * weight
        a[idx], z[idx] = an, zn
        da, dz = an - centre[0], zn - centre[1]
        r, theta = np.hypot(da, dz), np.arctan2(np.abs(da), np.abs(dz))
        delta = r - radius(theta)
        b_tot = np.sqrt(bv["x"] ** 2 + bv["y"] ** 2 + bv["z"] ** 2) / abs(b0)
        while True:
            kk = k[idx]
            hit = (kk < len(lev)) & (delta <= lev[np.minimum(kk, len(lev) - 1)])
            if not hit.any():
                break
            j, kh = idx[hit], kk[hit]
            # the values at the level, between the last two points of the path
            w = np.clip((delta_prev[j] - lev[kh]) / np.maximum(delta_prev[j] - delta[hit], 1e-12), 0.0, 1.0)
            out["phi"][j, kh] = before[hit] + w * (phi[j] - before[hit])
            out["phi_null"][j, kh] = before_n[hit] + w * (phi_n[j] - before_n[hit])
            out["b"][j, kh] = b_tot[hit]
            foot = kh == k_zero
            theta_foot[j[foot]] = theta[hit][foot]
            k[j] += 1
        delta_prev[idx] = delta
        alive[idx] = (dz * sz[idx] > 0.5) & (delta < (START + 4.0) * d_i) & (r < r_outer + 1.0) & (k[idx] < len(lev))
    quadrant = np.repeat(np.arange(4), n0)
    return {"theta_foot": theta_foot, "quadrant": quadrant, "clipped": n_clip / max(n_all, 1), **out}


def bin_paths(res, n_units: int = 4):
    """Sums over the paths per (quadrant, bin of the foot angle, level): phi, phi_null, b and the number."""
    n_lev = res["phi"].shape[1]
    sums = np.zeros((4, n_units, N_BIN, n_lev))
    ok = np.isfinite(res["theta_foot"])
    ib = np.minimum((np.degrees(res["theta_foot"][ok]) / (90.0 / N_BIN)).astype(int), N_BIN - 1)
    cell = res["quadrant"][ok] * N_BIN + ib
    for j, key in enumerate(("phi", "phi_null", "b")):
        v = res[key][ok]
        good = np.isfinite(v)
        for lv in range(n_lev):
            sums[j, :, :, lv] = np.bincount(cell[good[:, lv]], weights=v[good[:, lv], lv],
                                            minlength=n_units * N_BIN).reshape(n_units, N_BIN)
            if j == 0:
                sums[3, :, :, lv] = np.bincount(cell[good[:, lv]], minlength=n_units * N_BIN).reshape(n_units, N_BIN)
    return sums


def reduce(sums, min_paths: int = 8):
    """From sums [output, quantity, unit, bin, level] to the profiles and phi at the maximum of the mean field.

    Returns a dict: phi_prof, b_prof [output, bin, level]; k_max [output, bin]; phi, se (scatter of the units over
    the root of their number), phi_null [output, bin]; phi_unit [output, unit, bin]; n_paths [output, bin]."""
    num = sums[:, 3]
    tot = num.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        prof = {key: np.where(tot >= min_paths, sums[:, j].sum(axis=1) / tot, np.nan)
                for j, key in enumerate(("phi", "phi_null", "b"))}
        unit = np.where(num >= 2, sums[:, 0] / num, np.nan)
    b_s = prof["b"].copy()
    b_s[..., 1:-1] = (prof["b"][..., :-2] + prof["b"][..., 1:-1] + prof["b"][..., 2:]) / 3.0
    k_max = np.argmax(np.where(np.isfinite(b_s), b_s, -np.inf), axis=-1)
    seen = np.isfinite(b_s).any(axis=-1)
    pick = np.take_along_axis
    out = {"phi_prof": prof["phi"], "b_prof": prof["b"], "k_max": np.where(seen, k_max, -1),
           "n_paths": pick(tot, k_max[..., None], -1)[..., 0]}
    for key in ("phi", "phi_null", "b"):
        out[key if key != "b" else "b_max"] = np.where(seen, pick(prof[key], k_max[..., None], -1)[..., 0], np.nan)
    u = pick(unit, np.broadcast_to(k_max[:, None, :, None], unit.shape[:3] + (1,)), -1)[..., 0]
    out["phi_unit"] = np.where(seen[:, None, :], u, np.nan)
    n_u = np.isfinite(out["phi_unit"]).sum(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out["se"] = np.where(n_u >= 4, np.nanstd(out["phi_unit"], axis=1, ddof=1) / np.sqrt(np.maximum(n_u, 1)), np.nan)
    return out


def windowed(red, window: int = WINDOW):
    """Running mean over the outputs, per unit, and its error: phi_use, err [output, bin], noise [bin] (rms of
    the same mean of the null run), measured (err below REL_ERR_MAX of phi_use, phi_use > 0)."""
    n_t = red["phi"].shape[0]
    half = window // 2
    use, se, null = (np.full(red["phi"].shape, np.nan) for _ in range(3))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for i in range(n_t):
            sl = slice(max(i - half, 0), min(i + half + 1, n_t))
            u = np.nanmean(red["phi_unit"][sl], axis=0)  # [unit, bin]
            n_u = np.isfinite(u).sum(axis=0)
            use[i] = np.where(n_u >= 4, np.nanmean(u, axis=0), np.nan)
            se[i] = np.where(n_u >= 4, np.nanstd(u, axis=0, ddof=1) / np.sqrt(np.maximum(n_u, 1)), np.nan)
            null[i] = np.nanmean(red["phi_null"][sl], axis=0)
        noise = np.sqrt(np.nanmean(null ** 2, axis=0))
    noise = np.where(np.isfinite(noise), noise, 0.0)
    err = np.maximum(se, noise[None, :])
    return {"phi_use": use, "err": err, "noise": noise, "null_use": null,
            "measured": np.isfinite(use) & np.isfinite(err) & (use > 0.0) & (err < REL_ERR_MAX * np.abs(use))}


def bin_mean(red, win, flag):
    """Mean over the outputs with a detected front, per bin of the foot angle: phi_bin, err_bin, measured_bin,
    outputs_bin. The error is the scatter of the units, or the noise of the windowed null mean reduced by the
    root of the number of windows if that is larger."""
    use = (np.asarray(flag) >= 1) & np.isfinite(red["phi"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        u = np.nanmean(np.where(use[:, None, :], red["phi_unit"], np.nan), axis=0)  # [unit, bin]
        n_u = np.isfinite(u).sum(axis=0)
        mean = np.where(n_u >= 4, np.nanmean(u, axis=0), np.nan)
        se = np.where(n_u >= 4, np.nanstd(u, axis=0, ddof=1) / np.sqrt(np.maximum(n_u, 1)), np.nan)
    n_t = use.sum(axis=0)
    err = np.maximum(se, win["noise"] / np.sqrt(np.maximum(n_t / WINDOW, 1.0)))
    return {"phi_bin": mean, "err_bin": err, "outputs_bin": n_t,
            "measured_bin": np.isfinite(mean) & np.isfinite(err) & (mean > 0.0) & (err < REL_ERR_MAX * np.abs(mean))}


class PhiMap:
    """phi(theta, t) for the count model, from analysis/phi.npz: the windowed mean of the nearest output and bin
    where it is measured (source 2); else the mean of the bin over the outputs where that is measured (source 1);
    else 0 (source 0)."""

    def __init__(self, t, theta, phi_use, measured, phi_bin, measured_bin):
        self.t, self.theta = np.asarray(t, float), np.asarray(theta, float)
        self.phi, self.measured = np.asarray(phi_use, float), np.asarray(measured, bool)
        self.phi_bin, self.measured_bin = np.asarray(phi_bin, float), np.asarray(measured_bin, bool)

    @classmethod
    def load(cls, run_dir):
        path = Path(run_dir) / "analysis" / "phi.npz"
        if not path.exists():
            return None
        z = np.load(path)
        return cls(z["t"], z["theta"], z["phi_use"], z["measured"], z["phi_bin"], z["measured_bin"])

    def __call__(self, theta, t):
        theta, t = np.broadcast_arrays(np.asarray(theta, float), np.asarray(t, float))
        it = np.abs(t[..., None] - self.t).argmin(axis=-1)
        ib = np.abs(theta[..., None] - self.theta).argmin(axis=-1)
        fine, coarse = self.measured[it, ib], self.measured_bin[ib]
        phi = np.where(fine, self.phi[it, ib], np.where(coarse, self.phi_bin[ib], 0.0))
        return phi, np.where(fine, 2, np.where(coarse, 1, 0))


def cmd_phi(args) -> int:
    run_dir = Path(args.run)
    run = census._open(run_dir)
    null = census._open(Path(args.null)) if args.null else None
    model = census.FrontModel.from_run(run, run_dir)
    d_i, b0 = model.d_i, model.b0
    r_outer = 0.5 * float(run.L[2]) - 8.0
    planes = [(run.series(name), null.series(name) if null else None, across, ax)
              for name, across, ax in PLANES if name in run.names()]
    t = np.asarray(planes[0][0].t, float)
    if null is not None and not np.allclose(t, np.asarray(planes[0][1].t, float)):
        print("the outputs of the run and of the null run are not at the same times")
        return 1
    lev = levels(d_i)
    sums = np.zeros((len(t), 4, 4 * len(planes), N_BIN, len(lev)))
    clipped = []
    for i, ti in enumerate(t):
        if not model.t_range[0] <= ti <= model.t_range[1]:
            continue
        for p, (ser, ser0, across, ax) in enumerate(planes):
            fields = {f + c: (ser.read(f + c, i), ser.plane_axes(f + c)) for f in "EB" for c in "xyz"}
            e_null = {"E" + c: (ser0.read("E" + c, i), ser0.plane_axes("E" + c)) for c in "xyz"} if ser0 else None
            res = trace_plane(fields, (model.centre[ax], model.centre[2]),
                              lambda th, ti=ti: model.front(th, np.full(np.shape(th), ti))[0], d_i, b0, across,
                              e_null=e_null, r_outer=r_outer)
            sums[i, :, 4 * p:4 * p + 4] = bin_paths(res)
            clipped.append(res["clipped"])
    red = reduce(sums)
    win = windowed(red)
    theta = np.radians((np.arange(N_BIN) + 0.5) * 90.0 / N_BIN)
    tt, th = np.meshgrid(t, theta, indexing="ij")
    flag = model.flag(th, tt)
    temp = np.array([model.temp_upstream(x) for x in t])
    f = model.front(th, tt)
    v_n = f[1] * f[0] / np.hypot(f[0], f[2])
    mean = bin_mean(red, win, flag)
    out = summarize(t, theta, lev, red, win, flag, temp, v_n, float(run.config.get("mass_ratio", 16.0)),
                    null is not None)
    out["mean over the outputs with a detected front"] = [
        {"theta_deg": float(np.degrees(th)), "outputs": int(mean["outputs_bin"][j]),
         "phi": None if not np.isfinite(mean["phi_bin"][j]) else float(mean["phi_bin"][j]),
         "err": None if not np.isfinite(mean["err_bin"][j]) else float(mean["err_bin"][j]),
         "measured": bool(mean["measured_bin"][j])} for j, th in enumerate(theta)]
    out.update({"run": str(run_dir), "null": str(args.null) if args.null else None,
                "commit": run.info.get("git_commit"), "clipped steps": float(np.mean(clipped)) if clipped else None})
    np.savez_compressed(run_dir / "analysis" / "phi.npz", t=t, theta=theta, levels=lev, flag=flag, **red, **win,
                        **mean)
    (run_dir / "analysis" / "phi.json").write_text(json.dumps(out, indent=1))
    report(out)
    print(f"wrote {run_dir / 'analysis' / 'phi.json'} and phi.npz")
    return 0


def _med(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else None


def summarize(t, theta, lev, red, win, flag, temp, v_n, mass_ratio: float, has_null: bool) -> dict:
    rows = []
    for dom, fmin in (("detected front", 1), ("valid domain", 2)):
        for j, th in enumerate(theta):
            s = (flag[:, j] >= fmin) & np.isfinite(red["phi"][:, j])
            row = {"domain": dom, "theta_deg": float(np.degrees(th)), "outputs": int(s.sum())}
            if s.any():
                with np.errstate(divide="ignore", invalid="ignore"):
                    ion = win["phi_use"][s, j] / (0.5 * mass_ratio * v_n[s, j] ** 2)
                row.update({"phi": _med(red["phi"][s, j]), "se": _med(red["se"][s, j]),
                            "phi windowed": _med(win["phi_use"][s, j]), "err windowed": _med(win["err"][s, j]),
                            "measured": int(win["measured"][s, j].sum()),
                            "noise of the null run": float(win["noise"][j]),
                            "phi over T_e": _med(win["phi_use"][s, j] / temp[s]),
                            "phi over the ion ram energy": _med(ion), "b_max on the paths": _med(red["b_max"][s, j]),
                            "level of b_max (d_e behind the front)": _med(-lev[red["k_max"][s, j]]),
                            "paths": _med(red["n_paths"][s, j])})
            rows.append(row)
    seen = (flag >= 1) & np.isfinite(red["phi"])
    return {"start_d_i": START, "depth_d_i": DEPTH, "window_outputs": WINDOW, "null run used": has_null,
            "definition": "phi = - sum of (E.B) / B_pol ds along the in-plane field line from START d_i ahead of the "
                          "density front to the level of the largest mean |B|; error: scatter of the quadrants over "
                          "the root of their number, or the rms of the same mean of the null run if larger",
            "rows": rows, "bins on the detected front": int(seen.sum()),
            "measured (error below 50 %)": int((win["measured"] & seen).sum()),
            "measured, valid domain": int((win["measured"] & seen & (flag >= 2)).sum()),
            "bins of the valid domain": int((seen & (flag >= 2)).sum()),
            "significantly negative": int((seen & (win["phi_use"] < -2.0 * win["err"])).sum())}


def report(out) -> None:
    print(f"phi of {out.get('run', 'synthetic')}: rise of the potential to the largest mean |B| "
          f"(m_e c^2 / e; window of {out['window_outputs']} outputs)")
    for r in out["rows"]:
        if r["outputs"]:
            ion = r["phi over the ion ram energy"]
            print(f"  {r['domain']:14s} theta {r['theta_deg']:5.1f}  outputs {r['outputs']:2d}  phi {r['phi']:8.4f} "
                  f"(se {r['se'] if r['se'] is not None else float('nan'):6.4f})  windowed "
                  f"{r['phi windowed']:8.4f} +- {r['err windowed']:6.4f}  measured {r['measured']:2d}  "
                  f"phi/T_e {r['phi over T_e']:6.2f}  over ion ram {ion if ion is not None else float('nan'):6.3f}  "
                  f"b_max {r['b_max on the paths']:5.2f}  null noise {r['noise of the null run']:6.4f}")
    for r in out.get("mean over the outputs with a detected front", []):
        if r["phi"] is not None:
            print(f"  mean over {r['outputs']:2d} outputs  theta {r['theta_deg']:5.1f}  phi {r['phi']:8.4f} +- "
                  f"{r['err']:6.4f}  " + ("measured" if r["measured"] else "not measured"))
    print(f"windowed: measured in {out['measured (error below 50 %)']} of {out['bins on the detected front']} bins of "
          f"the detected front, {out['measured, valid domain']} of {out['bins of the valid domain']} of the valid "
          f"domain; significantly negative in {out['significantly negative']}")


# ------------------------------------------------------------------ checks on synthetic data
def selfcheck(check) -> None:
    n, dx, d_i, b0, r_f = 400, 0.5, 4.0, 0.3, 60.0
    ax = (np.arange(n) + 0.5) * dx
    c = 0.5 * n * dx
    aa, zz = np.meshgrid(ax - c, ax - c, indexing="ij")
    r = np.hypot(aa, zz)

    def radius(th):
        return np.full(np.shape(th), r_f)

    def pack(**comp):
        return {k: (np.broadcast_to(np.asarray(v, float), (n, n)).copy(), (ax, ax)) for k, v in comp.items()}

    # 1. a potential that depends on the radius, a tilted and sheared field with an out-of-plane part: the rise
    # along any field line is the difference of the potential
    pot0, width = 0.08, 3.0

    def pot(rr):
        return 0.5 * pot0 * (1.0 - np.tanh((rr - r_f) / width))

    dpot = -0.5 * pot0 / width / np.cosh((r - r_f) / width) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        ex, ez = np.where(r > 0, -dpot * aa / r, 0.0), np.where(r > 0, -dpot * zz / r, 0.0)
    fields = pack(Ex=ex, Ey=0.0, Ez=ez, Bx=0.15 * b0 * np.tanh(aa * zz / 400.0), By=0.5 * b0,
                  Bz=b0 * (1.0 + 0.3 * pot(r) / pot0))
    res = trace_plane(fields, (c, c), radius, d_i, b0, "x")
    lev = levels(d_i)
    want = pot(r_f + lev) - pot(r_f + START * d_i)
    err = np.nanmax(np.abs(res["phi"] - want[None, :]))
    reached = np.isfinite(res["theta_foot"])
    red = reduce(bin_paths(res)[None])
    deep = red["phi"][0, :12]  # the field lines that meet the front beyond 60 degrees do not go deep
    ok1 = err < 1e-3 and reached.sum() > 500 and np.max(np.abs(deep - pot0)) < 0.01 * pot0 \
        and red["phi"][0, -1] < 0.8 * pot0
    check("phi, potential that depends on the radius", ok1,
          f"largest error {err:.1e} for a rise of {pot0}; {int(reached.sum())} paths reach the front; at the "
          f"largest |B| below 60 degrees: {deep.min():.4f} to {deep.max():.4f}; at 87.5 degrees, where the "
          f"field lines graze the front: {red['phi'][0, -1]:.4f}")

    # 2. the out-of-plane part: E_y B_y along a straight field line; the sign follows the hemisphere
    e_y, c_y = 0.01, 0.5
    fields = pack(Ex=0.0, Ey=e_y, Ez=0.0, Bx=0.0, By=c_y * b0, Bz=b0)
    res = trace_plane(fields, (c, c), radius, d_i, b0, "x")
    k0 = int(round(START * d_i / LEVEL_STEP))
    th0 = np.radians(np.arange(0.5 * DTHETA0, 90.0, DTHETA0))
    rho = np.tile((r_f + START * d_i) * np.sin(th0), 4)
    sz = np.repeat([1.0, 1.0, -1.0, -1.0], len(th0))
    hit = np.isfinite(res["theta_foot"])
    length = np.sqrt((r_f + START * d_i) ** 2 - rho[hit] ** 2) - np.sqrt(r_f ** 2 - rho[hit] ** 2)
    want2 = sz[hit] * e_y * c_y * length
    err2 = np.max(np.abs(res["phi"][hit, k0] - want2))
    foot_err = np.max(np.abs(np.sin(res["theta_foot"][hit]) - rho[hit] / r_f))
    ok2 = err2 < 0.02 * np.max(np.abs(want2)) and foot_err < 0.01 and np.any(want2 < 0)
    check("phi, out-of-plane part and the angle of the foot", ok2,
          f"largest error {err2:.1e} for values up to {np.max(np.abs(want2)):.3f}; sine of the foot angle within "
          f"{foot_err:.1e}")

    # 3. noise: the E of a null run along the paths of test 1 gives no rise beyond its own scatter, and a rise
    # with noise is called measured only where the error is small
    rng = np.random.default_rng(5)
    fields = pack(Ex=ex, Ey=0.0, Ez=ez, Bx=0.0, By=0.0, Bz=b0 * (1.0 + 0.3 * pot(r) / pot0))
    sums = []
    for _ in range(6):
        noise = {k: (rng.normal(0.0, 0.03, (n, n)), (ax, ax)) for k in ("Ex", "Ey", "Ez")}
        noisy = dict(fields)
        for k in noise:
            noisy[k] = (fields[k][0] + rng.normal(0.0, 0.03, (n, n)), (ax, ax))
        sums.append(bin_paths(trace_plane(noisy, (c, c), radius, d_i, b0, "x", e_null=noise)))
    red = reduce(np.array(sums))
    win = windowed(red)
    seen = np.isfinite(win["phi_use"])
    seen[:, 12:] = False
    pull = np.abs(win["null_use"][seen]) / np.maximum(win["err"][seen], 1e-12)
    dev = np.abs(win["phi_use"][seen] - pot0) / win["err"][seen]
    mean = bin_mean(red, win, np.ones(win["phi_use"].shape, int))
    pm = PhiMap(np.arange(6.0), np.radians((np.arange(N_BIN) + 0.5) * 5.0), win["phi_use"],
                win["measured"] & (np.arange(6)[:, None] < 3), mean["phi_bin"], mean["measured_bin"])
    got, src = pm(np.radians([12.0, 12.0, 88.0]), np.array([0.2, 4.6, 4.6]))
    ok3 = np.all(pull < 4.0) and np.mean(dev < 3.0) > 0.9 and win["measured"][seen].mean() > 0.8 \
        and np.all(mean["measured_bin"][:12]) \
        and np.all(np.abs(mean["phi_bin"][:12] - pot0) < 3.0 * mean["err_bin"][:12]) \
        and src.tolist()[:2] == [2, 1] and got[0] == win["phi_use"][0, 2] and got[1] == mean["phi_bin"][2]
    check("phi, noise of a null run", ok3,
          f"null within {pull.max():.1f} errors; the rise within 3 errors in {100 * np.mean(dev < 3.0):.0f} % of "
          f"{int(seen.sum())} bins (median error {np.median(win['err'][seen]):.4f}); measured in "
          f"{100 * win['measured'][seen].mean():.0f} %; mean over the outputs within "
          f"{np.max(np.abs(mean['phi_bin'][:12] - pot0) / mean['err_bin'][:12]):.1f} errors")
