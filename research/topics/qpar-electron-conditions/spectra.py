"""Shock position, downstream spectra and upstream monitors of one wall run (README, measured quantities 1 and 2).

    python research/topics/qpar-electron-conditions/spectra.py <run dir> [--json <file>] [--floor <N1 run dir>]

CPU only (reads the output with ``pypic.analysis``). It applies the definitions that README.md ("Approach") fixes;
where the README leaves a detail open, the choice made here is listed below and in notebook.md. The defaults are
V1's values; a short pilot needs ``--blocks``, ``--early``, ``--region``, ``--q-zone`` and ``--edge`` of its own.

Units: lengths in d_i, times in 1 / Omega_ci, momenta in m_i v_pt (v_pt = v0, the speed of the flow in the box).
The box is the downstream frame, and the histograms ``spectrum_<species>`` hold |u| = |p| / (m c) in that frame.

- ``shock``: x_sh(t) is the largest x where the ion density, averaged over ``--smooth`` d_i, is above ``--level``
  (in units of the loaded density). The same with a tenth of the smoothing length (``fine``) and with 1.5 times
  the level (``high``) is given beside it.
- ``early`` (V1's early stop): the slope of x_sh(t) over the outputs of ``--early`` (t = 10-20), against v0 / 3;
  the compression is the ion density averaged over 0.2-0.8 x_sh and over those outputs. Criteria: the speed within
  +-10 % of v0 / 3, the compression within 3.6-4.4.
- ``blocks``: the spectrum of each species in D_fix = [x_sh - 15 M_A d_i, x_sh - 2.5 M_A d_i] (``--region``),
  summed over the outputs of each block of ``--blocks`` (t0, t1, length; V1: 100, 150, 10). The region moves with
  x_sh of each output; a bin of x enters with the fraction of its width that lies in the region. An output is used
  only if the far edge of the region is at least ``--clear`` M_A d_i (10) from the wall. Per block and species:
  - p_th = (momentum of the peak of p^4 f) / sqrt(2); the peak is the vertex of a parabola in ln p through the
    highest bin of p * counts and its two neighbours. p_0 = 5 p_th.
  - eta_nt = the share of the particles with p >= p_0 (a bin enters by the part of its width in ln p above p_0).
  - Electrons: phi = 4 pi p_ref^3 f(p_ref) / n at p_ref = 6 sqrt(m_e / m_i) m_i v_pt, taken as
    (share of the particles within +-0.1 dex of p_ref) / ln(10^0.2), i.e. p^3 f averaged over ln p. ``phi_valid``:
    every bin of that band holds >= 100 particles and, with ``--floor``, phi >= 5 times the same number of the N1
    run (whole box, the outputs of the same block; an empty band counts as one particle).
    eta_fixed = the share with p >= 3 sqrt(m_e / m_i) m_i v_pt.
    q = minus the slope of ln f against ln p over [0.6, 2.4] m_i v_pt (``--fit``), f = counts / p^3 per bin,
    weighted by the counting error; its standard error is that of the fit. Criterion in V1 (last block):
    |q - 4| <= 0.3.
  - ``d_lnp``: counts / (all counts) per bin, for figures.
- ``series``: p_th, eta_nt, phi, eta_fixed per output of the blocks (the region of that output alone).
- ``sigma_pt``: the standard deviation of the block values about a straight line through them (needs >= 4 blocks).
- ``sensitivity`` (last block): the same numbers in [x_sh - 15, x_sh] M_A d_i and in the two halves of D_fix;
  ``spread`` = the largest |value / value in D_fix - 1| of the three. Criterion: above 0.3 the spread is stated as
  the error of the definition.
- ``q_monitor``: Q = (<B_y> / B0y) / (<n_i> / n0) in the zone from ``--q-zone`` (5 to 200) d_i inside the absorbing
  layer, per block of ``--q-block`` (5); the means are over the zone and the outputs of the block. Criterion:
  |Q - 1| <= 0.05 in every block. The distance of the shock from the zone is given.
- ``edge`` (in place of N6 in V1, part b): dB_perp = the rms of B_y - <B_y> and of B_z over the last ``--edge-zone``
  (200) d_i before the absorbing layer, in units of B0, at the outputs nearest to ``--edge`` (t = 100, 125, 150).
  With ``--floor``: the same over the whole box of the N1 run at the nearest output; criterion <= 2 times that.
  Beside it the rms within windows of the N1 box length (the mean square over the windows), which leaves out the
  wavelengths that the N1 box cannot hold.

Not in this file: p_max (its definition is not yet fixed), N12 and the escape record of N6 (a), which need tracked
particles and a diagnostic at the face.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from pypic.analysis.run import open_run

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pilot_check import line, shock_position  # noqa: E402

DEX = 0.1  # half width of the band of phi


class Spectrum:
    """Histogram (x, u) of one species as particles per bin: ``cnt`` (n_t, n_x, n_u), ``t`` in 1 / Omega_ci,
    ``xe`` in d_i, ``p`` and ``pe`` (bin centres and edges) in m_i v_pt."""

    def __init__(self, run, name: str):
        info = run.meta["planner"]
        h = run.series(f"spectrum_{name}")
        (_, ny, nz), dx = run.info["n_cells"], run.info["dx"]
        per_weight = info["per_column"] / (2 * ny * nz) / (dx[0] * dx[1] * dx[2])
        mass = next(s["mass"] for s in run.info["species"] if s["name"] == name)
        self.cnt = np.asarray(h[name], dtype=np.float64) * per_weight
        self.t = np.asarray(h.t) * info["omega_ci"]
        self.xe = np.asarray(h.edges[0]) / info["d_i"]
        mi = next(s["mass"] for s in run.info["species"] if s["name"] == "i")
        self.pe = np.asarray(h.edges[1]) * mass / (mi * info["v0"])
        self.p = np.sqrt(self.pe[1:] * self.pe[:-1])
        self.dlnp = float(np.log(self.pe[1] / self.pe[0]))

    def region(self, j: int, a: float, b: float):
        """Particles per bin of p in a <= x <= b (d_i) at the output j."""
        w = np.clip((np.minimum(self.xe[1:], b) - np.maximum(self.xe[:-1], a)) / np.diff(self.xe), 0.0, 1.0)
        return w @ self.cnt[j]


def above(c, pe, p: float) -> float:
    """Particles with momentum >= p; the bin that holds p enters by the part of its width in ln p above p."""
    if p <= pe[0]:
        return float(c.sum())
    if p >= pe[-1]:
        return 0.0
    k = int(np.searchsorted(pe, p, side="right")) - 1
    return float(c[k + 1:].sum() + c[k] * np.log(pe[k + 1] / p) / np.log(pe[k + 1] / pe[k]))


def peak_p4f(c, p) -> float:
    """Momentum of the peak of p^4 f (p * counts for bins of equal width in ln p)."""
    y = c * p
    k = int(np.argmax(y))
    if k == 0 or k == len(y) - 1:
        return float(p[k])
    a, b, d = y[k - 1], y[k], y[k + 1]
    shift = 0.5 * (a - d) / (a - 2.0 * b + d) if a - 2.0 * b + d else 0.0
    return float(p[k] * np.exp(shift * np.log(p[k + 1] / p[k])))


def slope(c, p, lo: float, hi: float):
    """(q, standard error, bins, smallest count) of f ~ p^-q over lo <= p <= hi; f = counts / p^3 per bin."""
    sel = (p >= lo) & (p <= hi) & (c > 0)
    if sel.sum() < 4:
        return None
    co, cov = np.polyfit(np.log(p[sel]), np.log(c[sel] / p[sel] ** 3), 1, w=np.sqrt(c[sel]), cov="unscaled")
    return {"q": float(-co[0]), "stderr": float(np.sqrt(cov[0, 0])), "bins": int(sel.sum()),
            "min_count": float(c[sel].min())}


def numbers(c, sp: Spectrum, electron: bool, root: float, fit, floor=None) -> dict:
    """The numbers of one spectrum ``c`` (particles per bin of p). ``root`` = sqrt(m_e / m_i)."""
    n = float(c.sum())
    out = {"particles": n}
    if n <= 0:
        return out
    p_th = peak_p4f(c, sp.p) / np.sqrt(2.0)
    n_nt = above(c, sp.pe, 5.0 * p_th)
    out |= {"p_th": p_th, "p_0": 5.0 * p_th, "eta_nt": n_nt / n, "particles_above_p_0": n_nt}
    if not electron:
        return out
    p_ref = 6.0 * root
    lo, hi = p_ref * 10.0 ** -DEX, p_ref * 10.0 ** DEX
    band = above(c, sp.pe, lo) - above(c, sp.pe, hi)
    inside = (sp.pe[1:] > lo) & (sp.pe[:-1] < hi)
    width = 2.0 * DEX * np.log(10.0)
    out |= {"p_ref": p_ref, "phi": band / n / width, "phi_band_particles": band,
            "phi_min_bin": float(c[inside].min()) if inside.any() else 0.0,
            "eta_fixed": above(c, sp.pe, 3.0 * root) / n, "slope": slope(c, sp.p, *fit)}
    valid = out["phi_min_bin"] >= 100.0
    if floor is not None:
        nf = float(floor.sum())
        band_f = above(floor, sp.pe, lo) - above(floor, sp.pe, hi)
        out["phi_floor"] = max(band_f, 1.0) / nf / width
        out["phi_floor_band_particles"] = band_f
        out["eta_nt_floor"] = above(floor, sp.pe, 5.0 * p_th) / nf
        valid = valid and out["phi"] >= 5.0 * out["phi_floor"]
    out["phi_valid"] = bool(valid)
    return out


def detrended(t, y) -> float:
    c = np.polyfit(t, y, 1)
    return float(np.sqrt(((y - np.polyval(c, t)) ** 2).sum() / (len(t) - 2)))


def rms_perp(by, bz, b0: float, chunk: int | None = None):
    """rms of B_y - <B_y> and of B_z over one lineout, in units of b0; with ``chunk`` the mean square within
    windows of that many samples."""
    if chunk and chunk < len(by):
        k = len(by) // chunk
        by, bz = by[:k * chunk].reshape(k, chunk), bz[:k * chunk].reshape(k, chunk)
        var = ((by - by.mean(axis=1, keepdims=True)) ** 2).mean() + (bz * bz).mean()
    else:
        var = ((by - by.mean()) ** 2).mean() + (bz * bz).mean()
    return float(np.sqrt(var) / b0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--json")
    ap.add_argument("--floor", help="the N1 run of the same plasma")
    ap.add_argument("--smooth", type=float, default=5.0, help="d_i")
    ap.add_argument("--level", type=float, default=2.0)
    ap.add_argument("--early", type=float, nargs=2, default=(10.0, 20.0))
    ap.add_argument("--blocks", type=float, nargs=3, default=(100.0, 150.0, 10.0), help="t0 t1 length")
    ap.add_argument("--region", type=float, nargs=2, default=(15.0, 2.5), help="M_A d_i behind the shock")
    ap.add_argument("--clear", type=float, default=10.0, help="M_A d_i between the wall and the region")
    ap.add_argument("--fit", type=float, nargs=2, default=(0.6, 2.4), help="m_i v_pt")
    ap.add_argument("--q-zone", type=float, nargs=2, default=(5.0, 200.0), help="d_i inside the absorbing layer")
    ap.add_argument("--q-block", type=float, default=5.0)
    ap.add_argument("--edge", type=float, nargs="*", default=(100.0, 125.0, 150.0))
    ap.add_argument("--edge-zone", type=float, default=200.0, help="d_i")
    a = ap.parse_args()

    run = open_run(a.run)
    info = run.meta["planner"]
    if not info["wall"]:
        raise SystemExit("not a wall run")
    d_i, w_ci, v0 = info["d_i"], info["omega_ci"], info["v0"]
    m_a = info["shock"]["mach_r4"]
    mi = next(s["mass"] for s in run.info["species"] if s["name"] == "i")
    root = 1.0 / np.sqrt(mi)
    cfg = run.info["config"]
    out = {"run": a.run, "commit": run.info.get("git_commit"), "mach_alfven": m_a, "mass_ratio": mi,
           "M_A_d_i_in_d_i": m_a}

    # shock position and the early-stop numbers
    prof = run.series("profile")
    x, t = prof.x("n:i") / d_i, prof.t * w_ci
    dx = x[1] - x[0]
    n0 = float(prof["n:i"][0][(x > 0.1 * x[-1]) & (x < 0.9 * x[-1])].mean())
    n_i = prof["n:i"] / n0
    xs = shock_position(x, n_i, max(1, int(round(a.smooth / dx))), a.level)
    xs_fine = shock_position(x, n_i, max(1, int(round(0.1 * a.smooth / dx))), a.level)
    xs_high = shock_position(x, n_i, max(1, int(round(a.smooth / dx))), 1.5 * a.level)
    pick = np.unique(np.linspace(0, len(t) - 1, 16).astype(int))
    expected = v0 / 3.0
    to_c = d_i * w_ci  # (d_i per 1 / Omega_ci) in c
    out["shock"] = {"smooth_di": a.smooth, "level": a.level, "omega_ci_t": t[pick].tolist(),
                    "x_shock_di": xs[pick].tolist(), "x_shock_di_fine": xs_fine[pick].tolist(),
                    "x_shock_di_high": xs_high[pick].tolist(),
                    "largest_difference_di": {"fine": float(np.abs(xs - xs_fine).max()),
                                              "high": float(np.abs(xs - xs_high).max())},
                    "mhd": {"r": info["shock"]["r_mhd"], "speed_c": v0 / (info["shock"]["r_mhd"] - 1.0)}}
    sel = (t >= a.early[0]) & (t <= a.early[1] + 1e-6)
    if sel.sum() >= 4 and t[-1] >= a.early[1] - 0.5 * (t[1] - t[0]):
        s, e = line(t[sel], xs[sel])
        sf, _ = line(t[sel], xs_fine[sel])
        sh, _ = line(t[sel], xs_high[sel])
        comp = [float(n_i[j][(x >= 0.2 * xs[j]) & (x <= 0.8 * xs[j])].mean()) for j in np.nonzero(sel)[0]]
        out["early"] = {"window": list(a.early), "outputs": int(sel.sum()), "speed_c": s * to_c,
                        "speed_stderr_c": e * to_c, "speed_c_fine": sf * to_c, "speed_c_high": sh * to_c,
                        "expected_c": expected,
                        "speed_over_expected": s * to_c / expected, "compression": float(np.mean(comp)),
                        "compression_range": [float(min(comp)), float(max(comp))],
                        "pass_speed": bool(abs(s * to_c / expected - 1.0) <= 0.10),
                        "pass_compression": bool(3.6 <= np.mean(comp) <= 4.4)}

    # spectra in the region behind the shock
    t0, t1, length = a.blocks
    edges = np.arange(t0, t1 + 1e-9, length)
    far, near = a.region[0] * m_a, a.region[1] * m_a
    regions = {"D_fix": (far, near), "to_shock": (far, 0.0), "far_half": (far, 0.5 * (far + near)),
               "near_half": (0.5 * (far + near), near)}
    floor_run = open_run(a.floor) if a.floor else None
    blocks, series, last = {}, {}, {}
    for name in ("e", "i"):
        if f"spectrum_{name}" not in run.names("histogram"):
            continue
        sp = Spectrum(run, name)
        fl = Spectrum(floor_run, name) if floor_run else None
        x_h = np.interp(sp.t, t, xs)
        ok = x_h - far >= a.clear * m_a - 1e-9
        rows, per_output = [], []
        for k in range(len(edges) - 1):
            hi_in = k == len(edges) - 2
            inb = (sp.t >= edges[k] - 1e-6) & ((sp.t <= edges[k + 1] + 1e-6) if hi_in else (sp.t < edges[k + 1] - 1e-6))
            use = np.nonzero(inb & ok)[0]
            row = {"t": [float(edges[k]), float(edges[k + 1])], "outputs": int(len(use)),
                   "outputs_too_close_to_wall": int((inb & ~ok).sum())}
            if len(use):
                floor = None
                if fl is not None:
                    fb = (fl.t >= edges[k] - 1e-6) & (fl.t <= edges[k + 1] + 1e-6)
                    floor = fl.cnt[fb].sum(axis=(0, 1)) if fb.any() else None
                reg = {r: sum(sp.region(j, x_h[j] - lo, x_h[j] - hi) for j in use) for r, (lo, hi) in regions.items()}
                row |= numbers(reg["D_fix"], sp, name == "e", root, a.fit, floor)
                row["x_shock_di"] = [float(x_h[use[0]]), float(x_h[use[-1]])]
                row["d_lnp"] = (reg["D_fix"] / reg["D_fix"].sum()).tolist()
                last = last | {name: {r: numbers(c, sp, name == "e", root, a.fit, floor) for r, c in reg.items()}}
                for j in use:
                    one = numbers(sp.region(j, x_h[j] - far, x_h[j] - near), sp, name == "e", root, a.fit)
                    per_output.append({"t": float(sp.t[j])} | {q: one[q] for q in ("p_th", "eta_nt", "phi", "eta_fixed")
                                                             if q in one})
            rows.append(row)
        blocks[name] = {"p": sp.p.tolist(), "rows": rows}
        series[name] = per_output
    out["blocks"] = blocks
    out["series"] = series

    keys = ("p_th", "eta_nt", "phi", "eta_fixed", "q")

    def value(rec, q):
        if q == "q":
            return rec["slope"]["q"] if rec.get("slope") else None
        return rec.get(q)

    sigma, sens = {}, {}
    for name, b in blocks.items():
        rows = [r for r in b["rows"] if r["outputs"]]
        tm = np.array([0.5 * (r["t"][0] + r["t"][1]) for r in rows])
        sigma[name] = {"blocks": len(rows)}
        for q in keys:
            v = [value(r, q) for r in rows]
            if len(rows) >= 4 and all(u is not None for u in v):
                v = np.array(v)
                sigma[name][q] = {"values": v.tolist(), "sigma": detrended(tm, v),
                                  "sigma_over_last": detrended(tm, v) / v[-1] if v[-1] else None}
        if name in last:
            ref = last[name]["D_fix"]
            sens[name] = {}
            for q in keys:
                if value(ref, q) is None:
                    continue
                v = {r: value(rec, q) for r, rec in last[name].items()}
                rel = [abs(u / v["D_fix"] - 1.0) for r, u in v.items()
                       if r != "D_fix" and u is not None and v["D_fix"] and np.isfinite(v["D_fix"])]
                sens[name][q] = v | {"spread": max(rel) if rel else None,
                                     "above_30_percent": bool(rel and max(rel) > 0.30)}
    out["sigma_pt"] = sigma
    out["sensitivity"] = sens
    if "e" in last and last["e"]["D_fix"].get("slope"):
        q = last["e"]["D_fix"]["slope"]["q"]
        out["slope_last_block"] = last["e"]["D_fix"]["slope"] | {"pass": bool(abs(q - 4.0) <= 0.3)}

    # upstream monitors
    xb = prof.x("By") / d_i
    layer = (run.info["L"][0] - cfg["npml"] * run.info["dx"][0]) / d_i
    by0, b0 = info["by"], info["b"]
    if by0:
        zq = (xb >= layer - a.q_zone[1]) & (xb <= layer - a.q_zone[0])
        zn = (x >= layer - a.q_zone[1]) & (x <= layer - a.q_zone[0])
        rows = []
        for k in range(int(np.floor((t[-1] + 1e-6) / a.q_block))):
            sel = (t >= k * a.q_block - 1e-6) & (t < (k + 1) * a.q_block - 1e-6)
            if not sel.any():
                continue
            b_mean, n_mean = float(prof["By"][sel][:, zq].mean() / by0), float(n_i[sel][:, zn].mean())
            rows.append({"t": [k * a.q_block, (k + 1) * a.q_block], "By": b_mean, "n": n_mean, "Q": b_mean / n_mean,
                         "shock_to_zone_di": float(xb[zq][0] - xs[sel].max())})
        qs = np.array([r["Q"] for r in rows])
        out["q_monitor"] = {"zone_di": [float(xb[zq][0]), float(xb[zq][-1])], "block": a.q_block, "rows": rows,
                            "largest_deviation": float(np.abs(qs - 1.0).max()) if len(qs) else None,
                            "pass": bool(len(qs) and np.abs(qs - 1.0).max() <= 0.05)}
    ze = (xb >= layer - a.edge_zone) & (xb <= layer)
    fp = floor_run.series("profile") if floor_run else None
    rows = []
    for te in a.edge:
        if te > t[-1] + 0.5 * (t[1] - t[0]):
            continue
        j = int(np.argmin(np.abs(t - te)))
        row = {"t": float(t[j]), "dB_perp": rms_perp(prof["By"][j][ze], prof["Bz"][j][ze], b0),
               "shock_to_zone_di": float(xb[ze][0] - xs[j])}
        if fp is not None:
            tf = fp.t * floor_run.meta["planner"]["omega_ci"]
            jf = int(np.argmin(np.abs(tf - te)))
            nf = len(fp["By"][jf])
            row |= {"n1_t": float(tf[jf]), "n1_dB_perp": rms_perp(fp["By"][jf], fp["Bz"][jf], b0),
                    "dB_perp_n1_windows": rms_perp(prof["By"][j][ze], prof["Bz"][j][ze], b0, nf)}
            row["ratio"] = row["dB_perp"] / row["n1_dB_perp"]
            row["pass"] = bool(row["ratio"] <= 2.0)
        rows.append(row)
    out["edge"] = {"zone_di": [float(xb[ze][0]), float(xb[ze][-1])], "rows": rows}

    text = json.dumps(out, indent=1)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")
    brief = {k: v for k, v in out.items() if k not in ("blocks", "series")}
    brief["blocks"] = {n: [{k: v for k, v in r.items() if k != "d_lnp"} for r in b["rows"]] for n, b in blocks.items()}
    print(json.dumps(brief, indent=1))


if __name__ == "__main__":
    main()
