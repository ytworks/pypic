"""N12 (wall exclusions) and the escape record of N6 (a), from the particle dumps of a run (README, null tests).

    python research/topics/qpar-electron-conditions/n12_n6.py <run dir> --json <file> [options]

CPU only (NumPy). The regions, blocks, momentum bins and numbers are the ones of spectra.py (imported, not copied).
Units as there: x in d_i, t in 1 / Omega_ci, momenta in m_i v_pt (box frame), particle numbers in macro-particles;
energies in m_e c^2 per macro-particle.

N12. For every spectra output of the blocks, the dump ``origin_e`` (electrons with p >= p_d,e) of the same step is
binned into the (x, p) bins of the histogram, and the spectrum of the excluded electrons is subtracted from the
histogram's spectrum of the region. The particle number and p_th of every variant are the ones of the histogram
without exclusion. Variants: ``1_L`` (L) at X_c, ``2_none``, ``3_L_lo`` / ``3_L_hi`` (L) at the two other X_c,
``4_E`` (E), ``5_LE`` (L) or (E).

- (L): a loaded electron (id below cells x ppc of the electron stream) has x_0 = (cell + 0.5) dx from its id alone;
  left out if x_0 <= X_c. Injected electrons are never left out.
- (E): left out if the id is in any dump ``wall_e`` taken at or before the step of the spectra output.

Details that the README leaves open, fixed here before a dump of V1 is read:

- A number enters the verdict only if the lowest momentum it uses is at or above p_d,e (eta_nt: p_0; phi: the lower
  edge of its band; eta_fixed: 3 sqrt(m_e / m_i); q: the lower limit of the fit). Otherwise it is flagged
  (``valid`` false) and left out of the verdict; only its value without exclusion is to be reported.
- The change of p_th under (L) is taken from the dump without threshold nearest in time to the middle of the block.
- sigma_pt of a number is taken from the block series of variant ``2_none`` (spectra.detrended on the block
  mid-times), with at least four blocks.
- "The dump reproduces the histogram counts above p_d,e in D_fix within 1 %": compared on the momentum bins wholly
  above p_d,e, per spectra output; the worst output is the result. The comparison with the bin that holds p_d,e
  counted by its log fraction is given beside it.
- Bins of a reduced spectrum that are negative (rounding of x at a bin edge) are not clipped; their number and the
  most negative value are given.
- A block in which a spectra output has no dump ``origin_e`` of its step gets no reduced numbers (only
  ``2_none``), and the verdict is not given from it.

N6 (a). Per species and evaluation time (the spectra output nearest to it): p_0 = 5 p_th from the D_fix spectrum
of that output; E_nt from the whole-box histogram (a bin at its geometric centre, the bin holding p_0 by its log
fraction); N+ and N from all dumps ``edge_<species>`` at or before that step with dt = the dump interval and
exact cuts in x and p per particle; U_far from the dump nearest to the output: all its particles with p >= p_0
beyond the inner edge of the slab S. "R against time at every dump" is the cumulative N and N+ at every dump
divided by E_nt of the evaluation time. Verdict as the README states it.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import spectra as sx  # noqa: E402
from pilot_check import shock_position  # noqa: E402

EXCL = ("1_L", "3_L_lo", "3_L_hi", "4_E", "5_LE")
QUANT = ("eta_nt", "phi", "eta_fixed", "q")
LIMIT = {"eta_nt": 0.30, "phi": 0.30, "eta_fixed": 0.30, "q": 0.15}  # range that is "the error of the definition"
TOL = 0.05  # N6 (a)


def kinetic(p, mass: float, scale: float):
    """Kinetic energy (m_e c^2) of a particle of momentum p (m_i v_pt); scale = mass / (m_i v_pt)."""
    u2 = (np.asarray(p, dtype=np.float64) / scale) ** 2
    return mass * u2 / (np.sqrt(1.0 + u2) + 1.0)


def record(dump, it: int, name: str, mass: float, scale: float, d_i: float, per_weight: float):
    """The particles of one dump output in the units above, or None if the output or the species is absent."""
    try:
        d = dump.get(it, name)
    except KeyError:
        return None
    u = [np.asarray(d[k], dtype=np.float64) for k in ("ux", "uy", "uz")]
    u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2]
    g = np.sqrt(1.0 + u2)
    return {"x": np.asarray(d["x"], dtype=np.float64) / d_i, "p": np.sqrt(u2) * scale, "vx": u[0] / g,
            "ek": mass * u2 / (g + 1.0), "cnt": np.asarray(d["w"], dtype=np.float64) * per_weight,
            "id": np.asarray(d["id"]).astype(np.int64)}


def origin_x0(ids, ppc: int, n_cells, dx: float):
    """(x_0, loaded) of the ids of the first stream: loaded id = (global cell, x slowest) * ppc + slot."""
    nx, ny, nz = (int(n) for n in n_cells)
    ids = np.asarray(ids, dtype=np.int64)
    loaded = ids < nx * ny * nz * ppc
    return ((ids // ppc) // (ny * nz) + 0.5) * dx, loaded


def binned(sp, rec, a: float, b: float, sel=None):
    """Momentum spectrum of the dump particles in a <= x <= b, with the bins and x weights of the histogram."""
    x, p, c = rec["x"], rec["p"], rec["cnt"]
    if sel is not None:
        x, p, c = x[sel], p[sel], c[sel]
    return sp.weights(a, b) @ np.histogram2d(x, p, bins=(sp.xe, sp.pe), weights=c)[0]


def energy_above(c, pe, p0: float, mass: float, scale: float) -> float:
    """Kinetic energy of the counts ``c`` (per bin of ``pe``) with momentum >= p0."""
    lo = np.maximum(pe[:-1], p0)
    part = np.clip(np.log(pe[1:] / lo) / np.log(pe[1:] / pe[:-1]), 0.0, 1.0)
    return float(np.sum(c * part * kinetic(np.sqrt(lo * pe[1:]), mass, scale)))


def value(rec: dict, q: str):
    if q == "q":
        return rec["slope"]["q"] if rec.get("slope") else None
    return rec.get(q)


def masks(ids, x0, loaded, wall_ids, xc) -> dict:
    """The excluded particles of each variant. ``xc`` = X_c (headline, lower, upper) in the units of x0."""
    by_origin = [loaded & (x0 <= c) for c in xc]
    hit = np.isin(ids, wall_ids)
    return {"1_L": by_origin[0], "3_L_lo": by_origin[1], "3_L_hi": by_origin[2], "4_E": hit,
            "5_LE": by_origin[0] | hit}


def n12_block(sp, use, x_h, regs: dict, fetch, wall, x0_of, xc, root: float, fit, p_th_ratio=None) -> dict:
    """The numbers of one block for every region and variant.

    ``use``: the spectra outputs of the block; ``x_h``: the shock position at every output; ``fetch(j)``: the dump
    record of output j or None; ``wall(j)``: the ids of (E) as of output j; ``x0_of(ids)`` -> (x_0 in d_i, loaded).
    ``p_th_ratio``: p_th with (L) over p_th without, from a dump without threshold; if given, p_0 of ``1_L`` is
    recomputed with it and the numbers with the histogram's p_th are kept as ``1_L_before_recompute``.
    """
    n = len(sp.p)
    full = {r: np.zeros(n) for r in regs}
    excl = {r: {v: np.zeros(n) for v in EXCL} for r in regs}
    missing = 0
    for j in use:
        for r, (far, near) in regs.items():
            full[r] += sp.region(j, x_h[j] - far, x_h[j] - near)
        rec = fetch(j)
        if rec is None:
            missing += 1
            continue
        x0, loaded = x0_of(rec["id"])
        m = masks(rec["id"], x0, loaded, wall(j), xc)
        for r, (far, near) in regs.items():
            for v, sel in m.items():
                excl[r][v] += binned(sp, rec, x_h[j] - far, x_h[j] - near, sel)
    out = {"outputs": int(len(use)), "missing_dump": missing, "regions": {}}
    for r in regs:
        ref = sx.numbers(full[r], sp, True, root, fit)
        row = {"2_none": ref}
        for v in EXCL if not missing else ():  # a block with an output without dump gets no reduced numbers
            c = full[r] - excl[r][v]
            row[v] = sx.numbers(c, sp, True, root, fit, ref=ref) | {
                "excluded_particles": float(excl[r][v].sum()), "negative_bins": int((c < 0).sum()),
                "most_negative": float(c.min())}
            if v == "1_L" and p_th_ratio is not None and "p_th" in ref:
                row["1_L_before_recompute"] = row[v]
                row[v] = sx.numbers(c, sp, True, root, fit, ref=ref | {"p_th": ref["p_th"] * p_th_ratio}) | {
                    "p_th_ratio": p_th_ratio}
        out["regions"][r] = row
    return out


def n12_verdict(blocks: list, tmid, region: str, p_d: float, root: float, fit) -> dict:
    """The verdict of N12 in the last block of ``blocks`` for one region."""
    last = blocks[-1]["regions"][region]
    base = last["2_none"]
    lowest = {"eta_nt": base.get("p_0"), "phi": 6.0 * root * 10.0 ** -sx.DEX, "eta_fixed": 3.0 * root,
              "q": fit[0]}
    out = {"p_d": p_d, "numbers": {}}
    for q in QUANT:
        ys = [value(b["regions"][region]["2_none"], q) for b in blocks]
        row = {"2_none": value(base, q), "lowest_momentum": lowest[q],
               "valid": bool(lowest[q] is not None and lowest[q] >= p_d * (1.0 - 1e-9)),
               "sigma_pt": (sx.detrended(np.asarray(tmid), np.asarray(ys, dtype=float))
                            if len(ys) >= 4 and all(y is not None for y in ys) else None)}
        row |= {v: value(last[v], q) if v in last else None for v in EXCL}
        judged = [row[v] for v in ("1_L", "4_E", "5_LE")]
        row["verdict"] = None
        if row["valid"] and row["sigma_pt"] is not None and row["2_none"] is not None and None not in judged:
            ok = all(abs(y - row["2_none"]) <= 2.0 * row["sigma_pt"] for y in judged)
            row["verdict"] = "no wall effect" if ok else "wall effect detected"
        have = [row[v] for v in ("2_none",) + EXCL if row[v] is not None]
        if row["2_none"] is not None and len(have) > 1:
            spread = max(have) - min(have)
            row["range"] = spread if q == "q" else (spread / row["2_none"] if row["2_none"] > 0 else None)
            row["range_above_limit"] = bool(row["range"] is not None and row["range"] > LIMIT[q])
        out["numbers"][q] = row
    verdicts = [r["verdict"] for r in out["numbers"].values() if r["verdict"]]
    out["verdict"] = None if not verdicts else (
        "wall effect detected" if "wall effect detected" in verdicts else "no wall effect")
    out["not_judged"] = [q for q, r in out["numbers"].items() if r["verdict"] is None]
    q1, q2 = out["numbers"]["q"]["1_L"], out["numbers"]["q"]["2_none"]
    if q1 is not None and q2 is not None:
        out["slope_not_decided"] = bool((abs(q1 - 4.0) <= 0.3) != (abs(q2 - 4.0) <= 0.3))
    return out


def dump_against_histogram(sp, outputs, x_h, reg, fetch, p_d: float) -> dict:
    """The counts of the dump against those of the histogram in the region ``reg`` = (far, near), per output."""
    whole = sp.pe[:-1] >= p_d * (1.0 - 1e-9)
    rows = []
    for j in outputs:
        rec = fetch(j)
        if rec is None:
            continue
        a, b = x_h[j] - reg[0], x_h[j] - reg[1]
        h, d = sp.region(j, a, b), binned(sp, rec, a, b)
        rows.append({"t": float(sp.t[j]), "histogram": float(h[whole].sum()), "dump": float(d[whole].sum()),
                     "histogram_with_partial_bin": sx.above(h, sp.pe, p_d), "dump_all": float(d.sum())})
    dev = [abs(r["dump"] / r["histogram"] - 1.0) for r in rows if r["histogram"] > 0]
    return {"outputs": rows, "compared": len(dev), "worst": max(dev) if dev else None,
            "pass": bool(dev and max(dev) <= 0.01)}


def sample_check(sp, rec, x_sh: float, reg, x0, loaded, xc: float, p_d: float) -> dict:
    """The checks of N12 on a dump without energy threshold: the share of the electrons below p_d in the region
    that (L) would leave out, and the change of p_th under (L)."""
    inside = (rec["x"] >= x_sh - reg[0]) & (rec["x"] <= x_sh - reg[1])
    left_out = loaded & (x0 <= xc)
    low = inside & (rec["p"] < p_d)
    n_low = float(rec["cnt"][low].sum())
    c_all = np.histogram(rec["p"][inside], bins=sp.pe, weights=rec["cnt"][inside])[0]
    keep = inside & ~left_out
    c_kept = np.histogram(rec["p"][keep], bins=sp.pe, weights=rec["cnt"][keep])[0]
    out = {"particles": float(rec["cnt"][inside].sum()), "below_p_d": n_low,
           "below_p_d_left_out": float(rec["cnt"][low & left_out].sum()),
           "left_out": float(rec["cnt"][inside & left_out].sum())}
    out["share_below_p_d_left_out"] = out["below_p_d_left_out"] / n_low if n_low > 0 else None
    if c_all.sum() > 0 and c_kept.sum() > 0:
        out["p_th"] = sx.peak_p4f(c_all, sp.p) / np.sqrt(2.0)
        out["p_th_ratio"] = sx.peak_p4f(c_kept, sp.p) / np.sqrt(2.0) / out["p_th"]
        out["p_th_changes"] = bool(abs(out["p_th_ratio"] - 1.0) > 0.02)
    return out


def id_check(rec, x0, loaded, x_sh: float, margin: float, drift: float, face: float) -> dict:
    """Mechanics of the ids (not judged): loaded electrons upstream of the shock against x_0 - v0 t (``drift`` =
    v0 t), and the injected ones against the distance the inflow has covered from the face."""
    up = loaded & (rec["x"] > x_sh + margin)
    res = rec["x"][up] - (x0[up] - drift)
    inj = ~loaded
    return {"loaded_upstream": int(up.sum()), "residual_mean": float(res.mean()) if up.any() else None,
            "residual_rms": float(res.std()) if up.any() else None, "injected": int(inj.sum()),
            "injected_x_min": float(rec["x"][inj].min()) if inj.any() else None, "inflow_front": face - drift}


def n6_sums(recs: list, dt_ls: float, s_in: float, s_out: float, p0: float) -> list:
    """After each dump of ``recs`` (in time order) the cumulative N+ and N: energy (n_plus, n) and number
    (k_plus, k). ``dt_ls`` = the dump interval over the length of S (both in engine units; v_x in c)."""
    acc = np.zeros(4)
    rows = []
    for rec in recs:
        sel = (rec["x"] >= s_in) & (rec["x"] < s_out) & (rec["p"] >= p0)
        f = dt_ls * rec["cnt"][sel] * rec["vx"][sel]
        e, outward = rec["ek"][sel], rec["vx"][sel] > 0
        acc += [np.sum((f * e)[outward]), np.sum(f * e), np.sum(f[outward]), np.sum(f)]
        rows.append(acc.copy())
    return rows


def flux_spectrum(recs: list, dt_ls: float, s_in: float, s_out: float, pe):
    """Time-integrated energy flux through S per momentum bin: (signed, outward only)."""
    out = np.zeros((2, len(pe) - 1))
    for rec in recs:
        sel = (rec["x"] >= s_in) & (rec["x"] < s_out)
        f = dt_ls * rec["cnt"] * rec["vx"] * rec["ek"]
        out[0] += np.histogram(rec["p"][sel], bins=pe, weights=f[sel])[0]
        sel &= rec["vx"] > 0
        out[1] += np.histogram(rec["p"][sel], bins=pe, weights=f[sel])[0]
    return out


def n6_verdict(n_plus: float, n: float, u_far: float, e_nt: float, tol: float = TOL):
    if e_nt <= 0:
        return None
    if n_plus / e_nt <= tol or n / e_nt <= tol:
        return "pass"
    return "fail" if (n - u_far) / e_nt > tol else "not decided"


def n6_species(sp, j: int, c_fix, recs: list, steps, dt_ls: float, s_in: float, s_out: float, mass: float,
               scale: float, p_d: float) -> dict:
    """N6 (a) of one species at the spectra output j. ``c_fix``: the spectrum of D_fix at j; ``recs``, ``steps``:
    all edge dumps of the species in time order."""
    out = {"t": float(sp.t[j]), "particles_d_fix": float(c_fix.sum()), "thresholds": {}}
    use = [k for k, s in enumerate(steps) if s <= sp.step[j]]
    if c_fix.sum() <= 0 or not use:
        return out | {"verdict": None, "dumps": len(use)}
    p_th = sx.peak_p4f(c_fix, sp.p) / np.sqrt(2.0)
    box = sp.cnt[j].sum(axis=0)
    near = min(range(len(steps)), key=lambda k: abs(int(steps[k]) - int(sp.step[j])))
    far = recs[near]
    spec = flux_spectrum([recs[k] for k in use], dt_ls, s_in, s_out, sp.pe)
    out |= {"p_th": p_th, "dumps": len(use), "u_far_dump_steps_from_output": int(steps[near]) - int(sp.step[j]),
            "flux_spectrum_signed": spec[0].tolist(), "flux_spectrum_outward": spec[1].tolist()}
    for label, f in (("p_0", 1.0), ("p_0_half", 0.5), ("p_0_double", 2.0)):
        p0 = 5.0 * p_th * f
        e_nt, k_nt = energy_above(box, sp.pe, p0, mass, scale), sx.above(box, sp.pe, p0)
        cum = n6_sums([recs[k] for k in use], dt_ls, s_in, s_out, p0)
        n_plus, n, k_plus, k = (float(v) for v in cum[-1])
        sel = (far["x"] >= s_in) & (far["p"] >= p0)
        u_far = float(np.sum(far["cnt"][sel] * far["ek"][sel]))
        row = {"p_0": p0, "above_p_d": bool(p0 >= p_d * (1.0 - 1e-9)), "e_nt": e_nt, "particles_nt": k_nt,
               "n_plus": n_plus, "n": n, "u_far": u_far, "verdict": n6_verdict(n_plus, n, u_far, e_nt)}
        if e_nt > 0:
            row |= {"r_upper": n_plus / e_nt, "r": n / e_nt, "r_lower": (n - u_far) / e_nt,
                    "number_fraction": k / k_nt, "number_fraction_outward": k_plus / k_nt,
                    "r_at_dumps": [float(c[1] / e_nt) for c in cum],
                    "r_upper_at_dumps": [float(c[0] / e_nt) for c in cum]}
        out["thresholds"][label] = row
    out["verdict"] = out["thresholds"]["p_0"]["verdict"]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--json")
    ap.add_argument("--smooth", type=float, default=5.0, help="d_i")
    ap.add_argument("--level", type=float, default=2.0)
    ap.add_argument("--blocks", type=float, nargs=3, default=(100.0, 150.0, 10.0), help="t0 t1 length")
    ap.add_argument("--region", type=float, nargs=2, default=(15.0, 2.5), help="M_A d_i behind the shock")
    ap.add_argument("--clear", type=float, default=10.0, help="M_A d_i between the wall and the region")
    ap.add_argument("--fit", type=float, nargs=2, default=(0.6, 2.4), help="m_i v_pt")
    ap.add_argument("--xc", type=float, nargs=3, default=(20.0, 10.0, 40.0), help="X_c in M_A d_i: (1), lower, upper")
    ap.add_argument("--slab", type=float, nargs=2, default=(5.0, 200.0), help="S, d_i inside the absorbing layer")
    ap.add_argument("--edge", type=float, nargs="*", default=(100.0, 125.0, 150.0), help="times of N6 (a)")
    ap.add_argument("--margin", type=float, default=50.0, help="d_i ahead of the shock, for the id check")
    a = ap.parse_args()

    from pypic.analysis.run import open_run

    run = open_run(a.run)
    info = run.meta["planner"]
    d_i, w_ci, v0, m_a = info["d_i"], info["omega_ci"], info["v0"], info["shock"]["mach_r4"]
    if run.info["species"][0]["name"] != "e":
        raise SystemExit("x_0 from the id needs the electrons as the first stream")
    mass = {s["name"]: s["mass"] for s in run.info["species"]}
    scale = {n: mass[n] / (mass["i"] * v0) for n in ("e", "i")}
    root = 1.0 / np.sqrt(mass["i"])
    cfg, dt = run.info["config"], run.info["dt"]
    (_, ny, nz), dx = run.info["n_cells"], run.info["dx"]
    ppc = info["per_column"] // (2 * ny * nz)
    per_weight = info["per_column"] / (2 * ny * nz) / (dx[0] * dx[1] * dx[2])
    dumps, lx = info["dumps"], run.info["L"][0] / d_i
    p_d = dumps["p_d"]
    kinds = run.names("particles")

    def rec_of(series, it, name):
        return record(series, it, name, mass[name], scale[name], d_i, per_weight)

    def x0_of(ids):
        x0, loaded = origin_x0(ids, ppc, run.info["n_cells"], dx[0])
        return x0 / d_i, loaded

    _, x, t, n_i = sx.ion_density(run)
    xs = shock_position(x, n_i, max(1, int(round(a.smooth / (x[1] - x[0])))), a.level)
    sp = {n: sx.Spectrum(run, n) for n in ("e", "i")}
    spe = sp["e"]
    x_h = np.interp(spe.t, t, xs)
    far, near = a.region[0] * m_a, a.region[1] * m_a
    regs = sx.regions(far, near)
    clear = x_h - far >= a.clear * m_a - 1e-9
    xc = tuple(c * m_a for c in a.xc)
    out = {"run": a.run, "commit": run.info.get("git_commit"), "mach_alfven": m_a, "p_d": p_d,
           "region_di_behind_shock": [far, near], "x_c_di": list(xc), "macro_particles_per_weight": per_weight}

    # dumps without threshold: the checks of N12 and the mechanics of the ids
    samples = []
    for name in sorted(n for n in kinds if n.startswith("sample_e")):
        series = run.particles(name)
        for it in series.iterations:
            rec = rec_of(series, it, "e")
            if rec is None:
                continue
            ts = it * dt * w_ci
            x_s = float(np.interp(ts, t, xs))
            x0, loaded = x0_of(rec["id"])
            samples.append({"dump": name, "t": ts, "x_shock": x_s, "particles_in_dump": int(len(rec["x"])),
                            "check": sample_check(spe, rec, x_s, regs["D_fix"], x0, loaded, xc[0], p_d["e"]),
                            "ids": id_check(rec, x0, loaded, x_s, a.margin, v0 * it * dt / d_i, lx)})
    out["samples"] = samples

    # N12
    n12 = {"blocks": []}
    if "origin_e" in kinds and "wall_e" in kinds:
        origin, wall = run.particles("origin_e"), run.particles("wall_e")
        cache, seen = {}, {"k": 0, "ids": np.empty(0, dtype=np.int64)}
        w_its = list(wall.iterations)

        def fetch(j):
            it = int(spe.step[j])
            if it not in cache:
                cache.clear()
                cache[it] = rec_of(origin, it, "e")
            return cache[it]

        def wall_ids(j):  # called with increasing j
            while seen["k"] < len(w_its) and w_its[seen["k"]] <= spe.step[j]:
                rec = rec_of(wall, w_its[seen["k"]], "e")
                if rec is not None:
                    seen["ids"] = np.union1d(seen["ids"], rec["id"])
                seen["k"] += 1
            return seen["ids"]

        edges = np.arange(a.blocks[0], a.blocks[1] + 1e-9, a.blocks[2])
        tmid, used = [], []
        x_lo, x_hi = dumps["origin"]["x_di"]
        for k, (use, skipped) in enumerate(sx.block_outputs(spe.t, edges, clear)):
            if not len(use):
                continue
            mid = 0.5 * (edges[k] + edges[k + 1])
            with_ratio = [s for s in samples if "p_th_ratio" in s["check"]]
            ratio, source = None, None
            if with_ratio:
                s = min(with_ratio, key=lambda s: abs(s["t"] - mid))
                source = {"dump": s["dump"], "t": s["t"], "p_th_ratio": s["check"]["p_th_ratio"]}
                ratio = s["check"]["p_th_ratio"] if s["check"]["p_th_changes"] else None
            b = n12_block(spe, use, x_h, regs, fetch, wall_ids, x0_of, xc, root, a.fit, ratio)
            outside = sum(1 for j in use if x_h[j] - far < x_lo - 1e-9 or x_h[j] - near > x_hi + 1e-9)
            n12["blocks"].append({"t": [float(edges[k]), float(edges[k + 1])], "outputs_not_clear": skipped,
                                  "outputs_with_region_outside_dump": outside, "p_th_check": source,
                                  "wall_ids": int(len(seen["ids"]))} | b)
            tmid.append(mid)
            used += list(use)
        if n12["blocks"]:
            n12["verdict"] = {r: n12_verdict(n12["blocks"], tmid, r, p_d["e"], root, a.fit) for r in regs}
            n12["dump_against_histogram"] = dump_against_histogram(spe, used, x_h, regs["D_fix"], fetch, p_d["e"])
    out["n12"] = n12

    # N6 (a)
    layer = (run.info["L"][0] - cfg["npml"] * dx[0]) / d_i
    s_in, s_out = layer - a.slab[1], layer - a.slab[0]
    every = dumps["edge"]["every"]
    dt_ls = every * dt / ((s_out - s_in) * d_i)
    sc = run.scalars
    n6 = {"slab_di": [s_in, s_out], "dump_zone_di": dumps["edge"]["x_di"], "dump_interval_steps": every,
          "dump_interval_over_slab": dt_ls, "times": []}
    edge = {}
    for name in ("e", "i"):
        if f"edge_{name}" in kinds:
            series = run.particles(f"edge_{name}")
            steps = [it for it in series.iterations]
            recs = [rec_of(series, it, name) for it in steps]
            keep = [k for k, r in enumerate(recs) if r is not None]
            edge[name] = ([recs[k] for k in keep], [steps[k] for k in keep])
            n6[f"dumps_{name}"] = len(keep)
            n6[f"interval_regular_{name}"] = bool(np.all(np.diff(edge[name][1]) == every))
    half = 0.5 * (spe.t[1] - spe.t[0]) if len(spe.t) > 1 else 0.0
    for te in a.edge:
        if te > spe.t[-1] + half or len(edge) < 2:
            continue
        j = int(np.argmin(np.abs(spe.t - te)))
        row = {"t": float(spe.t[j]), "x_shock": float(x_h[j]), "region_clear_of_wall": bool(clear[j]),
               "engine_lost_kinetic": float(np.interp(spe.step[j], sc.step, sc["lost_kinetic"])) * per_weight}
        for name in ("e", "i"):
            c_fix = sp[name].region(j, x_h[j] - far, x_h[j] - near)
            row[name] = n6_species(sp[name], j, c_fix, *edge[name], dt_ls, s_in, s_out, mass[name], scale[name],
                                   p_d[name])
        verdicts = [row[n]["verdict"] for n in ("e", "i")]
        row["verdict"] = (None if None in verdicts else "pass" if verdicts == ["pass", "pass"]
                          else "fail" if "fail" in verdicts else "not decided")
        n6["times"].append(row)
    not_pass = [r["t"] for r in n6["times"] if r["verdict"] != "pass"]
    n6["first_time_not_pass"] = not_pass[0] if not_pass else None
    n6["pass"] = bool(n6["times"] and not not_pass)
    out["n6a"] = n6

    text = json.dumps(out, indent=1)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")
    else:
        print(text)


if __name__ == "__main__":
    main()
