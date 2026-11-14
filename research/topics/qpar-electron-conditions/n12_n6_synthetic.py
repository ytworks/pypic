"""The estimators of n12_n6.py on synthetic particles (a check of the code, not a null test of the README).

    python research/topics/qpar-electron-conditions/n12_n6_synthetic.py --json <file> [--seed 11]

CPU only (NumPy). Geometry, thresholds and cadences are the ones of V1 (configs/v1.yaml); the momentum bins are
240 logarithmic bins from 1.333e-4 to 100 m_i v_pt for both species.

1. ``n12``: a fixed set of electrons with known origin (id), momentum and step of first wall contact; the positions
   are drawn anew for each of 20 spectra outputs in t = 100-150. The histogram holds all of them, the dump the ones
   with p >= p_d,e. A "wall" population (loaded at x_0 <= 800 d_i, harder spectrum, partly below p_d,e, 70 % of it
   in a wall dump at some step) is added in the case ``wall``; the case ``null`` has none.
   Reference: the numbers of the histogram of the electrons that the exclusion keeps, built from the known origin
   and contact step of each electron (not from the dump), with the particle number and p_th of the full histogram.
   Fixed before the first run:
   - every number (eta_nt, phi, eta_fixed, q) of every block, region and variant equals its reference within
     1e-9 (relative; absolute for q);
   - x_0 and "loaded" from the id equal the known ones exactly;
   - the dump reproduces the histogram above p_d,e within 1e-9;
   - the verdict in D_fix is "no wall effect" with all four numbers judged in ``null`` and "wall effect detected"
     in ``wall``.
   Given without criterion: the share of the electrons below p_d,e that (L) leaves out, from every 101st id and
   from all electrons.
2. ``n6``: per species, particles are emitted at a constant mean rate 10 d_i before the slab S and move outward at
   constant v_x; 30 % turn back at a random place between the inner edge of S and the face and leave S again, the
   others are removed at the face. Dumps of the zone every 25,000 steps. The path of every particle is known, so
   the time integrals behind N+ and N, and the energy removed at the face E_esc, are known exactly.
   Fixed before the first run:
   - N and N+ of n12_n6.n6_species equal the exact time integrals within 1 %;
   - E_nt from the histogram equals the summed energy of the particles above p_0 within 1 %;
   - N - U_far <= E_esc <= N <= N+ holds for the estimates with a slack of 1 % of N;
   - n6_verdict returns pass / pass / not decided / fail on four constructed cases.
   Given without criterion: N and N+ against the integrals to half a dump interval after the output (the sum
   gives the dump at the output its full weight).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import n12_n6 as nn  # noqa: E402
import spectra as sx  # noqa: E402

NX, NY, NZ, PPC, DX = 270000, 2, 2, 250, 0.1  # V1: cells, particles per cell and species, dx in d_e
D_I, M_A, ROOT, V0 = 10.0, 20.0, 0.1, 0.1
DT, W_CI = 0.045, 1.0 / 1500.0
P_D = {"e": 0.1, "i": 1.5}
FIT = (0.6, 2.4)
N_IDS = NX * NY * NZ * PPC
PER = NY * NZ * PPC  # ids per x cell


def spectrum(xe, t, step):
    """A spectra.Spectrum with V1's electron bins of p and no counts."""
    sp = sx.Spectrum.__new__(sx.Spectrum)
    sp.xe, sp.pe = np.asarray(xe, dtype=float), np.logspace(np.log10(1.333e-4), 2.0, 241)
    sp.p = np.sqrt(sp.pe[1:] * sp.pe[:-1])
    sp.dlnp = float(np.log(sp.pe[1] / sp.pe[0]))
    sp.t, sp.step = np.asarray(t, dtype=float), np.asarray(step)
    sp.cnt = np.zeros((len(sp.t), len(sp.xe) - 1, len(sp.p)))
    return sp


def thermal(rng, n: int, p_th: float):
    """|p| of a Maxwellian f ~ exp(-p^2 / p_th^2)."""
    return (p_th / np.sqrt(2.0)) * np.sqrt((rng.standard_normal((3, n)) ** 2).sum(axis=0))


def power(rng, n: int, a: float, b: float, s: float):
    """p with dN/dp ~ p^-s between a and b (s != 1)."""
    e = 1.0 - s
    return (a ** e + rng.random(n) * (b ** e - a ** e)) ** (1.0 / e)


def loaded_ids(rng, n: int, cell_lo: int, cell_hi: int):
    """Distinct ids of loaded electrons with the x cell in [cell_lo, cell_hi), and their x_0 (d_i)."""
    k = (rng.choice((cell_hi - cell_lo) * PER, n, replace=False) + cell_lo * PER).astype(np.int64)
    return k, ((k // PER) + 0.5) * DX / D_I


def deviation(a, b, absolute: bool) -> float:
    if a is None and b is None:
        return 0.0
    if a is None or b is None:
        return float("inf")
    return abs(a - b) if absolute or b == 0 else abs(a / b - 1.0)


def n12_case(seed: int, wall_population: bool) -> dict:
    rng = np.random.default_rng([seed, 1, int(wall_population)])
    steps = 3_400_000 + 80_000 * np.arange(20)
    t = steps * DT * W_CI
    sp = spectrum(np.arange(0.0, 900.1, 5.0), t, steps)
    x_h = 5.06 * t  # the shock of V1 (0.03374 c), d_i
    xc = (20.0 * M_A, 10.0 * M_A, 40.0 * M_A)
    regs = sx.regions(15.0 * M_A, 2.5 * M_A)
    # the tails are small enough that the peak of p^4 f stays at the thermal population
    n_th, n_tail, n_wall = 1_500_000, 100_000, (30_000 if wall_population else 0)
    n_main, n = n_th + n_tail, n_th + n_tail + n_wall
    p = np.concatenate([thermal(rng, n_th, 0.03), power(rng, n_tail, 0.12, 20.0, 2.0),
                        power(rng, n_wall, 0.05, 20.0, 1.5)])
    ids, x0, loaded = np.empty(n, dtype=np.int64), np.full(n, np.inf), np.zeros(n, dtype=bool)
    m = np.zeros(n, dtype=bool)
    m[:n_main] = rng.permutation(n_main) < n_main // 2  # main population: half loaded beyond 800 d_i ...
    ids[m], x0[m] = loaded_ids(rng, int(m.sum()), 80000, NX)
    loaded[m] = True
    inj = ~m
    inj[n_main:] = False  # ... half injected
    ids[inj] = N_IDS + rng.choice(10 ** 12, int(inj.sum()), replace=False)
    band = rng.choice(3, n_wall, p=(0.4, 0.3, 0.3))  # wall population: x_0 <= 200, 200-400, 400-800 d_i
    for b, (lo, hi) in enumerate(((0, 20000), (20000, 40000), (40000, 80000))):
        sel = np.zeros(n, dtype=bool)
        sel[n_main:] = band == b
        ids[sel], x0[sel] = loaded_ids(rng, int(sel.sum()), lo, hi)
        loaded[sel] = True
    seen = np.full(n, -1, dtype=np.int64)  # step of the first wall dump that holds the id
    if wall_population:
        tagged = np.zeros(n, dtype=bool)
        tagged[n_main:] = rng.random(n_wall) < 0.7
        tagged[n_th:n_main] = rng.random(n_tail) < 0.01
        early = rng.random(n) < 0.5
        w_steps = np.arange(0, steps[-1] + 1, 20000)
        for part, pool in ((tagged & early, w_steps < steps[0]), (tagged & ~early, w_steps >= steps[0])):
            seen[part] = rng.choice(w_steps[pool], int(part.sum()))
    in_dump = p >= P_D["e"]
    one = np.ones(n)
    cache = {}

    def x_of(j):
        if j not in cache:
            cache.clear()
            cache[j] = np.random.default_rng([seed, 2, int(wall_population), int(j)]).random(n) * 900.0
        return cache[j]

    for j in range(len(steps)):
        sp.cnt[j] = np.histogram2d(x_of(j), p, bins=(sp.xe, sp.pe))[0]

    def fetch(j):
        return {"x": x_of(j)[in_dump], "p": p[in_dump], "cnt": one[in_dump], "id": ids[in_dump]}

    def wall(j):
        return np.unique(ids[(seen >= 0) & (seen <= steps[j])])

    def x0_of(i):
        x, is_loaded = nn.origin_x0(i, PPC, (NX, NY, NZ), DX)
        return x / D_I, is_loaded

    x0_id, loaded_id = x0_of(ids)
    out = {"particles": n, "in_dump": int(in_dump.sum()),
           "origin_exact": bool(np.array_equal(loaded_id, loaded) and np.array_equal(x0_id[loaded], x0[loaded]))}

    edges = np.arange(100.0, 150.1, 10.0)
    blocks, tmid, used, worst = [], [], [], {q: 0.0 for q in nn.QUANT}
    for k, (use, _) in enumerate(sx.block_outputs(t, edges, np.ones(len(t), dtype=bool))):
        blk = nn.n12_block(sp, use, x_h, regs, fetch, wall, x0_of, xc, ROOT, FIT)
        exact = {r: {v: np.zeros(len(sp.p)) for v in nn.EXCL} for r in regs}
        for j in use:
            x, e = x_of(j), (seen >= 0) & (seen <= steps[j])
            left = [loaded & (x0 <= c) for c in xc]
            drop = {"1_L": left[0], "3_L_lo": left[1], "3_L_hi": left[2], "4_E": e, "5_LE": left[0] | e}
            for v, d in drop.items():
                h = np.histogram2d(x[~d], p[~d], bins=(sp.xe, sp.pe))[0]
                for r, (far, near) in regs.items():
                    exact[r][v] += sp.weights(x_h[j] - far, x_h[j] - near) @ h
        for r in regs:
            ref = blk["regions"][r]["2_none"]
            for v in nn.EXCL:
                want = sx.numbers(exact[r][v], sp, True, ROOT, FIT, ref=ref)
                for q in nn.QUANT:
                    worst[q] = max(worst[q], deviation(nn.value(blk["regions"][r][v], q), nn.value(want, q), q == "q"))
        blocks.append(blk)
        tmid.append(0.5 * (edges[k] + edges[k + 1]))
        used += list(use)
    verdict = nn.n12_verdict(blocks, tmid, "D_fix", P_D["e"], ROOT, FIT)
    against = nn.dump_against_histogram(sp, used, x_h, regs["D_fix"], fetch, P_D["e"])
    # the dump without threshold: every 101st id, at the last output
    j = len(steps) - 1
    s = ids % 101 == 0
    x = x_of(j)
    chk = nn.sample_check(sp, {"x": x[s], "p": p[s], "cnt": one[s], "id": ids[s]}, x_h[j], regs["D_fix"], x0[s],
                          loaded[s], xc[0], P_D["e"])
    low = (x >= x_h[j] - regs["D_fix"][0]) & (x <= x_h[j] - regs["D_fix"][1]) & ~in_dump
    out |= {"largest_deviation": worst, "verdict": verdict["verdict"], "not_judged": verdict["not_judged"],
            "numbers_d_fix": {q: {k: r.get(k) for k in ("2_none", "1_L", "4_E", "5_LE", "sigma_pt", "verdict", "range")}
                              for q, r in verdict["numbers"].items()},
            "dump_against_histogram_worst": against["worst"],
            "share_below_p_d_left_out": {"every_101st_id": chk["share_below_p_d_left_out"],
                                         "all": float((low & loaded & (x0 <= xc[0])).sum() / low.sum())}}
    want = "wall effect detected" if wall_population else "no wall effect"
    out["pass"] = bool(max(worst.values()) <= 1e-9 and out["origin_exact"] and against["worst"] is not None
                       and against["worst"] <= 1e-9 and verdict["verdict"] == want
                       and (wall_population or not verdict["not_judged"]))
    return out


def n6_case(seed: int, name: str) -> dict:
    rng = np.random.default_rng([seed, 3, "ei".index(name)])
    mass = {"e": 1.0, "i": 100.0}[name]
    scale = mass / (100.0 * V0)
    layer, face = 2698.08, 2700.0  # inner edge of the absorbing layer and the face, d_i
    s_in, s_out = layer - 200.0, layer - 5.0
    every, n_steps = 25000, 5_000_000
    dt_ls = every * DT / ((s_out - s_in) * D_I)
    p_th, lo, hi = {"e": (0.03, 0.08, 20.0), "i": (0.42, 1.2, 30.0)}[name]
    n_emit, n_tail, n_bulk = 300_000, 100_000, 1_000_000
    p = power(rng, n_emit, lo, hi, 2.0)
    u = p / scale
    vx = u / np.sqrt(1.0 + u * u) * rng.uniform(0.2, 1.0, n_emit)
    ek = nn.kinetic(p, mass, scale)
    t0 = rng.random(n_emit) * n_steps * DT
    x_e = s_in - 10.0
    turn = np.where(rng.random(n_emit) < 0.3, rng.uniform(s_in, face, n_emit), face)
    reflected = turn < face
    static = np.concatenate([thermal(rng, n_bulk, p_th), power(rng, n_tail, lo, hi, 2.0)])
    one = np.ones(n_emit)

    def state(step):
        """Position, sign of v_x, "in the box", and the outward and signed path in S over L_S, at a step."""
        d = vx * (step * DT - t0) / D_I
        out = np.minimum(x_e + d, turn)
        back = np.where(reflected, np.maximum(x_e + d - turn, 0.0), 0.0)
        alive = (d >= 0.0) & (reflected | (x_e + d < face))
        path_out = np.clip(out, s_in, s_out) - s_in
        path_in = np.where(back > 0.0, np.clip(turn, s_in, s_out) - np.clip(turn - back, s_in, s_out), 0.0)
        return (out - back, np.where(back > 0.0, -1.0, 1.0), alive, path_out / (s_out - s_in),
                (path_out - path_in) / (s_out - s_in))

    steps = np.arange(0, n_steps + 1, every)
    recs = []
    for s in steps:
        x, sign, alive, _, _ = state(s)
        m = alive & (x >= s_in) & (p >= P_D[name])
        recs.append({"x": x[m], "p": p[m], "vx": sign[m] * vx[m], "ek": ek[m], "cnt": one[m]})
    outputs = np.array([3_300_000, 4_200_000, 5_000_000])
    sp = spectrum([0.0, face], outputs * DT * W_CI, outputs)
    c_fix = np.histogram(static, bins=sp.pe)[0].astype(float)
    rows, worst = [], {"n": 0.0, "n_plus": 0.0, "e_nt": 0.0, "order": 0.0}
    for j, s in enumerate(outputs):
        x, sign, alive, f_out, f_net = state(s)
        gone = (vx * (s * DT - t0) / D_I + x_e >= face) & ~reflected
        late = state(s + every // 2)
        sp.cnt[j, 0] = c_fix + np.histogram(p[alive], bins=sp.pe)[0]
        got = nn.n6_species(sp, j, c_fix, recs, steps, dt_ls, s_in, s_out, mass, scale, P_D[name])
        row = {"t": float(sp.t[j]), "p_th": got["p_th"], "thresholds": {}}
        for label in ("p_0", "p_0_double"):
            g = got["thresholds"][label]
            a = p >= g["p_0"]
            exact = {"n": float(ek[a] @ f_net[a]), "n_plus": float(ek[a] @ f_out[a]),
                     "e_esc": float(ek[a & gone].sum()),
                     "e_nt": float(ek[a & alive].sum() + nn.kinetic(static[static >= g["p_0"]], mass, scale).sum())}
            dev = {k: g[k] / exact[k] - 1.0 for k in ("n", "n_plus", "e_nt")}
            order = max(g["n"] - g["u_far"] - exact["e_esc"], exact["e_esc"] - g["n"], g["n"] - g["n_plus"]) / g["n"]
            row["thresholds"][label] = {
                "p_0": g["p_0"], "estimate": {k: g[k] for k in ("n", "n_plus", "u_far", "e_nt")}, "exact": exact,
                "deviation": dev, "order_violation_over_n": float(order),
                "deviation_from_integral_to_half_interval_later": {
                    "n": g["n"] / float(ek[a] @ late[4][a]) - 1.0,
                    "n_plus": g["n_plus"] / float(ek[a] @ late[3][a]) - 1.0}}
            for k in ("n", "n_plus", "e_nt"):
                worst[k] = max(worst[k], abs(dev[k]))
            worst["order"] = max(worst["order"], float(order))
        rows.append(row)
    return {"emitted": n_emit, "reflected": int(reflected.sum()), "dumps": len(steps), "dump_interval_over_slab": dt_ls,
            "largest_deviation": worst, "outputs": rows,
            "pass": bool(max(worst["n"], worst["n_plus"], worst["e_nt"]) <= 0.01 and worst["order"] <= 0.01)}


def verdict_cases() -> dict:
    """(n_plus, n, u_far, e_nt) -> the expected verdict of N6 (a) at a tolerance of 5 %."""
    cases = {"pass": (0.04, 0.03, 0.01, 1.0), "pass_signed": (0.08, 0.04, 0.01, 1.0),
             "not decided": (0.09, 0.07, 0.03, 1.0), "fail": (0.12, 0.10, 0.03, 1.0)}
    got = {k: nn.n6_verdict(*v) for k, v in cases.items()}
    return {"got": got, "pass": all(got[k] == k.split("_")[0] for k in cases)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    def git(*c):
        return subprocess.run(("git", "-C", here) + c, capture_output=True, text=True).stdout.strip()

    t0 = time.time()
    out = {"commit": git("rev-parse", "HEAD"), "clean_tree": git("status", "--short") == "", "seed": a.seed,
           "n12": {"null": n12_case(a.seed, False), "wall": n12_case(a.seed, True)},
           "n6": {n: n6_case(a.seed, n) for n in ("e", "i")}, "n6_verdict": verdict_cases()}
    out["n12_pass"] = bool(all(c["pass"] for c in out["n12"].values()))
    out["n6_pass"] = bool(all(c["pass"] for c in out["n6"].values()) and out["n6_verdict"]["pass"])
    out["seconds"] = time.time() - t0
    text = json.dumps(out, indent=1)
    print(text)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
