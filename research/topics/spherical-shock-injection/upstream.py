"""N3c of exploratory-protocol.md Sect. 5: the power of the field fluctuations upstream of the front in a run
against its null run at the same place and time, from the meridional slices.

The protocol does not fix the scale of the fluctuations, so three measures are taken and the note is carried
if any of them exceeds the ratio 2: the fields with their box mean over 1 d_i removed; the fields about their
mean along the front (cells of 5 degrees by one grid step of the distance to the front, per quadrant); and the
fields about their initial values, which holds the foot of the front too. A box mean wider than the distance
to the ramp is not used: it brings the ramp itself into the difference.
"""
from __future__ import annotations

import json
from pathlib import Path

import census
import numpy as np
from scipy.ndimage import uniform_filter

SHELLS = ((1.0, 4.0), (4.0, 12.0))  # distance ahead of the density front, in d_i
SCALE_NAMES = ("below 1 d_i", "about the mean along the front", "about the initial field")
DTHETA = np.radians(5.0)
DOMAINS = ("detected front", "valid domain")
PLANES = (("slice_y", 0), ("slice_x", 1))
RATIO_MAX = 2.0


def highpass(a, n: int):
    a = np.asarray(a, float)
    return a - uniform_filter(a, size=n, mode="nearest") if n > 1 else a


def along_front(a, delta, cell, dx: float, d_max: float):
    """a minus its mean over the points of the same cell along the front (cell: an integer per point that
    names the quadrant and the angle bin) and the same grid step of delta; NaN outside 0 <= delta < d_max."""
    ok = np.isfinite(delta) & (delta >= 0.0) & (delta < d_max)
    n_d = int(np.ceil(d_max / dx))
    key = (cell[ok] * n_d + np.minimum((delta[ok] / dx).astype(int), n_d - 1)).astype(np.int64)
    num = np.bincount(key, minlength=int(key.max()) + 1 if key.size else 1)
    mean = np.bincount(key, weights=a[ok], minlength=len(num)) / np.maximum(num, 1)
    out = np.full(np.shape(a), np.nan)
    out[ok] = a[ok] - mean[key]
    return out


def shell_sums(a_run, a_null, delta, flag, cell, d_i: float, dx: float):
    """Sums of the squared fluctuation [run or null, measure, shell, domain] and the number of points
    [shell, domain]. delta: distance ahead of the front on the grid of the field; flag: that of the front
    at the angle of the point (1 detected, 2 or more valid); cell: see along_front."""
    sums, counts = np.zeros((2, len(SCALE_NAMES), len(SHELLS), 2)), np.zeros((len(SHELLS), 2))
    masks = [[(delta >= lo * d_i) & (delta < hi * d_i) & (flag >= f) for f in (1, 2)] for lo, hi in SHELLS]
    for s, shell in enumerate(masks):
        counts[s] = [m.sum() for m in shell]
    d_max = SHELLS[-1][1] * d_i
    for k, a in enumerate((a_run, a_null)):
        a = np.asarray(a, float)
        for j, res in enumerate((highpass(a, int(round(d_i / dx))), along_front(a, delta, cell, dx, d_max), a)):
            sq = res ** 2
            for s, shell in enumerate(masks):
                sums[k, j, s] = [sq[m].sum() for m in shell]
    return sums, counts


def front_cell(aa, zz, theta):
    """Integer per point: quadrant of the plane and angle bin of DTHETA."""
    quad = (aa < 0).astype(int) + 2 * (zz < 0).astype(int)
    return quad * 64 + np.minimum((theta / DTHETA).astype(int), 63)


def cmd_n3c(args) -> int:
    run_dir, null_dir = Path(args.run), Path(args.null)
    run, null = census._open(run_dir), census._open(null_dir)
    model = census.FrontModel.from_run(run, run_dir)
    d_i, b0, dx = model.d_i, model.b0, model.dx[0]
    r_outer = 0.5 * float(run.L[2]) - 8.0
    sers = [(run.series(name), null.series(name), ax) for name, ax in PLANES if name in run.names()]
    t = np.asarray(sers[0][0].t, float)
    if not np.allclose(t, np.asarray(sers[0][1].t, float)):
        print("the outputs of the run and of the null run are not at the same times")
        return 1
    # [output, field (E, B), run or null, measure, shell, domain]
    sums = np.zeros((len(t), 2, 2, len(SCALE_NAMES), len(SHELLS), 2))
    counts = np.zeros((len(t), len(SHELLS), 2))
    for i, ti in enumerate(t):
        if not model.t_range[0] <= ti <= model.t_range[1]:
            continue
        for ser, ser0, ax in sers:
            for f, name in enumerate("EB"):
                for comp in "xyz":
                    key = name + comp
                    a, z = ser.plane_axes(key)
                    aa, zz = np.meshgrid(np.asarray(a) - model.centre[ax], np.asarray(z) - model.centre[2],
                                         indexing="ij")
                    r, theta = np.hypot(aa, zz), np.arctan2(np.abs(aa), np.abs(zz))
                    delta = np.where(r <= r_outer, r - model.front(theta, np.full(r.shape, ti))[0], np.nan)
                    flag = model.flag(theta, ti)
                    off = b0 if key == "Bz" else 0.0
                    s, c = shell_sums(ser.read(key, i).astype(float) - off, ser0.read(key, i).astype(float) - off,
                                      delta, flag, front_cell(aa, zz, theta), d_i, dx)
                    sums[i, f] += s
                    if key == "Bz":
                        counts[i] += c
    out = summarize(t, sums, counts, b0)
    out.update({"run": str(run_dir), "null": str(null_dir), "commit": run.info.get("git_commit"),
                "commit_null": null.info.get("git_commit")})
    (run_dir / "analysis" / "n3c.json").write_text(json.dumps(out, indent=1))
    report(out)
    print(f"wrote {run_dir / 'analysis' / 'n3c.json'}")
    return 0


def summarize(t, sums, counts, b0: float) -> dict:
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sums[:, :, 0] / sums[:, :, 1]  # [output, field, scale, shell, domain]
        rms_null = np.sqrt(sums[:, :, 1] / counts[:, None, None]) / b0
    rows, worst = [], dict.fromkeys(DOMAINS, 0.0)
    for s, (lo, hi) in enumerate(SHELLS):
        for d, dom in enumerate(DOMAINS):
            use = counts[:, s, d] >= 1000
            for f, field in enumerate("EB"):
                for j, scale in enumerate(SCALE_NAMES):
                    r = ratio[use, f, j, s, d]
                    row = {"shell_d_i": [lo, hi], "domain": dom, "field": field, "scale": scale,
                           "outputs": int(use.sum())}
                    if use.any():
                        k = int(np.argmax(r))
                        row.update({"median": float(np.median(r)), "max": float(r[k]), "t_max": float(t[use][k]),
                                    "above 2": int((r > RATIO_MAX).sum()),
                                    "rms_null_over_B0": float(np.median(rms_null[use, f, j, s, d]))})
                        worst[dom] = max(worst[dom], row["max"])
                    rows.append(row)
    return {"shells_d_i": SHELLS, "scales": SCALE_NAMES, "domains": DOMAINS, "times": np.asarray(t).tolist(),
            "points": counts.tolist(), "ratio": np.where(np.isfinite(ratio), ratio, -1.0).tolist(), "rows": rows,
            "largest ratio": worst, "note on N3 and N5": {d: bool(w > RATIO_MAX) for d, w in worst.items()}}


def report(out) -> None:
    print(f"N3c of {out.get('run', 'synthetic')}: fluctuation power ahead of the front, run / null run")
    for r in out["rows"]:
        if r["outputs"]:
            print(f"  {r['shell_d_i'][0]:4.0f}-{r['shell_d_i'][1]:2.0f} d_i  {r['domain']:14s} {r['field']}  "
                  f"{r['scale']:23s} median {r['median']:8.2f}  max {r['max']:9.2f} (t = {r['t_max']:5.1f})  "
                  f"above 2 in {r['above 2']:2d} of {r['outputs']:2d} outputs; null rms {r['rms_null_over_B0']:.4f} B0")
    for dom, w in out["largest ratio"].items():
        print(f"{dom}: largest ratio {w:.2f}: " + ("the floor of the null run does not stand for the upstream of "
              "this run; N3 and N5 carry the note" if out["note on N3 and N5"][dom] else "no note on N3 and N5"))


# ------------------------------------------------------------------ check on synthetic data
def selfcheck(check) -> None:
    rng = np.random.default_rng(11)
    n, d_i, dx = 384, 4.0, 0.5
    x = (np.arange(n) + 0.5) * dx - 0.5 * n * dx
    xx, zz = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(xx, zz)
    delta = r - 30.0
    theta = np.arctan2(np.abs(xx), np.abs(zz))
    flag = np.where(theta > np.radians(45.0), 2, 1)
    null = rng.normal(0.0, 1.0, (n, n))
    # the run: the same noise level, three times the power in the near shell, and a smooth ramp everywhere
    extra = rng.normal(0.0, np.sqrt(2.0), (n, n)) * ((delta >= 4.0) & (delta < 16.0))
    a_run = rng.normal(0.0, 1.0, (n, n)) + extra + 0.2 * xx
    sums, counts = shell_sums(a_run, null, delta, flag, front_cell(xx, zz, theta), d_i, dx)
    ratio = sums[0] / sums[1]  # [scale, shell, domain]
    area = np.pi * (46.0 ** 2 - 34.0 ** 2) / dx ** 2
    ok = np.all(np.abs(ratio[:2, 0] - 3.0) < 0.2) and np.all(np.abs(ratio[:2, 1] - 1.0) < 0.1) \
        and np.all(ratio[2] > 5.0) and abs(counts[0, 0] / area - 1.0) < 0.02 \
        and abs(counts[0, 1] / counts[0, 0] - 0.5) < 0.02
    check("N3c, fluctuation power ahead of the front", ok,
          f"ratio {np.round(ratio[:2, 0, 0], 2).tolist()} in the shell with three times the power, "
          f"{np.round(ratio[:2, 1, 0], 2).tolist()} outside it; the ramp shows only about the initial field "
          f"({ratio[2, 1, 0]:.0f}); {counts[0, 0]:.0f} points in the shell")
