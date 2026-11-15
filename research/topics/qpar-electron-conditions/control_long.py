"""K1 and K2 of the long flow-through controls (notebook, 2026-11-15, "controls to 150"), fixed before the runs.

    python control_long.py <run dir> [--other <run dir>] [--json out.json]

- ``k1``: in the zone (``--zone``, 30-37 d_i), the means over blocks of 5 / Omega_ci from t = 5 to the end of
  <B_y> / B0y, of <n_i> / n0 and of Q = (<B_y> / B0y) / (<n_i> / n0). n0: the first output between 10 % and 90 %
  of the box (as in ``spectra.py``). Criterion: every block of B_y and of Q within 1 +- 0.01. ``band`` names the
  outcome by the largest deviation of either: "pass" (<= 1 %), "owner" (1-5 %), "no V1" (> 5 %).
- ``k2``: the means of <B_y> / B0y over blocks of 6 / Omega_ci from t = 4 to 40 and the block that holds the
  largest one. The two earlier controls both had it at t = 22-28.
- Reported, no threshold: ``slope`` (straight line over t = 5 to the end; the standard error of independent
  outputs, and the factor by which the advection of the flow across the zone raises it: correlation
  1 - lag / crossing time between outputs); ``lag`` (the correlation between the residual series of the zone and
  of a far zone, ``--far``, 2-9 d_i, at lag 0 and at the advection time between them, after boxcar means over
  1, 15 and 30 outputs: a part that the flow carries peaks at the advection lag, a part that is uniform in the
  box at lag 0).
- ``--other``: a run with the same outputs (another layer depth or seed): this run's k1 block means of B_y minus
  the other's, and the correlation of the two residual series.
"""
import argparse
import json

import numpy as np

from pypic.analysis.run import open_run


def load(path):
    run = open_run(path)
    info = run.meta["planner"]
    prof = run.series("profile")
    x = prof.x("By") / info["d_i"]
    xn = prof.x("n:i") / info["d_i"]
    n0 = float(prof["n:i"][0][(xn > 0.1 * xn[-1]) & (xn < 0.9 * xn[-1])].mean())
    return {"t": prof.t * info["omega_ci"], "x": x, "xn": xn, "by": prof["By"] / info["by"],
            "n": prof["n:i"] / n0, "info": info}


def zone_mean(d, key, zone):
    x = d["x"] if key == "by" else d["xn"]
    return d[key][:, (x >= zone[0]) & (x <= zone[1])].mean(axis=1)


def blocks(t, y, start, stop, length):
    rows = []
    for a in np.arange(start, stop - 1e-6, length):
        m = (t >= a - 1e-6) & (t < a + length - 1e-6)
        if m.any():
            rows.append((float(a), float(a + length), int(m.sum()), float(y[m].mean())))
    return rows


def residual(t, y, window):
    sel = (t >= window[0] - 1e-6) & (t <= window[1] + 1e-6)
    tt, yy = t[sel], y[sel]
    fit, cov = np.polyfit(tt, yy, 1, cov=True)
    return tt, yy - np.polyval(fit, tt), float(fit[0]), float(np.sqrt(cov[0, 0]))


def corr(a, b, lag):
    return float(np.corrcoef(a[:-lag], b[lag:])[0, 1]) if lag else float(np.corrcoef(a, b)[0, 1])


def k1(d, zone, limit=0.01, owner=0.05):
    t, end = d["t"], float(d["t"][-1])
    by, n = zone_mean(d, "by", zone), zone_mean(d, "n", zone)
    rows = [{"t": [a, b], "outputs": m, "By": v, "n": w[3], "Q": v / w[3]}
            for (a, b, m, v), w in zip(blocks(t, by, 5.0, end, 5.0), blocks(t, n, 5.0, end, 5.0))]
    dev = {k: float(max(abs(r[k] - 1.0) for r in rows)) for k in ("By", "Q")}
    worst = max(dev.values())
    return {"zone_di": list(zone), "block": 5.0, "blocks": len(rows), "rows": rows, "largest_deviation": dev,
            "limit": limit, "pass": bool(worst <= limit),
            "band": "pass" if worst <= limit else "owner" if worst <= owner else "no V1"}


def k2(d, zone):
    rows = blocks(d["t"], zone_mean(d, "by", zone), 4.0, 40.0, 6.0)
    top = max(rows, key=lambda r: r[3])
    return {"blocks": [{"t": [a, b], "outputs": m, "By_minus_1": v - 1.0} for a, b, m, v in rows],
            "largest": [top[0], top[1]], "largest_is_22_28": bool(abs(top[0] - 22.0) < 1e-6)}


def crossing(d, length_di):
    """Time (1 / Omega_ci) in which the flow covers ``length_di``."""
    info = d["info"]
    return length_di * info["d_i"] / info["v0"] * info["omega_ci"]


def slope(d, zone):
    t = d["t"]
    tt, res, rate, err = residual(t, zone_mean(d, "by", zone), (5.0, float(t[-1])))
    step, tau = float(tt[1] - tt[0]), crossing(d, zone[1] - zone[0])
    rho = [1.0 - k * step / tau for k in range(1, int(tau / step) + 1)]
    factor = float(np.sqrt(1.0 + 2.0 * sum(r for r in rho if r > 0)))
    return {"window": [5.0, float(t[-1])], "outputs": int(len(tt)), "output_step": step, "per_omega_ci": rate,
            "stderr_independent": err, "advection_factor": factor, "stderr_advection": err * factor,
            "rms_residual": float(res.std()), "autocorrelation": {str(k): corr(res, res, k) for k in (1, 2, 3)},
            "autocorrelation_advection": {str(k + 1): max(0.0, r) for k, r in enumerate(rho[:3])}}


def lag(d, zone, far):
    t = d["t"]
    window = (5.0, float(t[-1]))
    tt, ra, _, _ = residual(t, zone_mean(d, "by", zone), window)
    _, rb, _, _ = residual(t, zone_mean(d, "by", far), window)
    shift = int(round(crossing(d, zone[0] - far[0]) / float(tt[1] - tt[0])))
    rows = []
    for n in (1, 15, 30):
        k = np.ones(n) / n
        sa, sb = np.convolve(ra, k, "valid"), np.convolve(rb, k, "valid")
        rows.append({"boxcar_outputs": n, "rms": [float(sa.std()), float(sb.std())],
                     "correlation_lag_0": corr(sa, sb, 0), "correlation_advection_lag": corr(sa, sb, shift)})
    return {"zone_di": list(zone), "far_zone_di": list(far), "advection_lag_outputs": shift, "rows": rows}


def other(d, o, zone):
    t, end = d["t"], float(d["t"][-1])
    if len(t) != len(o["t"]) or np.abs(t - o["t"]).max() > 1e-6:
        raise SystemExit("the two runs do not have the same outputs")
    a = blocks(t, zone_mean(d, "by", zone), 5.0, end, 5.0)
    b = blocks(t, zone_mean(o, "by", zone), 5.0, end, 5.0)
    diff = [p[3] - q[3] for p, q in zip(a, b)]
    _, ra, _, _ = residual(t, zone_mean(d, "by", zone), (5.0, end))
    _, rb, _, _ = residual(t, zone_mean(o, "by", zone), (5.0, end))
    return {"block_means_this_minus_other": diff, "largest": float(np.abs(diff).max()),
            "mean": float(np.mean(diff)), "rms": float(np.std(diff)),
            "correlation_of_residuals": corr(ra, rb, 0), "rms_residual": [float(ra.std()), float(rb.std())]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run")
    ap.add_argument("--zone", type=float, nargs=2, default=(30.0, 37.0), help="d_i")
    ap.add_argument("--far", type=float, nargs=2, default=(2.0, 9.0), help="d_i, the far zone of the lag")
    ap.add_argument("--other", help="a run with the same outputs (another layer depth or seed)")
    ap.add_argument("--json")
    a = ap.parse_args()
    d = load(a.run)
    out = {"run": str(a.run), "omega_ci_t_end": float(d["t"][-1]), "k1": k1(d, a.zone), "k2": k2(d, a.zone),
           "slope": slope(d, a.zone), "lag": lag(d, a.zone, a.far)}
    if a.other:
        out["other"] = {"run": str(a.other)} | other(d, load(a.other), a.zone)
    text = json.dumps(out, indent=1)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
