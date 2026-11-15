"""Diagnostic look at a flow-through control's <B_y>/B_0y series in the face zone.

Not a criterion: the pass criteria of the control are those of `pilot_check.py`, fixed before the runs. This
script was written after the first control had been read, to describe the series whose two halves disagreed:
block means, the autocorrelation of the residuals about the straight-line fit, and the largest spectral lines.

    python control_look.py <run dir> [--pair <other run dir>] [--json out.json]

``--pair`` compares the series with that of a second control of another seed: the correlation of the two
residual series. A change of the field that the boundaries cause at fixed times is common to both runs; a
fluctuation that depends on the seed is not. The standard error of the correlation under "no common part" uses
the effective number of outputs N (1 - r_a r_b) / (1 + r_a r_b) with the two lag-1 autocorrelations (Bartlett).
"""
import argparse
import json

import numpy as np

from pypic.analysis.run import open_run


def series(path, zone, window):
    run = open_run(path)
    info = run.meta["planner"]
    prof = run.series("profile")
    t = prof.t * info["omega_ci"]
    x = prof.x("By") / info["d_i"]
    by = prof["By"][:, (x >= zone[0]) & (x <= zone[1])].mean(axis=1) / info["by"]
    sel = (t >= window[0]) & (t <= window[1])
    tt, y = t[sel], by[sel]
    return tt, y, y - np.polyval(np.polyfit(tt, y, 1), tt)


def pair(a, b, zone, window, block):
    ta, ya, ra = series(a, zone, window)
    tb, yb, rb = series(b, zone, window)
    if len(ta) != len(tb) or np.abs(ta - tb).max() > 1e-6:
        raise SystemExit("the two runs do not have the same outputs in the window")
    r1a, r1b = (float(np.corrcoef(r[:-1], r[1:])[0, 1]) for r in (ra, rb))
    n_eff = len(ra) * (1 - r1a * r1b) / (1 + r1a * r1b)
    edges = np.arange(window[0], window[1], block)
    ma, mb = ([float(y[(ta >= e) & (ta < e + block + 1e-9)].mean() - 1.0) for e in edges] for y in (ya, yb))
    return {
        "other": str(b), "outputs": int(len(ra)),
        "correlation_of_residuals": float(np.corrcoef(ra, rb)[0, 1]),
        "effective_outputs": float(n_eff), "stderr_if_no_common_part": float(1.0 / np.sqrt(n_eff)),
        "block_means_minus_1": {"this": ma, "other": mb, "correlation": float(np.corrcoef(ma, mb)[0, 1])},
        "rms_of_half_sum_residual": float((0.5 * (ra + rb)).std()),
        "rms_of_half_difference_residual": float((0.5 * (ra - rb)).std()),
    }


def look(path, zone, window, block, lags):
    tt, y, res = series(path, zone, window)
    step = float(tt[1] - tt[0])
    auto = {str(lag): float(np.corrcoef(res[:-lag], res[lag:])[0, 1]) for lag in lags if lag < len(res) - 2}
    r1 = auto["1"]
    blocks = []
    for a in np.arange(window[0], window[1], block):
        m = (tt >= a) & (tt < a + block + 1e-9)
        if m.sum():
            blocks.append({"from": float(a), "to": float(a + block), "outputs": int(m.sum()),
                           "mean_minus_1": float(y[m].mean() - 1.0)})
    power = np.abs(np.fft.rfft(res * np.hanning(len(res)))) ** 2
    top = [int(k) for k in np.argsort(power)[::-1][:4] if k]
    span = float(tt[-1] - tt[0])
    return {
        "run": str(path), "zone_di": list(zone), "window": list(window), "outputs": int(len(tt)),
        "output_step": step, "rms_residual": float(res.std()),
        "autocorrelation_by_lag": auto,
        "ar1_error_inflation": float(np.sqrt((1 + r1) / (1 - r1))),
        "blocks": blocks,
        "spectral_lines": [{"period": span / k, "share": float(power[k] / power.sum())} for k in top],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run")
    ap.add_argument("--zone", type=float, nargs=2, default=(30.0, 37.0), help="d_i")
    ap.add_argument("--window", type=float, nargs=2, default=(4.0, 40.0), help="1 / Omega_ci")
    ap.add_argument("--block", type=float, default=6.0, help="length of the block means, 1 / Omega_ci")
    ap.add_argument("--lags", type=int, nargs="+", default=(1, 2, 5, 10, 20, 40), help="in outputs")
    ap.add_argument("--pair", help="a second control (another seed) to compare the series with")
    ap.add_argument("--json")
    a = ap.parse_args()
    out = look(a.run, a.zone, a.window, a.block, a.lags)
    if a.pair:
        out["pair"] = pair(a.run, a.pair, a.zone, a.window, a.block)
    text = json.dumps(out, indent=1)
    if a.json:
        with open(a.json, "w") as f:
            f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
