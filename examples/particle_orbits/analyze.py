"""Compare tracked test-particle orbits with the exact relativistic solutions (theory.py).

python examples/particle_orbits/analyze.py output/particle_orbits/<case> [--refined <run with dt_scale=0.5>]
verdict: relative position error <= 2e-2 at the base step and, with --refined, second-order
convergence: error(dt) / error(dt/2) in [3, 5].
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import theory as T  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402

Q, M = -1.0, 1.0


def exact(info, t, x0, u0):
    case = info["case"]
    E0, B0 = info["E0"], info["B0"]
    if case == "gyration":
        return T.gyration(t, x0, u0, Q, M, B0[2])
    if case == "exb":
        return T.exb(t, x0, Q, M, E0[1], B0[2])
    if case == "hyperbolic":
        return T.hyperbolic(t, x0, Q, M, E0[0], axis=0)
    if case == "e_gt_b":
        return T.e_gt_b(t, x0, Q, M, E0[1], B0[2])
    return T.plane_wave(t, x0, Q, M, info["E_wave"], info["omega"])


def orbit_error(run_dir):
    run = open_run(run_dir)
    info = run.meta["planner"]
    tr = run.tracks("tracks", "electron")
    t = tr.t
    L = np.asarray(run.L)
    ids = tr.ids.astype(int)
    errs, orbits = [], []
    for col, pid in enumerate(ids):
        p = info["particles"][pid]
        X, _ = exact(info, t, p[:3], p[3:6])
        xs = np.stack([tr[c][:, col] for c in "xyz"])
        # unwrap the periodic coordinates onto the exact (unbounded) orbit
        xs = xs + L[:, None] * np.round((X - xs) / L[:, None])
        scale = max(float(np.max(np.abs(X - X[:, :1]))), 1e-3)
        errs.append(float(np.max(np.abs(xs - X))) / scale)
        orbits.append((xs, X))
    return run, info, t, np.asarray(errs), orbits


def analyze(run_dir: str, show_plot: bool = True, refined: str | None = None) -> dict:
    run, info, t, errs, orbits = orbit_error(run_dir)
    res = {"case": info["case"], "dt": info["dt"], "max_rel_position_error": float(errs.max()),
           "per_particle": errs.tolist()}
    ok = res["max_rel_position_error"] <= 2e-2
    if info["case"] == "plane_wave":
        tr = run.tracks("tracks", "electron")
        p0 = info["particles"][int(tr.ids[0])]
        x = orbits[0][0][0]
        res["drift_measured"] = float((x[-1] - p0[0]) / t[-1])
        res["drift_theory_node"] = info["drift_theory"]
    if refined:
        _, info2, _, errs2, _ = orbit_error(refined)
        ratio = float(errs.max() / max(errs2.max(), 1e-300))
        res.update(refined_dt=info2["dt"], refined_max_error=float(errs2.max()), convergence_ratio=ratio)
        ok = ok and 3.0 <= ratio <= 5.0
    res["passed"] = bool(ok)
    res["gauss_max"] = 0.0  # test particles deposit nothing: no Gauss check applies
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(1, 2, figsize=(11, 4.2))
        for k, (xs, X) in enumerate(orbits[:6]):
            ax[0].plot(X[0], X[1], **{**P.THEORY, "lw": 0.9}, label="exact" if k == 0 else None)
            ax[0].plot(xs[0], xs[1], color=P.color(k), lw=1.0, alpha=0.8, label=f"PIC #{k}")
        ax[0].set_xlabel("x [c/ωpe]")
        ax[0].set_ylabel("y [c/ωpe]")
        ax[0].legend(ncol=2)
        ax[0].set_title(f"{info['case']}: orbits")
        for k, (xs, X) in enumerate(orbits[:6]):
            ax[1].semilogy(t, np.max(np.abs(xs - X), axis=0) + 1e-16, color=P.color(k))
        ax[1].set_xlabel(r"$\omega t$")
        ax[1].set_title("|x_PIC − x_exact|")
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "orbits.png"))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--refined")
    a = ap.parse_args()
    print(json.dumps(analyze(a.run_dir, refined=a.refined), indent=2))
