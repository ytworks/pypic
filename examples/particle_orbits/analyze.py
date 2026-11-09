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


def summary(results: dict[str, dict], path: str) -> str:
    """Second-order convergence: maximum relative position error at dt and dt/2 for every case, with
    dt^2 guide lines (a ratio of 4 per halving)."""
    P.style()
    fig, ax = P.plt.subplots(figsize=(6.4, 4.2))
    for i, (name, r) in enumerate(results.items()):
        e1, e2 = r["max_rel_position_error"], r.get("refined_max_error")
        if e2 is None:
            continue
        ax.loglog([1.0, 0.5], [e1, e2], "o-", ms=7, color=P.color(i),
                  label=f"{name} (ratio {r['convergence_ratio']:.2f})")
        g = np.array([1.15, 0.43])
        ax.loglog(g, e1 * g ** 2, **{**P.THEORY, "lw": 0.8}, label=(r"$\propto \Delta t^2$" if i == 0 else None))
    ax.set_xlabel(r"$\Delta t / \Delta t_0$")
    ax.set_ylabel("max |x_PIC − x_exact| / orbit size")
    ax.set_xticks([0.5, 1.0], ["1/2", "1"])
    ax.xaxis.set_minor_formatter(P.plt.NullFormatter())
    ax.legend(fontsize=7)
    fig.tight_layout()
    return P.save(fig, path, dpi=110)


def movie(runs: dict[str, str], out: str, fps: int = 12, n_frames: int = 90, dpi: int = 90,
          n_particles: int = 3) -> str:
    """One panel per case: the PIC orbit (solid, growing) traced over the exact orbit (dashed).
    The plane-wave case is drawn in the frame drifting with the exact mean velocity (figure-8)."""
    from matplotlib import animation

    P.style()
    cases = []
    for name, d in runs.items():
        run, info, t, _, orbits = orbit_error(d)
        if info["case"] == "gyration":  # the three largest Larmor radii (u_perp = 0.1 ... 10)
            pick = list(range(len(orbits)))[-n_particles:]
        elif info["case"] == "plane_wave":  # the particle the drift theory refers to (starts where A = 0)
            pick = [0]
        else:
            pick = list(range(len(orbits)))[:n_particles]
        drift = info.get("drift_theory", 0.0) if info["case"] == "plane_wave" else 0.0
        curves = []
        for k in pick:
            xs, X = orbits[k]
            if info["case"] == "hyperbolic":  # 1-D motion: show x(t)
                a, b, A, B = t, xs[0], t, X[0]
            else:
                a, b, A, B = xs[0] - drift * t, xs[1], X[0] - drift * t, X[1]
            # centre each orbit on its own start so several particles share one panel
            if info["case"] != "hyperbolic":
                a, A = a - A[0], A - A[0]
                b, B = b - B[0], B - B[0]
            curves.append((a, b, A, B))
        cases.append((name, info, t, curves))
    nc = len(cases)
    fig, axes = P.plt.subplots(1, nc, figsize=(3.3 * nc, 3.4), squeeze=False)
    lines = []
    for ax, (name, info, t, curves) in zip(axes[0], cases):
        for k, (a, b, A, B) in enumerate(curves):
            ax.plot(A, B, **{**P.THEORY, "lw": 0.9}, label="exact" if k == 0 else None)
            ln, = ax.plot([], [], color=P.color(k), lw=1.6, label="PIC" if k == 0 else None)
            pt, = ax.plot([], [], "o", ms=5, color=P.color(k))
            lines.append((ln, pt, a, b, t))
        ttl = {"plane_wave": "plane wave (drift frame)", "hyperbolic": "hyperbolic: x(t)",
               "exb": "E×B drift", "e_gt_b": "E > B", "gyration": "gyration"}.get(info["case"], name)
        ax.set_title(ttl)
        ax.set_xlabel("ω t" if info["case"] == "hyperbolic" else "x − x₀")
        ax.set_ylabel("x" if info["case"] == "hyperbolic" else "y − y₀")
        if info["case"] == "gyration":
            ax.set_aspect("equal", adjustable="datalim")
        ax.legend(loc="best", fontsize=7)
    fig.tight_layout()
    t_all = np.linspace(0.0, 1.0, n_frames)

    def update(f):
        art = []
        for ln, pt, a, b, t in lines:
            j = max(1, int(round(t_all[f] * (len(t) - 1))))
            ln.set_data(a[: j + 1], b[: j + 1])
            pt.set_data([a[j]], [b[j]])
            art += [ln, pt]
        return art

    ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=False)
    ani.save(out, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    P.plt.close(fig)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--refined")
    a = ap.parse_args()
    print(json.dumps(analyze(a.run_dir, refined=a.refined), indent=2))
