"""Analyse a two-stream run: growth rate vs linear theory, energy/conservation history, phase space.

python examples/two_stream/analyze.py output/ei            (reads history.npz + meta.json)
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def linear_window(t, amp, above_noise=2.0, sat_drop=1.0):
    """Index range of the linear phase: amplitude from `above_noise` e-folds above the initial
    noise up to `sat_drop` e-folds below the first saturation peak."""
    la = np.log(np.maximum(np.abs(amp), 1e-300))
    n0 = max(3, len(la) // 20)
    noise = np.median(la[:n0])
    isat = int(np.argmax(la))
    lo, hi = noise + above_noise, la[isat] - sat_drop
    if hi <= lo:
        lo, hi = noise + 1.0, la[isat] - 0.5
    above = np.where(la[: isat + 1] >= lo)[0]
    i0 = int(above[0]) if len(above) else 0
    below = np.where(la[i0: isat + 1] <= hi)[0]
    i1 = i0 + int(below[-1]) if len(below) else isat
    return i0, i1


def fit_growth(t, z, min_pts=8, omega_r=0.0):
    """Growth rate of one Fourier mode over its linear phase.

    omega_r == 0 (purely growing): straight-line fit of log|z|.
    omega_r != 0 (oscillating): the two counter-propagating growing waves make |z| beat with
    period pi/omega_r, so the complex coefficients are fitted with
        z(t) = exp(gamma t) (a exp(-i w t) + b exp(i w t))
    by variable projection (a, b linear; gamma, w nonlinear, started from the log fit and the
    theoretical omega_r), using relative residuals."""
    t = np.asarray(t, dtype=float)
    z = np.asarray(z)
    i0, i1 = linear_window(t, z)
    sel = np.arange(i0, i1 + 1)
    if len(sel) < min_pts:
        return float("nan"), None
    ts, zs = t[sel], z[sel]
    g = float(np.polyfit(ts, np.log(np.abs(zs)), 1)[0])
    if np.iscomplexobj(z) and abs(omega_r) > 1e-3:
        from scipy.optimize import minimize

        tt = ts - ts[0]
        wgt = 1.0 / np.abs(zs)

        def resid(pv):
            gg, ww = pv
            M = np.exp(gg * tt)[:, None] * np.stack([np.exp(-1j * ww * tt), np.exp(1j * ww * tt)], 1)
            coef, *_ = np.linalg.lstsq(M * wgt[:, None], zs * wgt, rcond=None)
            return float(np.sum(np.abs((M @ coef - zs) * wgt) ** 2))

        best = min((minimize(resid, [g, w0], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-12})
                    for w0 in (abs(omega_r), 0.8 * abs(omega_r), 1.25 * abs(omega_r))), key=lambda r: r.fun)
        g = float(best.x[0])
    a = np.abs(zs)
    c0 = np.exp(np.mean(np.log(a)) - g * np.mean(ts))
    return g, (ts, c0 * np.exp(g * ts))


def analyze(run_dir: str, show_plot: bool = True):
    h = np.load(os.path.join(run_dir, "history.npz"), allow_pickle=True)
    meta_path = os.path.join(run_dir, "meta.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    lin = meta.get("linear", {})
    nwl = meta.get("params", {}).get("n_wavelengths", 4)
    t = h["t"]
    modes = h["modes"]
    mt = h["mode_t"]
    res = {}
    if modes.size:
        m = min(nwl, modes.shape[1]) - 1  # mode index nwl corresponds to k_max
        g_fit, line = fit_growth(mt, modes[:, m], omega_r=lin.get("omega_r", 0.0))
        gth = lin.get("gamma_max", float("nan"))
        res.update(gamma_fit_kmax=g_fit, gamma_theory=gth, relative_error=(g_fit - gth) / gth)
    e0 = h["budget"][0]
    res["energy_error_max"] = float(np.max(np.abs(h["budget"] - e0)) / abs(e0))
    g = h["gauss"][np.isfinite(h["gauss"])]
    b = h["divB"][np.isfinite(h["divB"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    res["divB_max"] = float(b.max()) if b.size else float("nan")

    if show_plot:
        species = [k[len("kinetic_"):] for k in h.files if k.startswith("kinetic_")]
        phases = [k for k in h.files if k.startswith("phase_") and k not in ("phase_t", "phase_umax")]
        fig = plt.figure(figsize=(14, 9))
        ax1 = fig.add_subplot(2, 3, 1)
        for s in species:
            ax1.plot(t, h[f"kinetic_{s}"], label=f"K {s}")
        ax1.plot(t, h["electric"], label="E field")
        ax1.plot(t, h["magnetic"], label="B field")
        ax1.set_yscale("log")
        ax1.set_xlabel(r"$\omega_{pe} t$")
        ax1.set_title("energy")
        ax1.legend(fontsize=7)
        ax2 = fig.add_subplot(2, 3, 2)
        ax2.semilogy(t, np.abs(h["budget"] - e0) / abs(e0) + 1e-17)
        ax2.set_title(r"|$\Delta E_{tot}$| / $E_0$")
        ax2.set_xlabel(r"$\omega_{pe} t$")
        ax3 = fig.add_subplot(2, 3, 3)
        fin = np.isfinite(h["gauss"])
        ax3.semilogy(t[fin], h["gauss"][fin] + 1e-20, "o-", label=r"max|$\nabla\cdot E-\rho$|")
        ax3.semilogy(t[fin], h["divB"][fin] + 1e-20, "s-", label=r"max|$\nabla\cdot B$|")
        ax3.legend(fontsize=8)
        ax3.set_title("conservation")
        ax4 = fig.add_subplot(2, 3, 4)
        if modes.size:
            for j in range(min(8, modes.shape[1])):
                ax4.semilogy(mt, np.abs(modes[:, j]), lw=2.5 if j == m else 0.8, label=f"m={j + 1}")
            if line is not None:
                ax4.semilogy(line[0], line[1] * 2, "k--",
                             label=f"fit $\\gamma$={g_fit:.3f} (theory {res['gamma_theory']:.3f})")
            ax4.legend(fontsize=7)
        ax4.set_xlabel(r"$\omega_{pe} t$")
        ax4.set_title(r"$|E_x(k_m)|$, $k_m = 2\pi m/L_x$")
        el = [p for p in phases if "electron" in p] or phases
        if el:
            ph = h[el[0]]
            pt = h["phase_t"]
            umax = float(h["phase_umax"])
            Lx = float(h["Lx"])
            for i, idx in enumerate([len(ph) // 2, len(ph) - 1]):
                ax = fig.add_subplot(2, 3, 5 + i)
                ax.imshow(ph[idx].T, origin="lower", aspect="auto", extent=[0, Lx, -umax, umax], cmap="magma")
                ax.set_title(f"{el[0][6:]} (x, u_x) at t={pt[idx]:.1f}")
                ax.set_xlabel("x [c/omega_pe]")
                ax.set_ylabel("u_x / c")
        fig.suptitle(os.path.basename(os.path.normpath(run_dir)))
        fig.tight_layout()
        out = os.path.join(run_dir, "summary.png")
        fig.savefig(out, dpi=110)
        res["plot"] = out
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    a = ap.parse_args()
    print(json.dumps(analyze(a.run_dir), indent=2))
