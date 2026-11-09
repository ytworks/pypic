"""Analyse a plasma-expansion run against Mora's isothermal self-similar solution.

python examples/plasma_expansion/analyze.py output/plasma_expansion/isothermal [--movie]
gate (proton mass config): the self-similar ion velocity v = c_s + x / t of the quasi-neutral
      rarefaction (-0.8 < x / c_s t < 3) within 10 % (median), and E_max(t) within 20 % (median over
      2 <= tau <= tau_end) of 2 Z T_e [ln(tau + sqrt(tau^2 + 1))]^2.
Collisionless PIC electrons are not the ideal isothermal (Boltzmann, infinite heat bath) electrons
of the model: the expanding region is cooler than the slab, and the measured velocities (~8 %) and
E_max (~15 %) lie below the isothermal solution consistently with a ~15 % lower effective T_e there.
The deviation does not change with resolution, particle number or vacuum length, and shrinks with the
mass ratio (the model needs v_ion << v_te), so the criteria are set above this physical offset.
reported: front position, ion spectrum at the end, electron temperature.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import theory as T  # noqa: E402

from pypic.analysis import open_run  # noqa: E402
from pypic.analysis import plot as P  # noqa: E402
from pypic.analysis.fields import front_position  # noqa: E402
from pypic.analysis.particles import fit_temperature, spectrum  # noqa: E402


def analyze(run_dir: str, show_plot: bool = True) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    Te, mi, Z = info["Te"], info["mass_ratio"], info["Z"]
    sp = run.series("ion_spectrum")
    edges = sp.edges[0]
    t = sp.t
    h = sp["ion"]
    tot = h.sum(axis=1)
    emax = []
    for i in range(len(t)):
        nz = np.nonzero(h[i] > 1e-9 * max(tot[i], 1e-300))[0]
        emax.append(edges[nz[-1] + 1] if nz.size else np.nan)
    emax = np.asarray(emax)
    tau = T.tau(t, Te, mi, Z)
    eth = T.max_energy(t, Te, mi, Z)
    sel = tau >= 2.0
    rel = (emax[sel] - eth[sel]) / eth[sel]
    lo = run.series("lineout")
    x = lo.x("n:ion")
    ni = lo["n:ion"]
    xf = np.array([front_position(x, ni[i], 1e-3 / Z) for i in range(len(lo.t))]) - info["edge"]
    xf_th = T.front_position(lo.t, Te, mi, Z)
    fsel = T.tau(lo.t, Te, mi, Z) >= 2.0
    E, f = spectrum(edges, h[-1] / info["A_perp"])
    f_th = T.ion_spectrum(E, t[-1], Te, mi, Z)
    win = (E > 0.1 * eth[-1]) & (E < 0.7 * eth[-1]) & (f > 0) & (f_th > 0)
    # self-similar ion velocity profile v = c_s (1 + xi), xi = (x - edge) / (c_s t), at the last output
    ph = run.series("ion_phase")
    pex, peu = ph.edges
    hp = ph.read("ion", -1)
    xc = 0.5 * (pex[1:] + pex[:-1])
    uc = 0.5 * (peu[1:] + peu[:-1])
    wsum = hp.sum(axis=1)
    vmean = (hp * uc[None, :]).sum(axis=1) / np.maximum(wsum, 1e-300)
    cs = info["cs"]
    xi = (xc - info["edge"]) / (cs * ph.t[-1])
    vsel = (xi > -0.8) & (xi < 3.0) & (wsum > 0)
    v_err = np.abs(vmean[vsel] - cs * (1 + xi[vsel])) / (cs * (1 + xi[vsel]))
    es = run.series("electron_spectrum")
    Ee, fe = spectrum(es.edges[0], es["electron"][-1])
    theta_end, _ = fit_temperature(1.0 + Ee, fe)
    res = {"tau_end": float(tau[-1]), "E_max_final": float(emax[-1]), "E_max_theory_final": float(eth[-1]),
           "E_max_median_rel_error": float(np.nanmedian(np.abs(rel))),
           "E_max_max_rel_error": float(np.nanmax(np.abs(rel))),
           "front_median_rel_error": float(np.nanmedian(np.abs((xf[fsel] - xf_th[fsel]) / xf_th[fsel]))),
           "spectrum_rms_log_ratio": float(np.sqrt(np.mean(np.log(f[win] / f_th[win]) ** 2))) if win.any() else None,
           "electron_temperature_final_over_initial": theta_end / Te}
    res["velocity_profile_median_rel_error"] = float(np.median(v_err)) if v_err.size else float("nan")
    if info["params"]["slab_lambda_D"] < 50:
        # thin foil: the electrons cool (adiabatic expansion), so the isothermal law is an upper bound
        res["regime"] = "thin foil: isothermal E_max is an upper bound"
        res["passed"] = bool(np.all(emax[sel] <= 1.05 * eth[sel]))
    else:
        res["regime"] = "isothermal (thick slab)"
        res["passed"] = bool(res["velocity_profile_median_rel_error"] <= 0.10 and res["E_max_median_rel_error"] <= 0.20)
    s = run.scalars
    g = s["gauss"][np.isfinite(s["gauss"])]
    res["gauss_max"] = float(g.max()) if g.size else float("nan")
    if show_plot:
        P.style()
        fig, ax = P.plt.subplots(1, 3, figsize=(15, 4.2))
        ax[0].plot(tau, emax / Te, "o", ms=3, color=P.color(0), label="PIC (max ion energy)")
        ax[0].plot(tau, eth / Te, **P.THEORY, label="Mora 2003")
        ax[0].set_xlabel(r"$\tau = \omega_{pi} t / \sqrt{2e}$")
        ax[0].set_ylabel(r"$E_{max} / Z T_e$")
        ax[0].legend()
        tl = T.tau(lo.t, Te, mi, Z)
        ax[1].plot(tl, xf / math.sqrt(Te), "o", ms=3, color=P.color(1), label="PIC (ion front)")
        ax[1].plot(tl, xf_th / math.sqrt(Te), **P.THEORY, label="Mora 2003")
        ax[1].set_xlabel(r"$\tau$")
        ax[1].set_ylabel(r"$x_f / \lambda_{D0}$")
        ax[1].legend()
        ok = f > 0
        ax[2].loglog(E[ok] / Te, f[ok], color=P.color(0), label=f"PIC, τ={tau[-1]:.1f}")
        okt = f_th > 0
        ax[2].loglog(E[okt] / Te, f_th[okt], **P.THEORY, label="Mora 2003")
        ax[2].set_xlabel(r"$E / Z T_e$")
        ax[2].set_ylabel("dN/dE per area")
        ax[2].legend()
        fig.tight_layout()
        res["plot"] = P.save(fig, os.path.join(run_dir, "expansion.png"))
    return res


def movie(run_dir: str) -> str:
    from pypic.analysis import Movie

    run = open_run(run_dir)
    m = Movie(run, title="plasma expansion into vacuum")
    m.phase("ion_phase", "ion").lineout("lineout", ["n:ion", "n:electron"]).lineout("lineout", ["Ex"])
    m.spectrum("ion_spectrum", "ion")
    return m.save(os.path.join(run_dir, "expansion.mp4"), fps=15, n_frames=100)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    r = analyze(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2))
