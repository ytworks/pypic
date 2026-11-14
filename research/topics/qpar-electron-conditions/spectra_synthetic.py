"""The estimators of spectra.py on a known distribution (a check of the code, not a null test of the README).

    python research/topics/qpar-electron-conditions/spectra_synthetic.py

CPU only (NumPy). The distribution is a Maxwellian f ~ exp(-p^2 / p_th^2) plus a tail ~ p^-q that sets in near
3 p_th and is cut off at 30 m_i v_pt. It is put into the bins of V1's electron histogram (240 bins in ln p), once
without noise and once as Poisson counts of 3e8 particles. Printed: each number of ``spectra.numbers`` and its
deviation from the value of the distribution itself (computed on bins 400 times finer). phi is by its definition a
mean over +-0.1 dex; for f ~ p^-4 that mean is 0.9 % above the value at p_ref.
"""
from __future__ import annotations

import os
import sys
import types

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import spectra as sx  # noqa: E402

ROOT = 0.1  # sqrt(m_e / m_i) of V1
FINE = 400


def main():
    pe = np.logspace(np.log10(1.333e-4), 2.0, 241)  # V1's electron bins in m_i v_pt
    sp = types.SimpleNamespace(pe=pe, p=np.sqrt(pe[1:] * pe[:-1]), dlnp=np.log(pe[1] / pe[0]))
    for p_th, q, amp in ((0.058, 4.0, 3e-3), (0.058, 4.3, 3e-3), (0.04, 3.8, 1e-3)):
        fine = np.logspace(np.log10(pe[0]), np.log10(pe[-1]), 240 * FINE + 1)
        pc = np.sqrt(fine[1:] * fine[:-1])
        f = np.exp(-(pc / p_th) ** 2) / (np.pi ** 1.5 * p_th ** 3)
        f += (amp / (4 * np.pi * (5 * p_th) ** 3) * (pc / (5 * p_th)) ** -q
              * (1 - np.exp(-(pc / (3 * p_th)) ** 8)) * np.exp(-pc / 30.0))
        dn = 4 * np.pi * pc ** 3 * f * np.log(fine[1] / fine[0])
        n = dn.sum()
        exact = {"p_th": pc[np.argmax(pc ** 4 * f)] / np.sqrt(2), "eta_fixed": dn[pc >= 3 * ROOT].sum() / n}
        exact["eta_nt"] = dn[pc >= 5 * exact["p_th"]].sum() / n
        k = np.argmin(np.abs(pc - 6 * ROOT))
        exact["phi"] = 4 * np.pi * pc[k] ** 3 * f[k] / n
        sel = (pc >= 0.6) & (pc <= 2.4)
        exact["q"] = -np.polyfit(np.log(pc[sel]), np.log(f[sel]), 1)[0]
        for total in (None, 3e8):
            c = dn.reshape(240, FINE).sum(axis=1) / n * (total or 1e12)
            if total:
                c = np.random.default_rng(3).poisson(c).astype(float)
            got = sx.numbers(c, sp, True, ROOT, (0.6, 2.4))
            got["q"] = got["slope"]["q"]
            print(f"p_th={p_th} q={q} amp={amp} particles={'no noise' if not total else f'{total:.0e}'} "
                  f"q_stderr={got['slope']['stderr']:.4f} "
                  + " ".join(f"{k}: {got[k]:.5g} ({got[k] / v - 1:+.2%})" for k, v in exact.items()))


if __name__ == "__main__":
    main()
