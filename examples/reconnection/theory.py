"""Relativistic pair Harris sheet: equilibrium (analytic) and published reconnection scalings
(empirical, not analytic).

Normalisation: background density n_b per species = 1 (omega_pe = 1, c/omega_pe = 1), magnetisation
of the cold background sigma = B0^2 / (n_b,tot m c^2) = B0^2 / 2 (eps0 = 1 units).
Harris sheet of half-thickness delta and peak density n_cs per species, drifting at +-beta_d along z:
    B_x(y) = B0 tanh(y / delta),  n(y) = n_cs sech^2(y / delta)
    Ampere:            n_cs beta_d = B0 / (2 delta)
    pressure balance:  B0^2 / 2 = 2 n_cs' theta_cs,  n_cs' = n_cs / Gamma_d (rest-frame density)
Empirical results (relativistic pair reconnection, 2-D):
    reconnection rate  v_in / v_A ~ 0.1 (e.g. Werner et al. 2018; Sironi & Spitkovsky 2014)
    v_A = c sqrt(sigma / (1 + sigma))
    power-law index    p ~ 1.9 + 0.7 / sqrt(sigma) for dN/dgamma (Werner et al., ApJL 816, L8 (2016));
                       Sironi & Spitkovsky (2014): p ~ 4, 3, 2, 1.5 for sigma = 1, 3, 10, 50
"""
from __future__ import annotations

import math


def harris(sigma: float, n_cs_over_nb: float, delta: float, n_b: float = 1.0) -> dict:
    B0 = math.sqrt(2.0 * n_b * sigma)
    n_cs = n_cs_over_nb * n_b
    beta_d = B0 / (2.0 * delta * n_cs)
    if beta_d >= 1.0:
        raise ValueError(f"Harris sheet needs beta_d = {beta_d:.3f} < 1: thicker sheet or denser")
    gd = 1.0 / math.sqrt(1.0 - beta_d ** 2)
    theta = B0 ** 2 * gd / (4.0 * n_cs)
    return {"B0": B0, "n_cs": n_cs, "beta_d": beta_d, "u_d": gd * beta_d, "Gamma_d": gd, "theta_cs": theta,
            "v_A": math.sqrt(sigma / (1.0 + sigma))}


def werner_index(sigma: float) -> float:
    return 1.9 + 0.7 / math.sqrt(sigma)
