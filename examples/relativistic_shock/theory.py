"""Jump conditions of a relativistic shock launched by a cold flow hitting a reflecting wall.

In the wall (= downstream, simulation) frame a cold upstream of lab density n1 moves at -beta0
(Lorentz factor gamma0) into the wall; the shocked plasma is at rest and the shock moves at
+beta_sh. Particle, energy and momentum conservation across the shock give

    n2 beta_sh = n1 (beta0 + beta_sh),   <gamma>_2 = gamma0,   P2 = n2 T2 = n1 gamma0 beta0 (beta0 + beta_sh) m
    -> T2 = gamma0 beta0 beta_sh m,   beta_sh = theta2 / (gamma0 beta0),   n2 / n1 = 1 + beta0 / beta_sh,

with theta2 from the equation of state of the downstream (Synge) gas: <gamma>(theta2) = gamma0,
    3-D:  <gamma> = K1(1/theta)/K2(1/theta) + 3 theta
    2-D:  <gamma> = (1 + 2 theta + 2 theta^2) / (1 + theta)      (momenta confined to a plane)
For a constant adiabatic index Gamma (P = (Gamma - 1)(e - n m)):
    beta_sh = (Gamma - 1)(gamma0 - 1) / (gamma0 beta0)   (Spitkovsky 2008: Gamma = 3/2 in 2-D -> n2/n1 = 3.13
                                                          and beta_sh = 0.47 c for gamma0 = 15).
Magnetised flow (jump_conditions_magnetised): the upstream carries a field B1 perpendicular to the
flow and the motional field E = -v x B; sigma = B1^2 / (gamma0 n1 m) with n1 m the lab mass
density of all species (the ratio of the magnetic to the particle enthalpy flux, frame
independent). Flux freezing gives B2 / B1 = n2 / n1 = r = 1 + beta0 / beta_sh; with the Maxwell
stress and the Poynting flux, momentum and energy conservation in units of n1 gamma0 m read

    beta0^2 + sigma (1 + beta0^2) / 2 + beta_sh beta0 (1 + sigma) = r theta2 / gamma0 + sigma r^2 / 2
    beta0 (1 + sigma) + beta_sh (1 + sigma (1 + beta0^2) / 2)     = beta_sh (r <gamma>_2 / gamma0 + sigma r^2 / 2)

(two equations for beta_sh and theta2). sigma = 0 gives the relations above; gamma0 -> infinity in
3-D gives Kennel & Coroniti (1984): the downstream four-velocity in the shock frame obeys
    u2^2 = [8 sigma^2 + 10 sigma + 1 + sqrt(64 sigma^2 (sigma + 1)^2 + 20 sigma (sigma + 1) + 1)] / (16 (sigma + 1)).
Particle acceleration (known, not analytic): Weibel-mediated unmagnetised pair shocks develop a
non-thermal tail dN/dgamma ~ gamma^-p with p ~ 2.4 (2-D, Spitkovsky 2008b; Sironi et al. 2013),
the theoretical isotropic-diffusion index for ultra-relativistic shocks is p ~ 2.2-2.3
(Kirk et al. 2000; Achterberg et al. 2001), and gamma_max grows ~ t^(1/2) (Sironi et al. 2013).
"""
from __future__ import annotations

import math

from scipy.optimize import brentq
from scipy.special import kve


def mean_gamma(theta: float, dims: int = 3) -> float:
    if dims == 3:
        return kve(1, 1 / theta) / kve(2, 1 / theta) + 3 * theta
    if dims == 2:
        return (1 + 2 * theta + 2 * theta * theta) / (1 + theta)
    raise ValueError("dims must be 2 or 3")


def jump_conditions(gamma0: float, dims: int = 3) -> dict:
    """Exact (Synge gas) downstream temperature, shock speed and density compression (lab frame)."""
    beta0 = math.sqrt(1 - 1 / gamma0 ** 2)
    theta2 = brentq(lambda th: mean_gamma(th, dims) - gamma0, 1e-8, 10 * gamma0 + 10)
    beta_sh = theta2 / (gamma0 * beta0)
    return {"theta2": theta2, "beta_sh": beta_sh, "compression": 1 + beta0 / beta_sh, "beta0": beta0,
            "gamma0": gamma0, "dims": dims}


def jump_conditions_gamma(gamma0: float, Gamma: float) -> dict:
    """Constant adiabatic index Gamma (e.g. 4/3 for 3-D, 3/2 for 2-D ultra-relativistic gases)."""
    beta0 = math.sqrt(1 - 1 / gamma0 ** 2)
    beta_sh = (Gamma - 1) * (gamma0 - 1) / (gamma0 * beta0)
    return {"beta_sh": beta_sh, "compression": 1 + beta0 / beta_sh, "theta2": (Gamma - 1) * (gamma0 - 1)}


def jump_conditions_magnetised(gamma0: float, sigma: float, dims: int = 3) -> dict:
    """Shock speed, compression (= B2 / B1) and downstream temperature of a cold flow with a
    perpendicular field of magnetisation sigma (Synge gas downstream, lab = downstream frame)."""
    beta0 = math.sqrt(1 - 1 / gamma0 ** 2)

    def pressure(bs):  # theta2 r / gamma0 from the momentum equation
        r = 1 + beta0 / bs
        return beta0 ** 2 + 0.5 * sigma * (1 + beta0 ** 2) + bs * beta0 * (1 + sigma) - 0.5 * sigma * r * r

    def energy(bs):
        r = 1 + beta0 / bs
        theta = gamma0 * pressure(bs) / r
        return (bs * (r * mean_gamma(max(theta, 1e-8), dims) / gamma0 + 0.5 * sigma * r * r) - beta0 * (1 + sigma)
                - bs * (1 + 0.5 * sigma * (1 + beta0 ** 2)))

    # the downstream pressure is positive above lo; a slower shock would need a colder gas than T = 0
    lo = brentq(pressure, 1e-12, 1.0) * (1 + 1e-12) if sigma > 0 else 1e-12
    if energy(lo) >= 0 or energy(1.0) <= 0:
        raise ValueError(f"no shock solution for gamma0={gamma0:g}, sigma={sigma:g}")
    beta_sh = brentq(energy, lo, 1.0, xtol=1e-14)
    r = 1 + beta0 / beta_sh
    return {"theta2": gamma0 * pressure(beta_sh) / r, "beta_sh": beta_sh, "compression": r, "beta0": beta0,
            "gamma0": gamma0, "dims": dims, "sigma": sigma}
