"""Fields of the drivers of a plasma wake.

Units: c = eps0 = 1, lengths in c/omega_pe, densities in n0, charges in e.

Gaussian bunch. A rigid bunch of the charge Q with the density

    rho = Q / ((2 pi)^(3/2) sigma_x sigma_r^2) exp(-xi^2 / (2 sigma_x^2) - r^2 / (2 sigma_r^2)),   xi = x - beta t,

has in vacuum the potential phi with (k_perp^2 + k_x^2 / gamma^2) phi = rho, the fields
E_r = -d phi / dr, B_theta = beta E_r. In its rest frame it is a Gaussian at rest with the length
gamma sigma_x. With 1 / |x| = (2 / sqrt(pi)) int_0^inf exp(-s^2 x^2) ds the convolution with a
Gaussian is closed, and in the mid-plane xi = 0 (derived here, no paper was consulted)

    E_r(r) = gamma Q / (pi^(3/2)) r int_0^inf s^2 exp(-s^2 r^2 / (1 + 2 s^2 sigma_r^2))
                                        / ((1 + 2 s^2 sigma_r^2)^2 sqrt(1 + 2 s^2 gamma^2 sigma_x^2)) ds.

Limits: sigma -> 0 and gamma = 1 give Q / (4 pi r^2); gamma sigma_x >> r, sigma_r gives the line
charge, lambda / (2 pi r) (1 - exp(-r^2 / (2 sigma_r^2))) with lambda = Q / (sqrt(2 pi) sigma_x).

Start of the bunch (PIC-D21). The run loads the bunch on top of a partner of the opposite charge at
rest, with E = B = 0, and the bunch moves from t = 0 on. Mode by mode the transverse potential
obeys a'' + w^2 a = J_T exp(-i W t), a(0) = a'(0) = 0, with w = |k| and W = k_x beta; the
longitudinal field follows the charge at once. In the mid-plane of the bunch, relative to the
steady field of the same mode,

    B_theta:  1 - [(1 + W/w) / 2] exp(-i (w - W) t) - [(1 - W/w) / 2] exp(i (w + W) t)
    E_r:      the transverse part beta^2 k_x^2 / k^2 of the steady field with W/w -> w/W in the
              bracket, the longitudinal part unchanged, and the field of the partner, which
              is -(1 / k^2) exp(i W t) against 1 / (k_perp^2 + k_x^2 / gamma^2).

The part of the start that travels with the bunch (small k_perp) loses its phase at the rate
k_perp^2 / (2 k_x) + k_x / (2 gamma^2): the field is there after a time of the order
2 gamma^2 sigma_x, not gamma r.

Linear wake (cold plasma of the density n0 = 1, ions at rest; E in m c omega_pe / e). The density
answers the charge of a driver and the ponderomotive force of a pulse point by point,

    (d^2/dt^2 + 1) dn = n_b                 (electron bunch of the density n_b),
    (d^2/dt^2 + 1) phi = <a^2> / 2          (pulse; <a^2> = a_env^2 / 2 for linear polarisation),

so behind the driver every point oscillates at omega_pe and the pattern has the wavenumber
k_w = 1 / v, v the speed of the driver (the group velocity of the pulse, beta of the bunch).

Pulse with the envelope a0 exp(-(t / tau)^2) exp(-r^2 / w^2) of the field (the `gaussian` envelope
of the antenna, tau = duration): behind it E = -grad phi with

    phi = A exp(-2 r^2 / w^2) sin(t - k_w x),    A = a0^2 (tau / 4) sqrt(pi / 2) exp(-tau^2 / 8),
    E_x = k_w A exp(-2 r^2 / w^2) cos(t - k_w x),    E_r = (4 r / w^2) A exp(-2 r^2 / w^2) sin(t - k_w x),

largest for tau = 2. a0 and w are the local ones of the beam (`beam_at`): the wake of a pulse that
diffracts follows a0(x)^2 and w(x). a0 is the amplitude inside the plasma: a pulse that enters
through a ramp much longer than its wavelength keeps its energy flux, so a^2 = a0_vacuum^2 / eta
with the refractive index eta = sqrt(1 - 1 / omega^2) (`amplitude_in_plasma`). On the Yee grid the pulse travels with
v_g = (sin(k dx) / dx) / (sin(omega dt) / dt), sin^2(omega dt / 2) / dt^2 = sin^2(k dx / 2) / dx^2
+ 1 / 4 (the leapfrog of the particles gives the plasma term exactly; the shape of the particles
lowers it by a factor of the order 1 - (k dx)^2 / 6, which is not included).

Rigid Gaussian bunch of the charge q per particle, peak density n_b: with the fields
exp(i k (x - beta t)) in the dielectric 1 - 1 / omega^2 the longitudinal field is
-i k rho (1 / eps - beta^2) / (k_perp^2 + 1 + k^2 / gamma^2), and the pole of 1 / eps at
k = k_w = 1 / beta leaves behind the bunch

    E_x = -q n_b sqrt(2 pi) sigma_x exp(-k_w^2 sigma_x^2 / 2) R(r) cos(k_w xi),
    R(r) = k_w^2 int_0^inf r' exp(-r'^2 / (2 sigma_r^2)) I_0(k_w r_<) K_0(k_w r_>) dr',

R(0) = s exp(s) E_1(s) with s = (k_w sigma_r)^2 / 2. Limits: a wide bunch has
R = exp(-r^2 / (2 sigma_r^2)) and the wake of a sheet, whose field jumps by its charge per area;
a thin one has R = (k_w sigma_r)^2 K_0(k_w r).

Blowout (`blowout`; no gate, the numbers a run is reported against). A pulse with a0 >~ 2 and the
matched waist w = 2 sqrt(a0) expels the electrons from a nearly spherical cavity (Lu et al.,
Phys. Rev. ST Accel. Beams 10, 061301 (2007), Eqs. 1, 4, 6 and the text around them):
radius R = 2 sqrt(a0); on the axis E_x = (x - x_c) / 2 through the centre x_c of the cavity, so
its largest value is sqrt(a0); the wake follows the pulse at v_g minus the speed 1 / omega0^2 at
which the front of the pulse is etched away, 1 - 3 / (2 omega0^2); an electron stays in phase over
the length (2 / 3) omega0^2 R and gains (2 / 3) omega0^2 a0 there.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import constants
from scipy.integrate import quad
from scipy.special import i0e, j1, k0e


def bunch_charge(density: float, sigma_x: float, sigma_r: float) -> float:
    """Charge (in e n0 (c/omega_pe)^3, positive number) of a Gaussian bunch with the peak density."""
    return density * (2.0 * math.pi) ** 1.5 * sigma_x * sigma_r ** 2


def bunch_field(r, sigma_x: float, sigma_r: float, gamma: float, charge: float = 1.0):
    """E_r in the mid-plane of the steady Gaussian bunch at the radii r (B_theta = beta E_r)."""
    def integrand(s, rr):
        a = 1.0 + 2.0 * (s * sigma_r) ** 2
        return s * s * math.exp(-(s * rr) ** 2 / a) / (a * a * math.sqrt(1.0 + 2.0 * (s * gamma * sigma_x) ** 2))

    out = [gamma * charge / math.pi ** 1.5 * rr * quad(integrand, 0.0, math.inf, args=(rr,), epsabs=0.0,
                                                    epsrel=1e-10, limit=400)[0] for rr in np.atleast_1d(r)]
    return np.asarray(out)


def bunch_startup(r, t: float, sigma_x: float, sigma_r: float, gamma: float, nk: int = 1500, nx: int = 1201):
    """(E_r, B_theta) in the mid-plane at the time t after the start on the partner, each relative
    to its steady value, at the radii r. Trapezoid rule over k_perp <= 8 / sigma_r and
    |k_x| <= 7 / sigma_x."""
    beta = math.sqrt(1.0 - 1.0 / gamma ** 2)
    kp = np.linspace(1e-4 / sigma_r, 8.0 / sigma_r, nk)[:, None]
    kx = np.linspace(-7.0 / sigma_x, 7.0 / sigma_x, nx)[None, :]
    k2 = kp ** 2 + kx ** 2
    w, W = np.sqrt(k2), kx * beta
    steady = 1.0 / (kp ** 2 + (kx / gamma) ** 2)
    near, far = np.exp(-1j * (w - W) * t), np.exp(1j * (w + W) * t)
    b = steady * (1.0 - 0.5 * (1.0 + W / w) * near - 0.5 * (1.0 - W / w) * far)
    # beta^2 kx^2 / k^2 times the bracket with w / W, written without the division by W
    et = steady * beta / k2 * (0.5 * (W * kx + w * kx) * near + 0.5 * (W * kx - w * kx) * far)
    e = steady - et - np.exp(1j * W * t) / k2
    rho = np.exp(-0.5 * (kp * sigma_r) ** 2 - 0.5 * (kx * sigma_x) ** 2)

    def mid(f, rr):
        return np.trapezoid(np.trapezoid(kp * kp * j1(kp * rr) * rho * f, kx[0], axis=1), kp[:, 0]).real

    out = [(mid(e, rr) / mid(steady, rr), mid(b, rr) / mid(steady, rr)) for rr in np.atleast_1d(r)]
    return np.asarray(out).T


def beam_at(d: float, a0: float, waist: float, rayleigh: float):
    """(a0, waist) of a Gaussian beam at the distance d from its focus."""
    f = math.sqrt(1.0 + (d / rayleigh) ** 2)
    return a0 / f, waist * f


def yee_group_velocity(omega: float, dx: float, dt: float, omega_p: float = 1.0) -> float:
    """Group velocity along x of a wave of the frequency omega on the Yee grid in a cold plasma."""
    s2 = (dx / dt * math.sin(0.5 * omega * dt)) ** 2 - (0.5 * omega_p * dx) ** 2
    if not 0.0 < s2 < 1.0:
        raise ValueError(f"the grid carries no wave at omega = {omega} in a plasma with omega_p = {omega_p}")
    k = 2.0 / dx * math.asin(math.sqrt(s2))
    return (math.sin(k * dx) / dx) / (math.sin(omega * dt) / dt)


def laser_wake(a0: float, duration: float, waist: float, r: float = 0.0, speed: float = 1.0):
    """Amplitudes (E_x, E_r) at the radius r of the linear wake behind a linearly polarised pulse
    with the field envelope a0 exp(-(t / duration)^2) exp(-r^2 / waist^2) and the group velocity
    ``speed``, a0 inside the plasma: E = -grad phi, E_x = phi / speed and E_r = (4 r / waist^2) phi,
    a quarter of a period apart."""
    phi = a0 * a0 * 0.25 * duration * math.sqrt(0.5 * math.pi) * math.exp(-duration ** 2 / 8.0)
    phi *= math.exp(-2.0 * (r / waist) ** 2)
    return phi / speed, phi * 4.0 * r / waist ** 2


def amplitude_in_plasma(a0: float, omega: float, density: float = 1.0) -> float:
    """a0 of a pulse inside a plasma of the density ``density`` that it entered from vacuum through
    a long ramp: the energy flux a^2 omega k is kept."""
    return a0 * (1.0 - density / omega ** 2) ** -0.25


def blowout(a0: float, omega: float) -> dict:
    """Scalings of the blowout regime for a pulse of the matched waist (Lu et al. 2007)."""
    radius = 2.0 * math.sqrt(a0)
    return {"radius": radius, "slope": 0.5, "Ex_max": math.sqrt(a0), "wake_speed": 1.0 - 1.5 / omega ** 2,
            "dephasing_length": (2.0 / 3.0) * omega ** 2 * radius, "energy_gain": (2.0 / 3.0) * omega ** 2 * a0}


def electrons_per_unit(omega: float, wavelength: float = 0.8e-6) -> float:
    """Electrons in n0 (c / omega_pe)^3 for a laser of the given vacuum wavelength in metres with
    omega0 = omega omega_pe."""
    n0 = constants.epsilon_0 * constants.m_e * (2.0 * math.pi * constants.c / wavelength) ** 2 \
        / (constants.e ** 2 * omega ** 2)
    return n0 * (omega * wavelength / (2.0 * math.pi)) ** 3


def beam_wake_profile(r, sigma_r: float, kw: float = 1.0):
    """R(r): the wake of a bunch with the radial density exp(-r^2 / (2 sigma_r^2)) relative to the
    wake of a sheet; ``kw`` = 1 / beta."""
    def inner(x, rr):  # I_0(kw x) K_0(kw rr), x <= rr
        return x * math.exp(-0.5 * (x / sigma_r) ** 2 + kw * (x - rr)) * i0e(kw * x) * k0e(kw * rr)

    def outer(x, rr):  # I_0(kw rr) K_0(kw x), x >= rr
        return x * math.exp(-0.5 * (x / sigma_r) ** 2 + kw * (rr - x)) * i0e(kw * rr) * k0e(kw * x)

    out = []
    for rr in np.atleast_1d(r):
        a = quad(inner, 0.0, rr, args=(rr,), epsabs=0.0, epsrel=1e-10, limit=400)[0] if rr > 0.0 else 0.0
        b = quad(outer, rr, rr + 40.0 * sigma_r, args=(rr,), epsabs=0.0, epsrel=1e-10, limit=400)[0]
        out.append(kw * kw * (a + b))
    return np.asarray(out)


def beam_wake(density: float, sigma_x: float, sigma_r: float, gamma: float, r=0.0):
    """Amplitude of E_x at the radii r of the linear wake behind a rigid Gaussian bunch of the peak
    density ``density``; the pattern is cos(k_w xi) with k_w = 1 / beta and has the sign of -q
    right behind the bunch."""
    kw = 1.0 / math.sqrt(1.0 - 1.0 / gamma ** 2)
    return (density * math.sqrt(2.0 * math.pi) * sigma_x * math.exp(-0.5 * (kw * sigma_x) ** 2)
            * beam_wake_profile(r, sigma_r, kw))
