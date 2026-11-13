"""Minimal encounter model of multiple connection (review of 2026-11-13, Major 4).

Run: .venv/bin/python research/topics/shock-alfvenic-turbulence/design/model.py > .../design/model.md

Picture (2-D, B0 along y in the plane, slab delta B_x(y) convected into a perpendicular shock at U1):
- A field line is x = x0 - U1 t + xi(y). With y periodic it is a periodic curve, so once it touches the front it meets
  it at least twice. The upstream arc between two intersections A and B shrinks to zero over t_conn = xi_pp / U1.
- An electron reflected at A (adiabatic SDA in the de Hoffmann-Teller frame, HTF, of A) runs along the upstream arc
  and must meet B before the arc vanishes, because A and B approach each other while the electron moves towards B.
  In a periodic slab the return is therefore certain; nature has it only where field lines come back.
- B has the opposite slope, so its HTF moves along the front opposite to A's. Seen from B's HTF the electron arrives
  with |v_par| larger by about 2 U1 / cos(theta) (frames differ by 2 U1 tan(theta) along y, projected on the field).
- At each encounter: reflection if the HTF pitch angle lies outside the loss cone, sin^2(alpha) > 1 / r_B;
  otherwise transmission downstream. v_perp is conserved upstream (magnetic moment).
- A transmitted electron keeps its HTF energy and its magnetic moment: v_perp2^2 = r_B v_perp^2. Its energy is
  counted in the downstream plasma frame, which moves along B2 in the HTF at U2 / cos(theta2), with
  tan(theta2) = r_B tan(theta) (normal component of B continuous, tangential component times r_B).

Neglected (named, to be checked against the runs): cross-shock potential; pitch-angle scattering (SSDA, whistlers);
ripples and corrugation; the finite thickness and time of each encounter; theta of A and B taken equal in |theta|.
The local-theta null of the Approach is the same model with one encounter: a reflected electron escapes upstream.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq


def compression(M_A, beta_tot, gamma=5.0 / 3.0):
    """Density (= B) compression of a perpendicular MHD shock. Units rho1 = B1 = 1 (v_A = 1), p1 = beta_tot / 2."""
    U1, p1 = M_A, beta_tot / 2.0

    def p2_of(r):  # momentum flux
        return U1 * U1 + p1 + 0.5 - U1 * U1 / r - 0.5 * r * r

    def energy(r):
        U2, p2 = U1 / r, p2_of(r)
        e1 = 0.5 * U1**3 + gamma / (gamma - 1) * p1 * U1 + U1
        e2 = 0.5 * U1 * U2 * U2 + gamma / (gamma - 1) * p2 * U2 + r * r * U2
        return e1 - e2

    return brentq(energy, 1.0001, 3.999)


def encounters(theta_rms_db, M_A, mr, beta_e, r_B, n=400000, n_max=20, seed=3):
    """Monte Carlo of the encounter chain. Returns, for single (null) and multiple connection:
    fractions reflected at the first encounter, the distribution of the number of encounters, and the fraction of
    transmitted electrons with HTF energy above thresholds (in units of T_e upstream)."""
    rng = np.random.default_rng(seed)
    vte = math.sqrt(beta_e * mr / 2.0)  # sqrt(T/m) in v_A
    bx = rng.normal(0.0, theta_rms_db, n)
    cos_t = np.clip(np.abs(bx) / np.sqrt(1.0 + bx * bx), 1e-6, 1.0)
    V = M_A / cos_t  # HTF inflow speed along the field
    w = rng.normal(0.0, vte, (n, 3))  # thermal velocity; w[:,0] along B
    vpar = -V + w[:, 0]  # towards the shock if negative
    vperp2 = w[:, 1] ** 2 + w[:, 2] ** 2
    incoming = vpar < 0
    loss = 1.0 / r_B

    def reflects(vp):
        return vperp2 / (vp * vp + vperp2) > loss

    U2 = M_A / r_B
    sin_t = np.sqrt(1.0 - cos_t * cos_t)
    cos2 = cos_t / np.sqrt(cos_t * cos_t + (r_B * sin_t) ** 2)
    V2 = U2 / cos2  # downstream flow speed along B2 in the HTF

    def e_down(vp, sel):
        """Kinetic energy (m = 1) in the downstream plasma frame of an electron transmitted with HTF |v_par| = vp."""
        vperp2_2 = r_B * vperp2[sel]
        vpar2 = np.sqrt(np.maximum(vp[sel] ** 2 + vperp2[sel] - vperp2_2, 0.0))
        return 0.5 * ((vpar2 - V2[sel]) ** 2 + vperp2_2)

    n_enc = np.zeros(n, dtype=int)
    energy = np.full(n, np.nan)  # HTF energy at transmission, in m v^2 / 2 units
    active = incoming.copy()
    vp = np.abs(vpar)
    refl_first = incoming & reflects(vp)
    for _ in range(n_max):
        n_enc[active] += 1
        ref = active & reflects(vp)
        trans = active & ~ref
        energy[trans] = e_down(vp, trans)
        active = ref
        vp = np.where(ref, vp + 2.0 * M_A / cos_t, vp)  # arrival speed seen from the next intersection's HTF
        if not active.any():
            break
    # null: one encounter; reflected electrons escape upstream (no energy recorded downstream)
    e_null = np.full(n, np.nan)
    tr0 = incoming & ~refl_first
    e_null[tr0] = e_down(np.abs(vpar), tr0)
    Te = vte * vte  # T_e = m v_te^2 with m = 1
    out = {"f_in": incoming.mean(), "f_ref1": refl_first[incoming].mean()}
    for name, e in (("null", e_null), ("multi", energy)):
        ok = ~np.isnan(e)
        for thr in (3, 10, 30):
            out[f"{name}_gt{thr}"] = float(np.mean(e[ok] / Te > thr) * ok.mean()) if ok.any() else 0.0
    counts = np.bincount(n_enc[incoming], minlength=5)[1:5] / incoming.sum()
    out["n_enc"] = counts
    return out


def main():
    beta_e = beta_i = 1.0
    M_A = 4.5
    r_B = compression(M_A, beta_e + beta_i)
    print("# Minimal encounter model (generated by `design/model.py`)\n")
    print(f"M_A = {M_A}, beta_e = beta_i = {beta_e}: perpendicular MHD compression r = B2/B1 = {r_B:.3f}; "
          f"loss cone sin^2(alpha) < 1/r = {1 / r_B:.3f}. Assumptions and neglected effects: module docstring.\n")
    print("Energies are those of electrons transmitted downstream, in the downstream plasma frame, in units of the "
          "upstream T_e (fraction of all incoming electrons above the threshold). "
          "'null': one encounter, reflected electrons lost upstream (the local-θ_Bn null). "
          "'multi': encounters repeat at the next intersection while the electron reflects (periodic slab).\n")
    print("| δB_x/B0 | m_i/m_e | reflected at 1st | N=1 | N=2 | N=3 | N≥4 | >3 T_e null / multi | "
          ">10 T_e null / multi | >30 T_e null / multi |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for db in (0.15, 0.3, 0.6, 1.0):
        for mr in (100, 400, 1836):
            o = encounters(db, M_A, mr, beta_e, r_B)
            c = o["n_enc"]
            ge4 = max(0.0, 1.0 - c[:3].sum())
            print(f"| {db} | {mr} | {o['f_ref1']:.3f} | {c[0]:.3f} | {c[1]:.3f} | {c[2]:.4f} | {ge4:.4f} | "
                  f"{o['null_gt3']:.2e} / {o['multi_gt3']:.2e} | {o['null_gt10']:.2e} / {o['multi_gt10']:.2e} | "
                  f"{o['null_gt30']:.2e} / {o['multi_gt30']:.2e} |")
    print()
    print("Reading:")
    print("- 'null' loses its reflected electrons upstream; 'multi' brings them back with a larger v_par. Compare the")
    print("  columns above 10 and 30 T_e: there the excess is a gain in energy, not only a gain in number.")
    print("- The model assumes the return; it has no Ξ in it. In a periodic slab the return is certain once the")
    print("  upstream arc lies outside the shock layer (Ξ ≳ 1, not modelled). Above that the excess should not")
    print("  depend on L_y.")
    print("  It depends on δB_x/B0 (local θ) and on m_i/m_e.")
    print("- Without scattering the chain stops after a few encounters: each return raises v_par while v_perp is")
    print("  fixed, so the electron falls into the loss cone. A power-law tail needs pitch-angle scattering.")


if __name__ == "__main__":
    main()
