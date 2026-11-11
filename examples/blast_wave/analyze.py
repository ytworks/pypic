"""Analyse a blast-wave run.

python examples/blast_wave/analyze.py output/blast_wave/n1

Null test (a run without debris; criteria in docs/plans/2026-11-11-blast-wave.md, Phase 0): the
ambient plasma at rest must stay as loaded. Samples are the 3-D output; regions are the interior and,
for each face, the two block layers at the face and the next two. With a body, the blocks within two
dipole radii of its centre or of the axis through it along B0 are left out. Limits: Gauss law
<= 1e-10 and div B <= 1e-12; density of each species within 2.5 % of step 0; T_e <= 1.05 of step 0;
mean velocity < 0.02 v_te; with a body, its potential against the plasma at two dipole radii
<= 5 T_e / e in both slice planes. The charge of the body and its change over the last third are
reported without a limit (the plan, Phase 0, Amendment).
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from pypic.analysis import open_run

GAUSS_MAX, DIVB_MAX = 1e-10, 1e-12
DENSITY_TOL, TE_MAX, VELOCITY_MAX = 0.025, 1.05, 0.02  # against step 0; velocity in v_te
POTENTIAL_MAX = 5.0  # in T_e / e
SLAB, INNER = 2, 4  # block layers of a face slab; the interior starts behind two slabs
EXCLUDE = 2.0  # dipole radii around the body and its axis along B0 that the regions leave out
LAST = 2.0 / 3.0  # the last third of the run starts here


def _finite_max(a) -> float:
    """Largest finite value; nan when nothing was measured (which fails a limit)."""
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    return float(a.max()) if a.size else float("nan")


def regions(nb) -> dict:
    """Block indices per axis of the interior and of the two slabs at each face."""
    full = [np.arange(n) for n in nb]
    reg = {"interior": tuple(np.arange(INNER, n - INNER) for n in nb)}
    for a, ax in enumerate("xyz"):
        for side in "-+":
            for k, name in enumerate(("face", "next")):
                lo = k * SLAB
                sel = list(full)
                sel[a] = np.arange(lo, lo + SLAB) if side == "-" else np.arange(nb[a] - lo - SLAB, nb[a] - lo)
                reg[f"{ax}{side} {name}"] = tuple(sel)
    return reg


def region_history(run, info) -> dict:
    """Density, temperature (both against step 0) and mean speed (in v_te) of each species in each region
    at every 3-D output. T = (2/3) (sum ekin / sum n - m |sum flux / sum n|^2 / 2)."""
    f = run.fields("coarse")
    its = f.iterations
    first = f.get(its[0], "n:e")
    keep = np.ones(first.data.shape, dtype=bool)
    if info["body"]:
        x, y, z = (c - c0 for c, c0 in zip(np.meshgrid(*first.coords, indexing="ij"), info["center"]))
        r_ex = EXCLUDE * info["dipole_radius"]
        keep = (x * x + y * y >= r_ex * r_ex) & (x * x + y * y + z * z >= r_ex * r_ex)
    reg = {name: np.ix_(*sel) for name, sel in regions(first.data.shape).items()}
    vte = info["params"]["vte"]
    out = {"regions": list(reg), "steps": [int(it) for it in its],
           "blocks": {name: int(keep[ix].sum()) for name, ix in reg.items()}}
    for sp in ("e", "i"):
        m = run.species_info[sp]["mass"]
        n, t, v = (np.zeros((len(reg), len(its))) for _ in range(3))
        for j, it in enumerate(its):
            d = {q: f.get(it, f"{q}:{sp}").data.astype(float) * keep for q in ("n", "fluxx", "fluxy", "fluxz", "ekin")}
            for i, ix in enumerate(reg.values()):
                with np.errstate(divide="ignore", invalid="ignore"):  # a region without a block gives nan
                    n[i, j] = d["n"][ix].sum()
                    u2 = sum(d[q][ix].sum() ** 2 for q in ("fluxx", "fluxy", "fluxz")) / n[i, j] ** 2
                    t[i, j] = (2.0 / 3.0) * (d["ekin"][ix].sum() / n[i, j] - 0.5 * m * u2)
                    v[i, j] = np.sqrt(u2) / vte
        with np.errstate(divide="ignore", invalid="ignore"):
            out[sp] = {"density": n / n[:, :1], "temperature": t / t[:, :1], "speed": v}
    return out


def body_potential(run, info, name) -> float:
    """Potential of the body against the plasma at two dipole radii, in T_e / e, from the slice ``name``:
    the integral of the ring mean of E_r from the surface of the body outwards, E the mean of the
    slice samples of the last third of the run. Each in-plane component is projected at its own
    sample positions. nan when a ring holds no sample."""
    s = run.series(name)
    a, b = (i for i in range(3) if "xyz"[i] != s.layout["axis"])
    ca, cb = info["center"][a], info["center"][b]
    step = s.step
    sel = slice(int(np.nonzero(step >= LAST * step[-1])[0][0]), None)
    r0, r1 = info["body_radius"], EXCLUDE * info["dipole_radius"]
    edges = np.linspace(r0, r1, max(1, int(round((r1 - r0) / run.dx[a]))) + 1)
    e_r = np.zeros(len(edges) - 1)
    for k, ax in enumerate((a, b)):
        comp = "E" + "xyz"[ax]
        pa, pb = s.plane_axes(comp)
        da, db = np.meshgrid(pa - ca, pb - cb, indexing="ij")
        r = np.hypot(da, db)
        e = s.read(comp, sel).astype(float).mean(axis=0)
        ring = np.digitize(r, edges) - 1
        for j in range(len(e_r)):
            m = ring == j
            e_r[j] += (e[m] * (da, db)[k][m] / r[m]).mean() if m.any() else np.nan
    return float(np.sum(e_r * np.diff(edges)) / info["params"]["vte"] ** 2)


def null_test(run_dir: str) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    if info["debris"]:
        raise ValueError("the null test needs a run without debris")
    s = run.scalars
    h = region_history(run, info)
    names = h["regions"]

    def worst(a):  # the largest value over regions and samples and where it is; a nan wins
        a = np.where(np.isfinite(a), a, np.inf)
        i, j = np.unravel_index(int(np.argmax(a)), a.shape)
        return {"value": float(a[i, j]) if np.isfinite(a[i, j]) else float("nan"), "region": names[i],
                "step": h["steps"][j]}

    res = {"gauss_max": _finite_max(s["gauss"]), "divB_max": _finite_max(s["divB"]),
           "blocks": h["blocks"], "samples": len(h["steps"]),
           "density": {sp: worst(np.abs(h[sp]["density"] - 1.0)) for sp in "ei"},
           "temperature": {sp: worst(h[sp]["temperature"]) for sp in "ei"},
           "temperature_min": {sp: float(np.min(h[sp]["temperature"])) for sp in "ei"},
           "speed": {sp: worst(h[sp]["speed"]) for sp in "ei"},
           "interior_end": {sp: {q: float(h[sp][q][0, -1]) for q in ("density", "temperature", "speed")}
                            for sp in "ei"},
           "fill_max": float(np.max(s["fill"])), "particles_start": float(s["n_particles"][0]),
           "particles_end": float(s["n_particles"][-1])}
    stats = run.meta.get("run")
    if stats:
        res.update(ms_per_step=1e3 * stats["wall_time_s"] / stats["steps"],
                   ns_per_particle_step=stats["ns_per_particle_step"])
    checks = {
        "gauss": bool(res["gauss_max"] <= GAUSS_MAX),  # False for nan: no measurement, no pass
        "divB": bool(res["divB_max"] <= DIVB_MAX),
        "density": bool(all(res["density"][sp]["value"] <= DENSITY_TOL for sp in "ei")),
        "temperature": bool(res["temperature"]["e"]["value"] <= TE_MAX),
        "velocity": bool(all(res["speed"][sp]["value"] < VELOCITY_MAX for sp in "ei")),
    }
    if info["body"]:
        q, step = s["absorbed_charge"].astype(float), s["step"]
        i0 = int(np.nonzero(step >= LAST * step[-1])[0][0])
        q_max = float(np.max(np.abs(q)))
        res.update(potential={n: body_potential(run, info, n) for n in ("slice_z", "slice_y")},
                   body_charge=float(q[-1]), body_charge_max=q_max, body_charge_change=float(abs(q[-1] - q[i0])),
                   n_absorbed=int(s["n_absorbed"][-1]))
        checks["potential"] = bool(all(abs(p) <= POTENTIAL_MAX for p in res["potential"].values()))
    res.update(checks=checks, passed=all(checks.values()))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    a = ap.parse_args()
    print(json.dumps(null_test(a.run_dir), indent=2, default=float))
