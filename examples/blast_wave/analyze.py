"""Analyse a blast-wave run.

python examples/blast_wave/analyze.py output/blast_wave/gate [--movie]

With debris: the gate of docs/plans/2026-11-11-blast-wave.md, Phase 3 with its amendment. Measured in
the plane z = centre (perpendicular to B0), B_z of the total field.
- Front r_f on eight rays (the axes and the diagonals): the outermost radius at which B_z / B0,
  averaged over one d_i along the ray, exceeds 1.25; 0 on a ray without one. Window: the slice
  samples with the mean r_f between 0.3 and 0.7 R_B.
- gate: Gauss law <= 1e-10 and div B <= 1e-12; the speed of the front (fit over the window) above
  the fast magnetosonic speed; the largest ring mean of B_z / B0 in the window above 1.5; the mean
  r_f of the axes within 5 % of that of the diagonals at every sample of the window; the largest
  extent D98 of the debris ions (Eq. 3 of Le et al. 2021 on the ring mean of their density) within
  25 % of Eq. 4 with D = 3, and not at the last sample; the ambient densities and |B| in the 3-D
  blocks ahead of the front (more than 0.2 R_B beyond the largest r_f so far, outside the flux tube
  of the body) within 5 % of step 0 at every 3-D output. A quantity without a measurement is nan
  and fails its criterion.
- reported: the magnetic cavity, the compression of a perpendicular MHD shock, the distance of the
  front and the debris, the extent of the debris along and across B0, the potential and the charge
  of the body, the energy budget, the fill.

Without debris: the null test (Phase 0): the ambient plasma at rest must stay as loaded. Samples are
the 3-D output; regions are the interior and, for each face, the two block layers at the face and
the next two. With a body, the blocks within two dipole radii of its centre or of the axis through
it along B0 are left out. Limits: Gauss law <= 1e-10 and div B <= 1e-12; density of each species
within 2.5 % of step 0; T_e <= 1.05 of step 0; mean velocity < 0.02 v_te; with a body, its potential
against the plasma at two dipole radii <= 5 T_e / e in both slice planes. The charge of the body and
its change over the last third are reported without a limit (the plan, Phase 0, Amendment).
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
from pypic.analysis.fields import front_position  # noqa: E402
from pypic.setup import mhd  # noqa: E402

GAUSS_MAX, DIVB_MAX = 1e-10, 1e-12
# gate with debris
FRONT_LEVEL, CAVITY_LEVEL, COMPRESSION_MIN = 1.25, 0.5, 1.5  # B_z / B0
WINDOW, AHEAD = (0.3, 0.7), 0.2  # mean front radius of the window; margin ahead of the front; in R_B
AMBIENT_TOL, SYMMETRY_TOL = 0.05, 0.05
EXTENT_FRACTION, EXTENT_TOL = 0.98, 0.25  # Eq. 3 and the band around Eq. 4 of Le et al. 2021
DEBRIS_EDGE = 0.5  # debris-ion density of the reported edge of the debris, in the ambient density
ANGLES = np.arange(8) * np.pi / 4.0  # rays in the plane: even = axes, odd = diagonals
# null test
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


def _sample_spacing(a, dx) -> float:
    """Distance between the samples a of a slice along one axis: a whole number of cells of the width
    dx (the block size the slice was written with). A ring this wide around any point of the plane
    holds a sample; a narrower one can be empty."""
    return max(1, int(round((a[1] - a[0]) / dx))) * dx


def body_potential(run, info, name) -> float:
    """Potential of the body against the plasma at two dipole radii, in T_e / e, from the slice ``name``:
    the integral of the ring mean of E_r from the surface of the body outwards, E the mean of the
    slice samples of the last third of the run. Each in-plane component is projected at its own
    sample positions. The rings are at least as wide as the samples are apart."""
    s = run.series(name)
    a, b = (i for i in range(3) if "xyz"[i] != s.layout["axis"])
    ca, cb = info["center"][a], info["center"][b]
    step = s.step
    sel = slice(int(np.nonzero(step >= LAST * step[-1])[0][0]), None)
    r0, r1 = info["body_radius"], EXCLUDE * info["dipole_radius"]
    dr = _sample_spacing(s.plane_axes("E" + "xyz"[a])[0], run.dx[a])
    edges = np.linspace(r0, r1, max(1, int((r1 - r0) / dr)) + 1)
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


def _running_mean(a, w):
    return np.convolve(np.asarray(a, dtype=float), np.ones(w) / w, "valid")


def bilinear(data, a, b, pa, pb):
    """data[i, j] on the evenly spaced axes a, b, interpolated to the points (pa, pb) of the plane."""
    fa = np.clip((np.asarray(pa, dtype=float) - a[0]) / (a[1] - a[0]), 0.0, len(a) - 1.0)
    fb = np.clip((np.asarray(pb, dtype=float) - b[0]) / (b[1] - b[0]), 0.0, len(b) - 1.0)
    i, j = np.minimum(fa.astype(int), len(a) - 2), np.minimum(fb.astype(int), len(b) - 2)
    ta, tb = fa - i, fb - j
    return ((1 - ta) * (1 - tb) * data[i, j] + ta * (1 - tb) * data[i + 1, j]
            + (1 - ta) * tb * data[i, j + 1] + ta * tb * data[i + 1, j + 1])


def _plane_radius(a, b, center) -> float:
    """Largest radius around center that the plane holds on the axes."""
    return float(min(center[0] - a[0], a[-1] - center[0], center[1] - b[0], b[-1] - center[1]))


def ray_fronts(field, a, b, center, r0, dr, window, level) -> np.ndarray:
    """Front radius on the eight rays from center: the outermost radius at which the running mean of
    ``field`` over ``window`` samples (one every dr, from r0 to the radius the plane holds) reaches
    ``level``; 0 on a ray without one."""
    r = np.arange(r0, _plane_radius(a, b, center) + 1e-9, dr)
    rs = _running_mean(r, window)
    out = np.zeros(len(ANGLES))
    for k, th in enumerate(ANGLES):
        line = bilinear(field, a, b, center[0] + r * np.cos(th), center[1] + r * np.sin(th))
        rf = front_position(rs, _running_mean(line, window), level)
        out[k] = rf if np.isfinite(rf) else 0.0
    return out


def ring_mean(field, a, b, center, dr):
    """Mean of field[i, j] over the rings of width dr around center, out to the radius the plane holds
    on the axes: (ring centres, means); nan for a ring without a sample."""
    r = np.hypot(*np.meshgrid(a - center[0], b - center[1], indexing="ij"))
    n = int(_plane_radius(a, b, center) / dr)
    k = (r / dr).astype(int)
    m = k < n
    count = np.bincount(k[m], minlength=n)
    total = np.bincount(k[m], weights=np.asarray(field, dtype=float)[m], minlength=n)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (np.arange(n) + 0.5) * dr, total / count


def extent(r, profile, fraction=EXTENT_FRACTION) -> float:
    """Radius that holds ``fraction`` of the line integral of the profile from the centre (Eq. 3 of Le et
    al. 2021); the profile is given on rings of equal width with the centres r. nan without a
    measurement: a ring without a value, or no density at all."""
    dr = r[1] - r[0]
    c = np.concatenate(([0.0], np.cumsum(np.asarray(profile, dtype=float)))) * dr
    if not np.isfinite(c[-1]) or c[-1] <= 0.0:
        return float("nan")
    return float(np.interp(fraction * c[-1], c, np.concatenate((r - 0.5 * dr, [r[-1] + 0.5 * dr]))))


def plane_history(run, info) -> dict:
    """Front, compression, cavity and debris extent at every sample of the slice z = centre. The rings
    and the samples on the rays are as far apart as the samples of the slice."""
    sz = run.series("slice_z")
    c, b0 = info["center"][:2], info["b0"]
    r0 = info["body_radius"] if info["body"] else 0.0
    a, b = sz.plane_axes("Bz")
    na, nb = sz.plane_axes("n:di")
    dr = _sample_spacing(a, run.dx[0])
    w = max(1, int(round(info["d_i"] / dr)))
    step = np.asarray(sz.step)
    fronts = np.zeros((len(step), len(ANGLES)))
    peak, cavity, d98, edge = (np.full(len(step), np.nan) for _ in range(4))
    for j in range(len(step)):
        bz = np.asarray(sz.read("Bz", j), dtype=float) / b0
        fronts[j] = ray_fronts(bz, a, b, c, r0, dr, w, FRONT_LEVEL)
        r, ring = ring_mean(bz, a, b, c, dr)
        peak[j] = _finite_max(ring[r >= r0])
        low = np.nonzero((r >= r0) & (r < fronts[j].mean()) & (ring < CAVITY_LEVEL))[0]
        cavity[j] = r[low[-1]] if low.size else np.nan
        rn, nd = ring_mean(sz.read("n:di", j), na, nb, c, dr)
        d98[j] = extent(rn, nd)
        edge[j] = front_position(rn, nd, DEBRIS_EDGE)
    return {"step": step, "t": step * run.dt, "fronts": fronts, "front": fronts.mean(axis=1), "peak": peak,
            "cavity": cavity, "extent": d98, "edge": edge}


def window_measures(t, fronts, peak, r_b) -> dict:
    """Speed of the front, compression and asymmetry over the window (the samples with the mean front
    radius between WINDOW[0] and WINDOW[1] R_B); nan where the window holds too few samples."""
    front = fronts.mean(axis=1)
    win = np.nonzero((front >= WINDOW[0] * r_b) & (front <= WINDOW[1] * r_b))[0]
    nan = float("nan")
    out = {"window": win, "front_speed": nan, "compression": nan, "asymmetry": nan}
    if win.size:
        f = fronts[win]
        out["compression"] = _finite_max(peak[win])
        out["asymmetry"] = float(np.max(np.abs(f[:, ::2].mean(axis=1) - f[:, 1::2].mean(axis=1)) / front[win]))
    if win.size >= 2:
        out["front_speed"] = float(np.polyfit(t[win], front[win], 1)[0])
    return out


def ahead_of_front(run, info, step, front_max) -> dict:
    """Ambient densities and |B| in the 3-D blocks ahead of the front, against step 0 on the same
    blocks, at every 3-D output. Ahead: farther from the centre than the largest front radius of the
    slice samples up to that output plus AHEAD R_B and, with a body, EXCLUDE dipole radii or more
    from the axis through it along B0. |B| is the magnitude of the block means of the total field.
    nan for an output without such a block."""
    f = run.fields("coarse")
    its = list(f.iterations)
    first = f.get(its[0], "n:e")
    x, y, z = (c - c0 for c, c0 in zip(np.meshgrid(*first.coords, indexing="ij"), info["center"]))
    r = np.sqrt(x * x + y * y + z * z)
    tube = np.ones(r.shape, dtype=bool)
    if info["body"]:
        tube = x * x + y * y >= (EXCLUDE * info["dipole_radius"]) ** 2

    def sample(it):
        d = {q: f.get(it, q).data.astype(float) for q in ("n:e", "n:i")}
        d["B"] = np.sqrt(sum(f.get(it, q).data.astype(float) ** 2 for q in ("Bx", "By", "Bz")))
        return d

    ref = sample(its[0])
    out = {"steps": [int(it) for it in its], "blocks": [], "radius": [], "n:e": [], "n:i": [], "B": []}
    for it in its:
        r_min = float(np.max(front_max[step <= it])) + AHEAD * info["r_b"]
        keep = tube & (r > r_min)
        d = sample(it)
        out["blocks"].append(int(keep.sum()))
        out["radius"].append(r_min)
        for q in ref:
            out[q].append(float(d[q][keep].sum() / ref[q][keep].sum()) if keep.any() else float("nan"))
    # np.max keeps a nan: an output without a block ahead of the front is not "no deviation"
    out["deviation"] = float(np.max(np.abs(np.array([out[q] for q in ref]) - 1.0)))
    return out


def gate_checks(res: dict) -> dict:
    """The criteria of the gate on the measured values; a nan fails its criterion."""
    return {
        "gauss": bool(res["gauss_max"] <= GAUSS_MAX),
        "divB": bool(res["divB_max"] <= DIVB_MAX),
        "ambient": bool(res["ambient_deviation"] <= AMBIENT_TOL),
        "front_speed": bool(res["front_speed"] > res["fast_speed"]),
        "compression": bool(res["compression"] > COMPRESSION_MIN),
        "symmetry": bool(res["asymmetry"] <= SYMMETRY_TOL),
        "debris_extent": bool(abs(res["debris_extent"] / res["debris_extent_le"] - 1.0) <= EXTENT_TOL
                              and res["debris_extent_step"] < res["last_step"]),
    }


def analyze(run_dir: str, show_plot: bool = True, fig_dir: str | None = None) -> dict:
    run = open_run(run_dir)
    info = run.meta["planner"]
    if not info["debris"]:
        raise ValueError("the gate needs a run with debris (the null test is for a run without)")
    s = run.scalars
    r_b, r_m, m_a = info["r_b"], info["r_m"], info["mach_alfven"]
    h = plane_history(run, info)
    wm = window_measures(h["t"], h["fronts"], h["peak"], r_b)
    win = wm["window"]
    ahead = ahead_of_front(run, info, h["step"], h["fronts"].max(axis=1))
    j_ext = int(np.argmax(np.where(np.isfinite(h["extent"]), h["extent"], -np.inf)))
    res = {
        "gauss_max": _finite_max(s["gauss"]), "divB_max": _finite_max(s["divB"]),
        "r_b": r_b, "r_m": r_m, "fast_speed": info["v_fast"], "last_step": int(h["step"][-1]),
        "window_samples": int(win.size), "front_speed": wm["front_speed"], "compression": wm["compression"],
        "asymmetry": wm["asymmetry"],
        # np.max keeps a nan: a sample without a measured extent fails the criterion
        "debris_extent": float(np.max(h["extent"])), "debris_extent_step": int(h["step"][j_ext]),
        "debris_extent_le": T.debris_extent(r_m, m_a), "debris_extent_start": float(h["extent"][0]),
        "ambient_deviation": ahead["deviation"],
        "ahead": {k: ahead[k] for k in ("steps", "blocks", "n:e", "n:i", "B")},
        # reported
        "front_mach_fast": wm["front_speed"] / info["v_fast"],
        "debris_extent_ratio": float(np.max(h["extent"])) / T.debris_extent(r_m, m_a),
        "front_end": float(h["front"][-1]), "front_max": float(h["fronts"].max()),
        "cavity_max": _finite_max(h["cavity"]),
        "cavity_static": 2.0 ** (1.0 / 3.0) * info["dipole_radius"] if info["body"] else float("nan"),
        "compression_mhd": mhd.compression(m_a, info["beta"]),
        "compression_mhd_front": mhd.compression(wm["front_speed"] / info["v_a"], info["beta"]),
    }
    if win.size:
        j = int(win[-1])
        res.update(front_debris_distance=float(h["front"][j] - h["edge"][j]), front_debris_step=int(h["step"][j]))
    f = run.fields("coarse")
    last = f.get(f.iterations[-1], "n:di")
    x, y, z = (c - c0 for c, c0 in zip(np.meshgrid(*last.coords, indexing="ij"), info["center"]))
    n = last.data.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        res.update(debris_rms_along=float(np.sqrt((n * z * z).sum() / n.sum())),
                   debris_rms_across=float(np.sqrt((n * (x * x + y * y)).sum() / (2.0 * n.sum()))))
    k0 = float(s["kinetic_de"][0] + s["kinetic_di"][0])
    res.update(
        energy_budget_change=float((s["budget"][-1] - s["budget"][0]) / k0),
        debris_ion_energy_left=float(s["kinetic_di"][-1] / s["kinetic_di"][0]),
        ambient_ion_energy_gain=float((s["kinetic_i"][-1] - s["kinetic_i"][0]) / k0),
        fill_max=float(np.max(s["fill"])), particles_start=float(s["n_particles"][0]),
        particles_end=float(s["n_particles"][-1]))
    if info["body"]:
        q = s["absorbed_charge"].astype(float)
        res.update(potential={name: body_potential(run, info, name) for name in ("slice_z", "slice_y")},
                   body_charge=float(q[-1]), body_charge_max=float(np.max(np.abs(q))),
                   n_absorbed=int(s["n_absorbed"][-1]))
    stats = run.meta.get("run")
    if stats:
        res.update(ms_per_step=1e3 * stats["wall_time_s"] / stats["steps"],
                   ns_per_particle_step=stats["ns_per_particle_step"])
    checks = gate_checks(res)
    res.update(checks=checks, passed=all(checks.values()))
    if show_plot:
        res["plot"] = _figure(run, info, h, res, int(win[-1]) if win.size else len(h["step"]) - 1,
                              os.path.join(fig_dir or run_dir, "blast_wave.png"))
    return res


def _mark_body(ax, info, a, b):
    if info["body"]:
        ax.add_patch(P.plt.Circle((info["center"][a], info["center"][b]), info["body_radius"],
                                  fc=P.INK["muted"], ec="none"))
    ax.set_aspect("equal")


def _figure(run, info, h, res, j, path):
    """B_z in the plane across B0 and the debris ions in a plane along B0 at slice sample j; radii
    against time."""
    P.style()
    c, b0 = info["center"], info["b0"]
    fig, ax = P.plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    sz, sy = run.series("slice_z"), run.series("slice_y")
    t_j = h["t"][j]
    bz = np.array(sz.read("Bz", j), dtype=float) / b0
    a, b = sz.plane_axes("Bz")
    if info["body"]:  # the dipole field in there would set the colour scale
        bz[np.hypot(a[:, None] - c[0], b[None, :] - c[1]) < info["body_radius"]] = np.nan
    P.field_map(ax[0], bz, sz.extent("Bz"), signed=False, vlim=(0.0, max(2.0, 1.1 * res["compression"])),
                label="B_z / B0", xlabel="x [c/ωpe]", ylabel="y [c/ωpe]")
    rf = h["fronts"][j]
    ax[0].plot(c[0] + rf * np.cos(ANGLES), c[1] + rf * np.sin(ANGLES), "o", ms=5, mfc="none",
               color=P.INK["primary"], label="front on the eight rays")
    ax[0].legend(loc="upper right", fontsize=7)
    ax[0].set_title(f"B_z, plane z = centre (across B0), t = {t_j:.0f}")
    _mark_body(ax[0], info, 0, 1)
    P.field_map(ax[1], np.asarray(sy.read("n:di", j), dtype=float), sy.extent("n:di"), signed=False,
                label="debris ions, n / n0", xlabel="x [c/ωpe]", ylabel="z [c/ωpe]")
    ax[1].add_patch(P.plt.Circle((c[0], c[2]), res["debris_extent_le"], fc="none", ec=P.INK["primary"], ls="--",
                                 lw=0.8))
    ax[1].set_title(f"debris-ion density, plane y = centre (B0 upwards), t = {t_j:.0f}\n"
                    "dashed: Eq. 4 of Le et al. 2021 with D = 3")
    _mark_body(ax[1], info, 0, 2)
    g, t = ax[2], h["t"]
    le = res["debris_extent_le"]
    g.axhspan((1 - EXTENT_TOL) * le, (1 + EXTENT_TOL) * le, color=P.INK["grid"], alpha=0.7, lw=0,
              label="gate: Eq. 4 of Le et al. 2021 (D = 3) ± 25 %")
    g.plot(t, h["fronts"][:, ::2].mean(axis=1), color=P.color(0), label="front, axes")
    g.plot(t, h["fronts"][:, 1::2].mean(axis=1), color=P.color(0), ls="--", label="front, diagonals")
    g.plot(t, h["extent"], color=P.color(1), label="debris ions, D98")
    g.plot(t, h["cavity"], color=P.color(2), label="magnetic cavity (B_z < 0.5 B0)")
    g.plot(t, info["shell_radius"] + info["v_d"] * t, **P.THEORY, label="free expansion at V_d")
    win = np.nonzero((h["front"] >= WINDOW[0] * res["r_b"]) & (h["front"] <= WINDOW[1] * res["r_b"]))[0]
    if win.size:
        t0, r0 = t[win[0]], h["front"][win[0]]
        g.plot(t, r0 + info["v_fast"] * (t - t0), color=P.INK["primary"], ls=":", lw=1.0,
               label="fast magnetosonic speed")
        g.axvspan(t[win[0]], t[win[-1]], color=P.color(0), alpha=0.08, lw=0, label="window of the front criteria")
    g.set_xlim(0.0, t[-1])
    g.set_ylim(0.0, 1.05 * max(res["r_b"], float(np.max(h["fronts"]))))
    g.set_xlabel("t [1/ωpe]")
    g.set_ylabel("radius in the plane z = centre [c/ωpe]")
    g.set_title("front, debris and cavity")
    g.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7)
    fig.tight_layout()
    return P.save(fig, path)


def movie(run_dir: str, out: str | None = None, fps: int = 12, n_frames: int = 100, dpi: int = 110) -> str:
    """B_z and the ambient ions in the plane across B0, the debris ions in a plane along it."""
    from pypic.analysis import Movie
    from pypic.analysis.animate import SlicePanel

    run = open_run(run_dir)
    info = run.meta["planner"]
    kw = dict(signed=False, blank=(info["center"], info["body_radius"]) if info["body"] else None)
    m = Movie(run, title="blast wave around a magnetised body")
    bz = dict(scale=1.0 / info["b0"], vlim=(0.0, 3.0), label="B_z / B0 (z = centre)")
    m.add(SlicePanel(run, "slice_z", "Bz", **bz, **kw))
    m.add(SlicePanel(run, "slice_z", "n:i", vlim=(0.0, 3.0), label="ambient ions, n / n0 (z = centre)", **kw))
    if info["debris"]:
        m.add(SlicePanel(run, "slice_y", "n:di", log=True, label="debris ions, n / n0 (y = centre)", **kw))
    return m.save(out or os.path.join(run_dir, "blast_wave.mp4"), fps=fps, n_frames=n_frames, dpi=dpi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--movie", action="store_true")
    a = ap.parse_args()
    with_debris = open_run(a.run_dir).meta["planner"]["debris"]
    r = analyze(a.run_dir) if with_debris else null_test(a.run_dir)
    if a.movie:
        r["movie"] = movie(a.run_dir)
    print(json.dumps(r, indent=2, default=float))
