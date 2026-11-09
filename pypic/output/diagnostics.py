"""GPU computation of the output diagnostics (OUT-D2). Every ``compute`` is collective: all
ranks call it at the same step; the result (NumPy) is returned on rank 0 only (None elsewhere).
"""
from __future__ import annotations

import math

import cupy as cp
import numpy as np

from ..backend import launch
from ..parallel.comm import gather_blocks, gather_columns
from ..parallel.decomp import Decomposition
from .h5 import round_mantissa
from .names import MOMENT_KIND, STAGGER, expand_fields, record_of, stagger_of  # noqa: F401
from .schedule import region_cells


class FieldSource:
    """Physical-domain arrays of the requested fields at the current state (computed once per
    output step and shared by all diagnostics of that step)."""

    def __init__(self, sim):
        self.sim = sim
        self.cache = {}

    def clear(self):
        self.cache = {}

    def get(self, comp: str) -> cp.ndarray:
        """Full local array (with ghosts). Deposited quantities are collective (ghost reduction)."""
        if comp in self.cache:
            return self.cache[comp]
        sim, g = self.sim, self.sim.grid
        if comp in STAGGER or comp in ("Jx", "Jy", "Jz"):
            arr = getattr(g, comp)
        elif comp == "rho":
            arr = cp.zeros(g.shape)
            for sp in sim.species:
                if not getattr(sp, "test", False):
                    launch(sim.k_rho, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.q), arr))
            g.reduce_ghosts((arr,))
        else:
            base, name = comp.split(":", 1)
            sp = sim.species_by_name(name)
            arr = cp.zeros(g.shape)
            launch(sim.k_moment, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), np.int32(MOMENT_KIND[base]),
                                        np.float64(sp.m), arr))
            g.reduce_ghosts((arr,))
        self.cache[comp] = arr
        return arr

    def phys(self, comp: str) -> cp.ndarray:
        return self.get(comp)[self.sim.grid.physical_slice()]


def _finish(arr: cp.ndarray, d) -> np.ndarray:
    """Cast (and optionally round the mantissa) on the GPU, then copy to the host."""
    a = arr.astype(cp.float32 if d.dtype == "float32" else cp.float64)
    if d.keep_bits is not None:
        a = round_mantissa(a, d.keep_bits, xp=cp)
    return cp.asnumpy(a)


def _finish_complex(arr: cp.ndarray, d) -> np.ndarray:
    re = _finish(arr.real, d)
    im = _finish(arr.imag, d)
    return (re + 1j * im).astype(np.complex64 if d.dtype == "float32" else np.complex128)


class Diag:
    """Base class: one configured diagnostic bound to a simulation."""

    def __init__(self, sim, d, resolved):
        self.sim, self.d, self.r = sim, d, resolved
        self.root = sim.rank == 0
        self.accum = None
        self.n_acc = 0
        self.t_acc = 0.0

    def time_offset(self, comp: str) -> float:
        """Time level of a quantity relative to the output time (in 1/omega_pe)."""
        dt = self.sim.cfg.dt
        if comp.startswith("J"):
            return -0.5 * dt  # current of the last step: t - dt/2
        if (comp.startswith("flux") or comp.startswith("ekin")) and not self.sim.implicit:
            return -0.5 * dt  # explicit leapfrog: u at t - dt/2
        return 0.0

    def species_list(self):
        names = self.d.species or tuple(s.name for s in self.sim.species)
        return [self.sim.species_by_name(n) for n in names]


# ---------------------------------------------------------------------- blocks
def _rank_phys(sim, r):
    """(phys_lo, n_phys) per axis of rank r (deterministic from the decomposition)."""
    if r == sim.rank:
        axes = sim.decomp.axes
    else:
        cfg = sim.cfg
        axes = Decomposition.build(cfg.n_cells, sim.decomp.dims, r, cfg.boundaries, cfg.npml).axes
    return [(a.phys_lo, a.n_phys) for a in axes]


def _coarse_layout(lo, n, i0, i1, c):
    """Intersection of a rank's cells [lo, lo+n) with the region [i0, i1) and the coarse cell
    range it touches: (a_lo, a_hi, c_lo, c_hi) (a_hi <= a_lo when empty)."""
    a_lo, a_hi = max(lo, i0), min(lo + n, i1)
    if a_hi <= a_lo:
        return a_lo, a_lo, 0, 0
    return a_lo, a_hi, (a_lo - i0) // c, -(-(a_hi - i0) // c)


def coarse_partial(block, layouts, phys_lo, region, coarsen, axes=(0, 1, 2)):
    """Block sums of the part of ``block`` (local physical array) inside the region, on the
    coarse cells it touches. Works for any alignment of ranks, region and coarsening."""
    sl, pads, shape = [], [], []
    for a in axes:
        a_lo, a_hi, c_lo, c_hi = layouts[a]
        c = coarsen[a]
        sl.append(slice(a_lo - phys_lo[a], a_hi - phys_lo[a]))
        before = (a_lo - region[a][0]) - c_lo * c
        after = (c_hi - c_lo) * c - before - (a_hi - a_lo)
        pads.append((before, after))
        shape.append(c_hi - c_lo)
    sub = block[tuple(sl)]
    if sub.size == 0:
        return cp.zeros(shape)
    if any(p != (0, 0) for p in pads):
        sub = cp.pad(sub, pads)
    rs = []
    for k, a in enumerate(axes):
        rs += [shape[k], coarsen[a]]
    return sub.reshape(rs).sum(axis=tuple(range(1, 2 * len(axes), 2)))


class BlockGather:
    """Coarse-grained block average of a 3-D (or a plane of a) field over a region, assembled
    on rank 0 from per-rank partial sums."""

    def __init__(self, sim, region_idx, coarsen, axes=(0, 1, 2), active=None):
        self.sim, self.region, self.coarsen, self.axes = sim, region_idx, tuple(coarsen), tuple(axes)
        self.active = list(active) if active is not None else [True] * sim.comm.size
        self.layouts = []
        for r in range(sim.comm.size):
            ph = _rank_phys(sim, r)
            if self.active[r]:
                lay = [_coarse_layout(ph[a][0], ph[a][1], region_idx[a][0], region_idx[a][1], self.coarsen[a])
                       for a in range(3)]
            else:
                lay = [(0, 0, 0, 0)] * 3
            self.layouts.append(lay)
        self.shape = tuple(-(-(region_idx[a][1] - region_idx[a][0]) // self.coarsen[a]) for a in self.axes)
        # number of fine cells in every coarse cell (edge cells may be partial)
        cnt = None
        for a in self.axes:
            i0, i1 = region_idx[a]
            c = self.coarsen[a]
            v = np.array([min(c, i1 - (i0 + k * c)) for k in range(self.shape[self.axes.index(a)])], dtype=float)
            cnt = v if cnt is None else np.multiply.outer(cnt, v)
        self.count = cp.asarray(cnt)

    def __call__(self, block) -> cp.ndarray | None:
        """block: local physical array restricted to self.axes (3-D or a plane)."""
        sim = self.sim
        phys_lo = [a.phys_lo for a in sim.decomp.axes]
        if self.active[sim.rank]:
            part = coarse_partial(block, self.layouts[sim.rank], phys_lo, self.region, self.coarsen, self.axes)
        else:
            part = cp.zeros(tuple(0 for _ in self.axes))
        shapes = [tuple(lay[a][3] - lay[a][2] for a in self.axes) for lay in self.layouts]
        blocks = gather_blocks(sim.comm, part, shapes)
        if blocks is None:
            return None
        out = cp.zeros(self.shape)
        for lay, b in zip(self.layouts, blocks):
            if b.size == 0:
                continue
            out[tuple(slice(lay[a][2], lay[a][3]) for a in self.axes)] += b
        return out / self.count

    def positions(self, stagger):
        """openPMD 'position' (in coarse cells) of the block-averaged samples."""
        return tuple(((self.coarsen[a] - 1) / 2 + stagger[a]) / self.coarsen[a] for a in self.axes)


# ------------------------------------------------------------------- kinds
class FieldsDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        self.comps = expand_fields(d.fields)
        self.region = region_cells(sim.cfg, d.region)
        for a in range(3):
            if sim.cfg.n_cells[a] < d.coarsen[a]:
                raise ValueError(f"diagnostic {d.name!r}: coarsen larger than the grid")
        self.bg = BlockGather(sim, self.region, d.coarsen)

    def compute(self, src: FieldSource):
        meshes = {}
        for comp in self.comps:
            arr = self.bg(src.phys(comp))
            if not self.root:
                continue
            rec, c = record_of(comp)
            m = meshes.setdefault(rec, {"components": {}, "position": {}, "time_offset": self.time_offset(comp)})
            m["components"][c] = _finish(arr, self.d)
            m["position"][c] = self.bg.positions(stagger_of(comp))
        if not self.root:
            return None
        dx = self.sim.cfg.dx
        grid = {"spacing": tuple(dx[a] * self.d.coarsen[a] for a in range(3)),
                "offset": tuple(self.region[a][0] * dx[a] for a in range(3))}
        return {"meshes": meshes, "grid": grid}


class SliceDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        cfg = sim.cfg
        self.comps = expand_fields(d.fields)
        self.a = "xyz".index(d.axis)
        self.others = tuple(b for b in range(3) if b != self.a)
        pos = d.position if d.position is not None else 0.5 * cfg.n_cells[self.a] * cfg.dx[self.a]
        reg = [(0, cfg.n_cells[b]) for b in range(3)]
        self.plane, self.bg = {}, {}
        for comp in self.comps:
            i = int(round(pos / cfg.dx[self.a] - stagger_of(comp)[self.a]))
            i = min(max(i, 0), cfg.n_cells[self.a] - 1)
            self.plane[comp] = i
            active = []
            for r in range(sim.comm.size):
                lo, n = _rank_phys(sim, r)[self.a]
                active.append(lo <= i < lo + n)
            self.bg[comp] = BlockGather(sim, reg, d.coarsen, axes=self.others, active=active)

    def plane_coords(self):
        cfg = self.sim.cfg
        return {c: (self.plane[c] + stagger_of(c)[self.a]) * cfg.dx[self.a] for c in self.comps}

    def compute(self, src: FieldSource):
        sim = self.sim
        ax = sim.decomp.axes[self.a]
        out = {}
        for comp in self.comps:
            bg = self.bg[comp]
            block = src.phys(comp)  # deposited quantities are collective: every rank computes them
            plane = None
            if bg.active[sim.rank]:
                sl = [slice(None)] * 3
                sl[self.a] = self.plane[comp] - ax.phys_lo
                plane = block[tuple(sl)]
            res = bg(plane)
            if self.root:
                out[comp] = _finish(res, self.d)
        return out if self.root else None


class LineoutDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        self.comps = expand_fields(d.fields)
        self.a = "xyz".index(d.axis)

    def _local(self, src):
        sim = self.sim
        ax = sim.decomp.axes[self.a]
        N = sim.cfg.n_cells[self.a]
        others = tuple(b for b in range(3) if b != self.a)
        vals = cp.zeros((len(self.comps), N))
        for k, comp in enumerate(self.comps):
            vals[k, ax.phys_lo: ax.phys_lo + ax.n_phys] = cp.sum(src.phys(comp), axis=others)
        return vals

    def accumulate(self, src):
        v = self._local(src)
        self.accum = v if self.accum is None else self.accum + v
        self.n_acc += 1
        self.t_acc += self.sim.time

    def compute(self, src: FieldSource):
        sim = self.sim
        if self.r.average and self.accum is not None:
            v, n, t = self.accum, self.n_acc, self.t_acc / self.n_acc
            self.accum, self.n_acc, self.t_acc = None, 0, 0.0
        else:
            v, n, t = self._local(src), 1, sim.time
        ntr = math.prod(sim.cfg.n_cells) // sim.cfg.n_cells[self.a]
        v = sim.comm.allreduce(v) / (n * ntr)
        if not self.root:
            return None
        return {"t_mean": t, "data": {c: _finish(v[k], self.d) for k, c in enumerate(self.comps)}}


class ModesDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        if not d.modes:
            raise ValueError(f"diagnostic {d.name!r}: modes list is empty")
        self.comps = expand_fields(d.fields)
        cfg = sim.cfg
        g = sim.grid
        self.L = g.L
        self.N = math.prod(cfg.n_cells)
        self.modes = [tuple(int(m) for m in mm) for mm in d.modes]
        # modes grouped by (my, mz): sum over y, z with a phase vector (plain sum for 0, 0), then one
        # matrix product with the x phases of every mode of the group (local physical cells only)
        self.groups = {}
        for comp in self.comps:
            st = stagger_of(comp)
            xs = []
            for a in range(3):
                ax = sim.decomp.axes[a]
                xs.append((np.arange(ax.phys_lo, ax.phys_lo + ax.n_phys) + st[a]) * cfg.dx[a])
            groups = {}
            for j, (mx, my, mz) in enumerate(self.modes):
                groups.setdefault((my, mz), []).append((j, mx))
            gl = []
            for (my, mz), lst in groups.items():
                py = None if my == 0 else cp.asarray(np.exp(-2j * np.pi * my * xs[1] / self.L[1]))
                pz = None if mz == 0 else cp.asarray(np.exp(-2j * np.pi * mz * xs[2] / self.L[2]))
                Mx = cp.asarray(np.exp(-2j * np.pi * np.outer([mx for _, mx in lst], xs[0]) / self.L[0]))
                gl.append((cp.asarray([j for j, _ in lst]), py, pz, Mx))
            self.groups[comp] = gl

    def _local(self, src):
        out = cp.zeros((len(self.comps), len(self.modes)), dtype=cp.complex128)
        for k, comp in enumerate(self.comps):
            F = src.phys(comp)
            for idx, py, pz, Mx in self.groups[comp]:
                A = F.sum(axis=2) if pz is None else cp.tensordot(F, pz, axes=([2], [0]))
                A = A.sum(axis=1) if py is None else cp.tensordot(A, py, axes=([1], [0]))
                out[k, idx] = Mx @ A.astype(cp.complex128)
        return out

    def accumulate(self, src):
        v = self._local(src)
        self.accum = v if self.accum is None else self.accum + v
        self.n_acc += 1
        self.t_acc += self.sim.time

    def compute(self, src):
        sim = self.sim
        if self.r.average and self.accum is not None:
            v, n, t = self.accum, self.n_acc, self.t_acc / self.n_acc
            self.accum, self.n_acc, self.t_acc = None, 0, 0.0
        else:
            v, n, t = self._local(src), 1, sim.time
        ri = sim.comm.allreduce(cp.stack([v.real, v.imag]))
        v = (ri[0] + 1j * ri[1]) / (n * self.N)
        if not self.root:
            return None
        return {"t_mean": t, "data": {c: _finish_complex(v[k], self.d) for k, c in enumerate(self.comps)}}


class KSpectrumDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        self.comps = expand_fields(d.fields)
        cfg = sim.cfg
        self.bg = BlockGather(sim, [(0, n) for n in cfg.n_cells], d.coarsen)
        shp = self.bg.shape
        h = [cfg.dx[a] * d.coarsen[a] for a in range(3)]
        k = [2 * np.pi * np.fft.fftfreq(shp[a], d=h[a]) for a in range(3)]
        k[2] = 2 * np.pi * np.fft.rfftfreq(shp[2], d=h[2])
        kk = np.sqrt(k[0][:, None, None] ** 2 + k[1][None, :, None] ** 2 + k[2][None, None, :] ** 2)
        self.dk = 2 * np.pi / max(sim.grid.L)
        nb = int(kk.max() / self.dk) + 1
        self.k_bins = (np.arange(nb) + 0.5) * self.dk
        if self.root:
            self.bin = cp.asarray(np.minimum((kk / self.dk).astype(np.int64), nb - 1).ravel())
            # rfft keeps half of the k_z modes: weight 2 for the duplicated ones
            wz = np.full(shp[2] // 2 + 1, 2.0)
            wz[0] = 1.0
            if shp[2] % 2 == 0:
                wz[-1] = 1.0
            self.wz = cp.asarray(np.broadcast_to(wz[None, None, :], kk.shape).ravel())
        self.nb = nb

    def compute(self, src):
        out = {}
        for comp in self.comps:
            arr = self.bg(src.phys(comp))
            if not self.root:
                continue
            F = cp.fft.rfftn(arr) / arr.size
            p = (cp.abs(F) ** 2).ravel() * self.wz
            spec = cp.bincount(self.bin, weights=p, minlength=self.nb) / self.dk
            out[comp] = _finish(spec, self.d)
        return out if self.root else None


QUANTITIES = ("x", "y", "z", "ux", "uy", "uz", "u", "gamma", "ekin")


def default_umax(cfg) -> float:
    """|u| range covering all beams (drift + 4 thermal widths), doubled."""
    umax = 0.0
    for sc in cfg.species:
        for b in sc.beams:
            umax = max(umax, math.sqrt(sum(x * x for x in b.drift)) + 4.0 * math.sqrt(max(b.theta, 0.0)) + 1e-3)
    return 2.0 * umax if umax > 0 else 1.0


class HistogramDiag(Diag):
    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        cfg = sim.cfg
        nd = len(d.quantities)
        if nd not in (1, 2, 3) or len(d.bins) != nd:
            raise ValueError(f"diagnostic {d.name!r}: give 1-3 quantities and as many bin counts")
        for q in d.quantities:
            if q not in QUANTITIES:
                raise ValueError(f"diagnostic {d.name!r}: quantity {q!r} not in {QUANTITIES}")
        ranges = list(d.ranges) + [None] * (nd - len(d.ranges))
        logs = list(d.log) + [False] * (nd - len(d.log))
        umax = default_umax(cfg)
        self.edges = []
        for q, nb, rg, lg in zip(d.quantities, d.bins, ranges, logs):
            if rg is None:
                if q in ("x", "y", "z"):
                    rg = (0.0, sim.grid.L["xyz".index(q)])
                elif q in ("ux", "uy", "uz"):
                    rg = (-umax, umax)
                elif q == "u":
                    rg = (1e-4 * umax, umax) if lg else (0.0, umax)
                else:  # gamma / ekin: logarithmic decades by default
                    rg = (1e-6, 1e3) if q == "ekin" else (1.0, 1e3)
                    lg = True if d.log == () else lg
            lo, hi = float(rg[0]), float(rg[1])
            if lg:
                if lo <= 0:
                    raise ValueError(f"diagnostic {d.name!r}: log bins need a positive lower bound for {q}")
                e = np.geomspace(lo, hi, int(nb) + 1)
            else:
                e = np.linspace(lo, hi, int(nb) + 1)
            self.edges.append(e)
        self.spec = []  # (lo, hi, n, log) per axis: bins are uniform in q or log q
        for e in self.edges:
            lg = bool(e[0] > 0 and np.allclose(np.diff(np.log(e)), np.log(e[1] / e[0]), rtol=1e-9)
                      and not np.allclose(np.diff(e), e[1] - e[0], rtol=1e-9))
            self.spec.append((float(e[0]), float(e[-1]), len(e) - 1, lg))
        self.shape = tuple(n for _, _, n, _ in self.spec)
        self.region = None if d.region is None else d.region

    def _quantity(self, sp, q):
        P = sp.P[:, : sp.n]
        if q in ("x", "y", "z"):
            return P["xyz".index(q)]
        if q in ("ux", "uy", "uz"):
            return P[3 + "xyz".index(q[1])]
        u2 = P[3] ** 2 + P[4] ** 2 + P[5] ** 2
        if q == "u":
            return cp.sqrt(u2)
        g = cp.sqrt(1.0 + u2)
        return g if q == "gamma" else sp.m * (g - 1.0)

    def compute(self, src):
        sim = self.sim
        out = {}
        for sp in self.species_list():
            if sp.n:
                w = sp.P[6, : sp.n].copy()
                if self.region is not None:
                    for a, r in enumerate(self.region):
                        if r is not None:
                            x = sp.P[a, : sp.n]
                            w *= (x >= r[0]) & (x < r[1])
                idx = cp.zeros(sp.n, dtype=cp.int64)
                for q, (lo, hi, nb, lg) in zip(self.d.quantities, self.spec):
                    v = self._quantity(sp, q)
                    if lg:
                        f = (cp.log(cp.maximum(v, 1e-300)) - np.log(lo)) * (nb / (np.log(hi) - np.log(lo)))
                    else:
                        f = (v - lo) * (nb / (hi - lo))
                    i = cp.floor(f).astype(cp.int64)
                    i = cp.where(v == hi, nb - 1, i)  # the last bin includes its right edge
                    ok = (i >= 0) & (i < nb)
                    w = cp.where(ok, w, 0.0)
                    idx = idx * nb + cp.clip(i, 0, nb - 1)
                nbt = int(np.prod(self.shape))
                # spread the atomic adds over 64 private copies (few bins -> heavy contention)
                ncp = 64 if nbt <= 4096 else 1
                if ncp > 1:
                    idx = idx + nbt * (cp.arange(sp.n, dtype=cp.int64) & (ncp - 1))
                h = cp.bincount(idx, weights=w, minlength=nbt * ncp).reshape(ncp, nbt).sum(axis=0)
                h = h.reshape(self.shape)
            else:
                h = cp.zeros(self.shape)
            h = sim.comm.allreduce(cp.ascontiguousarray(h))
            if self.root:
                out[sp.name] = _finish(h, self.d)
        return out if self.root else None


def _select(sp, stride, select=None, region=None):
    P = sp.P[:, : sp.n]
    mask = cp.fmod(P[7], float(stride)) == 0 if stride > 1 else cp.ones(sp.n, dtype=cp.bool_)
    if select:
        for k, v in select.items():
            if k == "gamma_min":
                mask &= cp.sqrt(1.0 + P[3] ** 2 + P[4] ** 2 + P[5] ** 2) >= float(v)
            elif k == "ekin_min":
                mask &= sp.m * (cp.sqrt(1.0 + P[3] ** 2 + P[4] ** 2 + P[5] ** 2) - 1.0) >= float(v)
            else:
                raise ValueError(f"unknown particle selection {k!r} (gamma_min, ekin_min)")
    if region is not None:
        for a, r in enumerate(region):
            if r is not None:
                mask &= (P[a] >= r[0]) & (P[a] < r[1])
    return mask


class ParticlesDiag(Diag):
    def compute(self, src):
        sim = self.sim
        out = {}
        for sp in self.species_list():
            mask = _select(sp, self.d.stride, self.d.select, self.d.region)
            Q = sp.P[:, : sp.n][:, mask] if sp.n else cp.zeros((sp.P.shape[0], 0))
            Q = gather_columns(sim.comm, cp.ascontiguousarray(Q[:8]))
            if not self.root:
                continue
            order = cp.argsort(Q[7])
            Q = Q[:, order]
            dt = np.float32 if self.d.dtype == "float32" else np.float64
            rows = {k: _finish(Q[i], self.d).astype(dt) for i, k in enumerate(("x", "y", "z", "ux", "uy", "uz", "w"))}
            rows["id"] = cp.asnumpy(Q[7]).astype(np.uint64)
            rows.update(charge=sp.q, mass=sp.m, momentum_offset=0.0 if sim.implicit else -0.5 * sim.cfg.dt)
            out[sp.name] = rows
        return out if self.root else None


class TracksDiag(Diag):
    QUANT = ("x", "y", "z", "ux", "uy", "uz")
    FIELDQ = ("Ex", "Ey", "Ez", "Bx", "By", "Bz")

    def __init__(self, sim, d, resolved):
        super().__init__(sim, d, resolved)
        self.ids = {}  # species -> sorted tracked ids (rank 0), fixed at the first output

    def compute(self, src):
        sim = self.sim
        out = {}
        nq = 6 + (6 if self.d.with_fields else 0)
        for sp in self.species_list():
            mask = _select(sp, self.d.stride)
            sel = cp.nonzero(mask)[0] if sp.n else cp.zeros(0, dtype=cp.int64)
            m = int(sel.size)
            Q = cp.empty((nq + 1, m))
            if m:
                Q[:6] = sp.P[:6][:, sel]
                Q[nq] = sp.P[7][sel]
                if self.d.with_fields:
                    g = sim.grid
                    f = cp.empty((6, m))
                    launch(sim.k_gather, m, (sp.P, np.int64(sp.cap), sel.astype(cp.int64), np.int32(m),
                                             g.Ex, g.Ey, g.Ez, g.Bx, g.By, g.Bz, f))
                    Q[6:12] = f
            Q = gather_columns(sim.comm, Q)
            if not self.root:
                continue
            Qh = cp.asnumpy(Q)
            ids = Qh[nq].astype(np.uint64)
            if sp.name not in self.ids:
                self.ids[sp.name] = np.sort(ids)
            ref = self.ids[sp.name]
            col = np.searchsorted(ref, ids)
            ok = (col < ref.size)
            ok[ok] = ref[col[ok]] == ids[ok]
            names = self.QUANT + (self.FIELDQ if self.d.with_fields else ())
            dt = np.float32 if self.d.dtype == "float32" else np.float64
            rows = {}
            for i, k in enumerate(names):
                v = np.full(ref.size, np.nan, dtype=dt)
                v[col[ok]] = Qh[i][ok]
                if self.d.keep_bits is not None:
                    v = round_mantissa(v, self.d.keep_bits)
                rows[k] = v
            out[sp.name] = rows
        return out if self.root else None


KINDS = {"fields": FieldsDiag, "slice": SliceDiag, "lineout": LineoutDiag, "modes": ModesDiag,
         "kspectrum": KSpectrumDiag, "histogram": HistogramDiag, "particles": ParticlesDiag, "tracks": TracksDiag}
