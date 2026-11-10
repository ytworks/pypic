"""Particle species storage (SoA on GPU), loading, sorting and inter-rank migration."""
from __future__ import annotations

import cupy as cp
import numpy as np

from .backend import launch
from .config import Beam  # noqa: F401  (re-exported for backwards compatibility)

NROW = 8  # x y z ux uy uz w id   (id: integer-valued float64, PIC-D8)
ID_LIMIT = 2 ** 53  # float64 represents every integer below this exactly


# The kernels index particles with int (thread index, n); row offsets row * cap are 64-bit.
MAX_CAPACITY = 2**31 - 1


def check_capacity(cap: int) -> int:
    """Particles of one species on one GPU are limited by the kernels' int particle index."""
    cap = int(cap)
    if cap > MAX_CAPACITY:
        raise ValueError(f"particle capacity {cap} per species and GPU exceeds 2**31 - 1 "
                         "(the kernels' int particle index); use more GPUs")
    return cap


class Species:
    def __init__(self, name: str, charge: float, mass: float, capacity: int):
        self.name = name
        self.q = float(charge)
        self.m = float(mass)
        self.cap = check_capacity(capacity)
        self.P = cp.zeros((NROW, self.cap), dtype=cp.float64)
        self.n = 0
        self.Pn = None  # implicit scheme work array
        self.test = False  # test species: no charge deposit, not in the energy budget (PIC-D9)
        self.qdep = self.q

    # ---------------------------------------------------------------- storage
    def ensure_capacity(self, need: int):
        if need <= self.cap:
            return
        check_capacity(need)
        cap = min(int(need * 1.25) + 1024, MAX_CAPACITY)
        P = cp.zeros((self.P.shape[0], cap), dtype=cp.float64)
        P[:, : self.n] = self.P[:, : self.n]
        self.P, self.cap = P, cap
        if self.Pn is not None:
            Pn = cp.zeros((6, cap), dtype=cp.float64)
            Pn[:, : self.n] = self.Pn[:, : self.n]
            self.Pn = Pn

    @property
    def data(self):
        return self.P[:, : self.n]

    # ---------------------------------------------------------------- loading
    def load(self, grid, kernel, beam: Beam, ppc: int, seed: int, stream: int):
        """Append ppc particles per physical cell of this rank. Positions depend only on
        (seed, global cell, slot): beams loaded with equal seed share positions. Particle ids are
        ((stream - 1) * N_cells + global cell) * ppc + slot, independent of the rank layout."""
        ax = grid.d.axes
        nl = [a.n_phys for a in ax]
        ntot = nl[0] * nl[1] * nl[2] * ppc
        self.ensure_capacity(self.n + ntot)
        w = beam.density * grid.cell_volume / ppc
        gn = grid.d.n_cells
        id_base = (stream - 1) * gn[0] * gn[1] * gn[2] * ppc
        launch(kernel, ntot, (
            self.P, np.int64(self.cap), np.int64(self.n), np.int32(ppc),
            np.int32(nl[0]), np.int32(nl[1]), np.int32(nl[2]),
            np.int32(ax[0].phys_lo), np.int32(ax[1].phys_lo), np.int32(ax[2].phys_lo),
            np.int32(gn[0]), np.int32(gn[1]), np.int32(gn[2]),
            np.uint64(seed), np.uint64(stream), np.float64(w), np.float64(beam.theta),
            np.float64(beam.drift[0]), np.float64(beam.drift[1]), np.float64(beam.drift[2]), np.float64(id_base)))
        self.n += ntot

    def sort(self, kernel):
        if self.n == 0:
            return
        key = cp.empty(self.n, dtype=cp.int64)
        launch(kernel, self.n, (self.P, np.int32(self.n), np.int64(self.cap), key))
        order = cp.argsort(key)
        self.P[:, : self.n] = self.P[:, order]
        if self.Pn is not None:
            self.Pn[:, : self.n] = self.Pn[:, order]

    # --------------------------------------------------------- remove/migrate
    def flag_list(self):
        """Index list filled by the push kernels with particles to remove or migrate."""
        need = max(1024, self.cap // 8)
        if getattr(self, "_list", None) is None or self._list.size < need:
            self._list = cp.empty(need, dtype=cp.int32)
        return self._list

    @staticmethod
    def _leaving(Q, grid, a):
        """Masks of particles leaving the sub-domain through the low / high side of axis a. At a
        global non-periodic face there is no neighbour: a particle exactly on the face (x = L
        on a reflecting wall) stays with the edge rank instead of being dropped."""
        ax = grid.d.axes[a]
        x = Q[a]
        if ax.lo is not None:
            to_lo = x < ax.phys_lo * grid.dx[a]
        else:
            to_lo = cp.zeros(x.shape, dtype=cp.bool_)
        if ax.hi is not None:
            to_hi = x >= (ax.phys_lo + ax.n_phys) * grid.dx[a]
        else:
            to_hi = cp.zeros(x.shape, dtype=cp.bool_)
        return to_lo, to_hi

    def compact_and_migrate(self, grid, comm):
        """Remove dead particles (w == 0) and send particles that left the local sub-domain to
        neighbour ranks, axis by axis (diagonal moves are handled by successive passes)."""
        P = self.P[:, : self.n]
        keep = P[6] != 0
        P = P[:, keep]
        for a in range(3):
            ax = grid.d.axes[a]
            if ax.nproc == 1:
                continue
            to_lo, to_hi = self._leaving(P, grid, a)
            stay = ~(to_lo | to_hi)
            s_lo = cp.ascontiguousarray(P[:, to_lo])
            s_hi = cp.ascontiguousarray(P[:, to_hi])
            P = P[:, stay]
            # exchange counts (all neighbours exist on decomposed axes except non-periodic edges)
            cnt_send, cnt_recv = [], []
            c_lo = cp.asarray([s_lo.shape[1]], dtype=cp.int64)
            c_hi = cp.asarray([s_hi.shape[1]], dtype=cp.int64)
            r_hi = cp.zeros(1, dtype=cp.int64)
            r_lo = cp.zeros(1, dtype=cp.int64)
            if ax.lo is not None:
                cnt_send.append((c_lo, ax.lo))
            if ax.hi is not None:
                cnt_send.append((c_hi, ax.hi))
            if ax.hi is not None:
                cnt_recv.append((r_hi, ax.hi))
            if ax.lo is not None:
                cnt_recv.append((r_lo, ax.lo))
            comm.exchange(cnt_send, cnt_recv)
            n_hi, n_lo = int(r_hi[0]), int(r_lo[0])
            b_hi = cp.empty((P.shape[0], n_hi))
            b_lo = cp.empty((P.shape[0], n_lo))
            sends, recvs = [], []
            if ax.lo is not None:
                sends.append((s_lo, ax.lo))
            if ax.hi is not None:
                sends.append((s_hi, ax.hi))
            if ax.hi is not None:
                recvs.append((b_hi, ax.hi))
            if ax.lo is not None:
                recvs.append((b_lo, ax.lo))
            comm.exchange(sends, recvs)
            P = cp.concatenate([P, b_hi, b_lo], axis=1)
            # periodic wrap at the global boundary
            if grid.d.bc[a][0] == "periodic":
                L = grid.L[a]
                P[a] -= L * cp.floor(P[a] / L)
        n = P.shape[1]
        self.ensure_capacity(n)
        self.P[:, :n] = P
        self.n = n

    # ------------------------------------------------------------ diagnostics
    def count(self, comm) -> int:
        return int(comm.allreduce_scalar(float(self.n)))


def _order_by(dest, nkeys: int):
    """Permutation that groups columns by dest (0..nkeys-1), stable, with no host sync."""
    m = dest.size
    return cp.argsort(dest.astype(cp.int64) * m + cp.arange(m, dtype=cp.int64))


def migrate(species, grid, comm, flag):
    """Take the particles flagged by the push kernels (flag[s] of them for species s, indices in
    its flag list) out of every species, drop the dead ones (w == 0) and send the rest to the
    neighbour ranks, axis by axis, all species in one message group per axis (issue #7).

    Host syncs: one for the flag counts (with the global total), one for the hole counts, and one
    per decomposed axis for the message sizes. Every rank makes the same collective calls: if any
    rank's flag list overflowed, every rank takes all its particles out (the former compaction)."""
    S = len(species)
    caps = cp.asarray([sp.flag_list().size for sp in species], dtype=cp.int32)
    loc = flag[:S].astype(cp.float64)
    glob = cp.stack([loc.sum(), (flag[:S] > caps).sum().astype(cp.float64)])
    if comm.size > 1:
        glob = comm.allreduce(glob)
    h = cp.asnumpy(cp.concatenate([loc, glob]))
    if h[S] == 0:
        return
    full = h[S + 1] > 0
    # 1. take the flagged particles out; fill the holes below n - k from the tail
    Qs, idxs, meta = [], [], []
    for sp, k in zip(species, h[:S].astype(np.int64)):
        k = sp.n if full else int(k)
        idx = (cp.arange(sp.n, dtype=cp.int64) if full
               else cp.sort(sp.flag_list()[:k]).astype(cp.int64))
        Q = sp.P[:, idx]
        idxs.append((idx, k))
        Qs.append(Q)
        meta.append(cp.stack([cp.count_nonzero(idx < sp.n - k), cp.count_nonzero(Q[6] == 0)]))
    hm = cp.asnumpy(cp.stack(meta)) if meta else np.zeros((0, 2), np.int64)
    for i, (sp, (idx, k)) in enumerate(zip(species, idxs)):
        nh, ndead = int(hm[i, 0]), int(hm[i, 1])
        nk = sp.n - k
        if nh:
            tail = cp.zeros(k, dtype=cp.int8)
            tail[idx[nh:] - nk] = 1  # flagged positions of the tail [nk, n)
            movers = nk + _order_by(tail, 2)[:nh]
            sp.P[:, idx[:nh]] = sp.P[:, movers]
        sp.n = nk
        Q = Qs[i]
        Qs[i] = Q[:, _order_by((Q[6] == 0).astype(cp.int8), 2)[: Q.shape[1] - ndead]]
    # 2. axis-by-axis exchange (diagonal moves take successive axes)
    for a in range(3):
        ax = grid.d.axes[a]
        if ax.nproc == 1:
            continue
        lo_x = ax.phys_lo * grid.dx[a]
        hi_x = (ax.phys_lo + ax.n_phys) * grid.dx[a]
        cnts = []
        for i, Q in enumerate(Qs):
            x = Q[a]
            # 0: to lo, 1: to hi, 2: stays (at a global non-periodic face there is no neighbour:
            # a particle exactly on the face, x = L on a reflecting wall, stays with the edge rank)
            to_lo = (x < lo_x) if ax.lo is not None else cp.zeros(x.shape, cp.bool_)
            to_hi = (x >= hi_x) if ax.hi is not None else cp.zeros(x.shape, cp.bool_)
            dest = cp.where(to_lo, 0, cp.where(to_hi, 1, 2)).astype(cp.int8)
            Qs[i] = Q[:, _order_by(dest, 3)]
            # the three known classes counted on the device: no host sync, and an empty dest gives
            # zeros (CuPy 14.0's bincount takes max() first, which raises on an empty array)
            cnts.append(cp.stack([cp.count_nonzero(dest == v) for v in range(3)]).astype(cp.int64))
        c = cp.stack(cnts) if cnts else cp.zeros((0, 3), cp.int64)
        c_lo, c_hi = cp.ascontiguousarray(c[:, 0]), cp.ascontiguousarray(c[:, 1])
        r_hi, r_lo = cp.zeros(S, cp.int64), cp.zeros(S, cp.int64)
        sends = ([(c_lo, ax.lo)] if ax.lo is not None else []) + ([(c_hi, ax.hi)] if ax.hi is not None else [])
        recvs = ([(r_hi, ax.hi)] if ax.hi is not None else []) + ([(r_lo, ax.lo)] if ax.lo is not None else [])
        comm.exchange(sends, recvs)
        hc = cp.asnumpy(cp.concatenate([c.ravel(), r_hi, r_lo])).astype(np.int64)
        cs, rh, rl = hc[: 3 * S].reshape(S, 3), hc[3 * S: 4 * S], hc[4 * S:]
        s_lo, s_hi, b_hi, b_lo = [], [], [], []
        for i, Q in enumerate(Qs):
            n0, n1 = cs[i, 0], cs[i, 1]
            s_lo.append((Q[:, :n0], ax.lo))
            s_hi.append((Q[:, n0:n0 + n1], ax.hi))
            b_hi.append((cp.empty((Q.shape[0], int(rh[i]))), ax.hi))
            b_lo.append((cp.empty((Q.shape[0], int(rl[i]))), ax.lo))
        sends = (s_lo if ax.lo is not None else []) + (s_hi if ax.hi is not None else [])
        recvs = (b_hi if ax.hi is not None else []) + (b_lo if ax.lo is not None else [])
        comm.exchange(sends, recvs)
        for i, Q in enumerate(Qs):
            n01 = cs[i, 0] + cs[i, 1]
            Q = cp.concatenate([Q[:, n01:], b_hi[i][0], b_lo[i][0]], axis=1)
            if grid.d.bc[a][0] == "periodic":
                L = grid.L[a]
                Q[a] -= L * cp.floor(Q[a] / L)
            Qs[i] = Q
    # 3. append the particles that now belong to this rank
    for sp, Q in zip(species, Qs):
        m = Q.shape[1]
        if m:
            sp.ensure_capacity(sp.n + m)
            sp.P[:, sp.n: sp.n + m] = Q
            sp.n += m
