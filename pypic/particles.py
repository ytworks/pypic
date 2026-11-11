"""Particle species storage (SoA on GPU), loading, sorting and inter-rank migration."""
from __future__ import annotations

import cupy as cp
import numpy as np

from .backend import launch
from .config import MAX_CAPACITY, Beam  # noqa: F401  (Beam re-exported for backwards compatibility)

NROW = 8  # x y z ux uy uz w id   (id: integer-valued float64, PIC-D8)


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
        self._sorted = 0  # the largest n sorted so far
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
    def load(self, grid, kernel, beam: Beam, ppc: int, seed: int, stream: int, block=None, id_base=None):
        """Append ppc particles per cell of ``block`` = (first global cell, cells) per axis, a part of
        this rank's physical cells (default: all of them). Positions depend only on
        (seed, global cell, slot): beams loaded with equal seed share positions. Particle ids are
        ((stream - 1) * N_cells + global cell) * ppc + slot, independent of the rank layout;
        ``id_base`` replaces (stream - 1) * N_cells * ppc (planes of the moving window, PIC-D22)."""
        ax = grid.d.axes
        g0, nl = block if block is not None else ([a.phys_lo for a in ax], [a.n_phys for a in ax])
        ntot = nl[0] * nl[1] * nl[2] * ppc
        if ntot == 0:
            return
        self.ensure_capacity(self.n + ntot)
        w = beam.density * grid.cell_volume / ppc
        gn = grid.d.n_cells
        if id_base is None:
            id_base = (stream - 1) * gn[0] * gn[1] * gn[2] * ppc
        origin = beam.origin if beam.origin is not None else (0.0, 0.0, 0.0)
        launch(kernel, ntot, (
            self.P, np.int64(self.cap), np.int64(self.n), np.int32(ppc),
            np.int32(nl[0]), np.int32(nl[1]), np.int32(nl[2]),
            np.int32(g0[0]), np.int32(g0[1]), np.int32(g0[2]),
            np.int32(gn[0]), np.int32(gn[1]), np.int32(gn[2]),
            np.uint64(seed), np.uint64(stream), np.float64(w), np.float64(beam.theta),
            np.float64(beam.drift[0]), np.float64(beam.drift[1]), np.float64(beam.drift[2]),
            np.float64(beam.drift_radial), np.float64(origin[0]), np.float64(origin[1]), np.float64(origin[2]),
            np.float64(id_base)))
        self.n += ntot

    def sort(self, kernel):
        """Order the particles by cell. The rows are gathered one at a time into the key array,
        which is free once the order is known: no second copy of the particle arrays."""
        if self.n == 0:
            return
        key = cp.empty(self.n, dtype=cp.int64)
        launch(kernel, self.n, (self.P, np.int32(self.n), np.int64(self.cap), key))
        order = cp.argsort(key)
        buf = key.view(cp.float64)
        for A in (self.P, self.Pn):
            if A is not None:
                for r in range(A.shape[0]):
                    cp.take(A[r], order, out=buf)
                    A[r, : self.n] = buf
        # After a sort larger than any before, the temporaries (40 bytes per particle) go back to
        # the device, not only to CuPy's pool: a free block serves only requests that fit, so with a
        # growing n (injection) every sort would allocate new blocks and keep the old ones. With a
        # constant n the pool keeps the blocks and the next sort reuses them.
        del key, order, buf
        if self.n > self._sorted:
            self._sorted = self.n
            cp.get_default_memory_pool().free_all_blocks()

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


MIGRATE_KERNELS = ("mig_take", "mig_fill", "mig_classify", "mig_pack", "mig_append", "mig_rank32", "mig_rank64")
RANK_SORT_MAX = 16384  # keys; about where a radix sort overtakes the n * n compares of the counting kernel


def ascending(a, k_rank):
    """The keys of ``a`` in ascending order: a new array, or ``a`` itself sorted in place. The
    lists of a migration are short, and a radix sort costs several launches whatever the length;
    up to RANK_SORT_MAX keys one launch of the counting kernel gives the same order."""
    n = a.size
    if n < 2:
        return a
    if n > RANK_SORT_MAX:
        a.sort()
        return a
    out = cp.empty_like(a)
    launch(k_rank, 32 * n, (a, np.int32(n), out))
    return out


def migrate(species, grid, comm, flag, kernels, shift=0):
    """Take the particles flagged by the push kernels (flag[s] of them for species s, indices in
    its flag list; flag[S + s] counts the particles appended behind the n of species s, which
    are taken in first) out of every species, drop the dead ones (w == 0) and send the rest to the
    neighbour ranks, axis by axis, all species in one message group per axis (issue #7).
    ``kernels``: the functions MIGRATE_KERNELS of the particle module (issue #19). ``shift``: cells
    the moving window has gone along x (PIC-D22); the rank bounds along x move with it.

    Host syncs: one for the flag counts and one per decomposed axis for the message sizes. With
    more than one rank every rank makes the same collective calls, whatever it has flagged. A
    species whose flag list overflowed has all its particles taken out (the former compaction)."""
    k_take, k_fill, k_classify, k_pack, k_append, k_rank32, k_rank64 = kernels
    S = len(species)
    h = cp.asnumpy(flag)
    for sp, k in zip(species, h[S:]):  # particles the injection appended this step (PIC-D13)
        sp.n += int(k)
    h = h[:S]
    if comm.size == 1 and not h.any():
        return
    # 1. take the flagged particles out; fill the holes below n - k from the tail
    Qs = []  # per species: particle-major buffer, particles in it
    for sp, k in zip(species, h):
        lst = sp.flag_list()
        full = k > lst.size
        k = sp.n if full else int(k)
        nrow = np.int32(sp.P.shape[0])
        Q = cp.empty(k * int(nrow))
        if k:
            if not full:
                lst = ascending(lst[:k], k_rank32)  # the kernels fill the list in any order
            launch(k_take, k, (sp.P, np.int64(sp.cap), nrow, lst, np.int32(k), np.int32(full), Q))
            if not full:
                launch(k_fill, k, (sp.P, np.int64(sp.cap), nrow, lst, np.int32(k), np.int32(sp.n - k)))
            sp.n -= k
        Qs.append((Q, k))
    # 2. axis-by-axis exchange (diagonal moves take successive axes); without a decomposed axis
    #    one pass (axis -1) drops the dead particles
    for a in [a for a in range(3) if grid.d.axes[a].nproc > 1] or [-1]:
        ax = grid.d.axes[a] if a >= 0 else None
        has_lo, has_hi = (ax.lo is not None, ax.hi is not None) if ax else (False, False)
        off = shift if a == 0 else 0  # the same doubles as the kernel window_set gets
        x_lo = (ax.phys_lo + off) * grid.dx[a] if ax else 0.0
        x_hi = (ax.phys_lo + ax.n_phys + off) * grid.dx[a] if ax else 0.0
        wrap = grid.L[a] if ax and grid.d.bc[a][0] == "periodic" else 0.0
        # rows: to lo, to hi, stays, dead, from hi, from lo; one column per species
        C = cp.zeros((6, S), dtype=cp.int32)
        keys = []
        for i, (Q, m) in enumerate(Qs):
            key = cp.empty(m, dtype=cp.int64)
            if m:
                launch(k_classify, m, (Q, np.int32(species[i].P.shape[0]), np.int32(m), np.int32(a),
                                       np.int32(has_lo), np.int32(has_hi), np.float64(x_lo), np.float64(x_hi),
                                       key, C.ravel()[i:], np.int32(S)))
                key = ascending(key, k_rank64)
            keys.append(key)
        if ax:
            comm.exchange(([(C[0], ax.lo)] if has_lo else []) + ([(C[1], ax.hi)] if has_hi else []),
                          ([(C[4], ax.hi)] if has_hi else []) + ([(C[5], ax.lo)] if has_lo else []))
        hc = cp.asnumpy(C).astype(np.int64)
        s_lo, s_hi, r_hi, r_lo = [], [], [], []
        for i, (Q, m) in enumerate(Qs):
            n0, n1, n2, _, rh, rl = (int(v) for v in hc[:, i])
            nrow = species[i].P.shape[0]
            send = cp.empty((n0 + n1) * nrow)
            keep = cp.empty((n2 + rh + rl) * nrow)
            launch(k_pack, n0 + n1 + n2, (Q, np.int32(nrow), np.int32(m), np.int32(a), keys[i], np.int32(n0 + n1),
                                          np.int32(n2), np.float64(wrap), send, keep))
            if ax:
                s_lo.append((send[: n0 * nrow], ax.lo))
                s_hi.append((send[n0 * nrow:], ax.hi))
                r_hi.append((keep[n2 * nrow: (n2 + rh) * nrow], ax.hi))
                r_lo.append((keep[(n2 + rh) * nrow:], ax.lo))
            Qs[i] = (keep, n2 + rh + rl)
        if ax:
            comm.exchange((s_lo if has_lo else []) + (s_hi if has_hi else []),
                          (r_hi if has_hi else []) + (r_lo if has_lo else []))
    # 3. append the particles that now belong to this rank
    for sp, (Q, m) in zip(species, Qs):
        if m:
            sp.ensure_capacity(sp.n + m)
            launch(k_append, m, (sp.P, np.int64(sp.cap), np.int32(sp.P.shape[0]), np.int32(sp.n), Q, np.int32(m)))
            sp.n += m
