"""Particle species storage (SoA on GPU), loading, sorting and inter-rank migration."""
from __future__ import annotations

from contextlib import contextmanager

import cupy as cp
import numpy as np

from .backend import launch
from .config import MAX_CAPACITY, Beam  # noqa: F401  (Beam re-exported for backwards compatibility)

NROW = 8  # x y z ux uy uz w id   (id: integer-valued float64, PIC-D8)
TWO32 = 4294967296.0  # positions of the compact storage are kept in 2^-32 of a cell (PIC-D19)


def check_capacity(cap: int) -> int:
    """Particles of one species on one GPU are limited by the kernels' int particle index."""
    cap = int(cap)
    if cap > MAX_CAPACITY:
        raise ValueError(f"particle capacity {cap} per species and GPU exceeds 2**31 - 1 "
                         "(the kernels' int particle index); use more GPUs")
    return cap


class SortMemory:
    """The temporaries of the particle sort (`SORT_BYTES` per particle) in a memory pool of their own.

    In CuPy's default pool the argsort of a sort took new device memory while the pool held more
    free bytes than the request, and the pool kept both (measured: blast wave, every rank). In this
    pool the blocks of a sort meet only the next sort: one of fewer particles fits them, and one
    of more particles than any before returns them to the device first. So the pool holds the
    temporaries of the largest sort and nothing else. One process drives one device."""

    def __init__(self):
        self.pool = cp.cuda.MemoryPool()
        self.n = 0  # particles of the sort that the blocks were allocated for

    def release(self):
        self.pool.free_all_blocks()
        self.n = 0

    @contextmanager
    def __call__(self, n: int):
        """Allocate inside the block from the pool of the sort; nothing allocated may outlive it."""
        if n > self.n:
            self.pool.free_all_blocks()
            cp.get_default_memory_pool().free_all_blocks()  # what the set-up and the steps before left free
            self.n = n
        held = self.pool.total_bytes()
        with cp.cuda.using_allocator(self.pool.malloc):
            yield
        if self.pool.total_bytes() > held > 0:  # the kept blocks did not fit: the next sort starts again
            self.n = 0


sort_memory = SortMemory()


class Species:
    def __init__(self, name: str, charge: float, mass: float, capacity: int):
        self.name = name
        self.q = float(charge)
        self.m = float(mass)
        self.cap = check_capacity(capacity)
        self._allocate()
        self.n = 0
        self.Pn = None  # implicit scheme work array
        self.test = False  # test species: no charge deposit, not in the energy budget (PIC-D9)
        self.qdep = self.q

    # ---------------------------------------------------------------- storage
    def _allocate(self):
        self.P = cp.zeros((NROW, self.cap), dtype=cp.float64)

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

    # -------------------------------------------------------------- accessors
    # Code outside this module and the kernels reaches the particles only through these, in FP64
    # physical values, so the layout of the storage can change with the build (PIC-D19). A returned
    # array may be a view of the storage: read it, do not write into it.
    @property
    def store(self):
        """The storage as the arguments of a kernel's PSTORE parameters."""
        return (self.P,)

    def _range(self, lo, hi):
        return slice(lo, self.n if hi is None else hi)

    def pos(self, a, lo=0, hi=None):
        """Position along axis a of the particles lo..hi (default: all)."""
        return self.P[a, self._range(lo, hi)]

    def mom(self, a, lo=0, hi=None):
        """Momentum u = gamma v along axis a."""
        return self.P[3 + a, self._range(lo, hi)]

    def weight(self, lo=0, hi=None):
        return self.P[6, self._range(lo, hi)]

    def ids(self, lo=0, hi=None):
        return self.P[7, self._range(lo, hi)]

    def rows(self, sel=None):
        """(8, m) block x, y, z, ux, uy, uz, w, id of all particles, or of those of ``sel`` (a mask
        or an index array)."""
        P = self.P[:NROW, : self.n]
        return P if sel is None else P[:, sel]

    def scale_weight(self, f, lo=0, hi=None):
        self.P[6, self._range(lo, hi)] *= f

    def add_mom(self, a, du, lo=0, hi=None):
        self.P[3 + a, self._range(lo, hi)] += du

    def append(self, block):
        """Append the particles of an (8, m) block x, y, z, ux, uy, uz, w, id."""
        m = block.shape[1]
        self.ensure_capacity(self.n + m)
        self.P[:NROW, self.n: self.n + m] = cp.asarray(block)
        self.n += m

    # The remaining displacement of the set-up (density and displacement perturbations): three more
    # rows of the storage while it lasts, so that it migrates with its particle.
    @property
    def displacing(self):
        return self.P.shape[0] > NROW

    def set_displacement(self, disp):
        self.P = cp.concatenate([self.P, disp], axis=0)

    def displacement_max(self, a):
        return float(cp.max(cp.abs(self.P[NROW + a, : self.n])))

    def end_displacement(self):
        if self.displacing:
            self.P = cp.ascontiguousarray(self.P[:NROW])

    # What the migration needs to know of the layout: the doubles of a particle in its buffers, the
    # storage as the arguments of the kernels that take the displacement too, and the range of this
    # rank and the period along an axis in the coordinate of the buffers.
    @property
    def nrow(self):
        return self.P.shape[0]

    @property
    def disp_store(self):
        return (self.P,)

    def frame(self, grid, a):
        """(low limit, high limit, period) of axis a: the physical cells of this rank and the box."""
        ax = grid.d.axes[a]
        return ax.phys_lo * grid.dx[a], (ax.phys_lo + ax.n_phys) * grid.dx[a], grid.L[a]

    def snap(self, x, a, xp=np):
        """The position x along axis a in the coordinate of ``frame``, as the storage will hold it."""
        return x

    def _records(self):
        """All particles as rows of doubles, the first three in the coordinate of ``frame``."""
        return self.P[:, : self.n]

    def _set_records(self, Q):
        n = Q.shape[1]
        self.ensure_capacity(n)
        self.P[:, :n] = Q
        self.n = n

    # ---------------------------------------------------------------- loading
    def load(self, grid, kernel, beam: Beam, ppc: int, seed: int, stream: int, block=None):
        """Append ppc particles per cell of ``block`` = (first global cell, cells) per axis, a part of
        this rank's physical cells (default: all of them). Positions depend only on
        (seed, global cell, slot): beams loaded with equal seed share positions. Particle ids are
        ((stream - 1) * N_cells + global cell) * ppc + slot, independent of the rank layout."""
        ax = grid.d.axes
        g0, nl = block if block is not None else ([a.phys_lo for a in ax], [a.n_phys for a in ax])
        ntot = nl[0] * nl[1] * nl[2] * ppc
        if ntot == 0:
            return
        self.ensure_capacity(self.n + ntot)
        w = beam.density * grid.cell_volume / ppc
        gn = grid.d.n_cells
        id_base = (stream - 1) * gn[0] * gn[1] * gn[2] * ppc
        origin = beam.origin if beam.origin is not None else (0.0, 0.0, 0.0)
        launch(kernel, ntot, (
            *self.store, np.int64(self.cap), np.int64(self.n), np.int32(ppc),
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
        with sort_memory(self.n):
            key = cp.empty(self.n, dtype=cp.int64)
            launch(kernel, self.n, (*self.store, np.int32(self.n), np.int64(self.cap), key))
            order = cp.argsort(key)
            buf = key.view(cp.float64)
            for A in (self.P, self.Pn):
                if A is not None:
                    for r in range(A.shape[0]):
                        cp.take(A[r], order, out=buf)
                        A[r, : self.n] = buf
            del key, order, buf

    # --------------------------------------------------------- remove/migrate
    def flag_list(self):
        """Index list filled by the push kernels with particles to remove or migrate."""
        need = max(1024, self.cap // 8)
        if getattr(self, "_list", None) is None or self._list.size < need:
            self._list = cp.empty(need, dtype=cp.int32)
        return self._list

    def _leaving(self, Q, grid, a):
        """Masks of particles leaving the sub-domain through the low / high side of axis a. At a
        global non-periodic face there is no neighbour: a particle exactly on the face (x = L
        on a reflecting wall) stays with the edge rank instead of being dropped."""
        ax = grid.d.axes[a]
        x = Q[a]
        x_lo, x_hi, _ = self.frame(grid, a)
        if ax.lo is not None:
            to_lo = x < x_lo
        else:
            to_lo = cp.zeros(x.shape, dtype=cp.bool_)
        if ax.hi is not None:
            to_hi = x >= x_hi
        else:
            to_hi = cp.zeros(x.shape, dtype=cp.bool_)
        return to_lo, to_hi

    def compact_and_migrate(self, grid, comm):
        """Remove dead particles (w == 0) and send particles that left the local sub-domain to
        neighbour ranks, axis by axis (diagonal moves are handled by successive passes)."""
        P = self._records()
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
                L = self.frame(grid, a)[2]
                P[a] -= L * cp.floor(P[a] / L)
        self._set_records(P)

    # ------------------------------------------------------------ diagnostics
    def count(self, comm) -> int:
        return int(comm.allreduce_scalar(float(self.n)))


class CompactSpecies(Species):
    """The storage of ``precision: compact`` (PIC-D19), 56 bytes per particle: ``W`` weight and id
    (float64), ``C`` the cell in the arrays of this rank (int32, the flat index of the field
    arrays), ``F`` the position within the cell in 2^-32 of a cell (uint32 per axis), ``U`` the
    momentum (float64 per axis). The accessors return FP64 physical values; a position is
    (global cell + fraction) * dx with one rounding, the expression of the kernels. The set-up
    displacement has its own array ``Dsp`` while it lasts. The buffers of the migration hold the
    position as the global index coordinate x / dx, which is exact."""

    def __init__(self, name: str, charge: float, mass: float, capacity: int, grid):
        self.shape = tuple(grid.shape)
        self.origin = tuple(ax.cell_lo - grid.g for ax in grid.d.axes)  # global cell of local index 0
        self.h = tuple(grid.dx)
        self.cells = tuple(grid.d.n_cells)
        self.periodic = tuple(b[0] == "periodic" for b in grid.d.bc)
        super().__init__(name, charge, mass, capacity)

    # ---------------------------------------------------------------- storage
    def _allocate(self):
        self.W = cp.zeros((2, self.cap), dtype=cp.float64)
        self.C = cp.zeros(self.cap, dtype=cp.int32)
        self.F = cp.zeros((3, self.cap), dtype=cp.uint32)
        self.U = cp.zeros((3, self.cap), dtype=cp.float64)
        self.Dsp = None
        self._nodisp = cp.zeros(1, dtype=cp.float64)  # stands for Dsp where a kernel does not read it

    def _arrays(self):
        return [(k, A) for k, A in (("W", self.W), ("C", self.C), ("F", self.F), ("U", self.U), ("Dsp", self.Dsp))
                if A is not None]

    def ensure_capacity(self, need: int):
        if need <= self.cap:
            return
        check_capacity(need)
        cap = min(int(need * 1.25) + 1024, MAX_CAPACITY)
        for k, A in self._arrays():
            B = cp.zeros(A.shape[:-1] + (cap,), dtype=A.dtype)
            B[..., : self.n] = A[..., : self.n]
            setattr(self, k, B)
        self.cap = cap

    # -------------------------------------------------------------- accessors
    @property
    def store(self):
        return self.W, self.C, self.F, self.U

    def index_pos(self, a, lo=0, hi=None):
        """Position along axis a as the global index coordinate x / dx (exact)."""
        r = self._range(lo, hi)
        c, (_, ny, nz) = self.C[r], self.shape
        cell = c // (ny * nz) if a == 0 else (c // nz) % ny if a == 1 else c % nz
        return (cell + self.origin[a]) + self.F[a, r] * (1.0 / TWO32)

    def pos(self, a, lo=0, hi=None):
        return self.index_pos(a, lo, hi) * self.h[a]

    def mom(self, a, lo=0, hi=None):
        return self.U[a, self._range(lo, hi)]

    def weight(self, lo=0, hi=None):
        return self.W[0, self._range(lo, hi)]

    def ids(self, lo=0, hi=None):
        return self.W[1, self._range(lo, hi)]

    def rows(self, sel=None):
        R = [self.pos(a) for a in range(3)] + [self.mom(a) for a in range(3)] + [self.weight(), self.ids()]
        return cp.stack(R if sel is None else [r[sel] for r in R])

    def scale_weight(self, f, lo=0, hi=None):
        self.W[0, self._range(lo, hi)] *= f

    def add_mom(self, a, du, lo=0, hi=None):
        self.U[a, self._range(lo, hi)] += du

    def snap(self, x, a, xp=np):
        """The index coordinate of the position x along axis a, on the lattice of the storage. A
        position that rounds up to the end of the box is the first position of a periodic axis and
        the last one before a wall."""
        X = xp.rint(x / self.h[a] * TWO32) / TWO32
        top = float(self.cells[a])
        return xp.where(X >= top, X - top, X) if self.periodic[a] else xp.minimum(X, top - 1.0 / TWO32)

    def _set_index_pos(self, r, X):
        """Store index coordinates that are on the lattice already (three arrays) in the slots r."""
        flat = 0
        for a in range(3):
            g = cp.floor(X[a])
            c = g.astype(cp.int64) - self.origin[a]
            if c.size and not (int(c.min()) >= 0 and int(c.max()) < self.shape[a]):
                raise ValueError(f"species {self.name}: a particle lies outside the arrays of its rank (axis {a})")
            self.F[a, r] = ((X[a] - g) * TWO32).astype(cp.uint32)
            flat = flat * self.shape[a] + c
        self.C[r] = flat.astype(cp.int32)

    def append(self, block):
        block = cp.asarray(block, dtype=cp.float64)
        m = block.shape[1]
        self.ensure_capacity(self.n + m)
        r = slice(self.n, self.n + m)
        self._set_index_pos(r, [self.snap(block[a], a, cp) for a in range(3)])
        self.U[:, r] = block[3:6]
        self.W[:, r] = block[6:8]
        self.n += m

    @property
    def displacing(self):
        return self.Dsp is not None

    def set_displacement(self, disp):
        self.Dsp = cp.ascontiguousarray(disp, dtype=cp.float64)

    def displacement_max(self, a):
        return float(cp.max(cp.abs(self.Dsp[a, : self.n])))

    def end_displacement(self):
        self.Dsp = None

    @property
    def nrow(self):
        return NROW + (3 if self.displacing else 0)

    @property
    def disp_store(self):
        return (*self.store, self.Dsp if self.displacing else self._nodisp)

    def frame(self, grid, a):
        ax = grid.d.axes[a]
        return float(ax.phys_lo), float(ax.phys_lo + ax.n_phys), float(self.cells[a])

    def _records(self):
        R = [self.index_pos(a) for a in range(3)] + [self.mom(a) for a in range(3)] + [self.weight(), self.ids()]
        if self.displacing:
            R += [self.Dsp[a, : self.n] for a in range(3)]
        return cp.stack(R)

    def _set_records(self, Q):
        n = Q.shape[1]
        self.ensure_capacity(n)
        r = slice(0, n)
        self._set_index_pos(r, Q[:3])
        self.U[:, r] = Q[3:6]
        self.W[:, r] = Q[6:8]
        if self.displacing:
            self.Dsp[:, r] = Q[8:11]
        self.n = n

    def sort(self, kernel=None):
        """Order the particles by cell: the key is the stored cell, the rows are gathered one at a
        time through one buffer of 8 bytes per particle."""
        n = self.n
        if n == 0:
            return
        with sort_memory(n):
            order = cp.argsort(self.C[:n])
            buf = cp.empty(n, dtype=cp.float64)
            for _, A in self._arrays():
                out = buf.view(A.dtype)[:n]
                for row in A.reshape(-1, self.cap):
                    cp.take(row, order, out=out)
                    row[:n] = out
            del order, buf, out


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


def migrate(species, grid, comm, flag, kernels):
    """Take the particles flagged by the push kernels (flag[s] of them for species s, indices in
    its flag list; flag[S + s] counts the particles appended behind the n of species s, which
    are taken in first) out of every species, drop the dead ones (w == 0) and send the rest to the
    neighbour ranks, axis by axis, all species in one message group per axis (issue #7).
    ``kernels``: the functions MIGRATE_KERNELS of the particle module (issue #19).

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
        nrow = np.int32(sp.nrow)
        Q = cp.empty(k * int(nrow))
        if k:
            if not full:
                lst = ascending(lst[:k], k_rank32)  # the kernels fill the list in any order
            launch(k_take, k, (*sp.disp_store, np.int64(sp.cap), nrow, lst, np.int32(k), np.int32(full), Q))
            if not full:
                launch(k_fill, k, (*sp.disp_store, np.int64(sp.cap), nrow, lst, np.int32(k), np.int32(sp.n - k)))
            sp.n -= k
        Qs.append((Q, k))
    # 2. axis-by-axis exchange (diagonal moves take successive axes); without a decomposed axis
    #    one pass (axis -1) drops the dead particles
    for a in [a for a in range(3) if grid.d.axes[a].nproc > 1] or [-1]:
        ax = grid.d.axes[a] if a >= 0 else None
        has_lo, has_hi = (ax.lo is not None, ax.hi is not None) if ax else (False, False)
        periodic = bool(ax) and grid.d.bc[a][0] == "periodic"
        frames = [sp.frame(grid, a) if ax else (0.0, 0.0, 0.0) for sp in species]  # low, high, period
        # rows: to lo, to hi, stays, dead, from hi, from lo; one column per species
        C = cp.zeros((6, S), dtype=cp.int32)
        keys = []
        for i, (Q, m) in enumerate(Qs):
            key = cp.empty(m, dtype=cp.int64)
            if m:
                x_lo, x_hi, _ = frames[i]
                launch(k_classify, m, (Q, np.int32(species[i].nrow), np.int32(m), np.int32(a),
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
            nrow = species[i].nrow
            wrap = frames[i][2] if periodic else 0.0
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
            launch(k_append, m, (*sp.disp_store, np.int64(sp.cap), np.int32(sp.nrow), np.int32(sp.n), Q, np.int32(m)))
            sp.n += m
