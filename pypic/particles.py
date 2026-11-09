"""Particle species storage (SoA on GPU), loading, sorting and inter-rank migration."""
from __future__ import annotations

from dataclasses import dataclass

import cupy as cp
import numpy as np

from .backend import launch

NROW = 7  # x y z ux uy uz w


@dataclass
class Beam:
    """One drifting Maxwell-Juttner component of a species."""

    density: float  # in units of the reference density n0
    drift: tuple = (0.0, 0.0, 0.0)  # proper velocity u = gamma v / c
    theta: float = 0.0  # temperature T / (m c^2)


class Species:
    def __init__(self, name: str, charge: float, mass: float, capacity: int):
        self.name = name
        self.q = float(charge)
        self.m = float(mass)
        self.cap = int(capacity)
        self.P = cp.zeros((NROW, self.cap), dtype=cp.float64)
        self.n = 0
        self.Pn = None  # implicit scheme work array

    # ---------------------------------------------------------------- storage
    def ensure_capacity(self, need: int):
        if need <= self.cap:
            return
        cap = int(need * 1.25) + 1024
        P = cp.zeros((NROW, cap), dtype=cp.float64)
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
        (seed, global cell, slot): beams loaded with equal seed share positions."""
        ax = grid.d.axes
        nl = [a.n_phys for a in ax]
        ntot = nl[0] * nl[1] * nl[2] * ppc
        self.ensure_capacity(self.n + ntot)
        w = beam.density * grid.cell_volume / ppc
        launch(kernel, ntot, (
            self.P, np.int32(self.cap), np.int64(self.n), np.int32(ppc),
            np.int32(nl[0]), np.int32(nl[1]), np.int32(nl[2]),
            np.int32(ax[0].phys_lo), np.int32(ax[1].phys_lo), np.int32(ax[2].phys_lo),
            np.int32(grid.d.n_cells[0]), np.int32(grid.d.n_cells[1]),
            np.uint64(seed), np.uint64(stream), np.float64(w), np.float64(beam.theta),
            np.float64(beam.drift[0]), np.float64(beam.drift[1]), np.float64(beam.drift[2])))
        self.n += ntot

    def sort(self, kernel):
        if self.n == 0:
            return
        key = cp.empty(self.n, dtype=cp.int64)
        launch(kernel, self.n, (self.P, np.int32(self.n), np.int32(self.cap), key))
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

    def remove_and_migrate(self, grid, comm, k: int):
        """O(k) version of :meth:`compact_and_migrate`: the k flagged particles (indices in the
        flag list) are taken out, the holes are filled from the tail, and the departing particles
        are exchanged axis by axis. Falls back to the full compaction if the list overflowed."""
        lst = self.flag_list()
        if k > lst.size:
            return self.compact_and_migrate(grid, comm)
        n = self.n
        idx = cp.sort(lst[:k]).astype(cp.int64)
        Q = self.P[:, idx]  # departing / dead particles
        Q = Q[:, Q[6] != 0]  # drop removed (open boundary) particles
        nk = n - k
        # fill holes below nk with the surviving particles of the tail [nk, n)
        holes = idx[idx < nk]
        if holes.size:
            tail = cp.ones(k, dtype=cp.bool_)
            tail[(idx[idx >= nk] - nk)] = False
            movers = nk + cp.nonzero(tail)[0]
            self.P[:, holes] = self.P[:, movers]
        self.n = nk
        Q = self._exchange(Q, grid, comm)
        m = Q.shape[1]
        if m:
            self.ensure_capacity(self.n + m)
            self.P[:, self.n: self.n + m] = Q
            self.n += m

    def _exchange(self, Q, grid, comm):
        """Send particles in Q that lie outside the local sub-domain to neighbours, axis by axis;
        returns the particles that now belong to this rank (periodically wrapped)."""
        for a in range(3):
            ax = grid.d.axes[a]
            if ax.nproc == 1:
                continue
            lo = ax.phys_lo * grid.dx[a]
            hi = (ax.phys_lo + ax.n_phys) * grid.dx[a]
            x = Q[a]
            to_lo = x < lo
            to_hi = x >= hi
            s_lo = cp.ascontiguousarray(Q[:, to_lo])
            s_hi = cp.ascontiguousarray(Q[:, to_hi])
            Q = Q[:, ~(to_lo | to_hi)]
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
            b_hi = cp.empty((NROW, int(r_hi[0])))
            b_lo = cp.empty((NROW, int(r_lo[0])))
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
            Q = cp.concatenate([Q, b_hi, b_lo], axis=1)
            if grid.d.bc[a][0] == "periodic":
                L = grid.L[a]
                Q[a] -= L * cp.floor(Q[a] / L)
        return Q

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
            lo = ax.phys_lo * grid.dx[a]
            hi = (ax.phys_lo + ax.n_phys) * grid.dx[a]
            x = P[a]
            to_lo = x < lo
            to_hi = x >= hi
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
            b_hi = cp.empty((NROW, n_hi))
            b_lo = cp.empty((NROW, n_lo))
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
