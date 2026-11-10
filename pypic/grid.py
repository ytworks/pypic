"""Yee grid: field storage, halo exchange, boundary conditions, CPML and field diagnostics."""
from __future__ import annotations

import math

import cupy as cp
import numpy as np

from .backend import compile_module, launch
from .kernelsrc import FIELD_SOURCES, grid_defines
from .output.names import STAGGER  # noqa: F401  (staggering offsets of each component, re-exported)
from .parallel.decomp import Decomposition

FIELD_NAMES = ("Ex", "Ey", "Ez", "Bx", "By", "Bz")


class YeeGrid:
    def __init__(self, decomp: Decomposition, dx, dt: float, comm, ghost: int = 2, pml: bool = True):
        self.d = decomp
        self.comm = comm
        self.dx = tuple(float(v) for v in dx)
        self.dt = float(dt)
        self.g = int(ghost)
        g = self.g
        ax = decomp.axes
        self.shape = tuple(a.n_int + 2 * g for a in ax)
        self.L = tuple(n * h for n, h in zip(decomp.n_cells, self.dx))
        self.cell_volume = self.dx[0] * self.dx[1] * self.dx[2]
        self.use_pml = pml and any(a.pml_lo or a.pml_hi for a in ax)
        z = lambda: cp.zeros(self.shape, dtype=cp.float64)  # noqa: E731
        for n in FIELD_NAMES + ("Jx", "Jy", "Jz"):
            setattr(self, n, z())
        self.defines = grid_defines(decomp, self.dx, self.dt, self.g, self.use_pml)
        self.mod = compile_module(FIELD_SOURCES, self.defines)
        self.k_update_B = self.mod.get_function("update_B")
        self.k_update_E = self.mod.get_function("update_E")
        self.k_curl_E = self.mod.get_function("curl_E")
        self.k_curl_B = self.mod.get_function("curl_B")
        self.k_div_E = self.mod.get_function("div_E")
        self.k_div_B = self.mod.get_function("div_B")
        self.k_binomial = self.mod.get_function("binomial_pass3")
        self.k_copy_pass = self.mod.get_function("copy_pass3")
        self.n_update = int(np.prod([a.n_int for a in ax]))
        if self.use_pml:
            self._init_pml()
        self.boundary_flux = 0.0  # cumulative Poynting energy leaving the physical domain

    # ------------------------------------------------------------------ setup
    def global_positions(self, a: int, stagger: float) -> np.ndarray:
        """Global coordinate in cells of local array index i (+stagger) along axis a."""
        ax = self.d.axes[a]
        return np.arange(self.shape[a]) + ax.cell_lo - self.g + stagger

    def pml_sigma(self, a: int, stagger: float, order: int = 3, R: float = 1e-8) -> np.ndarray:
        """Graded conductivity profile along axis a (zero in the physical domain)."""
        x = self.global_positions(a, stagger)
        N = self.d.n_cells[a]
        npml = max(self.d.npml, 1)
        depth = np.zeros_like(x)
        if self.d.bc[a][0] == "open":
            depth = np.maximum(depth, -x / npml)
        if self.d.bc[a][1] == "open":
            depth = np.maximum(depth, (x - N) / npml)
        depth = np.clip(depth, 0.0, 1.0)
        smax = -(order + 1) * math.log(R) / (2.0 * npml * self.dx[a])
        return smax * depth ** order

    def _init_pml(self):
        self.psiE = [cp.zeros(self.shape) for _ in range(6)]
        self.psiB = [cp.zeros(self.shape) for _ in range(6)]
        self.pml_E, self.pml_B = [], []
        for a in range(3):
            sn = self.pml_sigma(a, 0.0)
            sh = self.pml_sigma(a, 0.5)
            bn = np.exp(-sn * self.dt)
            bh = np.exp(-sh * self.dt * 0.5)
            self.pml_E += [cp.asarray(bn), cp.asarray(bn - 1.0)]
            self.pml_B += [cp.asarray(bh), cp.asarray(bh - 1.0)]

    # ------------------------------------------------------------- halo comms
    def _slab(self, F, a, lo, hi):
        """Slab along spatial axis a; F may be a single field or a stacked (k, nx, ny, nz) array."""
        sl = [slice(None)] * F.ndim
        sl[a + F.ndim - 3] = slice(lo, hi)
        return F[tuple(sl)]

    @staticmethod
    def _as_list(arrays):
        # a stacked 4-D array is exchanged as one block (one message per face for all components)
        return [arrays] if isinstance(arrays, cp.ndarray) else list(arrays)

    def fill_ghosts(self, arrays, axes=(0, 1, 2)):
        g = self.g
        arrays = self._as_list(arrays)
        for a in axes:
            ax = self.d.axes[a]
            n = ax.n_int
            for F in arrays:
                sends, recvs = [], []
                if ax.lo is not None:
                    sends.append((self._slab(F, a, g, 2 * g), ax.lo))
                if ax.hi is not None:
                    sends.append((self._slab(F, a, n, n + g), ax.hi))
                if ax.hi is not None:
                    recvs.append((self._slab(F, a, n + g, n + 2 * g), ax.hi))
                if ax.lo is not None:
                    recvs.append((self._slab(F, a, 0, g), ax.lo))
                if sends or recvs:
                    self.comm.exchange(sends, recvs)

    def reduce_ghosts(self, arrays):
        """Add ghost-cell contributions (current/charge deposited outside) to the owners."""
        g = self.g
        arrays = self._as_list(arrays)
        for a in range(3):
            ax = self.d.axes[a]
            n = ax.n_int
            for F in arrays:
                sends, recvs, adds = [], [], []
                if ax.lo is not None:
                    sends.append((self._slab(F, a, 0, g).copy(), ax.lo))
                if ax.hi is not None:
                    sends.append((self._slab(F, a, n + g, n + 2 * g).copy(), ax.hi))
                if ax.hi is not None:
                    b = cp.empty_like(self._slab(F, a, n, n + g))
                    recvs.append((b, ax.hi))
                    adds.append((self._slab(F, a, n, n + g), b))
                if ax.lo is not None:
                    b = cp.empty_like(self._slab(F, a, g, 2 * g))
                    recvs.append((b, ax.lo))
                    adds.append((self._slab(F, a, g, 2 * g), b))
                if sends or recvs:
                    self.comm.exchange(sends, recvs)
                for dst, b in adds:
                    dst += b
                # ghost contributions have been handed over
                if ax.lo is not None:
                    self._slab(F, a, 0, g)[...] = 0
                if ax.hi is not None:
                    self._slab(F, a, n + g, n + 2 * g)[...] = 0

    # ------------------------------------------------------------- filtering
    def filter(self, arrays, staggers, passes: int, compensate: bool = False):
        """Binomial (1,2,1)/4 filter, `passes` times along every axis (+ compensator), in place on
        the owned cells (PIC-D10). Ghosts come from the neighbours (periodic / rank faces); at a
        conducting wall node-type data are mirrored about the face node (face row (3/4, 1/4)) and
        wall-normal staggered data have a zero ghost, which makes the filter symmetric and commute
        with the discrete divergence. Call after reduce_ghosts and before apply_bc_J; the hi face
        node of node-type data is filtered too (wall deposits live there)."""
        if passes <= 0:
            return
        g = self.g
        w = [(0.25, 0.5, 0.25)] * passes
        if compensate:
            w.append((-passes / 4.0, 1.0 + passes / 2.0, -passes / 4.0))
        arrays = list(arrays)
        assert len(arrays) <= 3 and all(F.flags.c_contiguous and F.shape == self.shape for F in arrays)
        tmp = self._filter_tmp(len(arrays))
        for a in range(3):
            ax = self.d.axes[a]
            lo_wall = self.d.bc[a][0] != "periodic" and ax.lo is None
            hi_wall = self.d.bc[a][1] != "periodic" and ax.hi is None
            f = g + ax.n_int  # hi face node index (ghost for periodic / rank faces)
            his = [f + 1 if (hi_wall and st[a] == 0) else f for st in staggers]
            # passes alternate between the arrays and the work buffers (ping-pong, issue #9)
            src, dst = arrays, tmp
            for wts in w:
                self.fill_ghosts(src, axes=(a,))
                for F, st in zip(src, staggers):
                    node = st[a] == 0
                    if lo_wall:
                        self._slab(F, a, g - 1, g)[...] = self._slab(F, a, g, g + 1) if node else 0.0
                    if hi_wall:
                        if node:
                            self._slab(F, a, f + 1, f + 2)[...] = self._slab(F, a, f, f + 1)
                        else:
                            self._slab(F, a, f, f + 1)[...] = 0.0
                self._binomial_pass(src, dst, a, g, his, wts)
                src, dst = dst, src
            if src is tmp:  # the last pass wrote the buffers: copy its range back
                self._copy_pass(arrays, tmp, a, g, his)
            else:  # the last pass read the buffers: copy the layers it wrote there (ghosts, walls)
                for F, T, st in zip(arrays, tmp, staggers):
                    node = st[a] == 0
                    lo_l = (g - 1, g) if lo_wall else ((0, g) if ax.lo is not None else None)
                    hi_l = ((f + 1, f + 2) if node else (f, f + 1)) if hi_wall else (
                        (f, f + g) if ax.hi is not None else None)
                    for lay in (lo_l, hi_l):
                        if lay:
                            self._slab(F, a, *lay)[...] = self._slab(T, a, *lay)

    def _pass_args(self, a, lo, his, n):
        hs = [np.int32(h) for h in his] + [np.int32(lo)] * (3 - len(his))
        count = (max(his) - lo) * math.prod(self.shape) // self.shape[a]
        return count, (np.int32(n), np.int32(a), np.int32(lo), *hs)

    def _binomial_pass(self, src, dst, a, lo, his, wts):
        """dst[lo:hi] = wl src[-1] + wc src[0] + wr src[+1] along axis a, every other index, for
        up to 3 arrays in one kernel (issue #9); dst outside [lo, hi) is left as it is."""
        count, args = self._pass_args(a, lo, his, len(src))
        if count > 0:
            S = list(src) + [src[0]] * (3 - len(src))
            D = list(dst) + [dst[0]] * (3 - len(dst))
            launch(self.k_binomial, count, (*S, *D, *args, *(np.float64(v) for v in wts)))

    def _copy_pass(self, arrays, tmp, a, lo, his):
        count, args = self._pass_args(a, lo, his, len(arrays))
        if count > 0:
            F = list(arrays) + [arrays[0]] * (3 - len(arrays))
            T = list(tmp) + [tmp[0]] * (3 - len(tmp))
            launch(self.k_copy_pass, count, (*F, *T, *args))

    def _filter_tmp(self, k):
        if getattr(self, "_ftmp", None) is None:  # allocated on the first filtered step only
            self._ftmp = [cp.empty(self.shape, dtype=cp.float64) for _ in range(3)]
        return self._ftmp[:k]

    # --------------------------------------------------- boundary conditions
    def _wall_planes(self):
        """Yield (axis, boundary node index, ghost slice) for each global non-periodic face."""
        g = self.g
        for a in range(3):
            ax = self.d.axes[a]
            if self.d.bc[a][0] != "periodic":
                if ax.lo is None:
                    yield a, g, slice(0, g)
                if ax.hi is None:
                    yield a, g + ax.n_int, slice(g + ax.n_int, None)

    def apply_bc_E(self, arrays=None):
        """Perfect conductor at global non-periodic faces (outer edge of the PML for open faces):
        tangential E = 0 on the face, ghost values beyond the face = 0."""
        E = arrays if arrays is not None else (self.Ex, self.Ey, self.Ez)
        for a, ib, gs in self._wall_planes():
            for c in range(3):
                sl = [slice(None)] * 3
                sl[a] = gs
                E[c][tuple(sl)] = 0
                if c != a:
                    sl[a] = ib
                    E[c][tuple(sl)] = 0

    def apply_bc_B(self, arrays=None):
        """Normal B = 0 on the conducting face, ghost values beyond the face = 0."""
        B = arrays if arrays is not None else (self.Bx, self.By, self.Bz)
        for a, ib, gs in self._wall_planes():
            for c in range(3):
                sl = [slice(None)] * 3
                sl[a] = gs
                B[c][tuple(sl)] = 0
                if c == a:
                    sl[a] = ib
                    B[c][tuple(sl)] = 0

    def apply_bc_J(self, arrays=None):
        J = arrays if arrays is not None else (self.Jx, self.Jy, self.Jz)
        for a, _ib, gs in self._wall_planes():
            for c in range(3):
                sl = [slice(None)] * 3
                sl[a] = gs
                J[c][tuple(sl)] = 0

    # ------------------------------------------------------------ updates
    def push_B(self, dtb: float):
        args = [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz, np.float64(dtb)]
        if self.use_pml:
            args += self.psiB + self.pml_B
        launch(self.k_update_B, self.n_update, tuple(args))
        self.apply_bc_B()
        self.fill_ghosts((self.Bx, self.By, self.Bz))

    def push_E(self, dte: float):
        args = [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz, self.Jx, self.Jy, self.Jz, np.float64(dte)]
        if self.use_pml:
            args += self.psiE + self.pml_E
        launch(self.k_update_E, self.n_update, tuple(args))
        self.apply_bc_E()
        self.fill_ghosts((self.Ex, self.Ey, self.Ez))

    def curl_E(self, E, out, s=1.0):
        launch(self.k_curl_E, self.n_update, (*E, *out, np.float64(s)))

    def curl_B(self, B, out, s=1.0):
        launch(self.k_curl_B, self.n_update, (*B, *out, np.float64(s)))

    def zero_J(self):
        self.Jx.fill(0)
        self.Jy.fill(0)
        self.Jz.fill(0)

    # ---------------------------------------------------------- diagnostics
    def physical_slice(self, include_hi_face=False):
        g = self.g
        sl = []
        for ax in self.d.axes:
            lo = g + ax.pml_lo
            sl.append(slice(lo, lo + ax.n_phys))
        return tuple(sl)

    def owned_slice(self):
        g = self.g
        return tuple(slice(g, g + a.n_int) for a in self.d.axes)

    def field_energy(self, E=None, B=None, region: str = "physical"):
        """(electric, magnetic) energy in the physical domain, or in all owned cells including
        absorbing layers (region="all"); local sums + allreduce."""
        E = E if E is not None else (self.Ex, self.Ey, self.Ez)
        B = B if B is not None else (self.Bx, self.By, self.Bz)
        s = self.physical_slice() if region == "physical" else self.owned_slice()
        ee = sum(float(cp.sum(F[s] ** 2)) for F in E)
        eb = sum(float(cp.sum(F[s] ** 2)) for F in B)
        v = cp.asarray([ee, eb]) * 0.5 * self.cell_volume
        v = self.comm.allreduce(v)
        return float(v[0]), float(v[1])

    def poynting_out(self):
        """Instantaneous electromagnetic power leaving the physical domain through open faces."""
        g = self.g
        P = 0.0
        E = (self.Ex, self.Ey, self.Ez)
        B = (self.Bx, self.By, self.Bz)
        s = self.physical_slice()
        for a in range(3):
            ax = self.d.axes[a]
            b, c = (a + 1) % 3, (a + 2) % 3
            dA = self.dx[b] * self.dx[c]
            for side, ok in ((0, self.d.bc[a][0] == "open" and ax.lo is None),
                             (1, self.d.bc[a][1] == "open" and ax.hi is None)):
                if not ok:
                    continue
                ib = g + ax.pml_lo + (ax.n_phys if side else 0)
                sl0 = list(s)
                sl1 = list(s)
                sl0[a] = ib
                sl1[a] = ib - 1
                # S_a = E_b B_c - E_c B_b at the face node plane; B averaged across the face
                Eb, Ec = E[b][tuple(sl0)], E[c][tuple(sl0)]
                Bc = 0.5 * (B[c][tuple(sl0)] + B[c][tuple(sl1)])
                Bb = 0.5 * (B[b][tuple(sl0)] + B[b][tuple(sl1)])
                flux = float(cp.sum(Eb * Bc - Ec * Bb)) * dA
                P += flux if side else -flux
        return self.comm.allreduce_scalar(P)

    def div_E(self, E=None):
        E = E if E is not None else (self.Ex, self.Ey, self.Ez)
        out = cp.zeros(self.shape)
        launch(self.k_div_E, self.n_update, (*E, out))
        return out

    def div_B(self):
        out = cp.zeros(self.shape)
        launch(self.k_div_B, self.n_update, (self.Bx, self.By, self.Bz, out))
        return out

    def gauss_slice(self):
        """Nodes where the discrete Gauss law must hold: physical domain without conducting faces."""
        g = self.g
        sl = []
        for a, ax in enumerate(self.d.axes):
            lo = g + ax.pml_lo
            hi = lo + ax.n_phys
            if self.d.bc[a][0] != "periodic" and ax.lo is None:
                lo += 1
            sl.append(slice(lo, hi))
        return tuple(sl)

    def divb_slice(self):
        g = self.g
        sl = []
        for a, ax in enumerate(self.d.axes):
            lo = g + ax.pml_lo
            hi = lo + ax.n_phys
            if self.d.bc[a][1] != "periodic" and ax.hi is None:
                hi -= 1
            sl.append(slice(lo, hi))
        return tuple(sl)
