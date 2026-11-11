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
        self.k_cg_matvec_end = self.mod.get_function("cg_matvec_end")
        self.k_cg_matvec = self.mod.get_function("cg_matvec")
        self.k_cg_scalars = self.mod.get_function("cg_scalars")
        self.k_cg_step = self.mod.get_function("cg_step")
        self.k_half_E = self.mod.get_function("half_E")
        self.k_half_B = self.mod.get_function("half_B")
        self.k_picard_change = self.mod.get_function("picard_change")
        self.k_div_E = self.mod.get_function("div_E")
        self.k_div_B = self.mod.get_function("div_B")
        self.k_marder_residual = self.mod.get_function("marder_residual")
        self.k_marder_apply = self.mod.get_function("marder_apply")
        self.k_binomial = self.mod.get_function("binomial_pass3")
        self.k_copy_pass = self.mod.get_function("copy_pass3")
        self.k_binomial_wrap = self.mod.get_function("binomial_wrap3")
        self.k_slab_copy = self.mod.get_function("slab_copy3")
        self.k_slab_zero = self.mod.get_function("slab_zero3")
        self.k_slab_pack = self.mod.get_function("slab_pack3")
        self.k_slab_wrap = self.mod.get_function("slab_wrap3")
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

    def _self_axis(self, ax):
        """The axis wraps onto this rank (periodic, not decomposed): ghosts are local copies."""
        r = self.comm.rank
        return ax.lo == r and ax.hi == r

    @staticmethod
    def _groups(arrays):
        """3-D arrays in groups of up to 3 (a stacked (k, nx, ny, nz) array counts as k arrays)."""
        flat = []
        for F in arrays:
            flat.extend(list(F) if F.ndim == 4 else [F])
        for F in flat:
            assert F.flags.c_contiguous
        return [flat[i:i + 3] for i in range(0, len(flat), 3)]

    def _slab_count(self, a, width, narr):
        return narr * width * (math.prod(self.shape) // self.shape[a])

    def _ptrs(self, grp):
        return tuple(grp) + (grp[0],) * (3 - len(grp))

    def _copy(self, grp, a, s0, d0, width, add):
        launch(self.k_slab_copy, self._slab_count(a, width, len(grp)),
               (*self._ptrs(grp), np.int32(len(grp)), np.int32(a), np.int32(s0), np.int32(d0), np.int32(width),
                np.int32(add)))

    def _zero(self, grp, a, s0, width):
        launch(self.k_slab_zero, self._slab_count(a, width, len(grp)),
               (*self._ptrs(grp), np.int32(len(grp)), np.int32(a), np.int32(s0), np.int32(width)))

    def _face_buf(self, key, n):
        """Persistent contiguous face buffer (all arrays' slabs of one face in one NCCL message);
        it stays referenced, as NCCL requires until groupEnd()."""
        bufs = self.__dict__.setdefault("_face_bufs", {})
        b = bufs.get(key)
        if b is None or b.size < n:
            b = bufs[key] = cp.empty(n, dtype=cp.float64)
        return b[:n]

    def _pack(self, grp, a, starts, width, buf, unpack=False, add=False):
        """Pack (or unpack) the slabs starting at ``starts`` (one or two faces) in one launch."""
        s = [*starts, -1]
        launch(self.k_slab_pack, buf.size,
               (*self._ptrs(grp), buf, np.int32(len(grp)), np.int32(a), np.int32(s[0]), np.int32(s[1]),
                np.int32(width), np.int32(unpack), np.int32(add)))

    def _wrap(self, grp, mask):
        """Ghost layers of the axes of ``mask`` (they wrap onto this rank) <- periodic images."""
        g, count, owned = self.g, 0, list(self.shape)
        for a in range(3):
            if mask >> a & 1:
                count += 2 * g * math.prod(owned) // owned[a]
                owned[a] -= 2 * g
        launch(self.k_slab_wrap, len(grp) * count, (*self._ptrs(grp), np.int32(len(grp)), np.int32(mask)))

    def _exchange_faces(self, grp, a, kind, send_lo, send_hi, recv_hi, recv_lo, add):
        """Pack the send slabs of every array (both faces in one launch), exchange one message per
        face, unpack (both faces in one launch unless they add into the same cells)."""
        ax, g = self.d.axes[a], self.g
        n = self._slab_count(a, g, len(grp))
        out = [(s, r) for s, r in ((send_lo, ax.lo), (send_hi, ax.hi)) if r is not None]
        inc = [(d, r) for d, r in ((recv_hi, ax.hi), (recv_lo, ax.lo)) if r is not None]
        sbuf = self._face_buf((kind, a, "s"), len(out) * n)
        rbuf = self._face_buf((kind, a, "r"), len(inc) * n)
        self._pack(grp, a, [s for s, _ in out], g, sbuf)
        self.comm.exchange([(sbuf[i * n:(i + 1) * n], r) for i, (_, r) in enumerate(out)],
                           [(rbuf[i * n:(i + 1) * n], r) for i, (_, r) in enumerate(inc)])
        starts = [d for d, _ in inc]
        if add and len(starts) == 2 and ax.n_int < 2 * g:  # [n, n+g) and [g, 2g) overlap: one by one
            for i, d in enumerate(starts):
                self._pack(grp, a, [d], g, rbuf[i * n:(i + 1) * n], unpack=True, add=True)
        else:
            self._pack(grp, a, starts, g, rbuf, unpack=True, add=add)

    def fill_ghosts(self, arrays, axes=(0, 1, 2)):
        """Copy the neighbours' boundary layers into the ghost layers. Fused kernels handle up to 3
        arrays per launch, and a face carries all of them in one message (issue #8). The axes that
        wrap onto this rank are filled last, in one launch, from the periodic images: the values the
        axis-by-axis sequence gives, edges and corners included."""
        g = self.g
        for grp in self._groups(self._as_list(arrays)):
            mask = 0
            for a in axes:
                ax = self.d.axes[a]
                if ax.lo is None and ax.hi is None:
                    continue
                if self._self_axis(ax):
                    mask |= 1 << a
                else:
                    self._exchange_faces(grp, a, "fill", g, ax.n_int, ax.n_int + g, 0, add=False)
            if mask:
                self._wrap(grp, mask)

    def reduce_ghosts(self, arrays):
        """Add ghost-cell contributions (current/charge deposited outside) to the owners, then zero
        the ghosts (same fused kernels and messages as fill_ghosts)."""
        g = self.g
        groups = self._groups(self._as_list(arrays))
        for a in range(3):
            ax = self.d.axes[a]
            n = ax.n_int
            if ax.lo is None and ax.hi is None:
                continue
            for grp in groups:
                if self._self_axis(ax):
                    self._copy(grp, a, 0, n, g, True)  # [n, n+g) += lo ghosts
                    self._copy(grp, a, n + g, g, g, True)  # [g, 2g) += hi ghosts
                else:
                    self._exchange_faces(grp, a, "reduce", 0, n + g, n, g, add=True)
                # ghost contributions have been handed over
                if ax.lo is not None:
                    self._zero(grp, a, 0, g)
                if ax.hi is not None:
                    self._zero(grp, a, n + g, g)

    # ------------------------------------------------------------- filtering
    def filter(self, arrays, staggers, passes: int, compensate: bool = False):
        """Binomial (1,2,1)/4 filter, `passes` times along every axis (+ compensator), in place on
        the owned cells (PIC-D10). Ghosts come from the neighbours (periodic / rank faces); at a
        conducting wall node-type data are mirrored about the face node (face row (3/4, 1/4)) and
        wall-normal staggered data have a zero ghost, which makes the filter symmetric and commute
        with the discrete divergence. Call after reduce_ghosts and before apply_bc_J; the hi face
        node of node-type data is filtered too (wall deposits live there). An open face is a wall
        of the same kind at the face of the box (PIC-D15): the absorbing layer behind it holds no
        deposit (the kernels cut the paths at the face) and is not filtered. Along an axis that
        wraps onto this rank the passes are applied as one composite stencil: the same filter,
        with the rounding of one sum instead of one per pass."""
        if passes <= 0:
            return
        g = self.g
        w = [(0.25, 0.5, 0.25)] * passes
        if compensate:
            w.append((-passes / 4.0, 1.0 + passes / 2.0, -passes / 4.0))
        arrays = list(arrays)
        assert len(arrays) <= 3 and all(F.flags.c_contiguous and F.shape == self.shape for F in arrays)
        tmp = self._filter_tmp(len(arrays))
        cur = arrays  # where the data are: the arrays, or the work buffers after a wrapped axis
        for a in range(3):
            ax = self.d.axes[a]
            if self._self_axis(ax):  # every pass in one launch, no ghosts needed (issue #37)
                nxt = tmp if cur is arrays else arrays
                self._wrapped_passes(cur, nxt, a, ax.n_int, w)
                cur = nxt
                continue
            if cur is tmp:
                self._copy_pass(arrays, tmp, a, 0, [self.shape[a]] * len(arrays))
                cur = arrays
            lo_wall = self.d.bc[a][0] != "periodic" and ax.lo is None
            hi_wall = self.d.bc[a][1] != "periodic" and ax.hi is None
            # the walls: the ends of the array, or the faces of the box in front of the layers
            lo = g + ax.pml_lo  # lo face node index (first owned index at a rank face)
            f = g + ax.n_int - ax.pml_hi  # hi face node index (ghost for periodic / rank faces)
            his = [f + 1 if (hi_wall and st[a] == 0) else f for st in staggers]
            # passes alternate between the arrays and the work buffers (ping-pong, issue #9)
            src, dst = arrays, tmp
            for wts in w:
                self.fill_ghosts(src, axes=(a,))
                for F, st in zip(src, staggers):
                    node = st[a] == 0
                    if lo_wall:
                        self._slab(F, a, lo - 1, lo)[...] = self._slab(F, a, lo, lo + 1) if node else 0.0
                    if hi_wall:
                        if node:
                            self._slab(F, a, f + 1, f + 2)[...] = self._slab(F, a, f, f + 1)
                        else:
                            self._slab(F, a, f, f + 1)[...] = 0.0
                self._binomial_pass(src, dst, a, lo, his, wts)
                src, dst = dst, src
            if src is tmp:  # the last pass wrote the buffers: copy its range back
                self._copy_pass(arrays, tmp, a, lo, his)
            else:  # the last pass read the buffers: copy the layers it wrote there (ghosts, walls)
                for F, T, st in zip(arrays, tmp, staggers):
                    node = st[a] == 0
                    lo_l = (lo - 1, lo) if lo_wall else ((0, g) if ax.lo is not None else None)
                    hi_l = ((f + 1, f + 2) if node else (f, f + 1)) if hi_wall else (
                        (f, f + g) if ax.hi is not None else None)
                    for lay in (lo_l, hi_l):
                        if lay:
                            self._slab(F, a, *lay)[...] = self._slab(T, a, *lay)
            # the mirror rows of node-type data lie in the layer of an open face: keep it empty
            for F, st in zip(arrays, staggers):
                if st[a] == 0 and lo_wall and ax.pml_lo:
                    self._slab(F, a, lo - 1, lo)[...] = 0.0
                if st[a] == 0 and hi_wall and ax.pml_hi:
                    self._slab(F, a, f + 1, f + 2)[...] = 0.0
        if cur is tmp:
            self._copy_pass(arrays, tmp, 0, 0, [self.shape[0]] * len(arrays))

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

    def _wrapped_passes(self, src, dst, a, n, w):
        """dst = all the passes `w` along the self-periodic axis a applied to src, at every index
        (the ghost layers of the axis as periodic images). The composite stencil is the convolution
        of the passes; its coefficients are exact (dyadic) for the binomial passes and compensator."""
        if getattr(self, "_fw", None) is None or self._fw[0] != w:
            c = np.ones(1)
            for wts in w:
                c = np.convolve(c, wts)
            self._fw = (w, cp.asarray(c[len(w):]))  # centre and one side: the stencil is symmetric
        coef = self._fw[1]
        S = list(src) + [src[0]] * (3 - len(src))
        D = list(dst) + [dst[0]] * (3 - len(dst))
        launch(self.k_binomial_wrap, math.prod(self.shape),
               (*S, *D, np.int32(len(src)), np.int32(a), np.int32(n), np.int32(len(w)), coef))

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
        """B -= dtb curl E on the update range and one ghost layer around it, from the exchanged
        ghosts of E (two layers): the neighbours' values of B, bit for bit, without an exchange.
        The outermost ghost layer of B is not updated; nothing reads it."""
        args = [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz, np.float64(dtb)]
        if self.use_pml:
            args += self.psiB + self.pml_B
        launch(self.k_update_B, math.prod(a.n_int + 2 for a in self.d.axes), tuple(args))
        self.apply_bc_B()

    def push_E(self, dte: float):
        args = [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz, self.Jx, self.Jy, self.Jz, np.float64(dte)]
        if self.use_pml:
            args += self.psiE + self.pml_E
        launch(self.k_update_E, self.n_update, tuple(args))
        self.apply_bc_E()
        self.fill_ghosts((self.Ex, self.Ey, self.Ez))

    def curl_E(self, E, out, s=1.0, low=False):
        """out = s curl E on the update range. ``low``: also on the ghost layer below it on every
        axis, computed from the ghosts of E (what curl_B of ``out`` reads)."""
        n = math.prod(a.n_int + 1 for a in self.d.axes) if low else self.n_update
        launch(self.k_curl_E, n, (*E, *out, np.float64(s), np.int32(low)))

    def curl_B(self, B, out, s=1.0):
        launch(self.k_curl_B, self.n_update, (*B, *out, np.float64(s)))

    def shift_x(self):
        """Moving window (PIC-D22): E, B and the memory of the PML go one plane down along x. The
        plane that enters at the upper end is zero, but for E_x of the new wall cell, which takes
        the value of the cell below it: div E = 0 at the node between them (the caller refuses
        charge there). B, and the memory of the PML, are valid afterwards with their ghosts; the
        caller ends with ``apply_bc_E()`` and ``fill_ghosts`` of E, after it has taken the current
        of the plane it loads out of E."""
        ax = self.d.axes[0]
        # The psi of a layer across the window lies in every rank of a row along x, and the plane
        # above comes from the neighbour (psiE is updated on the owned range only). The psi of the
        # layer at x- is zero above that layer: nothing to fetch, and the neighbour may have no
        # layer at all (the two sides of a face must make the same exchange)
        psi_x = self.use_pml and any(a.pml_lo or a.pml_hi for a in self.d.axes[1:])
        if psi_x:
            self.fill_ghosts(self.psiE, axes=(0,))
        if getattr(self, "_shift_tmp", None) is None:
            self._shift_tmp = cp.empty(self.shape)

        def down(arr):
            out = self._shift_tmp
            out[:-1] = arr[1:]
            out[-1] = 0
            self._shift_tmp = arr
            return out

        for name in FIELD_NAMES:
            setattr(self, name, down(getattr(self, name)))
        if self.use_pml:
            self.psiE = [down(p) for p in self.psiE]
            self.psiB = [down(p) for p in self.psiB]
        if ax.hi is None:
            w = self.g + ax.n_int
            self.Ex[w - 1] = self.Ex[w - 2]
        # the outermost ghost plane of B and of psiB is not kept up to date by push_B
        self.apply_bc_B()
        self.fill_ghosts((self.Bx, self.By, self.Bz), axes=(0,))
        if psi_x:
            self.fill_ghosts(self.psiB, axes=(0,))

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

    _sum_sq = cp.ReductionKernel("T x, T y, T z", "T s", "x * x + y * y + z * z", "a + b", "s = a", "0",
                                 "pypic_sum_sq")

    def field_energy_local(self, E=None, B=None, region: str = "physical"):
        """(electric, magnetic) energy of this rank in the physical domain, or in all owned cells
        including absorbing layers (region="all"), as a device array: no host synchronisation and
        no collective, the caller reduces it over the ranks."""
        E = E if E is not None else (self.Ex, self.Ey, self.Ez)
        B = B if B is not None else (self.Bx, self.By, self.Bz)
        s = self.physical_slice() if region == "physical" else self.owned_slice()
        # two stages: a full reduction to one scalar runs as a single thread block (issue #5)
        sums = [self._sum_sq(*(F[s] for F in V), axis=(1, 2)).sum() for V in (E, B)]
        return cp.stack(sums) * (0.5 * self.cell_volume)

    def field_energy(self, E=None, B=None, region: str = "physical"):
        """Global (electric, magnetic) energy: field_energy_local + allreduce."""
        v = cp.asnumpy(self.comm.allreduce(self.field_energy_local(E, B, region)))
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

    def gauss_slice(self, retreat: int = 0):
        """Nodes where the discrete Gauss law must hold: physical domain without conducting faces,
        and without the ``retreat`` nodes next to an open face that the filtered charge of the
        removed and entering particles and the Marder zone reach (PIC-D15, PIC-D16)."""
        g = self.g
        sl = []
        for a, ax in enumerate(self.d.axes):
            b, n = self.d.bc[a], self.d.n_cells[a]
            first = 0 if b[0] == "periodic" else 1 + (retreat if b[0] == "open" else 0)
            end = n - (retreat if b[1] == "open" else 0)  # the hi face node is not an owned node
            lo = max(ax.phys_lo, first)
            hi = max(lo, min(ax.phys_lo + ax.n_phys, end))
            off = g + ax.pml_lo - ax.phys_lo
            sl.append(slice(lo + off, hi + off))
        return tuple(sl)

    def marder_zone(self, rho, depth: int, number: float):
        """One Marder step in the zone of the open faces (PIC-D16): E += nu h^2 grad(div E - rho) on
        the edges whose two nodes are in the zone, with h^2 = 3 / sum(1/dx^2), so that the residual
        diffuses with the number ``number`` (<= 1/6) per step. The zone holds the absorbing layer,
        the face node and ``depth`` nodes inside, at every open face; the conducting walls take the
        residual (0 there). ``rho``: the charge density that div E must equal (ghosts reduced,
        filtered like the current); it is overwritten."""
        if getattr(self, "_zone", None) is None or self._zone[0] != depth:
            big = 1 << 30
            zb = np.empty(6, dtype=np.int32)
            for a, ax in enumerate(self.d.axes):
                off = self.g + ax.pml_lo - ax.phys_lo  # local index of the global node 0
                zb[a] = off + depth if self.d.bc[a][0] == "open" else -big
                zb[3 + a] = off + self.d.n_cells[a] - depth if self.d.bc[a][1] == "open" else big
            self._zone = (depth, cp.asarray(zb))
        zb = self._zone[1]
        launch(self.k_marder_residual, self.n_update, (self.Ex, self.Ey, self.Ez, rho, zb))
        self.fill_ghosts((rho,))
        for a, ib, gs in self._wall_planes():
            sl = [slice(None)] * 3
            sl[a] = gs
            rho[tuple(sl)] = 0
            sl[a] = ib
            rho[tuple(sl)] = 0
        h2 = 3.0 / sum(1.0 / d ** 2 for d in self.dx)
        launch(self.k_marder_apply, self.n_update, (
            self.Ex, self.Ey, self.Ez, rho, *(np.float64(number * h2 / d) for d in self.dx), zb))
        self.apply_bc_E()
        self.fill_ghosts((self.Ex, self.Ey, self.Ez))

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
