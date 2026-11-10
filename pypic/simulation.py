"""Main PIC loop: explicit (Yee + Buneman-Boris + Esirkepov) and implicit
(Crank-Nicolson + implicit midpoint, energy and charge conserving) schemes."""
from __future__ import annotations

import math
import time

import cupy as cp
import numpy as np

from .backend import compile_module, launch
from .config import SimConfig
from .grid import YeeGrid
from .kernelsrc import PARTICLE_SOURCES, particle_defines
from .parallel.comm import init_from_env
from .parallel.decomp import Decomposition, auto_dims
from .particles import ID_LIMIT, Species, migrate
from .profiles import invert_cos_density


class ImplicitStepError(RuntimeError):
    """An implicit step could not be completed; the state is left at the previous time level."""


class Simulation:
    def __init__(self, cfg: SimConfig, comm=None, verbose: bool = True):
        self.cfg = cfg
        self.comm = comm if comm is not None else init_from_env()
        self.rank = self.comm.rank
        self.verbose = verbose and self.rank == 0
        dims = cfg.dims or auto_dims(self.comm.size, cfg.n_cells)
        if int(np.prod(dims)) != self.comm.size:
            raise ValueError(f"dims {dims} do not match {self.comm.size} ranks")
        self.implicit = cfg.solver == "implicit"
        dmin = min(cfg.dx)
        ghost = 2 if not self.implicit else int(math.ceil(cfg.dt / dmin)) + 2
        if not self.implicit:
            cfl = cfg.dt * math.sqrt(sum(1.0 / h ** 2 for h in cfg.dx))
            if cfl >= 1.0:
                raise ValueError(f"explicit scheme violates the CFL condition (c dt sqrt(sum 1/dx^2) = {cfl:.3f} >= 1)")
        # Halo exchange and particle migration only reach the adjacent sub-domain: the ghost width
        # (and hence the largest particle displacement per step) must fit in every rank's cells.
        for a in range(3):
            n_loc = cfg.n_cells[a] // dims[a]
            if n_loc < ghost:
                raise ValueError(
                    f"axis {'xyz'[a]}: {n_loc} cells per rank < ghost width {ghost} "
                    f"(solver={cfg.solver}, c dt/dx={cfg.dt / cfg.dx[a]:.3g}); use more cells per rank, "
                    "fewer ranks along this axis or a smaller dt")
        if cfg.remove_mean_current and any(b != "periodic" for ax in cfg.boundaries for b in ax):
            raise ValueError("remove_mean_current needs periodic boundaries on every axis")
        self.n_filter = int(cfg.current_filter.passes)
        if self.n_filter and any(b == "open" for ax in cfg.boundaries for b in ax):
            raise ValueError("current_filter supports periodic and reflecting faces only")
        self.decomp = Decomposition.build(cfg.n_cells, dims, self.rank, cfg.boundaries, cfg.npml)
        self.grid = YeeGrid(self.decomp, cfg.dx, cfg.dt, self.comm, ghost=ghost, pml=not self.implicit)
        if cfg.precision not in ("fp64", "mixed"):
            raise ValueError("precision must be 'fp64' or 'mixed'")
        self.mixed = cfg.precision == "mixed"
        if self.mixed and cfg.solver == "implicit":
            raise ValueError("precision='mixed' is for the explicit solver only (PIC-D11)")
        self.pmod = compile_module(PARTICLE_SOURCES, particle_defines(self.grid.defines, mixed=self.mixed))
        self.k_push = self.pmod.get_function("push_deposit")
        self.k_rho = self.pmod.get_function("deposit_rho")
        self.k_kin = self.pmod.get_function("kinetic")
        self.k_cell = self.pmod.get_function("cell_index")
        self.k_load = self.pmod.get_function("load_particles")
        self.k_ipush = self.pmod.get_function("implicit_push")
        self.k_ifin = self.pmod.get_function("implicit_finalize")
        self.k_moment = self.pmod.get_function("deposit_moment")
        self.k_gather = self.pmod.get_function("gather_at_particles")
        self.k_displace = self.pmod.get_function("displace_deposit")
        self.k_kick = self.pmod.get_function("boris_kick")
        self.flag = cp.zeros(max(1, len(cfg.species)), dtype=cp.int32)  # per-species flag counters
        self.err = cp.zeros(2, dtype=cp.int32)
        self.itmax = cp.zeros(1, dtype=cp.int32)

        self.species: list[Species] = []
        self.sdiag = []
        self._load_species()
        for c, b in zip((self.grid.Bx, self.grid.By, self.grid.Bz), cfg.B0):
            c[...] = b
        self.grid.apply_bc_B()
        self.B0 = tuple(float(b) for b in cfg.B0)
        fields_set = self._set_initial_fields()
        moved = self._apply_perturbations()
        if any(b.profile is not None for sc in cfg.species for b in sc.beams):
            for sp in self.species:  # drop the particles a profile gave zero weight (collective)
                sp.compact_and_migrate(self.grid, self.comm)
        if not self.implicit and (moved or fields_set or any(self.B0)):
            self._initial_half_kick()
        if self.implicit:
            self._init_implicit()
        self.step_count = 0
        self.time = 0.0
        self.history = {k: [] for k in (
            "step", "t", "kinetic", "electric", "magnetic", "total", "lost_kinetic", "boundary_flux",
            "dissipated", "budget", "px", "py", "pz", "gauss", "divB", "picard", "n_particles")}
        self.kinetic_species = {s.name: [] for s in self.species}
        self.flux_out = 0.0
        self.dissipated = 0.0
        self._last = None
        n_tot = self.total_particles()  # collective: must run on every rank, not only the printing one
        if self.verbose:
            print(f"[pypic] {self.comm.size} GPU(s), rank grid {dims}, global cells {tuple(cfg.n_cells)}, "
                  f"solver={cfg.solver}, particles={n_tot:,}", flush=True)
        from .output.manager import OutputManager

        self.output = OutputManager(self)

    # ----------------------------------------------------------------- set-up
    def _load_species(self):
        """Load every beam (plasma species first, then test species, so adding a test species never
        changes the random numbers of the plasma), explicit test-particle lists, density profiles
        (weights) and check the charge neutrality slot by slot (all beams share positions)."""
        cfg = self.cfg
        n_beams = sum(len(sc.beams) for sc in cfg.species)
        n_glob = math.prod(cfg.n_cells)
        if n_beams * n_glob * cfg.ppc + sum(len(sc.particles) for sc in cfg.species) >= ID_LIMIT:
            raise ValueError("too many particles for exact float64 particle ids (n_beams * cells * ppc >= 2^53)")
        cells_local = int(np.prod([a.n_phys for a in self.decomp.axes]))
        order = [i for i, sc in enumerate(cfg.species) if not sc.test] + \
                [i for i, sc in enumerate(cfg.species) if sc.test]
        streams, stream = {}, 1
        for si in order:
            for bi in range(len(cfg.species[si].beams)):
                streams[(si, bi)] = stream
                stream += 1
        self._beam_ranges = {}
        list_id = float(n_beams * n_glob * cfg.ppc)  # explicit-list ids follow all beam ids
        for si, sc in enumerate(cfg.species):
            cap = int(1.3 * cells_local * cfg.ppc * max(1, len(sc.beams))) + 4096 + len(sc.particles)
            sp = Species(sc.name, sc.charge, sc.mass, cap)
            sp.test = bool(sc.test)
            sp.qdep = 0.0 if sp.test else sp.q
            for bi, b in enumerate(sc.beams):
                lo = sp.n
                sp.load(self.grid, self.k_load, b, cfg.ppc, cfg.seed, streams[(si, bi)])
                self._beam_ranges[(si, bi)] = (lo, sp.n)
                if b.profile is not None:
                    from .profiles import evaluate

                    a = "xyz".index(b.profile.get("axis", "x"))
                    f = evaluate(b.profile, sp.P[a, lo:sp.n], xp=cp)
                    if bool(cp.any(f < 0)):
                        raise ValueError(f"species {sc.name}: density profile is negative somewhere")
                    # negligible tails (sech^2, gauss) are cut: those particles get weight 0 and are
                    # removed; a fixed threshold keeps this decomposition independent and, for equal
                    # profiles, neutral
                    sp.P[6, lo:sp.n] *= cp.where(f < 1e-10, 0.0, f)
            if sc.particles:
                self._load_list(sp, sc.particles, list_id)
            list_id += len(sc.particles)
            self.species.append(sp)
            self.sdiag.append(cp.zeros(6))
        # neutrality: beams were loaded slot by slot on the same positions, so summing q w over the
        # beams of the plasma species position by position gives the net charge exactly
        plasma = [(si, bi) for si, sc in enumerate(cfg.species) if not sc.test for bi in range(len(sc.beams))]
        if plasma:
            n0 = {self._beam_ranges[k][1] - self._beam_ranges[k][0] for k in plasma}
            if len(n0) == 1:
                net = None
                for si, bi in plasma:
                    lo, hi = self._beam_ranges[(si, bi)]
                    v = self.species[si].q * self.species[si].P[6, lo:hi]
                    net = v if net is None else net + v
                scale = max(float(cp.max(cp.abs(net))) if net.size else 0.0, 0.0)
                wmax = max(float(cp.max(self.species[si].P[6, slice(*self._beam_ranges[(si, bi)])]))
                           if self._beam_ranges[(si, bi)][1] > self._beam_ranges[(si, bi)][0] else 0.0
                           for si, bi in plasma)
                bad = self.comm.allreduce_scalar(scale, "max") > 1e-12 * max(wmax, 1e-300)
                if bad and self.verbose:
                    print("[pypic] WARNING: the plasma is not charge neutral (net charge density != 0)", flush=True)

    def _load_list(self, sp, particles, id0: float):
        """Test particles given one by one: each rank keeps those inside its physical sub-domain."""
        arr = np.asarray(particles, dtype=float)
        L = self.grid.L
        if np.any(arr[:, :3] < 0) or np.any(arr[:, :3] >= np.asarray(L)):
            raise ValueError(f"species {sp.name}: test particles must lie inside [0, L)")
        mine = np.ones(len(arr), dtype=bool)
        for a, ax in enumerate(self.decomp.axes):
            lo, hi = ax.phys_lo * self.grid.dx[a], (ax.phys_lo + ax.n_phys) * self.grid.dx[a]
            mine &= (arr[:, a] >= lo) & (arr[:, a] < hi)
        idx = np.nonzero(mine)[0]
        m = len(idx)
        if not m:
            return
        sp.ensure_capacity(sp.n + m)
        blk = np.zeros((8, m))
        blk[:6] = arr[idx, :6].T
        blk[6] = arr[idx, 6] if arr.shape[1] == 7 else 1.0
        blk[7] = id0 + idx
        sp.P[:, sp.n: sp.n + m] = cp.asarray(blk)
        sp.n += m

    def _set_initial_fields(self) -> bool:
        """Uniform E0 and the initial_fields list (PIC-D9). Returns True if any field was set."""
        cfg, g = self.cfg, self.grid
        E0 = tuple(float(e) for e in cfg.E0)
        for c, e in enumerate(E0):
            if e != 0.0:
                for a in range(3):
                    if a != c and cfg.boundaries[a][0] != "periodic":
                        raise ValueError(f"E0 component {'xyz'[c]} is tangential to the conducting faces of "
                                         f"axis {'xyz'[a]}")
                getattr(g, "E" + "xyz"[c])[...] += e
        from .initial_fields import apply_initial_fields

        apply_initial_fields(self, cfg.initial_fields)
        g.apply_bc_E()
        g.apply_bc_B()
        g.fill_ghosts((g.Ex, g.Ey, g.Ez))
        g.fill_ghosts((g.Bx, g.By, g.Bz))
        return any(E0) or bool(cfg.initial_fields)

    def species_by_name(self, name: str) -> Species:
        for sp in self.species:
            if sp.name == name:
                return sp
        raise KeyError(f"no species {name!r} (have {[s.name for s in self.species]})")

    def _apply_perturbations(self) -> bool:
        """Beam perturbations (PIC-D9): velocity modes are added directly; position modes (density,
        displacement) move the particles in rounds of < 1 cell with the current of each move
        deposited by Esirkepov and migration in between, then E -= DT J. Since the loaded state is
        neutral (rho = 0, E = 0) and the deposit conserves charge, div E = rho holds exactly.
        Returns True if particles were moved (E is then non-zero)."""
        cfg, g = self.cfg, self.grid
        moving = False
        residual = 0.0  # largest |F(x) - x0| / L of the density inverses on this rank
        for si, (sc, sp) in enumerate(zip(cfg.species, self.species)):
            disp = None
            for bi, b in enumerate(sc.beams):
                lo, hi = self._beam_ranges[(si, bi)]
                for pt in b.perturbation:
                    a = "xyz".index(pt.get("axis", "x"))
                    if cfg.boundaries[a][0] != "periodic":
                        raise ValueError("perturbations need a periodic axis")
                    k = 2.0 * math.pi * int(pt.get("mode", 1)) / g.L[a]
                    ph = float(pt.get("phase", 0.0))
                    x0 = sp.P[a, lo:hi].copy()
                    if "velocity" in pt:
                        c = "xyz".index(pt.get("velocity_component", pt.get("axis", "x")))
                        sp.P[3 + c, lo:hi] += float(pt["velocity"]) * cp.sin(k * x0 + ph)
                    if "displacement" in pt or "density" in pt:
                        if disp is None:
                            disp = cp.zeros((3, sp.cap))
                        if "displacement" in pt:
                            disp[a, lo:hi] += float(pt["displacement"]) * cp.sin(k * x0 + ph)
                        else:  # exact inverse of the cumulative density (bracketed, converges for |alpha| < 1)
                            x, res = invert_cos_density(x0, float(pt["density"]), k, ph, cp)
                            residual = max(residual, res / g.L[a])
                            disp[a, lo:hi] += x - x0
            if disp is not None:
                sp.P = cp.concatenate([sp.P, disp], axis=0)  # rows 8..10: remaining displacement
                moving = True
        if not any(b.perturbation and any(("density" in p or "displacement" in p) for p in b.perturbation)
                   for sc in cfg.species for b in sc.beams):
            return False  # decided from the configuration: identical on every rank
        # every rank reaches this reduction (the branch above depends on the configuration only)
        residual = self.comm.allreduce_scalar(residual, "max")
        if residual > 1e-12:
            raise RuntimeError(f"density perturbation: the inverse of the cumulative density did not converge "
                               f"(max residual {residual:.3g} L)")
        h = [0.9 * d for d in g.dx]
        local = 0.0
        for sp in self.species:
            if sp.P.shape[0] > 8 and sp.n:
                local = max(local, max(float(cp.max(cp.abs(sp.P[8 + a, : sp.n]))) / h[a] for a in range(3)))
        rounds = int(math.ceil(self.comm.allreduce_scalar(local, "max")))
        g.zero_J()
        for _ in range(rounds):
            self.flag.fill(0)
            for i, sp in enumerate(self.species):
                if sp.P.shape[0] > 8:
                    lst = sp.flag_list()
                    launch(self.k_displace, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), np.float64(h[0]),
                                                   np.float64(h[1]), np.float64(h[2]), g.Jx, g.Jy, g.Jz,
                                                   np.float64(sp.qdep), self.flag[i:i + 1], lst, np.int32(lst.size)))
            self._migrate()
        g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
        self._filter_J()
        g.apply_bc_J()
        o = g.owned_slice()
        dt = cfg.dt
        for E, J in zip((g.Ex, g.Ey, g.Ez), (g.Jx, g.Jy, g.Jz)):
            E[o] -= dt * J[o]
        g.apply_bc_E()
        g.fill_ghosts((g.Ex, g.Ey, g.Ez))
        g.zero_J()
        for sp in self.species:
            if sp.P.shape[0] > 8:
                sp.P = cp.ascontiguousarray(sp.P[:8])
        return moving or rounds > 0

    def _initial_half_kick(self):
        """Explicit leapfrog: loaded momenta are at t = 0; with non-zero initial E or B push them
        back to t = -dt/2 by a Boris step of -dt/2 (positions unchanged)."""
        g = self.grid
        for sp in self.species:
            launch(self.k_kick, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), g.Ex, g.Ey, g.Ez, g.Bx, g.By, g.Bz,
                                       np.float64(sp.q), np.float64(sp.m), np.float64(-0.5 * self.cfg.dt)))

    def total_particles(self) -> int:
        return int(self.comm.allreduce_scalar(float(sum(s.n for s in self.species))))

    # ------------------------------------------------------------- explicit
    def _step_explicit(self, do_diag: bool):
        g = self.grid
        g.zero_J()
        self.flag.fill(0)
        if do_diag:
            for d in self.sdiag:
                d[:4] = 0
            self._fe = g.field_energy()
        fields = (g.Ex, g.Ey, g.Ez, g.Bx, g.By, g.Bz)
        if self.mixed:  # float32 copies for the gather (PIC-D11)
            if getattr(self, "_f32", None) is None:
                self._f32 = [cp.empty(g.shape, dtype=cp.float32) for _ in fields]
            for dst, src in zip(self._f32, fields):
                cp.copyto(dst, src, casting="same_kind")
            fields = self._f32
        for i, (sp, d) in enumerate(zip(self.species, self.sdiag)):
            lst = sp.flag_list()
            launch(self.k_push, sp.n, (
                sp.P, np.int32(sp.n), np.int64(sp.cap), *fields,
                g.Jx, g.Jy, g.Jz, np.float64(sp.q), np.float64(sp.m), np.float64(sp.qdep), d, self.flag[i:i + 1],
                lst, np.int32(lst.size), np.int32(do_diag)))
        g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
        self._filter_J()
        g.apply_bc_J()
        if self.cfg.remove_mean_current:
            self._remove_mean_current()
        self._migrate()
        dt = self.cfg.dt
        g.push_B(0.5 * dt)
        g.push_E(dt)
        g.push_B(0.5 * dt)
        if self._has_open:
            self.flux_out += dt * g.poynting_out()

    J_STAGGER = ((0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5))

    def _filter_J(self):
        """Binomial current filter (PIC-D10), between the ghost reduction and the wall BC."""
        if self.n_filter:
            g = self.grid
            g.filter((g.Jx, g.Jy, g.Jz), self.J_STAGGER, self.n_filter, self.cfg.current_filter.compensate)

    def _gather_E(self):
        """Electric field the implicit particles gather: E_h, or F E_h with the current filter F, so
        that sum E_h . (F J) = sum (F E_h) . J keeps the energy identity exact (F is symmetric)."""
        if not self.n_filter:
            return self.Eh
        g = self.grid
        G = [e.copy() for e in self.Eh]
        g.filter(G, self.J_STAGGER, self.n_filter, self.cfg.current_filter.compensate)
        g.fill_ghosts(G)
        return G

    def _remove_mean_current(self):
        """J -= <J> (box average, all ranks): div J is unchanged, so Gauss's law still holds; the
        mean E then stays exactly zero, so the removed current does no work (energy budget exact)."""
        g = self.grid
        o = g.owned_slice()
        Js = (g.Jx, g.Jy, g.Jz)
        mean = self.comm.allreduce(cp.stack([cp.sum(J[o]) for J in Js])) / math.prod(self.cfg.n_cells)
        for J, m in zip(Js, mean):
            J[o] -= m

    def _migrate(self):
        migrate(self.species, self.grid, self.comm, self.flag)

    @property
    def _has_open(self):
        return any(b == "open" for ax in self.cfg.boundaries for b in ax)

    # ------------------------------------------------------------- implicit
    def _init_implicit(self):
        g = self.grid
        for sp in self.species:
            sp.Pn = cp.zeros((6, sp.cap))
        shp = g.shape
        self.E1 = cp.zeros((3,) + shp)
        self.Eh = cp.zeros((3,) + shp)
        self.Bh = cp.zeros((3,) + shp)
        self.Eprev = None
        # absorbing (matched conductive) layer for open faces
        self.absorb = self._has_open
        if self.absorb:
            stag = {"Ex": (0.5, 0, 0), "Ey": (0, 0.5, 0), "Ez": (0, 0, 0.5),
                    "Bx": (0, 0.5, 0.5), "By": (0.5, 0, 0.5), "Bz": (0.5, 0.5, 0)}

            def sigma(name):
                s = stag[name]
                prof = [g.pml_sigma(a, s[a], order=2, R=1e-6) for a in range(3)]
                return cp.asarray(prof[0][:, None, None] + prof[1][None, :, None] + prof[2][None, None, :])

            self.sigE = [sigma(n) for n in ("Ex", "Ey", "Ez")]
            sigB = [sigma(n) for n in ("Bx", "By", "Bz")]
            self.sigB = sigB
            self.Dmat = [1.0 / (1.0 + 0.5 * self.cfg.dt * s) for s in sigB]
        # owned E unknowns for CG inner products
        self.own = tuple(slice(g.g, g.g + a.n_int) for a in self.decomp.axes)
        W = cp.zeros((3,) + shp)
        W[(slice(None),) + self.own] = 1.0
        g.apply_bc_E(W)  # conducting faces: tangential E is not an unknown
        self.Wmask = W
        # Jacobi preconditioner (diagonal of A)
        dt = self.cfg.dt
        h = g.dx
        diag = [1 + dt * dt / 4 * (2 / h[1] ** 2 + 2 / h[2] ** 2),
                1 + dt * dt / 4 * (2 / h[0] ** 2 + 2 / h[2] ** 2),
                1 + dt * dt / 4 * (2 / h[0] ** 2 + 2 / h[1] ** 2)]
        self.Adiag = cp.ones((3,) + shp)
        if self.absorb:
            for c in range(3):
                self.Adiag[c] += 0.5 * dt * self.sigE[c]
        self.Minv = cp.stack([1.0 / (d - 1.0 + self.Adiag[c]) for c, d in enumerate(diag)])

    _wdot = cp.ReductionKernel("T w, T a, T b", "T y", "w * a * b", "a + b", "y = a", "0", "pypic_wdot")

    def _dot(self, a, b):
        """Global inner product over owned, unconstrained E unknowns (stays on the device).
        Two stages (inner axes, then the rest): a full reduction to one scalar runs as a single
        thread block in CuPy's ReductionKernel, ~36x slower on a 1e6-cell grid (issue #5)."""
        return self.comm.allreduce(self._wdot(self.Wmask, a, b, axis=(2, 3)).sum().reshape(1)).reshape(())

    def _applyD(self, B):
        if self.absorb:
            for c in range(3):
                B[c] *= self.Dmat[c]

    def _matvec(self, v, out, tmp):
        """out = A v,  A = (1 + sigma dt/2) + dt^2/4 curl D curl  (v, out, tmp stacked (3, ...))."""
        g, dt = self.grid, self.cfg.dt
        g.apply_bc_E(v)
        g.fill_ghosts(v)
        g.curl_E(v, tmp)
        self._applyD(tmp)
        g.apply_bc_B(tmp)
        g.fill_ghosts(tmp)
        g.curl_B(tmp, out, dt * dt / 4.0)
        out += self.Adiag * v
        g.apply_bc_E(out)

    def _cg(self, x, b):
        """Jacobi-preconditioned CG for A x = b with A symmetric positive definite.
        x, b: stacked (3, nx, ny, nz) arrays. Scalars stay on the GPU; convergence is checked
        every few iterations to avoid host synchronisation."""
        icfg = self.cfg.implicit
        tmp = cp.empty_like(x)
        Ax = cp.empty_like(x)
        self._matvec(x, Ax, tmp)
        r = b - Ax
        self.grid.apply_bc_E(r)
        z = self.Minv * r
        p = z.copy()
        rz = self._dot(r, z)
        tol2 = (icfg.cg_tol ** 2) * float(self._dot(b, b))
        it = 0
        for it in range(icfg.cg_maxit):
            if it % 4 == 0 and float(self._dot(r, r)) <= tol2:
                break
            self._matvec(p, Ax, tmp)
            # The convergence test runs every few iterations only, so the residual can reach exactly
            # zero in between (e.g. a curl-free mode where A = I is solved in one step). Guard the
            # divisions on the device: a zero denominator means the step/update is zero.
            pAp = self._dot(p, Ax)
            alpha = cp.where(pAp != 0, rz / cp.where(pAp != 0, pAp, 1.0), 0.0)
            x += alpha * p
            r -= alpha * Ax
            cp.multiply(self.Minv, r, out=z)
            rz_new = self._dot(r, z)
            p *= cp.where(rz != 0, rz_new / cp.where(rz != 0, rz, 1.0), 0.0)
            p += z
            rz = rz_new
        g = self.grid
        g.apply_bc_E(x)
        g.fill_ghosts(x)
        return it

    def _step_implicit(self, do_diag: bool):
        g, cfg, icfg = self.grid, self.cfg, self.cfg.implicit
        dt = cfg.dt
        E0 = (g.Ex, g.Ey, g.Ez)
        B0 = (g.Bx, g.By, g.Bz)
        shp = g.shape
        # R0 = (1 - sigma dt/2) E0 + dt curlB( D (B0 - dt/4 curlE E0) )
        tB = cp.zeros((3,) + shp)
        g.curl_E(E0, tB, -dt / 4.0)
        for c in range(3):
            tB[c] += B0[c]
        self._applyD(tB)
        g.apply_bc_B(tB)
        g.fill_ghosts(tB)
        R0 = cp.zeros((3,) + shp)
        g.curl_B(tB, R0, dt)
        for c in range(3):
            R0[c] += (2.0 - self.Adiag[c]) * E0[c]  # (1 - sigma dt/2) E0
        # initial guess: linear extrapolation
        for c in range(3):
            if self.Eprev is not None:
                self.E1[c][...] = 2 * E0[c] - self.Eprev[c]
            else:
                self.E1[c][...] = E0[c]
        g.apply_bc_E(self.E1)
        g.fill_ghosts(self.E1)
        self.err.fill(0)
        rhs = cp.empty_like(R0)
        k = 0
        converged = False
        num, den = float("nan"), 1.0
        for k in range(icfg.picard_maxit):
            self._half_fields(E0, B0)
            G = self._gather_E()
            g.zero_J()
            self.itmax.fill(0)
            for sp in self.species:
                launch(self.k_ipush, sp.n, (
                    sp.P, sp.Pn, np.int32(sp.n), np.int64(sp.cap), *G, *self.Bh, g.Jx, g.Jy, g.Jz,
                    np.float64(sp.q), np.float64(sp.m), np.float64(sp.qdep), np.int32(1 if k == 0 else 0),
                    np.float64(icfg.particle_tol), np.int32(icfg.particle_maxit), self.err, self.itmax))
            g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
            self._filter_J()
            g.apply_bc_J()
            if cfg.remove_mean_current:
                self._remove_mean_current()
            for c, J in enumerate((g.Jx, g.Jy, g.Jz)):
                cp.subtract(R0[c], dt * J, out=rhs[c])
            Enew = self.E1.copy()
            self._cg(Enew, rhs)
            nd = cp.stack([cp.max(cp.abs(Enew - self.E1) * self.Wmask), cp.max(cp.abs(Enew) * self.Wmask)])
            num, den = (float(v) for v in self.comm.allreduce(nd, "max"))
            self.E1 = Enew
            converged = num <= icfg.picard_tol * max(den, 1e-30)
            if converged:
                break
        # Nothing has been committed yet (particles live in Pn, fields in E1), so raising here
        # leaves the simulation at step n. The error flags are per rank: reduce them so that every
        # rank takes the same branch (a rank-local raise would leave the others in a collective).
        errs = self.comm.allreduce(self.err, "max") if self.comm.size > 1 else self.err
        err = int(errs[0])
        if err == 2:
            raise ImplicitStepError("implicit step: particle displacement exceeds the ghost width (reduce dt)")
        if not converged:
            raise ImplicitStepError(
                f"implicit step {self.step_count}: Picard iteration did not converge in {icfg.picard_maxit} "
                f"iterations (relative change {num / max(den, 1e-300):.3e} > picard_tol={icfg.picard_tol:g}); "
                "the step was not committed. Reduce dt or raise implicit.picard_maxit.")
        if err == 3:
            raise ImplicitStepError("implicit step: energy closure impossible for a particle (gamma < 1); "
                                    "the step was not committed. Reduce dt.")
        self._picard_its = k + 1
        self.nonconverged = int(errs[1])
        # final half-step fields (consistent with converged E1)
        self._half_fields(E0, B0)
        if self.absorb:
            o = self.own
            loss = sum(float(cp.sum(sg[o] * e[o] * e[o])) for sg, e in zip(self.sigE, self.Eh))
            loss += sum(float(cp.sum(sg[o] * b[o] * b[o])) for sg, b in zip(self.sigB, self.Bh))
            self.dissipated += self.comm.allreduce_scalar(loss) * dt * g.cell_volume
        self.Eprev = [e.copy() for e in E0]
        for c in range(3):
            E0[c][...] = self.E1[c]
            B0[c][...] = 2.0 * self.Bh[c] - B0[c]
        g.apply_bc_E()
        g.apply_bc_B()
        g.fill_ghosts(E0)
        g.fill_ghosts(B0)
        self.flag.fill(0)
        for i, (sp, d) in enumerate(zip(self.species, self.sdiag)):
            lst = sp.flag_list()
            launch(self.k_ifin, sp.n, (sp.P, sp.Pn, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.q),
                                        np.float64(sp.m), np.float64(sp.qdep), d, self.flag[i:i + 1], lst,
                                        np.int32(lst.size)))
        self._migrate()

    def _half_fields(self, E0, B0):
        g, dt = self.grid, self.cfg.dt
        for c in range(3):
            self.Eh[c][...] = 0.5 * (E0[c] + self.E1[c])
        g.apply_bc_E(self.Eh)
        g.fill_ghosts(self.Eh)
        g.curl_E(self.Eh, self.Bh, -0.5 * dt)
        for c in range(3):
            self.Bh[c] += B0[c]
        self._applyD(self.Bh)
        g.apply_bc_B(self.Bh)
        g.fill_ghosts(self.Bh)

    # ------------------------------------------------------------ diagnostics
    def _kinetic_now(self):
        """Kinetic energy and momentum at the current integer time level (implicit scheme)."""
        out = []
        for sp, d in zip(self.species, self.sdiag):
            d[:4] = 0
            launch(self.k_kin, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.m), d))
            out.append(d[:4].copy())
        return out

    def gauss_error(self):
        g = self.grid
        rho = cp.zeros(g.shape)
        for sp in self.species:
            if not sp.test:
                launch(self.k_rho, sp.n, (sp.P, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.q), rho))
        g.reduce_ghosts((rho,))
        if self.n_filter:  # div E equals the identically filtered charge density
            g.filter((rho,), ((0, 0, 0),), self.n_filter, self.cfg.current_filter.compensate)
        div = g.div_E()
        s = g.gauss_slice()
        err = float(cp.max(cp.abs(div[s] - rho[s]))) if div[s].size else 0.0
        db = g.div_B()
        sb = g.divb_slice()
        errb = float(cp.max(cp.abs(db[sb]))) if db[sb].size else 0.0
        return self.comm.allreduce_scalar(err, "max"), self.comm.allreduce_scalar(errb, "max")

    def _record(self, kin_list, fe):
        h = self.history
        cfg = self.cfg
        # one allreduce and one host sync for all species: kinetic + momentum, lost energy, count
        S = len(self.species)
        vec = cp.concatenate([cp.stack(kin_list).ravel(), cp.stack([d[4] for d in self.sdiag]),
                              cp.asarray([float(sp.n) for sp in self.species])]) if S else cp.zeros(0)
        if self.comm.size > 1 and S:
            vec = self.comm.allreduce(vec)
        hv = cp.asnumpy(vec)
        kin, lost_s, n_s = hv[: 4 * S].reshape(S, 4), hv[4 * S: 5 * S], hv[5 * S:]
        # test species are recorded per species but are not part of the budget
        plasma = [k for k, sp in zip(kin, self.species) if not sp.test]
        ke = sum(float(k[0]) for k in plasma)
        mom = [sum(float(k[i]) for k in plasma) for i in (1, 2, 3)]
        lost = sum(float(v) for v, sp in zip(lost_s, self.species) if not sp.test)
        for sp, k in zip(self.species, kin):
            self.kinetic_species[sp.name].append(float(k[0]))
        tot = ke + fe[0] + fe[1]
        h["step"].append(self.step_count)
        h["t"].append(self.time)
        h["kinetic"].append(ke)
        h["electric"].append(fe[0])
        h["magnetic"].append(fe[1])
        h["total"].append(tot)
        h["lost_kinetic"].append(lost)
        h["boundary_flux"].append(self.flux_out)
        h["dissipated"].append(self.dissipated)
        h["budget"].append(tot + lost + self.flux_out + self.dissipated)
        h["px"].append(mom[0])
        h["py"].append(mom[1])
        h["pz"].append(mom[2])
        h["picard"].append(getattr(self, "_picard_its", 0))
        h["n_particles"].append(int(round(float(np.sum(n_s)))))
        if cfg.gauss_interval and self.step_count % cfg.gauss_interval == 0:
            ge, gb = self.gauss_error()
        else:
            ge, gb = float("nan"), float("nan")
        h["gauss"].append(ge)
        h["divB"].append(gb)

    # ------------------------------------------------------------------ run
    def step(self):
        cfg = self.cfg
        # file output of the state at time level n (collective; rank 0 writes)
        self.output.on_step()
        do_diag = cfg.diag_interval > 0 and self.step_count % cfg.diag_interval == 0
        if cfg.sort_interval and self.step_count % cfg.sort_interval == 0:
            for sp in self.species:
                sp.sort(self.k_cell)
        if self.implicit:
            if do_diag:
                kin = self._kinetic_now()
                # with an absorbing layer the field energy stored in the layer is part of the budget
                fe = self.grid.field_energy(region="all" if self.absorb else "physical")
                self._record(kin, fe)
            self._step_implicit(do_diag)
        else:
            self._step_explicit(do_diag)
            if do_diag:
                # kinetic energy centred at t^n (from the push), fields at t^n
                self._record([d[:4].copy() for d in self.sdiag], self._fe)
        self.step_count += 1
        self.time += cfg.dt

    def run(self, n_steps: int | None = None, progress: int = 0):
        n = self.cfg.n_steps if n_steps is None else n_steps
        t0 = time.perf_counter()
        for i in range(n):
            self.step()
            if progress and self.verbose and (i + 1) % progress == 0:
                cp.cuda.Device().synchronize()
                el = time.perf_counter() - t0
                h = self.history
                e0 = h["budget"][0] if h["budget"] else float("nan")
                msg = (f"step {self.step_count:7d} t={self.time:9.3f} "
                       f"{el / (i + 1) * 1e3:7.2f} ms/step")
                if h["budget"]:
                    msg += f" dE/E0={(h['budget'][-1] - e0) / e0:+.3e}"
                print(msg, flush=True)
        self.output.on_step()  # outputs due at the final time level
        self.output.flush()
        cp.cuda.Device().synchronize()
        self.wall_time = time.perf_counter() - t0
        self._run_steps = n
        # collective: every rank must reach this point (push_rate() itself is then local)
        self._n_total_end = self.total_particles()
        return self.history

    def push_rate(self) -> float:
        """GPU-ns per particle per step of the last run (wall time x number of GPUs / particles).
        Local (no communication), so it may be called on a single rank."""
        n = getattr(self, "_n_total_end", 0)
        steps = getattr(self, "_run_steps", self.cfg.n_steps)
        return self.wall_time / max(steps, 1) / max(n, 1) * 1e9 * self.comm.size

    def close(self):
        """Flush and close the output files (rank 0); the simulation can no longer write output."""
        self.output.close()
