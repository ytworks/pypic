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
from .particles import ID_LIMIT, Species


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
        self.decomp = Decomposition.build(cfg.n_cells, dims, self.rank, cfg.boundaries, cfg.npml)
        self.grid = YeeGrid(self.decomp, cfg.dx, cfg.dt, self.comm, ghost=ghost, pml=not self.implicit)
        self.pmod = compile_module(PARTICLE_SOURCES, particle_defines(self.grid.defines))
        self.k_push = self.pmod.get_function("push_deposit")
        self.k_rho = self.pmod.get_function("deposit_rho")
        self.k_kin = self.pmod.get_function("kinetic")
        self.k_cell = self.pmod.get_function("cell_index")
        self.k_load = self.pmod.get_function("load_particles")
        self.k_ipush = self.pmod.get_function("implicit_push")
        self.k_ifin = self.pmod.get_function("implicit_finalize")
        self.k_moment = self.pmod.get_function("deposit_moment")
        self.k_gather = self.pmod.get_function("gather_at_particles")
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
        cfg = self.cfg
        n_beams = sum(len(sc.beams) for sc in cfg.species)
        if n_beams * math.prod(cfg.n_cells) * cfg.ppc >= ID_LIMIT:
            raise ValueError("too many particles for exact float64 particle ids (n_beams * cells * ppc >= 2^53)")
        stream = 1
        cells_local = int(np.prod([a.n_phys for a in self.decomp.axes]))
        for sc in cfg.species:
            cap = int(1.3 * cells_local * cfg.ppc * max(1, len(sc.beams))) + 4096
            sp = Species(sc.name, sc.charge, sc.mass, cap)
            for b in sc.beams:
                sp.load(self.grid, self.k_load, b, cfg.ppc, cfg.seed, stream)
                stream += 1
            self.species.append(sp)
            self.sdiag.append(cp.zeros(6))
        # global charge neutrality check
        net = sum(sc.charge * b.density for sc in cfg.species for b in sc.beams)
        if abs(net) > 1e-12 and self.verbose:
            print(f"[pypic] WARNING: net charge density {net:+.3e} n0 e (plasma is not neutral)")

    def species_by_name(self, name: str) -> Species:
        for sp in self.species:
            if sp.name == name:
                return sp
        raise KeyError(f"no species {name!r} (have {[s.name for s in self.species]})")

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
        for i, (sp, d) in enumerate(zip(self.species, self.sdiag)):
            lst = sp.flag_list()
            launch(self.k_push, sp.n, (
                sp.P, np.int32(sp.n), np.int32(sp.cap), g.Ex, g.Ey, g.Ez, g.Bx, g.By, g.Bz,
                g.Jx, g.Jy, g.Jz, np.float64(sp.q), np.float64(sp.m), d, self.flag[i:i + 1],
                lst, np.int32(lst.size), np.int32(do_diag)))
        g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
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
        counts = cp.asnumpy(self.flag)  # one host sync per step
        if self.comm.size > 1:
            nflag = int(self.comm.allreduce_scalar(float(counts.sum())))
        else:
            nflag = int(counts.sum())
        if nflag > 0:
            for sp, k in zip(self.species, counts):
                sp.remove_and_migrate(self.grid, self.comm, int(k))

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
        """Global inner product over owned, unconstrained E unknowns (stays on the device)."""
        return self.comm.allreduce(self._wdot(self.Wmask, a, b).reshape(1)).reshape(())

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
            g.zero_J()
            self.itmax.fill(0)
            for sp in self.species:
                launch(self.k_ipush, sp.n, (
                    sp.P, sp.Pn, np.int32(sp.n), np.int32(sp.cap), *self.Eh, *self.Bh, g.Jx, g.Jy, g.Jz,
                    np.float64(sp.q), np.float64(sp.m), np.int32(1 if k == 0 else 0),
                    np.float64(icfg.particle_tol), np.int32(icfg.particle_maxit), self.err, self.itmax))
            g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
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
            launch(self.k_ifin, sp.n, (sp.P, sp.Pn, np.int32(sp.n), np.int32(sp.cap), np.float64(sp.q),
                                        np.float64(sp.m), d, self.flag[i:i + 1], lst, np.int32(lst.size)))
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
            launch(self.k_kin, sp.n, (sp.P, np.int32(sp.n), np.int32(sp.cap), np.float64(sp.m), d))
            out.append(d[:4].copy())
        return out

    def gauss_error(self):
        g = self.grid
        rho = cp.zeros(g.shape)
        for sp in self.species:
            launch(self.k_rho, sp.n, (sp.P, np.int32(sp.n), np.int32(sp.cap), np.float64(sp.q), rho))
        g.reduce_ghosts((rho,))
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
        kin = [self.comm.allreduce(k) for k in kin_list]
        ke = sum(float(k[0]) for k in kin)
        mom = [sum(float(k[i]) for k in kin) for i in (1, 2, 3)]
        lost = sum(float(d[4]) for d in self.sdiag)
        lost = self.comm.allreduce_scalar(lost)
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
        h["n_particles"].append(self.total_particles())
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
