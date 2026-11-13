"""Main PIC loop: explicit (Yee + Buneman-Boris + Esirkepov) and implicit
(Crank-Nicolson + implicit midpoint, energy and charge conserving) schemes."""
from __future__ import annotations

import itertools
import math
import time

import cupy as cp
import numpy as np

from . import antenna, body, checkpoint, injection, window
from .backend import compile_module, launch
from .config import FACES, ID_LIMIT, SLOT_BYTES, SimConfig
from .external import check_config, dipoles_B
from .grid import FIELD_NAMES, YeeGrid
from .kernelsrc import CG_BLOCK, PARTICLE_SOURCES, particle_defines
from .output.names import stagger_of
from .parallel.comm import init_from_env
from .parallel.decomp import Decomposition, auto_dims
from .particles import MIGRATE_KERNELS, CompactSpecies, Species, migrate
from .profiles import CUT, evaluate, invert_cos_density, support, support_interval

NDIAG = 12  # entries of the diagnostics buffer of a species (Simulation.sdiag)
# The shift of a moving window (PIC-D22) is a device global of the particle module: every
# simulation with a window compiles a module of its own (the number goes into the source).
_WINDOW_MODULES = itertools.count(1)


class ImplicitStepError(RuntimeError):
    """An implicit step could not be completed; the state is left at the previous time level."""


class Simulation:
    def __init__(self, cfg: SimConfig, comm=None, verbose: bool = True):
        self.cfg = cfg
        self.comm = comm if comm is not None else init_from_env()
        self.rank = self.comm.rank
        self.verbose = verbose and self.rank == 0
        checkpoint.check_config(cfg)
        # restart (PIC-D25): the checkpoint is read and its configuration checked before the set-up
        self._restart = checkpoint.open_restart(self) if cfg.restart else None
        self._restart_step = self._restart[1]["step"] if self._restart else None
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
        has_open = any(b == "open" for ax in cfg.boundaries for b in ax)
        if self.n_filter and has_open and self.implicit:
            raise ValueError("current_filter with an open face is for the explicit solver only (PIC-D15)")
        self.marder = float(cfg.marder_zone)
        if self.marder and (self.implicit or not has_open):
            raise ValueError("marder_zone needs the explicit solver and an open face (PIC-D16)")
        # Nodes from an open face that the charge left on the face nodes reaches through the filter
        # (PIC-D15); the Marder zone holds one node more, and the Gauss check starts behind them.
        reach = (self.n_filter + int(bool(cfg.current_filter.compensate))) if self.n_filter else 0
        self.zone_depth = reach + 1
        self.gauss_retreat = (reach + (1 if self.marder else 0)) if has_open else 0
        self.body = body.check_config(cfg, reach)
        antenna.check_config(cfg, self.gauss_retreat, reach)
        self.window = window.check_config(cfg, self.gauss_retreat, reach)
        self.decomp = Decomposition.build(cfg.n_cells, dims, self.rank, cfg.boundaries, cfg.npml)
        self.grid = YeeGrid(self.decomp, cfg.dx, cfg.dt, self.comm, ghost=ghost, pml=not self.implicit)
        if cfg.precision not in SLOT_BYTES:
            raise ValueError("precision must be 'fp64', 'mixed' or 'compact'")
        self.compact = cfg.precision == "compact"  # the storage of the particles (PIC-D19)
        self.mixed = cfg.precision != "fp64"  # the push gathers from float32 copies of the fields (PIC-D11)
        if self.mixed and cfg.solver == "implicit":
            raise ValueError(f"precision='{cfg.precision}' is for the explicit solver only (PIC-D11, PIC-D19)")
        self.injection = injection.check_config(cfg)
        # PIC-D26; checked even without an injecting face, which refuses a table that nothing carries in
        self.fluct = injection.check_fluctuation(cfg)
        self.pmod = compile_module(PARTICLE_SOURCES, particle_defines(
            self.grid.defines, mixed=self.mixed, paired_faces=injection.paired_faces(cfg), body=cfg.body,
            window=next(_WINDOW_MODULES) if self.window else 0, x_low=self.decomp.x_low(ghost),
            lattice=cfg.lattice, compact=self.compact))
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
        self.k_extsum = self.pmod.get_function("ext_sum")
        self.k_inject = self.pmod.get_function("inject_face")
        self.k_partners = self.pmod.get_function("inject_partners")
        self.k_mark_body = self.pmod.get_function("mark_body")
        self.k_flag_dead = self.pmod.get_function("flag_dead")
        self.k_migrate = tuple(self.pmod.get_function(name) for name in MIGRATE_KERNELS)
        # per species: [s] particles flagged for removal or migration, [S + s] particles injected
        self.flag = cp.zeros(2 * max(1, len(cfg.species)), dtype=cp.int32)
        self.err = cp.zeros(2, dtype=cp.int32)
        self.window_shift = 0  # cells that the moving window has gone (PIC-D22)
        if self.window:
            self.k_wset, self.k_wdrop, self.k_wmove = (
                self.pmod.get_function(name) for name in ("window_set", "window_drop", "window_move"))
            self._window_set()  # before any kernel that takes a position
            self._window_end = window.shift_after(cfg, cfg.n_steps)
            self._antenna_off = [window.antenna_off_shift(cfg, laser, self.gauss_retreat) for laser in cfg.lasers]
            # field energy that the shifts took out of this rank's physical cells (only the sum
            # over the ranks means something); per species, cumulative: [0] kinetic energy,
            # [1] charge and [2] number dropped at x-, [3] charged particles found next to the
            # x+ wall, [4] kinetic energy of the planes loaded
            self._window_fe = cp.zeros(1)
            self._window_J = None  # the current of the move of a warm plane (allocated at the first one)
            self._wd = [cp.zeros(5) for _ in cfg.species]

        self.species: list[Species] = []
        self.sdiag = []
        self._load_species()
        self._init_injection()
        for c, b in zip((self.grid.Bx, self.grid.By, self.grid.Bz), cfg.B0):
            c[...] = b
        self.grid.apply_bc_B()
        self.B0 = tuple(float(b) for b in cfg.B0)
        self._init_external()
        fields_set = self._set_initial_fields()
        moved = self._apply_perturbations()
        # Drop the particles with weight zero (collective; the choice depends on the configuration
        # only, so every rank makes the same calls). A species with a profile can lose most of what
        # it loaded: it is compacted, which copies only the particles it keeps. Of the other species
        # only the particles inside the body go: they are flagged and taken out without a copy of
        # the species (a flag list that overflows takes all particles of the species out).
        shaped = [any(b.profile is not None for b in sc.beams) for sc in cfg.species]
        for sp, sh in zip(self.species, shaped):
            if sh:
                sp.compact_and_migrate(self.grid, self.comm)
        if self.body and not all(shaped):
            self.flag.fill(0)
            for i, sp in enumerate(self.species):
                if not shaped[i]:
                    lst = sp.flag_list()
                    launch(self.k_flag_dead, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap), self.flag[i:i + 1],
                                                    lst, np.int32(lst.size)))
            self._migrate()
        # nodes at which Gauss's law is checked: all but those that hold the charge of the body
        self._gauss_keep = None
        if self.body:
            g = self.grid
            keep = body.gauss_mask(cp, [g.global_positions(a, 0.0) * g.dx[a] for a in range(3)], cfg.body, g.dx, reach)
            self._gauss_keep = cp.ascontiguousarray(keep[g.gauss_slice(self.gauss_retreat)])
        # body.relax: the edges inside the body, as flat indices into the field arrays (ghosts too,
        # so that a ghost keeps the value of the node it copies), and the factor of one step
        self._relax = None
        if self.body and cfg.body.relax:
            g = self.grid
            edges = body.relax_edges(cp, [g.global_positions(a, 0.0) * g.dx[a] for a in range(3)], cfg.body, g.dx)
            self._relax = ([cp.flatnonzero(m) for m in edges], math.exp(-cfg.body.relax * cfg.dt))
        if self.fluct is not None:
            self._tfsf_start()
        if not self.implicit and (moved or fields_set or any(self.B0) or self.external):
            if self.fluct is not None:
                self._tfsf_gather(0.0, +1.0)
            self._initial_half_kick()
            if self.fluct is not None:
                self._tfsf_gather(0.0, -1.0)
        if self.implicit:
            self._init_implicit()
        self.step_count = 0
        self.time = 0.0
        self.history = {k: [] for k in (
            "step", "t", "kinetic", "electric", "magnetic", "total", "lost_kinetic", "lost_charge",
            "injected_kinetic", "injected_charge", "n_injected", "absorbed_kinetic", "absorbed_charge",
            "n_absorbed", "boundary_flux", "dissipated", "antenna", "window_out", "budget",
            "px", "py", "pz", "gauss", "divB", "picard", "n_particles", "fill")}
        self.kinetic_species = {s.name: [] for s in self.species}
        self.flux_out = 0.0
        self.dissipated = 0.0
        # laser antennas (PIC-D20) on this rank's block, and the work of their currents on the
        # field in this rank's physical cells: up to now, and up to the last recorded time level
        index = [self.grid.global_positions(a, 0.0) for a in range(3)]
        self.antennas = [antenna.Antenna(cp, laser, cfg, index) for laser in cfg.lasers]
        self._antenna_work = cp.zeros(1)
        self._antenna_rec = cp.zeros(1)
        self._last = None
        n_tot = self.total_particles()  # collective: must run on every rank, not only the printing one
        if self.verbose:
            print(f"[pypic] {self.comm.size} GPU(s), rank grid {dims}, global cells {tuple(cfg.n_cells)}, "
                  f"solver={cfg.solver}, particles={n_tot:,}", flush=True)
        from .output.manager import OutputManager

        self.output = OutputManager(self)
        if self._restart:
            checkpoint.restore(self, *self._restart)
            if self.verbose:
                print(f"[pypic] continued from {self._restart[0]} (step {self.step_count})", flush=True)

    # ----------------------------------------------------------------- set-up
    def _beam_block(self, b):
        """Global cells ``(lo, hi)`` per axis that a beam loads: the box; for a radial profile the
        cells the ball of its support reaches (PIC-D17); for a list of profiles the cells in which
        every shape reaches the cut (PIC-D21). With a moving window the cells of the first load
        (PIC-D22). The same on every rank."""
        blk = self._box_block(b)
        if self.window:  # along x: what the first load takes of the cells that the beam loads during the run
            front = window.front_cell(self.cfg)
            blk = (tuple(min(v, front) for v in window.beam_cells(self.cfg, b)),) + blk[1:]
        return blk

    def _box_block(self, b):
        n = self.cfg.n_cells
        L, dx = self.grid.L, self.grid.dx
        if isinstance(b.profile, list):
            blk = [(0, n[a]) for a in range(3)]
            for p in b.profile:
                a = "xyz".index(p.get("axis", "x"))
                lo, hi = support_interval(p, L[a])
                lo = min(max(int(math.floor(lo / dx[a])), 0), n[a])
                blk[a] = (lo, max(min(int(math.ceil(hi / dx[a])), n[a]), lo))
            return tuple(blk)
        if b.profile is None or b.profile.get("axis", "x") != "r":
            return tuple((0, n[a]) for a in range(3))
        rs = support(b.profile, math.sqrt(sum(max(abs(o), abs(L[a] - o)) ** 2 for a, o in enumerate(b.origin))))
        blk = []
        for a, o in enumerate(b.origin):
            lo = min(max(int(math.floor((o - rs) / dx[a])), 0), n[a])
            blk.append((lo, max(min(int(math.ceil((o + rs) / dx[a])), n[a]), lo)))
        return tuple(blk)

    def _local_block(self, blk):
        """The part of the global block ``blk`` in this rank's physical cells: (first cell, cells)."""
        g0, nl = [], []
        for (lo, hi), ax in zip(blk, self.decomp.axes):
            g0.append(max(lo, ax.phys_lo))
            nl.append(max(min(hi, ax.phys_lo + ax.n_phys) - g0[-1], 0))
        return g0, nl if all(nl) else [0, 0, 0]

    def _load_species(self):
        """Load every beam (plasma species first, then test species, so adding a test species never
        changes the random numbers of the plasma), explicit test-particle lists, density profiles
        (weights) and check the charge neutrality slot by slot (beams that load the same cells
        share positions)."""
        cfg = self.cfg
        n_beams = sum(len(sc.beams) for sc in cfg.species)
        n_glob = math.prod(cfg.n_cells)
        self._n_ids = n_beams * n_glob * cfg.ppc + sum(len(sc.particles) for sc in cfg.species)
        if self._n_ids >= ID_LIMIT:
            raise ValueError("too many particles for exact float64 particle ids (n_beams * cells * ppc >= 2^53)")
        cells_local = int(np.prod([a.n_phys for a in self.decomp.axes]))
        box = tuple((0, n) for n in cfg.n_cells)
        front = None
        if self.window:  # the first load ends at the front; the planes behind it get ids of their own
            front = window.front_cell(cfg)
            box = ((0, front),) + box[1:]
            self._window_cells = {(si, bi): window.beam_cells(cfg, b)
                                  for si, sc in enumerate(cfg.species) for bi, b in enumerate(sc.beams)}
            planes = self._window_end * n_beams * cfg.n_cells[1] * cfg.n_cells[2] * cfg.ppc
            if self._n_ids + planes >= ID_LIMIT:
                raise ValueError("too many particles for exact float64 particle ids with the planes that the "
                                 "moving window loads (>= 2^53)")
        order = [i for i, sc in enumerate(cfg.species) if not sc.test] + \
                [i for i, sc in enumerate(cfg.species) if sc.test]
        streams, stream = {}, 1
        for si in order:
            for bi in range(len(cfg.species[si].beams)):
                streams[(si, bi)] = stream
                stream += 1
        self._streams = streams
        self._beam_ranges = {}
        blocks = self._blocks = {}
        list_id = float(n_beams * n_glob * cfg.ppc)  # explicit-list ids follow all beam ids
        for si, sc in enumerate(cfg.species):
            # room for what the beams load; of a block smaller than the box a rank holds at least
            # its even share, since the particles spread from the ranks that load them. A bunch (a
            # list of profiles) can pass through a rank as a whole: room for all of it there. A beam
            # that a moving window goes on loading fills the rank in the end
            cells = 0 if sc.beams else cells_local
            for bi, b in enumerate(sc.beams):
                blk = blocks[(si, bi)] = self._beam_block(b)
                mine = math.prod(self._local_block(blk)[1])
                whole = math.prod(hi - lo for lo, hi in blk)
                if blk == box or (front is not None and self._window_cells[(si, bi)][1] > front):
                    cells += cells_local
                elif isinstance(b.profile, list):
                    cells += max(mine, min(whole, cells_local))
                else:
                    cells += max(mine, -(-whole // self.comm.size))
            cap = int(cfg.capacity_factor * cells * cfg.ppc) + 4096 + len(sc.particles)
            if self.compact:
                sp = CompactSpecies(sc.name, sc.charge, sc.mass, cap, self.grid)
            else:
                sp = Species(sc.name, sc.charge, sc.mass, cap)
            sp.test = bool(sc.test)
            sp.qdep = 0.0 if sp.test else sp.q
            for bi, b in enumerate(sc.beams):
                lo = sp.n
                sp.load(self.grid, self.k_load, b, cfg.ppc, cfg.seed, streams[(si, bi)],
                        self._local_block(blocks[(si, bi)]))
                self._beam_ranges[(si, bi)] = (lo, sp.n)
                if b.profile and bool(self._shape_weights(sp, b, lo)):
                    raise ValueError(f"species {sc.name}: density profile is negative somewhere")
            if sc.particles:
                self._load_list(sp, sc.particles, list_id)
            list_id += len(sc.particles)
            if self.body:  # no particle starts inside the body: weight zero, removed after the set-up
                launch(self.k_mark_body, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap)))
            self.species.append(sp)
            # [0] kinetic energy, [1..3] momentum (per diagnostic step); cumulative: [4] kinetic energy and
            # [5] charge lost at open faces, [6] kinetic energy, [7] charge and [8] number injected,
            # [9] kinetic energy, [10] charge and [11] number absorbed by the body
            self.sdiag.append(cp.zeros(NDIAG))
        # neutrality: the beams that load the same cells were loaded slot by slot on the same
        # positions, so summing q w over them position by position gives their net charge exactly.
        # One reduction for all blocks: [largest |net|, largest weight] per block.
        groups = {}
        for (si, bi), blk in blocks.items():
            if not cfg.species[si].test:
                groups.setdefault(blk, []).append((si, bi))
        groups = sorted(groups.items())
        loc = np.zeros(2 * len(groups))
        for g, (_, beams) in enumerate(groups):
            net = None
            for si, bi in beams:
                lo, hi = self._beam_ranges[(si, bi)]
                w = self.species[si].weight(lo, hi)
                net = self.species[si].q * w if net is None else net + self.species[si].q * w
                if hi > lo:
                    loc[2 * g + 1] = max(loc[2 * g + 1], float(cp.max(w)))
            if net.size:
                loc[2 * g] = float(cp.max(cp.abs(net)))
        if groups:
            glob = cp.asnumpy(self.comm.allreduce(cp.asarray(loc), "max"))
            for g, (blk, _) in enumerate(groups):
                if not glob[2 * g] > 1e-12 * max(glob[2 * g + 1], 1e-300):
                    continue
                if blk != box:  # a charge left in a part of the box has a field that E = 0 does not hold
                    raise ValueError(f"the beams loaded into the cells {blk} are not charge neutral together "
                                     "(PIC-D17)")
                if self.verbose:
                    print("[pypic] WARNING: the plasma is not charge neutral (net charge density != 0)", flush=True)

    @staticmethod
    def _shape_weights(sp, b, lo):
        """Multiply the weights of the particles [lo, sp.n) of ``sp`` by the density profiles of the
        beam ``b``. Returns a device flag: a profile is negative somewhere."""
        bad = cp.zeros((), dtype=bool)
        for p in b.profile if isinstance(b.profile, list) else [b.profile] if b.profile else []:
            axis = p.get("axis", "x")
            if axis == "r":
                x = cp.sqrt(sum((sp.pos(a, lo) - o) ** 2 for a, o in enumerate(b.origin)))
            else:
                x = sp.pos("xyz".index(axis), lo)
            f = evaluate(p, x, xp=cp)
            bad = bad | cp.any(f < 0)
            # negligible tails (sech^2, gauss) are cut: those particles get weight 0 and are
            # removed; a fixed threshold keeps this decomposition independent and, for equal
            # profiles, neutral. In a product each shape is cut by itself
            sp.scale_weight(cp.where(f < CUT, 0.0, f), lo)
        return bad

    def _load_list(self, sp, particles, id0: float):
        """Test particles given one by one: each rank keeps those inside its physical sub-domain,
        decided from the position as the storage holds it."""
        arr = np.asarray(particles, dtype=float)
        L = self.grid.L
        if np.any(arr[:, :3] < 0) or np.any(arr[:, :3] >= np.asarray(L)):
            raise ValueError(f"species {sp.name}: test particles must lie inside [0, L)")
        mine = np.ones(len(arr), dtype=bool)
        for a in range(3):
            lo, hi, _ = sp.frame(self.grid, a)
            x = sp.snap(arr[:, a], a)
            mine &= (x >= lo) & (x < hi)
        idx = np.nonzero(mine)[0]
        m = len(idx)
        if not m:
            return
        blk = np.zeros((8, m))
        blk[:6] = arr[idx, :6].T
        blk[6] = arr[idx, 6] if arr.shape[1] == 7 else 1.0
        blk[7] = id0 + idx
        sp.append(blk)

    def _init_injection(self):
        """Injection at open faces (PIC-D13): what this rank launches every step.

        ``inj_launch``: (species, face, first cell, cells, paired) per block of ``face_slabs`` and
        species that fills it; ``inj_room``: the most particles a step can append per species;
        ``inj_fmask``: the face bits of the species for ``push_deposit``; ``inj_bt``: the beam
        tables of ``inject_face``; ``inj_rb``: the inflow tables of the neutraliser per face;
        ``bank`` / ``binfo``: what the columns of the paired faces of this rank owe;
        ``inj_pstate`` / ``inj_pface``: the partners waiting for their momentum (``inject_partners``)."""
        cfg, S = self.cfg, len(self.species)
        self.inj_launch, self.inj_room, self.inj_fmask = [], [0] * S, [0] * S
        self.inj_rb, self.inj_nrb = cp.zeros(6), 1
        # inflow fluctuation (PIC-D26): per species the mode rows, the face (-1: none)
        self.fl_face = -1 if self.fluct is None else self.fluct["face"]
        self.fl_rows = [cp.zeros(11) for _ in range(S)]
        self.fl_n = [0] * S
        if self.fluct is not None:
            for si, rows in self.fluct["rows"].items():
                self.fl_rows[si] = cp.asarray(rows, dtype=cp.float64)
                self.fl_n[si] = len(rows)
        self.bank, self.binfo = cp.zeros(1), cp.full(12, -1, dtype=cp.int32)
        self.bank0 = self.bank.copy()
        self.inj_neutral = None
        if not self.injection:
            return
        index = {sc.name: si for si, sc in enumerate(cfg.species)}
        names = [index[n] for n in injection.injected_species(cfg)]
        modes = {FACES.index(k): v for k, v in cfg.injection.faces.items()}
        paired = [f for f, v in modes.items() if v == "paired"]
        ne = self.inj_neutral = index[cfg.injection.neutralizer] if paired else None
        carriers = [si for si in names if si != ne]
        # ids (PIC-D8) of the injected particles follow those of the loaded ones: one block per
        # step, in it one sub-block of the extended grid (cells -1..N) per stream, then one per
        # stream for the partners born with its particles at a paired face
        ext = math.prod(n + 2 for n in cfg.n_cells) * cfg.ppc
        n_streams = len(self._streams)
        self._inj_block = injection.ids_per_step(cfg)  # = 2 * n_streams * ext
        self.inj_bt = {}
        for si in names:
            rows = []
            for bi, b in enumerate(cfg.species[si].beams):
                st = self._streams[(si, bi)]
                rows.append([b.density * self.grid.cell_volume / cfg.ppc, b.theta, *b.drift, st,
                             (st - 1) * ext, (n_streams + st - 1) * ext])
            self.inj_bt[si] = cp.asarray(rows, dtype=cp.float64)
        axes = self.decomp.axes
        if paired:
            beams = cfg.species[ne].beams
            rb = np.zeros((6, len(beams), 6))
            binfo = np.full(12, -1, dtype=np.int32)
            binfo[6:9] = [ax.phys_lo for ax in axes]
            binfo[9:12] = [ax.n_phys for ax in axes]
            n_bank = 0
            for f in paired:
                extra = self.fluct["du_n"].get(ne, 0.0) if self.fluct is not None and f == self.fl_face else 0.0
                rb[f] = injection.inflow_table(beams, f, extra)
                if (axes[f >> 1].hi if f & 1 else axes[f >> 1].lo) is None:
                    binfo[f] = n_bank
                    n_bank += math.prod(ax.n_phys for a, ax in enumerate(axes) if a != f >> 1)
                self.inj_fmask[ne] |= 1 << f
                for si in carriers:
                    self.inj_fmask[si] |= 1 << (8 + f)
            self.inj_rb, self.inj_nrb = cp.asarray(rb), len(beams)
            self.bank, self.binfo = cp.zeros(max(n_bank, 1)), cp.asarray(binfo)
            self.bank0 = self.bank.copy()  # the bank before a launch (inject_face)
        for face, g0, nl in injection.face_slabs(axes, cfg.n_cells, cfg.injection.faces):
            pair = modes[face] == "paired"
            for si in (carriers if pair else names):
                room = math.prod(nl) * cfg.ppc * len(cfg.species[si].beams)
                self.inj_room[si] += room
                if pair:
                    self.inj_room[ne] += room
                self.inj_launch.append((si, face, g0, nl, int(pair)))
        room = max(self.inj_room[ne], 1) if paired else 1
        self.inj_pstate, self.inj_pface = cp.zeros(room, dtype=cp.uint64), cp.zeros(room, dtype=cp.int32)

    def _inject(self, do_diag: bool):
        """Fill the cell layers outside the injecting faces (after the push of every species: the
        particles are appended behind those being pushed; migrate() adds them to the counts)."""
        cfg, g, S = self.cfg, self.grid, len(self.species)
        id_step = float(self._n_ids + self.step_count * self._inj_block)
        for si, face, g0, nl, pair in self.inj_launch:
            sp = self.species[si]
            pi = self.inj_neutral if pair else si  # the species of the partner (unused on a flux face)
            pt = self.species[pi]
            lst, lstp = sp.flag_list(), pt.flag_list()
            if pair:
                cp.copyto(self.bank0, self.bank)
            launch(self.k_inject, math.prod(nl) * cfg.ppc * len(self.inj_bt[si]), (
                *sp.store, np.int64(sp.cap), np.int32(sp.n), self.flag[S + si:S + si + 1], self.sdiag[si],
                self.flag[si:si + 1], lst, np.int32(lst.size), g.Jx, g.Jy, g.Jz,
                np.float64(sp.q), np.float64(sp.m), *(np.int32(v) for v in g0), *(np.int32(v) for v in nl),
                *(np.int32(v) for v in cfg.n_cells), np.int32(cfg.ppc), np.uint64(cfg.seed),
                np.int64(self.step_count), np.float64(id_step), self.inj_bt[si], np.int32(len(self.inj_bt[si])),
                np.int32(face), np.int32(pair), *pt.store, np.int64(pt.cap), np.int32(pt.n),
                self.flag[S + pi:S + pi + 1], self.sdiag[pi], self.flag[pi:pi + 1], lstp, np.int32(lstp.size),
                np.float64(pt.q), self.inj_pstate, self.inj_pface, self.bank, self.bank0, self.binfo,
                np.int32(do_diag), self.fl_rows[si], np.int32(self.fl_n[si] if face == self.fl_face else 0),
                np.float64(self.time + 0.5 * cfg.dt)))
        ne = self.inj_neutral
        if ne is not None:
            pt = self.species[ne]
            launch(self.k_partners, self.inj_room[ne], (
                *pt.store, np.int64(pt.cap), np.int32(pt.n), self.flag[S + ne:S + ne + 1], self.sdiag[ne],
                np.float64(pt.q), np.float64(pt.m), self.inj_rb, np.int32(self.inj_nrb), self.inj_pstate,
                self.inj_pface, np.int32(do_diag), self.fl_rows[ne], np.int32(self.fl_n[ne]),
                np.int32(self.fl_face), np.float64(self.time + 0.5 * cfg.dt)))

    def _init_external(self):
        """External static fields (PIC-D12): per component Ex..Bz the uniform value and, for B with
        dipoles, an array on the local grid (ghosts included) computed from global indices."""
        cfg, g = self.cfg, self.grid
        self.external = check_config(cfg)
        self._gather = None
        if not self.external:
            return
        ef = cfg.external_fields
        self.ext_uniform = (*ef.E, *ef.B)
        index = [cp.asarray(g.global_positions(a, 0.0), dtype=cp.float64) for a in range(3)]
        B = dipoles_B(cp, index, g.dx, ef.dipoles)
        self.ext_array = (None, None, None) + (tuple(cp.ascontiguousarray(b) for b in B) if B else (None,) * 3)

    def total_field(self, name: str):
        """Component Ex..Bz as the particles feel it, in FP64: the Yee array itself, or a new array
        if the component has an external part (PIC-D12)."""
        f = getattr(self.grid, name)
        if self.external:
            c = FIELD_NAMES.index(name)
            u, ext = self.ext_uniform[c], self.ext_array[c]
            if u != 0.0 or ext is not None:
                f = f + u  # the order of the sums is that of the kernel ext_sum
                if ext is not None:
                    f += ext
        return f

    def _gather_fields(self):
        """Arrays the explicit push gathers from: the Yee arrays, float32 copies of them (PIC-D11)
        and / or their sums with the external fields (PIC-D12), refreshed every step."""
        g = self.grid
        fields = [g.Ex, g.Ey, g.Ez, g.Bx, g.By, g.Bz]
        if self.external:
            has = [u != 0.0 or ext is not None for u, ext in zip(self.ext_uniform, self.ext_array)]
            if self._gather is None:
                dtype = cp.float32 if self.mixed else cp.float64
                self._gather = [cp.empty(g.shape, dtype=dtype) if (h or self.mixed) else None for h in has]
            for c, dst in enumerate(self._gather):
                if dst is not None:
                    src, ext = fields[c], self.ext_array[c]
                    launch(self.k_extsum, src.size, (src, ext if ext is not None else src, np.int32(ext is not None),
                                                     np.float64(self.ext_uniform[c]), dst, np.int64(src.size)))
                    fields[c] = dst
        elif self.mixed:
            if self._gather is None:
                self._gather = [cp.empty(g.shape, dtype=cp.float32) for _ in fields]
            for dst, src in zip(self._gather, fields):
                cp.copyto(dst, src, casting="same_kind")
            fields = self._gather
        return fields

    def _set_initial_fields(self) -> bool:
        """Uniform E0 and the initial_fields list (PIC-D9). Returns True if any field was set.
        E0 may be tangential to reflecting faces (PIC-D18, explicit solver): the wall condition
        applied below zeroes it on the wall, which is a conductor at rest in a plasma that
        carries E = -v x B."""
        cfg, g = self.cfg, self.grid
        E0 = tuple(float(e) for e in cfg.E0)
        for c, e in enumerate(E0):
            if e != 0.0:
                for a in range(3):
                    kinds = set(cfg.boundaries[a]) - {"periodic"}
                    # PIC-D26: a reflecting face opposite the open face that carries the inflow fluctuation,
                    # whose incident field holds E0 there
                    carried = self.fluct is not None and self.fluct["face"] >> 1 == a and \
                        cfg.boundaries[a][self.fluct["face"] & 1] == "open" and kinds <= {"reflecting", "open"}
                    if a != c and kinds and ((kinds != {"reflecting"} and not carried) or cfg.solver != "explicit"):
                        raise ValueError(f"E0 component {'xyz'[c]} is tangential to the faces of axis {'xyz'[a]}: "
                                         "allowed at reflecting faces with the explicit solver only (PIC-D18)")
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
        too_close = 0.0  # 1 if a shaped displacement moves a particle near a non-periodic face
        for si, (sc, sp) in enumerate(zip(cfg.species, self.species)):
            disp = None
            for bi, b in enumerate(sc.beams):
                lo, hi = self._beam_ranges[(si, bi)]
                for pt in b.perturbation:
                    if "shape" in pt:  # PIC-D24: amplitude x product of 1-D shapes at the loaded place
                        f = None
                        for ax_name, spec in pt["shape"].items():
                            v = evaluate(spec, sp.pos("xyz".index(ax_name), lo, hi), cp)
                            f = v if f is None else f * v
                        if "velocity" in pt:
                            sp.add_mom("xyz".index(pt["velocity_component"]), pt["velocity"] * f, lo, hi)
                        if "displacement" in pt:
                            if disp is None:
                                disp = cp.zeros((3, sp.cap))
                            disp["xyz".index(pt["displacement_component"]), lo:hi] += pt["displacement"] * f
                        continue
                    a = "xyz".index(pt.get("axis", "x"))
                    if cfg.boundaries[a][0] != "periodic":
                        raise ValueError("perturbations need a periodic axis")
                    k = 2.0 * math.pi * int(pt.get("mode", 1)) / g.L[a]
                    ph = float(pt.get("phase", 0.0))
                    x0 = sp.pos(a, lo, hi).copy()
                    if "velocity" in pt:
                        c = "xyz".index(pt.get("velocity_component", pt.get("axis", "x")))
                        sp.add_mom(c, float(pt["velocity"]) * cp.sin(k * x0 + ph), lo, hi)
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
                too_close = max(too_close, self._displacement_near_face(sp, disp))
                sp.set_displacement(disp)
                moving = True
        if not any(b.perturbation and any(("density" in p or "displacement" in p) for p in b.perturbation)
                   for sc in cfg.species for b in sc.beams):
            return False  # decided from the configuration: identical on every rank
        # every rank reaches this reduction (the branch above depends on the configuration only)
        residual = self.comm.allreduce_scalar(residual, "max")
        too_close = self.comm.allreduce_scalar(too_close, "max")
        if too_close > 0.0:
            raise ValueError("shaped displacement (PIC-D24): particles closer than one cell plus the largest "
                             "displacement to a non-periodic face of the displaced axis must not move; give the "
                             "shape exact zeros there (e.g. ramp or step)")
        if residual > 1e-12:
            raise RuntimeError(f"density perturbation: the inverse of the cumulative density did not converge "
                               f"(max residual {residual:.3g} L)")
        h = [0.9 * d for d in g.dx]
        local = 0.0
        for sp in self.species:
            if sp.displacing and sp.n:
                local = max(local, max(sp.displacement_max(a) / h[a] for a in range(3)))
        rounds = int(math.ceil(self.comm.allreduce_scalar(local, "max")))
        g.zero_J()
        for _ in range(rounds):
            self.flag.fill(0)
            for i, sp in enumerate(self.species):
                if sp.displacing:
                    lst = sp.flag_list()
                    launch(self.k_displace, sp.n, (*sp.disp_store, np.int32(sp.n), np.int64(sp.cap), np.float64(h[0]),
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
            sp.end_displacement()
        return moving or rounds > 0

    def _displacement_near_face(self, sp, disp) -> float:
        """1.0 if a particle of this rank with a displacement along a non-periodic axis c lies closer
        than one cell plus the largest such displacement (over all ranks) to a face of axis c."""
        cfg, g = self.cfg, self.grid
        flag = 0.0
        for c in range(3):
            if cfg.boundaries[c][0] == "periodic" and cfg.boundaries[c][1] == "periodic":
                continue
            d = cp.abs(disp[c, : sp.n])
            dmax = self.comm.allreduce_scalar(float(cp.max(d)) if sp.n else 0.0, "max")
            if dmax == 0.0 or not sp.n:
                continue
            x = sp.pos(c)
            margin = g.dx[c] + dmax
            near = (x < margin) | (x > g.L[c] - margin)
            if bool(cp.any(near & (d > 0.0))):
                flag = 1.0
        return flag

    def _initial_half_kick(self):
        """Explicit leapfrog: loaded momenta are at t = 0; with non-zero initial E or B push them
        back to t = -dt/2 by a Boris step of -dt/2 (positions unchanged)."""
        fields = [self.total_field(n) for n in FIELD_NAMES]
        for sp in self.species:
            launch(self.k_kick, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap), *fields,
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
            self._fe = g.field_energy_local()  # fields at t^n; reduced over the ranks in _record
        if self.fluct is not None:  # the push gathers the total field across the TF/SF face (PIC-D26)
            self._tfsf_gather(self.time, +1.0)
        fields = self._gather_fields()
        if self.injection:
            # decided from the configuration and the step: every rank raises together
            if self._n_ids + (self.step_count + 1) * self._inj_block >= ID_LIMIT:
                raise RuntimeError("the ids of the injected particles reach 2^53 at this step (PIC-D13)")
            for sp, room in zip(self.species, self.inj_room):
                sp.ensure_capacity(sp.n + room)  # before the kernels take the arrays
        seed, step = np.uint64(self.cfg.seed), np.int64(self.step_count)
        for i, (sp, d) in enumerate(zip(self.species, self.sdiag)):
            lst = sp.flag_list()
            launch(self.k_push, sp.n, (
                *sp.store, np.int32(sp.n), np.int64(sp.cap), *fields,
                g.Jx, g.Jy, g.Jz, np.float64(sp.q * self.cfg.dt / (2.0 * sp.m)), np.float64(sp.m),
                np.float64(sp.qdep), d, self.flag[i:i + 1],
                lst, np.int32(lst.size), np.int32(do_diag), seed, step, np.int32(self.inj_fmask[i]),
                self.inj_rb, np.int32(self.inj_nrb), self.bank, self.binfo, self.fl_rows[i], np.int32(self.fl_n[i]),
                np.int32(self.fl_face), np.float64(self.time + 0.5 * self.cfg.dt)))
        if self.fluct is not None:
            self._tfsf_gather(self.time, -1.0)
        if self.inj_launch:
            self._inject(do_diag)
        g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
        self._filter_J()
        g.apply_bc_J()
        if self.cfg.remove_mean_current:
            self._remove_mean_current()
        if self.antennas:
            driven = self._antenna_add(do_diag)
        # The field update reads J and the fields, the migration the particles: queueing the update
        # first lets the host prepare it while the push kernels run, as the migration starts with
        # a host synchronisation.
        dt = self.cfg.dt
        t0 = self.time
        g.push_B(0.5 * dt)
        if self.fluct is not None:
            self._tfsf_B(0.5 * dt, t0)
        g.push_E(dt)
        if self.fluct is not None:
            self._tfsf_E(dt, t0 + 0.5 * dt)
        if self.antennas:
            self._antenna_work_done(driven)
        if self._relax:
            for e, idx in zip((g.Ex, g.Ey, g.Ez), self._relax[0]):
                if idx.size:
                    flat = e.reshape(-1)
                    flat[idx] = flat[idx] * self._relax[1]
        g.push_B(0.5 * dt)
        if self.fluct is not None:
            self._tfsf_B(0.5 * dt, t0 + dt)
        self._migrate()
        if self.marder:
            if self.fluct is not None:  # the zone acts on the total field: add the incident one beyond the face
                self._tfsf_total(t0 + dt, +1.0)
            g.marder_zone(self._charge_density(zone=True), self.zone_depth, self.marder)
            if self.fluct is not None:
                self._tfsf_total(t0 + dt, -1.0)
        if self._has_open:
            if self.fluct is not None:  # the flux of the total field through the fluctuation face (PIC-D26)
                self._tfsf_half_B(t0 + dt, +1.0)
            self.flux_out += dt * g.poynting_out()
            if self.fluct is not None:
                self._tfsf_half_B(t0 + dt, -1.0)

    # ------------------------------------------- inflow fluctuation: total-field / scattered-field face (PIC-D26)
    def _incident(self, comp: str, coords, t: float):
        """Incident field component ``comp`` (Ex..Bz) at the physical points ``coords`` (three
        broadcastable CuPy arrays) and time t: the uniform E0 / B0 and the modes of the fluctuation."""
        c, kind = "xyz".index(comp[1]), comp[0]
        bg = float((self.cfg.E0 if kind == "E" else self.cfg.B0)[c])
        out = cp.full(cp.broadcast(*coords).shape, bg)
        for m in self.fluct["modes"]:
            C, S = m[kind]
            if C[c] == 0.0 and S[c] == 0.0:
                continue
            k = m["k"]
            th = k[0] * coords[0] + k[1] * coords[1] + k[2] * coords[2] - m["omega"] * t + m["phase"]
            out += C[c] * cp.cos(th) + S[c] * cp.sin(th)
        return out

    def _plane(self, target: str, at: float):
        """Local slice of field ``target`` on the plane of global position ``at`` (cells) along the
        normal of the fluctuation face (ghost copies included: every rank that holds the plane
        corrects it the same way, without an exchange), and the physical coordinates of its points;
        None if this rank has no such plane."""
        g, a = self.grid, self.fl_face >> 1
        st = stagger_of(target)
        pos = g.global_positions(a, st[a])
        hit = np.nonzero(np.abs(pos - at) < 1e-9)[0]
        if hit.size == 0:
            return None
        sl = [slice(None)] * 3
        sl[a] = slice(int(hit[0]), int(hit[0]) + 1)
        coords = []
        for b in range(3):
            shape = [1, 1, 1]
            x = (np.full(1, at) if b == a else g.global_positions(b, st[b])) * g.dx[b]
            shape[b] = -1
            coords.append(cp.asarray(x.reshape(shape)))
        return tuple(sl), coords

    def _tfsf_axes(self):
        a = self.fl_face >> 1
        side = self.fl_face & 1
        b, c = (a + 1) % 3, (a + 2) % 3
        N = float(self.cfg.n_cells[a])
        node, half = (N, N + 0.5) if side else (0.0, -0.5)
        sgn = -1.0 if side else 1.0
        return a, b, c, node, half, sgn

    def _tfsf_B(self, dtb: float, tE: float):
        """B on the first half plane outside the face (scattered) was updated with the total E of the
        face plane: subtract the incident part. Hi face: B_c -= dtb E_b,inc / d, B_b += dtb E_c,inc / d;
        lo face: the opposite signs."""
        a, b, c, node, half, sgn = self._tfsf_axes()
        g = self.grid
        d = g.dx[a]
        for tgt, src, fac in (("B" + "xyz"[c], "E" + "xyz"[b], sgn), ("B" + "xyz"[b], "E" + "xyz"[c], -sgn)):
            pl = self._plane(tgt, half)
            if pl is None:
                continue
            sl, xyz = pl
            xyz[a] = cp.full((1, 1, 1), node * g.dx[a])
            getattr(g, tgt)[sl] += (fac * dtb / d) * self._incident(src, xyz, tE)

    def _tfsf_E(self, dte: float, tB: float):
        """E on the face plane (total) was updated with the scattered B of the half plane outside: add
        the incident part. Hi face: E_b -= dte B_c,inc / d, E_c += dte B_b,inc / d; lo face: the
        opposite signs."""
        a, b, c, node, half, sgn = self._tfsf_axes()
        g = self.grid
        d = g.dx[a]
        for tgt, src, fac in (("E" + "xyz"[b], "B" + "xyz"[c], sgn), ("E" + "xyz"[c], "B" + "xyz"[b], -sgn)):
            pl = self._plane(tgt, node)
            if pl is None:
                continue
            sl, xyz = pl
            xyz[a] = cp.full((1, 1, 1), half * g.dx[a])
            getattr(g, tgt)[sl] += (fac * dte / d) * self._incident(src, xyz, tB)

    def _tfsf_gather(self, t: float, sign: float):
        """Add (sign = +1) or take away (-1) the incident field at time t on the half plane outside the
        fluctuation face for the components staggered across it (E_a, B_b, B_c): a particle within
        half a cell of the face interpolates between the face's two sides, and must see the total
        field on both (review of PR #50)."""
        a, b, c, node, half, sgn = self._tfsf_axes()
        for comp in ("E" + "xyz"[a], "B" + "xyz"[b], "B" + "xyz"[c]):
            pl = self._plane(comp, half)
            if pl is None:
                continue
            sl, xyz = pl
            getattr(self.grid, comp)[sl] += sign * self._incident(comp, xyz, t)

    def _tfsf_half_B(self, t: float, sign: float):
        """Add (sign = +1) or take away (-1) the incident B_b, B_c at time t on the half plane outside
        the fluctuation face, so that the Poynting flux through the face (which averages B across it)
        is that of the total field: the energy that the incident wave brings in enters the record."""
        a, b, c, node, half, sgn = self._tfsf_axes()
        for comp in ("B" + "xyz"[b], "B" + "xyz"[c]):
            pl = self._plane(comp, half)
            if pl is None:
                continue
            sl, xyz = pl
            getattr(self.grid, comp)[sl] += sign * self._incident(comp, xyz, t)

    def _tfsf_total(self, t: float, sign: float):
        """Add (sign = +1) or take away (-1) the incident E at time t on the points beyond the
        fluctuation face, where the stored E is the scattered field: the Marder zone (PIC-D16) then
        sees the divergence of the total field. Its change of E is linear, so the round trip is exact."""
        g = self.grid
        a = self.fl_face >> 1
        side = self.fl_face & 1
        N = float(self.cfg.n_cells[a])
        for comp in ("Ex", "Ey", "Ez"):
            st = stagger_of(comp)
            pos = g.global_positions(a, st[a])
            beyond = np.nonzero(pos > N + 1e-9)[0] if side else np.nonzero(pos < -1e-9)[0]
            if beyond.size == 0:
                continue
            i0, i1 = int(beyond[0]), int(beyond[-1]) + 1
            sl = [slice(None)] * 3
            sl[a] = slice(i0, i1)
            coords = []
            for b in range(3):
                x = g.global_positions(b, st[b])
                if b == a:
                    x = x[i0:i1]
                shape = [1, 1, 1]
                shape[b] = -1
                coords.append(cp.asarray((x * g.dx[b]).reshape(shape)))
            getattr(g, comp)[tuple(sl)] += sign * self._incident(comp, coords, t)

    def _tfsf_start(self):
        """The fields outside the fluctuation face hold the scattered field (total - incident). The
        initial state is taken to be the incident one out there, so they start at zero."""
        g = self.grid
        a = self.fl_face >> 1
        side = self.fl_face & 1
        N = float(self.cfg.n_cells[a])
        for comp in FIELD_NAMES:
            pos = g.global_positions(a, stagger_of(comp)[a])
            beyond = np.nonzero(pos > N + 1e-9)[0] if side else np.nonzero(pos < -1e-9)[0]
            if beyond.size:
                sl = [slice(None)] * 3
                sl[a] = slice(int(beyond[0]), int(beyond[-1]) + 1)
                getattr(g, comp)[tuple(sl)] = 0.0

    def _antenna_add(self, do_diag: bool):
        """Add the antenna currents at t^(n+1/2) to J (PIC-D20): after the filter, the wall
        condition and the mean-current removal, which are for the current of the particles.
        Returns, per piece of current in this rank's physical cells, what its work needs: the
        component, the index along x, the current and E^n there."""
        g = self.grid
        if do_diag:
            self._antenna_rec = self._antenna_work.copy()  # the fields of the record are at t^n
        J, E = (g.Jx, g.Jy, g.Jz), (g.Ex, g.Ey, g.Ez)
        sx, sy, sz = g.physical_slice()
        t = (self.step_count + 0.5) * self.cfg.dt
        driven = []
        for ant in self.antennas:
            for c, i, a in ant.pieces(t):
                J[c][i] += a
                if sx.start <= i < sx.stop:
                    driven.append((c, i, a[sy, sz], E[c][i, sy, sz].copy()))
        return driven

    def _antenna_work_done(self, driven):
        """Work of the antenna currents on the field in this step, -dt sum J . (E^n + E^(n+1)) / 2
        dV: the Yee update changes the electric energy by exactly this through the term -J."""
        g = self.grid
        E = (g.Ex, g.Ey, g.Ez)
        _sx, sy, sz = g.physical_slice()
        f = 0.5 * self.cfg.dt * g.cell_volume
        for c, i, a, e0 in driven:
            self._antenna_work -= f * cp.sum(a * (e0 + E[c][i, sy, sz]))

    J_STAGGER = ((0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5))

    def _filter_J(self, J=None):
        """Binomial current filter (PIC-D10), between the ghost reduction and the wall BC."""
        if self.n_filter:
            g = self.grid
            g.filter(J or (g.Jx, g.Jy, g.Jz), self.J_STAGGER, self.n_filter, self.cfg.current_filter.compensate)

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
        migrate(self.species, self.grid, self.comm, self.flag, self.k_migrate, self.window_shift)

    # ---------------------------------------------------------- moving window
    _WALL_ERROR = ("a charged particle came next to the x+ wall of the moving window, where the field of the cell "
                   "that enters is set for a vacuum: raise window.margin (PIC-D22)")

    def _window_set(self):
        """Hand the shift, the faces of the window along x and this rank's bounds along x to the
        kernels: the same doubles as migrate() computes."""
        s, n = self.window_shift, self.cfg.n_cells[0]
        ax, dx = self.decomp.axes[0], self.cfg.dx[0]
        launch(self.k_wset, 1, (np.float64(s), np.float64(s * dx), np.float64((s + n) * dx),
                                np.float64((ax.phys_lo + s) * dx), np.float64((ax.phys_lo + ax.n_phys + s) * dx)))

    def _window_advance(self, s: int):
        """Shift the window to ``s`` cells, one more than before (PIC-D22), after a step: the
        fields go one plane down; the particles below the new x- face are dropped; the laboratory
        plane that enters behind the vacuum margin is loaded as it was at the shift before and
        moved to the present in rounds of less than a cell, each followed by a migration, with the
        current of the move taken out of E (div E = rho stays exact). That current has arrays of its
        own: J stays the current of the step, one plane down like the fields (the outputs read it).
        The energy that leaves and enters is kept for the record. Every rank makes the same
        collective calls."""
        cfg, g = self.cfg, self.grid
        if s > self._window_end:
            raise RuntimeError(f"the moving window was set up for n_steps = {cfg.n_steps} steps: no plane to load "
                               "beyond them (PIC-D22)")
        fe0 = g.field_energy_local()
        g.shift_x()
        self.window_shift = s
        self._window_set()
        for a, off in zip(self.antennas, self._antenna_off):
            a.set_shift(s, s < off)
        dx = cfg.dx[0]
        xnear = (s + cfg.n_cells[0] - self.zone_depth) * dx  # one cell plus the reach of the filter from the wall
        self.flag.fill(0)
        for i, sp in enumerate(self.species):
            lst = sp.flag_list()
            launch(self.k_wdrop, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.m),
                                        np.float64(sp.qdep), np.float64(xnear), self._wd[i], self.flag[i:i + 1], lst,
                                        np.int32(lst.size)))
        self._migrate()
        n0 = [sp.n for sp in self.species]
        hot = self._window_load(s)
        if hot is not None:
            # the plane was the front for the time T: thermal particles have left it downwards
            T = window.steps_since_shift(cfg, self.step_count + 1) * cfg.dt
            rounds = (math.floor(T / min(cfg.dx)) + 1) if any(hot) else 1
            xw = (window.front_cell(cfg) + s) * dx
            Jm = (g.Jx, g.Jy, g.Jz)  # not written by a plane that does not move
            if any(hot):
                if self._window_J is None:
                    self._window_J = tuple(cp.empty(g.shape) for _ in range(3))
                Jm = self._window_J
                for J in Jm:
                    J.fill(0)
            for r in range(rounds):
                self.flag.fill(0)
                for i, sp in enumerate(self.species):
                    if r and not hot[i]:
                        continue
                    k = sp.n - n0[i]
                    if k > sp.flag_list().size:  # the list holds every particle of the plane: no particle
                        sp.ensure_capacity(8 * k)  # that was there before is taken out with them
                    lst = sp.flag_list()
                    launch(self.k_wmove, k, (
                        *sp.store, np.int32(n0[i]), np.int32(sp.n), np.int64(sp.cap),
                        np.float64(T / rounds if hot[i] else 0.0), np.float64(xw), *Jm,
                        np.float64(sp.m), np.float64(sp.qdep), np.int32(r == 0), self._wd[i], self.sdiag[i],
                        self.flag[i:i + 1], lst, np.int32(lst.size)))
                self._migrate()
            if any(hot):
                g.reduce_ghosts(Jm)
                self._filter_J(Jm)
                g.apply_bc_J(Jm)
                o = g.owned_slice()
                for E, J in zip((g.Ex, g.Ey, g.Ez), Jm):
                    E[o] -= cfg.dt * J[o]
        g.apply_bc_E()
        g.fill_ghosts((g.Ex, g.Ey, g.Ez))
        self._window_fe += cp.sum(fe0 - g.field_energy_local())

    def _window_load(self, s: int):
        """Load the laboratory plane that the shift to ``s`` brings into the window, on the rank
        that holds its cells, with the loader of the first load: the particles depend on (seed,
        laboratory cell, slot) only, and their ids follow those of the first load plane by plane.
        The beams of the plane must be neutral together (one reduction; every rank raises).
        Returns per species whether its part of the plane moves (a temperature or a drift), or
        None if no beam reaches the plane."""
        cfg = self.cfg
        front = window.front_cell(cfg)
        p = front + s - 1
        ax = self.decomp.axes[0]
        mine = ax.phys_lo <= front - 1 < ax.phys_lo + ax.n_phys
        per_plane = cfg.n_cells[1] * cfg.n_cells[2] * cfg.ppc
        hot = [False] * len(self.species)
        groups, bad, found = {}, cp.zeros((), dtype=bool), False
        for si, sc in enumerate(cfg.species):
            sp = self.species[si]
            for bi, b in enumerate(sc.beams):
                lo, hi = self._window_cells[(si, bi)]
                if not lo <= p < hi:
                    continue
                found = True
                hot[si] = hot[si] or b.theta != 0.0 or any(b.drift)
                a0 = sp.n
                blk = self._blocks[(si, bi)]
                if mine:
                    g0, nl = self._local_block(((front - 1, front),) + blk[1:])
                    stream = self._streams[(si, bi)]
                    # the loader adds (global cell * ppc + slot), which holds p * per_plane
                    base = self._n_ids + ((p - front) * len(self._streams) + stream - 1 - p) * per_plane
                    sp.load(self.grid, self.k_load, b, cfg.ppc, cfg.seed, stream, ([p] + g0[1:], nl), id_base=base)
                    bad = bad | self._shape_weights(sp, b, a0)
                if not sc.test:
                    groups.setdefault(blk[1:], []).append((si, a0, sp.n))
        if not found:
            return None
        groups = sorted(groups.items())
        loc = cp.zeros(2 * len(groups) + 1)  # per group [largest |net|, largest weight]; a negative profile
        loc[-1] = bad
        for gi, (_, beams) in enumerate(groups):
            net = None
            for si, a0, a1 in beams:
                if a1 > a0:
                    w = self.species[si].weight(a0, a1)
                    net = self.species[si].q * w if net is None else net + self.species[si].q * w
                    loc[2 * gi + 1] = cp.maximum(loc[2 * gi + 1], cp.max(w))
            if net is not None:
                loc[2 * gi] = cp.max(cp.abs(net))
        glob = cp.asnumpy(self.comm.allreduce(loc, "max"))
        if glob[-1]:
            raise ValueError(f"a density profile is negative in the plane {p} that the moving window loads")
        for gi in range(len(groups)):
            if glob[2 * gi] > 1e-12 * max(glob[2 * gi + 1], 1e-300):
                raise ValueError(f"the beams loaded into the plane {p} of the moving window are not charge neutral "
                                 "together (PIC-D22)")
        return hot

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
        self.cg_part = cp.zeros((-(-self.E1.size // CG_BLOCK), 3))  # inner products per thread block
        self.pc_part = cp.zeros((self.cg_part.shape[0], 2))  # Picard change per thread block
        # no wall on any rank (and then no absorbing layer): the CG matrix product is one kernel
        self.cg_fused = all(b == "periodic" for ax in self.cfg.boundaries for b in ax)

    def _applyD(self, B):
        if self.absorb:
            for c in range(3):
                B[c] *= self.Dmat[c]

    def _matvec(self, v, out, tmp, r):
        """out = A v,  A = (1 + sigma dt/2) + dt^2/4 curl D curl  (v, out, tmp, r stacked (3, ...)),
        and the inner products the CG takes with it: returns ((r, v), (out, v), (r, r)) over the
        owned, unconstrained E unknowns of all ranks, on the device (one reduction).
        One ghost exchange: curl E is also taken on the ghost layer that curl B reads."""
        g, dt = self.grid, self.cfg.dt
        if self.cg_fused:  # no wall, no absorbing layer on any rank: the whole product in one launch
            g.fill_ghosts(v)
            launch(g.k_cg_matvec, out.size,
                   (v, out, self.Adiag, r, self.Wmask, self.cg_part, np.float64(dt * dt / 4.0)), block=CG_BLOCK)
            return self.comm.allreduce(self.cg_part.sum(axis=0))
        g.apply_bc_E(v)
        g.fill_ghosts(v)
        g.curl_E(v, tmp, low=True)
        self._applyD(tmp)
        g.apply_bc_B(tmp)
        g.curl_B(tmp, out, dt * dt / 4.0)
        # the sums of each thread block, then of the blocks: a full reduction to one scalar by a
        # ReductionKernel runs as a single thread block (issue #5)
        launch(g.k_cg_matvec_end, out.size, (out, self.Adiag, v, r, self.Wmask, self.cg_part), block=CG_BLOCK)
        g.apply_bc_E(out)
        return self.comm.allreduce(self.cg_part.sum(axis=0))

    def _cg(self, x, b):
        """Jacobi-preconditioned CG for A x = b with A symmetric positive definite, in the form of
        Chronopoulos and Gear: the matrix product of the search direction follows from a
        recurrence, so that all inner products of an iteration are taken with its one matrix
        product (u = Minv r, w = A u) and reduced over the ranks together.
        x, b: stacked (3, nx, ny, nz) arrays. Scalars stay on the GPU; convergence is checked
        every few iterations to avoid host synchronisation."""
        icfg = self.cfg.implicit
        g = self.grid
        tmp = cp.empty_like(x)
        w = cp.empty_like(x)
        bb = self._matvec(x, w, tmp, b)[2]
        r = b - w
        g.apply_bc_E(r)
        u = self.Minv * r
        p = cp.zeros_like(x)
        s = cp.zeros_like(x)
        st = cp.zeros(3)
        sc = self._matvec(u, w, tmp, r)
        tol2 = (icfg.cg_tol ** 2) * float(bb)
        it = 0
        for it in range(icfg.cg_maxit):
            if it % 4 == 0 and float(sc[2]) <= tol2:
                break
            # The convergence test runs every few iterations only, so the residual can reach exactly
            # zero in between (e.g. a curl-free mode where A = I is solved in one step): cg_scalars
            # takes a zero denominator as a zero step.
            launch(g.k_cg_scalars, 1, (sc, st))
            launch(g.k_cg_step, x.size, (x, r, u, p, s, w, self.Minv, st))
            sc = self._matvec(u, w, tmp, r)
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
            for sp in self.species:
                launch(self.k_ipush, sp.n, (
                    sp.P, sp.Pn, np.int32(sp.n), np.int64(sp.cap), *G, *self.Bh, g.Jx, g.Jy, g.Jz,
                    np.float64(sp.q * dt / (2.0 * sp.m)), np.float64(sp.q / sp.m), np.float64(sp.qdep),
                    np.int32(1 if k == 0 else 0), np.float64(icfg.particle_tol), np.int32(icfg.particle_maxit),
                    self.err))
            g.reduce_ghosts((g.Jx, g.Jy, g.Jz))
            self._filter_J()
            g.apply_bc_J()
            if cfg.remove_mean_current:
                self._remove_mean_current()
            for c, J in enumerate((g.Jx, g.Jy, g.Jz)):
                cp.subtract(R0[c], dt * J, out=rhs[c])
            Enew = self.E1.copy()
            self._cg(Enew, rhs)
            num, den = self._picard_change(Enew)
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

    def _picard_change(self, Enew):
        """Largest |Enew - E1| and largest |Enew| over the E unknowns of all ranks (a NaN is kept)."""
        launch(self.grid.k_picard_change, Enew.size, (Enew, self.E1, self.Wmask, self.pc_part), block=CG_BLOCK)
        return tuple(float(v) for v in self.comm.allreduce(self.pc_part.max(axis=0), "max"))

    def _half_fields(self, E0, B0):
        g, dt = self.grid, self.cfg.dt
        # Eh on every point: the ghosts of E0 (filled when it was committed) and of E1 (filled with the
        # initial guess and at the end of each field solve) are those of their neighbours already
        launch(g.k_half_E, self.Eh[0].size, (*E0, self.E1, self.Eh))
        g.apply_bc_E(self.Eh)
        launch(g.k_half_B, g.n_update, (self.Eh, *B0, self.Bh, np.float64(-0.5 * dt)))
        self._applyD(self.Bh)
        g.apply_bc_B(self.Bh)
        g.fill_ghosts(self.Bh)

    # ------------------------------------------------------------ diagnostics
    def _kinetic_now(self):
        """Kinetic energy and momentum at the current integer time level (implicit scheme), into
        the diagnostics buffer of each species."""
        for sp, d in zip(self.species, self.sdiag):
            d[:4] = 0
            launch(self.k_kin, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.m), d))

    def _charge_density(self, zone: bool = False):
        """Charge density at the nodes that div E equals: deposited with the shape of the current
        deposit, ghosts reduced, and filtered like the current (a new array). ``zone``: right at the
        nodes of the Marder zone only, from the particles that the filter brings there."""
        g = self.grid
        rho = cp.zeros(g.shape)
        shell = np.float64(2 * self.zone_depth if zone else 0.0)  # zone nodes + reach of the filter + a cell
        for sp in self.species:
            if not sp.test:
                launch(self.k_rho, sp.n, (*sp.store, np.int32(sp.n), np.int64(sp.cap), np.float64(sp.q), shell, rho))
        g.reduce_ghosts((rho,))
        if self.n_filter:  # div E equals the identically filtered charge density
            g.filter((rho,), ((0, 0, 0),), self.n_filter, self.cfg.current_filter.compensate)
        return rho

    def gauss_error(self):
        g = self.grid
        rho = self._charge_density()
        div = g.div_E()
        s = g.gauss_slice(self.gauss_retreat)
        d = cp.abs(div[s] - rho[s])
        if self._gauss_keep is not None:
            d *= self._gauss_keep
        err = float(cp.max(d)) if d.size else 0.0
        db = g.div_B()
        sb = g.divb_slice()
        errb = float(cp.max(cp.abs(db[sb]))) if db[sb].size else 0.0
        return self.comm.allreduce_scalar(err, "max"), self.comm.allreduce_scalar(errb, "max")

    def _record(self, fe_local):
        """Append the diagnostics of this step to the history, from the diagnostics buffers of the
        species (kinetic energy, momentum, lost energy) and this rank's YeeGrid.field_energy_local."""
        h = self.history
        cfg = self.cfg
        # one allreduce and one host sync for everything
        S = len(self.species)
        # the fill of this rank's particle arrays, in its own slot: the sum over the ranks keeps
        # every rank's value and the host takes the largest
        fill = np.zeros(self.comm.size)
        fill[self.rank] = max((sp.n / sp.cap for sp in self.species), default=0.0)
        vec = cp.concatenate([*self.sdiag, cp.asarray([float(sp.n) for sp in self.species]), cp.asarray(fill),
                              fe_local, self._antenna_rec, *((self._window_fe, *self._wd) if self.window else ())])
        if self.comm.size > 1:
            vec = self.comm.allreduce(vec)
        hv = cp.asnumpy(vec)
        sd = hv[: NDIAG * S].reshape(S, NDIAG)
        kin, lost_s = sd[:, :4], sd[:, 4]
        hv = hv[NDIAG * S:]
        n_s, fill, fe = hv[:S], hv[S: S + self.comm.size], hv[S + self.comm.size:]
        # test species are recorded per species but are not part of the budget
        plasma = [k for k, sp in zip(kin, self.species) if not sp.test]
        ke = sum(float(k[0]) for k in plasma)
        mom = [sum(float(k[i]) for k in plasma) for i in (1, 2, 3)]
        lost = sum(float(v) for v, sp in zip(lost_s, self.species) if not sp.test)
        absorbed = sum(float(v) for v, sp in zip(sd[:, 9], self.species) if not sp.test)
        injected = float(np.sum(sd[:, 6]))  # test species are never injected and deposit no charge
        for sp, k in zip(self.species, kin):
            self.kinetic_species[sp.name].append(float(k[0]))
        tot = ke + float(fe[0]) + float(fe[1])
        # moving window: (field + kinetic energy) before - after its shifts, up to t^n
        wd = fe[4:].reshape(S, 5) if self.window else np.zeros((S, 5))
        window_out = float(fe[3]) if self.window else 0.0
        window_out += sum(float(w[0] - w[4]) for w, sp in zip(wd, self.species) if not sp.test)
        h["step"].append(self.step_count)
        h["t"].append(self.time)
        h["kinetic"].append(ke)
        h["electric"].append(float(fe[0]))
        h["magnetic"].append(float(fe[1]))
        h["total"].append(tot)
        h["lost_kinetic"].append(lost)
        h["lost_charge"].append(float(np.sum(sd[:, 5])))
        h["injected_kinetic"].append(injected)
        h["injected_charge"].append(float(np.sum(sd[:, 7])))
        h["n_injected"].append(int(round(float(np.sum(sd[:, 8])))))
        h["absorbed_kinetic"].append(absorbed)
        h["absorbed_charge"].append(float(np.sum(sd[:, 10])))
        h["n_absorbed"].append(int(round(float(np.sum(sd[:, 11])))))
        h["boundary_flux"].append(self.flux_out)
        h["dissipated"].append(self.dissipated)
        h["antenna"].append(float(fe[2]))  # work of the laser antennas on the field up to t^n
        h["window_out"].append(window_out)
        h["budget"].append(tot + lost + absorbed - injected + self.flux_out + self.dissipated - float(fe[2])
                           + window_out)
        h["px"].append(mom[0])
        h["py"].append(mom[1])
        h["pz"].append(mom[2])
        h["picard"].append(getattr(self, "_picard_its", 0))
        h["n_particles"].append(int(round(float(np.sum(n_s)))))
        h["fill"].append(float(np.max(fill)))
        if cfg.gauss_interval and self.step_count % cfg.gauss_interval == 0:
            ge, gb = self.gauss_error()
        else:
            ge, gb = float("nan"), float("nan")
        h["gauss"].append(ge)
        h["divB"].append(gb)
        if float(np.sum(wd[:, 3])) > 0.0:  # the reduced count: every rank raises
            raise RuntimeError(self._WALL_ERROR)

    def _output_step(self):
        """The outputs due at this time level. With a fluctuation face (PIC-D26) the gathers of the
        outputs (tracks with fields) see the total field across the face, as the push does; the
        incident field sits only on the half plane outside the physical domain meanwhile (review of
        PR #50)."""
        if self.fluct is not None:
            self._tfsf_gather(self.time, +1.0)
        try:
            self.output.on_step()
        finally:
            if self.fluct is not None:
                self._tfsf_gather(self.time, -1.0)

    # ------------------------------------------------------------------ run
    def step(self):
        cfg = self.cfg
        # file output of the state at time level n (collective; rank 0 writes)
        self._output_step()
        do_diag = cfg.diag_interval > 0 and self.step_count % cfg.diag_interval == 0
        if cfg.sort_interval and self.step_count % cfg.sort_interval == 0:
            for sp in self.species:
                sp.sort(self.k_cell)
        if self.implicit:
            if do_diag:
                self._kinetic_now()
                # with an absorbing layer the field energy stored in the layer is part of the budget
                fe = self.grid.field_energy_local(region="all" if self.absorb else "physical")
                self._record(fe)
            self._step_implicit(do_diag)
        else:
            self._step_explicit(do_diag)
            if do_diag:
                # kinetic energy centred at t^n (from the push), fields at t^n
                self._record(self._fe)
            if self.window:
                s = window.shift_after(cfg, self.step_count + 1)
                if s > self.window_shift:
                    self._window_advance(s)
        self.step_count += 1
        self.time += cfg.dt
        every = cfg.checkpoint.every
        if every and self.step_count % every == 0:
            checkpoint.save(self)

    def run(self, n_steps: int | None = None, progress: int = 0):
        """Run ``n_steps`` steps; by default cfg.n_steps, less the steps before the checkpoint that a
        continued run starts from (PIC-D25)."""
        n = max(self.cfg.n_steps - (self._restart_step or 0), 0) if n_steps is None else n_steps
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
                    de = h["budget"][-1] - e0
                    # a run that starts without energy (cold plasma at rest, no field) has no E0
                    msg += f" dE/E0={de / e0:+.3e}" if e0 else f" dE={de:+.3e}"
                print(msg, flush=True)
        self._output_step()  # outputs due at the final time level
        self.output.flush()
        cp.cuda.Device().synchronize()
        self.wall_time = time.perf_counter() - t0
        self._run_steps = n
        # collective: every rank must reach this point (push_rate() itself is then local)
        self._n_total_end = self.total_particles()
        if self.window and self._wd:  # collective, as above (the same on every rank: no species, no check)
            near = self.comm.allreduce(cp.stack([w[3] for w in self._wd]))
            if float(cp.sum(near)) > 0.0:
                raise RuntimeError(self._WALL_ERROR)
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
