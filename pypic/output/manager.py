"""Output manager: runs the resolved schedule inside the simulation loop (OUT-D2).

All ranks call :meth:`on_step` at the start of every step (state at time level n) and once
after the run; due diagnostics are computed collectively and written by rank 0 only. Every file
operation of rank 0 shares its outcome with all ranks (:meth:`OutputManager._root_io`): if it
fails, every rank raises :class:`OutputError` together, instead of the other ranks waiting in
the next collective.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import secrets
import shutil

import numpy as np

from .. import __version__, window
from ..provenance import commit_label, repo_root
from .diagnostics import KINDS, FieldSource
from .names import expand_fields
from .schedule import resolve

SCALAR_KEYS = ("kinetic", "electric", "magnetic", "total", "lost_kinetic", "lost_charge", "injected_kinetic",
               "injected_charge", "n_injected", "absorbed_kinetic", "absorbed_charge", "n_absorbed",
               "boundary_flux", "dissipated", "antenna", "window_out", "budget",
               "px", "py", "pz", "gauss", "divB", "picard", "n_particles", "fill")


class OutputError(RuntimeError):
    """A file operation of the output failed (raised on every rank)."""


def _git_commit() -> str | None:
    """Commit of the running code, "-dirty" if the tree differs from it (pypic.provenance)."""
    try:
        root = repo_root(os.path.dirname(os.path.abspath(__file__)))
        return commit_label(root) if root is not None else None
    except Exception:  # noqa: BLE001 - provenance is best effort
        return None


def _gpu_name() -> str | None:
    try:
        import cupy as cp

        return cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)["name"].decode()
    except Exception:  # noqa: BLE001 - provenance is best effort
        return None


def _bc_names(cfg, particle: bool):
    out = []
    for lo, hi in cfg.boundaries:
        for b in (lo, hi):
            out.append("absorbing" if (particle and b == "open") else b)
    return out


class OutputManager:
    def __init__(self, sim):
        self.sim = sim
        cfg = sim.cfg
        self.cfg = cfg.output
        self.enabled = self.cfg.dir is not None
        self.root = sim.rank == 0
        self.diags = []
        self.writers = {}
        self.last = {}
        self.acc_step = {}
        self.meta = {}
        self.history_seen = 0  # history rows seen by on_step (every rank: it decides a collective)
        self.scalar_rows = 0  # history rows taken by the scalars writer (rank 0)
        if not self.enabled:
            if self.cfg.diagnostics:
                raise ValueError("output.diagnostics are configured but output.dir is not set")
            return
        self.resolved = resolve(cfg)
        for d, r in zip(self.cfg.diagnostics, self.resolved):
            self.diags.append(KINDS[d.kind](sim, d, r))
        self.src = FieldSource(sim)
        from .h5 import Units

        self.units = Units(self.cfg.n0_m3)
        # provenance of the code that runs, taken once: run.json is rewritten later (metadata, run
        # statistics), when the working tree may have changed without affecting this run
        self.git_commit = _git_commit() if self.root else None
        # identity of this run, fixed at start-up: separate runs of one configuration (realizations)
        # stay apart in the recorded evidence (scripts/record_example.py)
        self.created = _dt.datetime.now().astimezone().isoformat(timespec="seconds")
        self.run_id = secrets.token_hex(4)
        self._claim_directory()

    # entries a run writes into its directory (removed by output.overwrite; nothing else is touched)
    RUN_ENTRIES = ("run.json", "scalars.h5", "reduced", "fields", "particles", "tracks")

    def _claim_directory(self):
        """Refuse a directory that already holds run output: files the new run does not rewrite
        (other dump iterations, other diagnostics) would be read as part of it. Decided on rank 0
        and reduced, so every rank raises together. output.overwrite removes those entries first."""
        d = self.cfg.dir
        present = [e for e in self.RUN_ENTRIES if os.path.lexists(os.path.join(d, e))] if self.root else []
        refuse = float(bool(present) and not self.cfg.overwrite)
        if self.sim.comm.allreduce_scalar(refuse, "max") > 0:
            found = f" ({', '.join(present)})" if present else ""
            raise FileExistsError(f"{d} already holds run output{found}: choose another output.dir, or set "
                                  "output.overwrite (--overwrite) to replace it")

        def prepare():
            os.makedirs(d, exist_ok=True)
            for e in present:
                p = os.path.join(d, e)
                if os.path.isdir(p) and not os.path.islink(p):
                    shutil.rmtree(p)
                else:
                    os.remove(p)
            self.write_run_json()

        self._root_io(f"preparing the output directory {d}", prepare)

    def _root_io(self, what: str, fn, *args):
        """Run a file operation on rank 0 and share its outcome (collective: every rank calls it at
        the same point). A failure raises OutputError on every rank before any further collective;
        rank 0 chains the original exception."""
        err = None
        if self.root:
            try:
                fn(*args)
            except Exception as exc:  # noqa: BLE001 - re-raised below, on every rank
                err = exc
        if self.sim.comm.allreduce_scalar(float(err is not None), "max") > 0:
            if err is not None:
                raise OutputError(f"{what} failed: {err!r}") from err
            raise OutputError(f"{what} failed on rank 0")

    # ------------------------------------------------------------------ metadata
    def edpic(self) -> dict:
        cfg = self.sim.cfg
        filt = getattr(cfg, "current_filter", None)
        passes = getattr(filt, "passes", 0) if filt is not None else 0
        d = {"fieldSolver": "Yee" if not self.sim.implicit else "other",
             "fieldBoundary": _bc_names(cfg, False), "particleBoundary": _bc_names(cfg, True),
             "currentSmoothing": "Binomial" if passes else "none", "chargeCorrection": "none",
             "particleShape": 1.0, "currentDeposition": "Esirkepov" if not self.sim.implicit else "other",
             "particlePush": "Boris" if not self.sim.implicit else "other",
             "particleInterpolation": "trilinear", "particleSmoothing": "none"}
        if passes:
            d["currentSmoothingParameters"] = f"period=1;numPasses={passes};compensator={int(filt.compensate)}"
        return d

    def time_centering(self) -> dict:
        imp = self.sim.implicit
        return {"position": 0.0, "momentum": 0.0 if imp else -0.5, "E": 0.0, "B": 0.0, "J": -0.5, "rho": 0.0,
                "n": 0.0, "flux": 0.0 if imp else -0.5, "ekin": 0.0 if imp else -0.5,
                "unit": "dt (time level of a quantity relative to the output time)"}

    def diag_layout(self, dg) -> dict:
        d, r = dg.d, dg.r
        dt = self.sim.cfg.dt
        e = {"name": d.name, "kind": d.kind, "every_steps": r.every, "dt_out": r.every * dt, "first_step": r.first,
             "last_step": r.last, "n_outputs": r.n_outputs, "bytes_per_output": r.bytes_per_output,
             "average": r.average, "rule": r.rule, "notes": r.notes, "dtype": d.dtype, "keep_bits": d.keep_bits}
        if d.kind in ("fields", "particles"):
            e["path"] = f"{d.kind}/{d.name}/{d.name}_%08d.h5"
        elif d.kind == "tracks":
            e["path"] = f"tracks/{d.name}_%s.h5"  # %s = species
        else:
            e["path"] = f"reduced/{d.name}.h5"
        if d.kind in ("fields", "slice", "lineout", "modes", "kspectrum"):
            e["components"] = expand_fields(d.fields)
        if d.kind == "slice":
            e["axis"] = d.axis
            e["plane_position"] = dg.plane_coords()
        if d.kind == "lineout":
            e["axis"] = d.axis
        if d.kind == "modes":
            e["modes"] = [list(m) for m in d.modes]
        if d.kind == "histogram":
            e["quantities"] = list(d.quantities)
            e["species"] = [s.name for s in dg.species_list()]
        if d.kind in ("particles", "tracks"):
            e["species"] = [s.name for s in dg.species_list()]
        return e

    def write_run_json(self):
        sim, cfg = self.sim, self.sim.cfg
        info = {
            "format": "pypic-output-1", "pypic_version": __version__, "git_commit": self.git_commit,
            "run_id": self.run_id, "created": self.created,
            "world_size": sim.comm.size, "gpu_name": _gpu_name(), "dims": list(sim.decomp.dims), "solver": cfg.solver,
            "dt": cfg.dt, "n_steps": cfg.n_steps, "n_cells": list(cfg.n_cells), "dx": list(cfg.dx),
            "L": list(sim.grid.L), "species": [{"name": s.name, "charge": s.q, "mass": s.m,
                                                "test": bool(getattr(s, "test", False))} for s in sim.species],
            "units": self.units.as_dict(), "time_centering": self.time_centering(), "edpic": self.edpic(),
            # E and B of every diagnostic are the fields the particles feel; the energies of the scalar
            # history are those of the Yee fields, which exclude the external fields (PIC-D12)
            "field_convention": {"E_B": "total", "field_energy": "yee",
                                 "external_fields": bool(sim.external)},
            # moving window (OUT-D5): the shift of an output, in cells, is the row ``window_shift`` of
            # scalars.h5 and of every reduced series, and is in gridGlobalOffset of a field dump
            "window": ({"start": cfg.window.start, "margin_cells": window.margin_cells(cfg),
                        "front_cell": window.front_cell(cfg), "window_shift_end": window.shift_after(cfg, cfg.n_steps)}
                       if sim.window else None),
            "timescales": dict(self.cfg.timescales), "codec": self.cfg.codec,
            "scalars": {"path": "scalars.h5", "every_steps": cfg.diag_interval} if self.cfg.scalars else None,
            "diagnostics": [self.diag_layout(dg) for dg in self.diags],
            "config": cfg.to_dict(), "meta": self.meta,
        }
        path = os.path.join(self.cfg.dir, "run.json")
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(info, f, indent=1, default=_json_default)
        os.replace(tmp, path)

    def update_meta(self, meta: dict):
        """Add experiment metadata (planner report, theory predictions, run statistics) to run.json.
        Collective: call it on every rank (rank 0 writes)."""
        self.meta.update(meta)
        if self.enabled:
            self._root_io("writing run.json", self.write_run_json)

    # ------------------------------------------------------------------ loop
    def on_step(self):
        if not self.enabled:
            return
        sim = self.sim
        step = sim.step_count
        self.src.clear()
        pending = []
        for dg in self.diags:
            r, name = dg.r, dg.d.name
            # time-averaged diagnostics sample every step of (S - every, S] for an output at S
            if r.average and r.first - r.every < step <= r.last and self.acc_step.get(name, -1) < step:
                dg.accumulate(self.src)
                self.acc_step[name] = step
            if r.due(step) and self.last.get(name, -1) < step:  # once per step, even if called twice
                pending.append((dg, dg.compute(self.src)))
                self.last[name] = step
        # whether to write is decided from the schedule and the history length, the same on every
        # rank; which history rows still need writing is rank 0's (rows are marked written only
        # once the scalars writer has taken them, so rows of a failed attempt are retried)
        n = len(sim.history["step"]) if self.cfg.scalars else 0
        new_rows, self.history_seen = n > self.history_seen, max(n, self.history_seen)
        if pending or new_rows:

            def write():
                for dg, data in pending:
                    self._write(dg, step, data)
                self._write_scalars()

            self._root_io(f"writing the output of step {step}", write)
        self.src.clear()

    def flush(self):
        """Write everything recorded so far (end of a run): pending scalar rows and partial chunks.
        Collective."""
        if not self.enabled:
            return

        def write():
            self._write_scalars()
            for w in self.writers.values():
                w.flush()

        self._root_io("flushing the output", write)

    def close(self):
        """Collective."""
        self.flush()
        if not self.enabled:
            return

        def close():
            for w in self.writers.values():
                w.close()

        try:
            self._root_io("closing the output files", close)
        finally:
            self.writers = {}

    # --------------------------------------------------------------- writing
    def _series(self, dg, attrs=None):
        from .h5 import SeriesWriter

        name = dg.d.name
        if name not in self.writers:
            a = {"diagnostic": name, "kind": dg.d.kind, "every_steps": dg.r.every, "dt": self.sim.cfg.dt}
            a.update(attrs or {})
            self.writers[name] = SeriesWriter(os.path.join(self.cfg.dir, "reduced", f"{name}.h5"), a,
                                              self.cfg.codec, self.cfg.level)
        return self.writers[name]

    def _row(self, data: dict) -> dict:
        """A row of a reduced series; with a moving window it carries the shift of its output, in
        cells (OUT-D5): the x of the row in the laboratory is the x of the series plus shift * dx."""
        if self.sim.window:
            return dict(data, window_shift=np.int64(self.sim.window_shift))
        return data

    def _write(self, dg, step, data):
        from . import h5

        sim, d = self.sim, dg.d
        t = sim.time
        meta = {"version": __version__}
        if d.kind == "fields":
            path = os.path.join(self.cfg.dir, "fields", d.name, f"{d.name}_{step:08d}.h5")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            h5.write_meshes(path, d.name, step, t, sim.cfg.dt, data["meshes"], data["grid"], self.units, meta,
                            self.cfg.codec, self.cfg.level, self.edpic())
        elif d.kind == "particles":
            path = os.path.join(self.cfg.dir, "particles", d.name, f"{d.name}_{step:08d}.h5")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            h5.write_particles(path, d.name, step, t, sim.cfg.dt, data, self.units, meta, self.cfg.codec,
                               self.cfg.level, self.edpic())
        elif d.kind == "tracks":
            for sp, rows in data.items():
                key = f"{d.name}/{sp}"
                if key not in self.writers:
                    ids = dg.ids[sp]
                    self.writers[key] = h5.TrackWriter(
                        os.path.join(self.cfg.dir, "tracks", f"{d.name}_{sp}.h5"),
                        {"diagnostic": d.name, "species": sp, "every_steps": dg.r.every, "dt": sim.cfg.dt,
                         "time_offsets": json.dumps({"position": 0.0, "momentum": self.time_centering()["momentum"],
                                                     "fields": 0.0})},
                        self.cfg.codec, self.cfg.level, ids)
                self.writers[key].append(step, t, rows)
        elif d.kind in ("lineout", "modes"):
            w = self._series(dg)
            w.append(step, data["t_mean"], self._row(data["data"]))
        elif d.kind == "slice":
            w = self._series(dg, {"axis": d.axis, "plane_position": json.dumps(dg.plane_coords())})
            w.append(step, t, self._row(data))
        elif d.kind == "kspectrum":
            w = self._series(dg, {"k": dg.k_bins.tolist()})
            w.append(step, t, self._row(data))
        elif d.kind == "histogram":
            w = self._series(dg, {"quantities": json.dumps(list(d.quantities)),
                                  "edges": json.dumps([e.tolist() for e in dg.edges])})
            w.append(step, t, self._row(data))

    def _write_scalars(self):
        """Hand the history rows not yet taken to the (chunk-buffered) scalars writer (rank 0).
        A row counts as written once the writer has taken it: if creating the file or its
        datasets fails, nothing is taken and the next write retries the same rows."""
        from .h5 import SeriesWriter

        if not self.cfg.scalars:
            return
        h = self.sim.history
        n = len(h["step"])
        if n <= self.scalar_rows:
            return
        if "__scalars__" not in self.writers:
            attrs = {"diagnostic": "scalars", "every_steps": self.sim.cfg.diag_interval}
            self.writers["__scalars__"] = SeriesWriter(os.path.join(self.cfg.dir, "scalars.h5"), attrs,
                                                       self.cfg.codec, self.cfg.level)
        w = self.writers["__scalars__"]
        sl = slice(self.scalar_rows, n)
        data = {k: np.asarray(h[k][sl], dtype=np.float64) for k in SCALAR_KEYS}
        for sp, v in self.sim.kinetic_species.items():
            data[f"kinetic_{sp}"] = np.asarray(v[sl], dtype=np.float64)
        if self.sim.window:  # the record of a step is taken before the shift at its end (PIC-D22)
            data["window_shift"] = np.array([window.shift_after(self.sim.cfg, int(s)) for s in h["step"][sl]],
                                            dtype=np.int64)
        try:
            w.buffer_rows(h["step"][sl], h["t"][sl], data)
        except Exception:
            if w.ds is None:  # no dataset yet: drop the writer, the next attempt recreates the file
                w.f.close()
                del self.writers["__scalars__"]
            raise
        self.scalar_rows = n
        w.write_due()


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)

