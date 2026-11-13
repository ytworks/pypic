"""Restart checkpoints (PIC-D25): the state of a run at the start of a step, written by every rank,
and its restore onto the same layout of ranks or onto another one.

A checkpoint is the directory ``<output.dir>/checkpoint/step-<n>/`` (n zero-padded to 8 digits):

- ``rank-<r>.npz`` (every rank): the field arrays with their ghosts (E, B, J, the CPML memory, the
  implicit ``Eprev``), the particle stores of every species as the build holds them (FP64 rows, or
  the compact cell, fraction, momentum, weight and id) with n and the capacity, the injection bank,
  the per-rank cumulative counters (``sdiag``, the moving-window and antenna sums), the running sums
  of the averaged diagnostics, the ids of the tracks (rank 0), and the layout of the rank.
- ``meta.json`` (rank 0): format, commit, the configuration, the layout, step and time, the scalar
  state (history, per-species kinetic energy, boundary flux, dissipated energy, Picard count, window
  shift), the per-rank counters reduced over the ranks, the particle numbers, and the state of the
  output manager.

The directory is written as ``step-<n>.tmp`` and renamed after every rank has written it, so a crash
never leaves a partial checkpoint that looks complete; ``.tmp`` directories are never read.

Restore onto the same layout copies the arrays back (the state is that of the run, bit for bit).
Onto another layout: the particles are read from every file in their global form and each rank
keeps those inside its domain (the test of the migration, so each particle is kept once); the owned
blocks of the fields are assembled from the blocks that overlap them (global index ranges) and the
ghosts filled by the exchange; the bank goes by global face column; the reduced counters go to rank
0 and zero elsewhere, so the sums over the ranks keep their totals.

The pure-Python part (configuration comparison, directories, meta.json) needs no CuPy.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import shutil

import numpy as np

FORMAT = "pypic-checkpoint-1"
ROOT = "checkpoint"
_STEP = re.compile(r"^step-(\d+)$")

# Configuration keys that may differ between the run that wrote a checkpoint and the run that continues
# from it: the run length, the layout of the ranks, the output cadences and the checkpoint block itself.
# "*" stands for every element of a list.
ALLOWED = ("n_steps", "dims", "checkpoint", "restart", "diag_interval", "gauss_interval", "output.dir",
           "output.overwrite", "output.budget_GB", "output.timescales",
           "output.diagnostics.*.every", "output.diagnostics.*.t_start", "output.diagnostics.*.t_end")


# ----------------------------------------------------------------- configuration (pure Python)
def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def jsonable(d):
    """``d`` as JSON gives it back: tuples become lists, NumPy scalars Python numbers."""
    return json.loads(json.dumps(d, default=_jsonable))


def _drop(d, path):
    head, _, rest = path.partition(".")
    if head == "*":
        for x in d if isinstance(d, list) else []:
            _drop(x, rest)
    elif isinstance(d, dict) and head in d:
        if rest:
            _drop(d[head], rest)
        else:
            del d[head]


def comparable(cfg_dict: dict) -> dict:
    """The configuration without the keys of ALLOWED."""
    d = jsonable(cfg_dict)
    for path in ALLOWED:
        _drop(d, path)
    return d


def _diff(a, b, path):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            out += _diff(a.get(k, "<missing>"), b.get(k, "<missing>"), f"{path}.{k}" if path else str(k))
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [p for i, (x, y) in enumerate(zip(a, b)) for p in _diff(x, y, f"{path}[{i}]")]
    return [] if a == b else [path or "<root>"]


def config_differences(saved: dict, new: dict) -> list[str]:
    """Paths of the configuration at which ``new`` differs from ``saved`` outside ALLOWED."""
    return _diff(comparable(saved), comparable(new), "")


def check_config(cfg) -> None:
    """A checkpoint goes into the output directory: refuse ``checkpoint.every`` without one."""
    if cfg.checkpoint.every and cfg.output.dir is None:
        raise ValueError("checkpoint.every needs output.dir: checkpoints are written to <output.dir>/checkpoint "
                         "(PIC-D25)")


def check_restart(cfg, meta: dict) -> None:
    diff = config_differences(meta["config"], cfg.to_dict())
    if diff:
        raise ValueError(f"restart from {cfg.restart}: the configuration differs from the stored one at "
                         f"{', '.join(diff)}; only {', '.join(ALLOWED)} may change (PIC-D25)")


# ------------------------------------------------------------------- directories (pure Python)
def root_of(cfg) -> str:
    return os.path.join(cfg.output.dir, ROOT)


def step_name(step: int) -> str:
    return f"step-{int(step):08d}"


def complete_steps(root: str) -> list[tuple[int, str]]:
    """The complete checkpoints in ``root`` as (step, path), oldest first. A ``.tmp`` directory (a
    checkpoint whose writing did not end) or one without meta.json is not a checkpoint."""
    out = []
    if os.path.isdir(root):
        for e in os.listdir(root):
            m = _STEP.match(e)
            p = os.path.join(root, e)
            if m and os.path.isfile(os.path.join(p, "meta.json")):
                out.append((int(m.group(1)), p))
    return sorted(out)


def resolve(path: str) -> str:
    """The checkpoint that ``restart`` names: a directory step-<n>, or a directory of checkpoints (its
    latest complete one)."""
    if os.path.isfile(os.path.join(path, "meta.json")) and _STEP.match(os.path.basename(os.path.normpath(path))):
        return path
    steps = complete_steps(path)
    if not steps:
        raise FileNotFoundError(f"restart: {path} is neither a complete checkpoint step-<n> nor a directory that "
                                "holds one (a step-<n>.tmp directory is a checkpoint whose writing did not end)")
    return steps[-1][1]


def write_meta(path: str, meta: dict) -> None:
    """meta.json of a checkpoint (rank 0); floats keep their bits (JSON repr), NaN included."""
    with open(os.path.join(path, "meta.json"), "w") as f:
        json.dump(meta, f, default=_jsonable)


def read_meta(path: str) -> dict:
    with open(os.path.join(path, "meta.json")) as f:
        meta = json.load(f)
    if meta.get("format") != FORMAT:
        raise ValueError(f"{path}: not a checkpoint of format {FORMAT} (found {meta.get('format')!r})")
    missing = [r for r in range(meta["layout"]["world"]) if not os.path.isfile(rank_file(path, r))]
    if missing:
        raise FileNotFoundError(f"{path}: the files of the ranks {missing} are missing")
    return meta


def rank_file(path: str, r: int) -> str:
    return os.path.join(path, f"rank-{r}.npz")


def prune(root: str, keep: int, current: int) -> None:
    """Keep the last ``keep`` complete checkpoints; remove the ``.tmp`` directories of earlier steps."""
    for _, p in complete_steps(root)[:-keep]:
        shutil.rmtree(p)
    for e in os.listdir(root):
        m = re.match(r"^step-(\d+)\.tmp$", e)
        if m and int(m.group(1)) < current:
            shutil.rmtree(os.path.join(root, e))


# ------------------------------------------------------------------------------- GPU side
def _counters(sim):
    """Per-rank cumulative sums (cupy arrays), in the order of meta['reduced']: sdiag per species,
    the window sums per species and of the field, the antenna work up to now and to the last record."""
    out = list(sim.sdiag)
    if sim.window:
        out += [*sim._wd, sim._window_fe]
    out += [sim._antenna_work, sim._antenna_rec]
    return out


def _rank_layout(sim):
    g, axes = sim.grid, sim.decomp.axes
    return np.array([*(a.cell_lo for a in axes), *(a.n_int for a in axes), *(a.phys_lo for a in axes),
                     *(a.n_phys for a in axes), g.g, *(a.coord for a in axes)], dtype=np.int64)


def _layout_parts(lay):
    lay = np.asarray(lay)
    return {"cell_lo": lay[0:3], "n_int": lay[3:6], "phys_lo": lay[6:9], "n_phys": lay[9:12], "ghost": int(lay[12]),
            "coord": lay[13:16]}


def _grid_arrays(sim):
    """name -> array of the field state (ghosts included)."""
    from .grid import FIELD_NAMES

    g = sim.grid
    out = {f"f.{n}": getattr(g, n) for n in FIELD_NAMES + ("Jx", "Jy", "Jz")}
    if g.use_pml:
        out.update({f"psiE.{i}": p for i, p in enumerate(g.psiE)})
        out.update({f"psiB.{i}": p for i, p in enumerate(g.psiB)})
    return out


def _accumulators(sim):
    """(index of the diagnostic, list of arrays) of the averaged diagnostics with a running sum."""
    out = []
    for k, dg in enumerate(sim.output.diags):
        acc = getattr(dg, "accum", None)
        if acc is not None:
            out.append((k, list(acc) if isinstance(acc, list) else [acc]))
    return out


def save(sim) -> str:
    """Write the checkpoint of the state at the present step (collective: every rank calls it)."""
    import cupy as cp

    cfg, comm = sim.cfg, sim.comm
    n = sim.step_count
    root = root_of(cfg)
    final = os.path.join(root, step_name(n))
    tmp = final + ".tmp"
    out = sim.output
    out.flush()  # rows of the reduced series up to this step are on disk: a restart cuts only later ones
    counters = _counters(sim)
    reduced = cp.concatenate([c.ravel() for c in counters] + [cp.asarray([float(sp.n) for sp in sim.species])])
    reduced = cp.asnumpy(comm.allreduce(reduced) if comm.size > 1 else reduced)
    S = len(sim.species)

    def prepare():
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)
        os.makedirs(tmp)

    out._root_io(f"preparing the checkpoint {tmp}", prepare)
    err = None
    try:
        arrays = {"layout": _rank_layout(sim)}
        arrays.update({k: cp.asnumpy(v) for k, v in _grid_arrays(sim).items()})
        if sim.implicit and sim.Eprev is not None:
            arrays["Eprev"] = cp.asnumpy(cp.stack(sim.Eprev))
        for i, sp in enumerate(sim.species):
            arrays[f"s{i}.nc"] = np.array([sp.n, sp.cap], dtype=np.int64)
            if sim.compact:
                arrays[f"s{i}.W"], arrays[f"s{i}.C"] = cp.asnumpy(sp.W[:, : sp.n]), cp.asnumpy(sp.C[: sp.n])
                arrays[f"s{i}.F"], arrays[f"s{i}.U"] = cp.asnumpy(sp.F[:, : sp.n]), cp.asnumpy(sp.U[:, : sp.n])
            else:
                arrays[f"s{i}.P"] = cp.asnumpy(sp.P[:, : sp.n])
        arrays["counters"] = cp.asnumpy(cp.concatenate([c.ravel() for c in counters]))
        arrays["bank"], arrays["binfo"] = cp.asnumpy(sim.bank), cp.asnumpy(sim.binfo)
        for k, acc in _accumulators(sim):
            for i, a in enumerate(acc):
                arrays[f"acc.{k}.{i}"] = cp.asnumpy(a)
        if sim.rank == 0:
            names = [sp.name for sp in sim.species]
            for k, dg in enumerate(out.diags):
                for name, ids in getattr(dg, "ids", {}).items():
                    arrays[f"tracks.{k}.{names.index(name)}"] = np.asarray(ids)
        np.savez(rank_file(tmp, sim.rank), **arrays)
    except Exception as exc:  # noqa: BLE001 - reduced below, raised on every rank
        err = exc
    if comm.allreduce_scalar(float(err is not None), "max") > 0:
        if err is not None:
            raise RuntimeError(f"writing the checkpoint {tmp} failed: {err!r}") from err
        raise RuntimeError(f"writing the checkpoint {tmp} failed on another rank")

    h = sim.history
    meta = {
        "format": FORMAT, "commit": out.git_commit if out.enabled else None,
        "written": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "config": jsonable(cfg.to_dict()),
        "layout": {"world": comm.size, "dims": list(sim.decomp.dims), "ghost": sim.grid.g},
        "step": n, "time": sim.time, "history": h, "kinetic_species": sim.kinetic_species,
        "flux_out": sim.flux_out, "dissipated": sim.dissipated, "picard_its": getattr(sim, "_picard_its", None),
        "nonconverged": getattr(sim, "nonconverged", None), "window_shift": sim.window_shift,
        "counter_sizes": [int(c.size) for c in counters], "reduced": reduced[:-S].tolist(),
        "particles": [int(round(v)) for v in reduced[-S:]] if S else [],
        "output": {"last": out.last, "acc_step": out.acc_step,
                   "acc": {str(k): {"n": out.diags[k].n_acc, "t": out.diags[k].t_acc,
                                    "list": isinstance(out.diags[k].accum, list)} for k, _ in _accumulators(sim)}},
    }

    def finish():
        write_meta(tmp, meta)
        if os.path.isdir(final):  # a run that went back to an earlier checkpoint reaches this step again
            shutil.rmtree(final)
        os.rename(tmp, final)
        prune(root, cfg.checkpoint.keep, n)

    out._root_io(f"completing the checkpoint {final}", finish)
    return final


def open_restart(sim):
    """Read and check the checkpoint that ``cfg.restart`` names, on every rank, before the set-up
    touches the output directory. Returns (path, meta); every rank raises together."""
    err, res = None, None
    try:
        path = resolve(sim.cfg.restart)
        meta = read_meta(path)
        check_restart(sim.cfg, meta)
        res = (path, meta)
    except Exception as exc:  # noqa: BLE001 - reduced below
        err = exc
    if sim.comm.allreduce_scalar(float(err is not None), "max") > 0:
        if err is not None:
            raise err
        raise RuntimeError(f"restart from {sim.cfg.restart}: reading the checkpoint failed on another rank")
    return res


def _overlap(lo_a, n_a, lo_b, n_b):
    lo, hi = max(lo_a, lo_b), min(lo_a + n_a, lo_b + n_b)
    return (lo, hi - lo) if hi > lo else None


def _index_positions(z, i, lay):
    """Index coordinates x / dx of the compact particles of species i in the file ``z`` (exact)."""
    shape = lay["n_int"] + 2 * lay["ghost"]
    origin = lay["cell_lo"] - lay["ghost"]
    c = z[f"s{i}.C"].astype(np.int64)
    f = z[f"s{i}.F"]
    cells = (c // (shape[1] * shape[2]), (c // shape[2]) % shape[1], c % shape[2])
    return np.stack([(cells[a] + origin[a]) + f[a] * (1.0 / 4294967296.0) for a in range(3)])


def _ghost_overlap(gidx, lo, n, ghost):
    """(local indices, source indices) of the cells ``gidx`` that the owned range [lo, lo + n) holds."""
    m = (gidx >= lo) & (gidx < lo + n)
    return np.nonzero(m)[0], gidx[m] - lo + ghost


def _global_records(sim, z, i, lay):
    """(8, m) rows of species i of a file: position in the coordinate of Species.frame, momentum,
    weight, id."""
    if sim.compact:
        W, U = z[f"s{i}.W"], z[f"s{i}.U"]
        return np.concatenate([_index_positions(z, i, lay), U, W])
    return z[f"s{i}.P"]


def _mine(sim, sp, Q):
    """The particles of Q inside this rank's domain: the bounds that the migration keeps, so that the
    domains of the ranks take every particle once."""
    keep = np.ones(Q.shape[1], dtype=bool)
    for a, ax in enumerate(sim.decomp.axes):
        if ax.nproc == 1:
            continue
        lo, hi, _ = sp.frame(sim.grid, a, sim.window_shift)
        if ax.lo is not None:
            keep &= Q[a] >= lo
        if ax.hi is not None:
            keep &= Q[a] < hi
    return keep


def restore(sim, path: str, meta: dict) -> None:
    """Set the state of ``sim`` (built from the same configuration) to the checkpoint (collective)."""
    import cupy as cp

    lay = meta["layout"]
    same = (lay["world"] == sim.comm.size and tuple(lay["dims"]) == tuple(sim.decomp.dims)
            and lay["ghost"] == sim.grid.g)
    if sim.window:  # first: the bounds of the ranks along x move with the shift (_mine)
        s = int(meta["window_shift"])
        sim.window_shift = s
        sim._window_set()
        for a, off in zip(sim.antennas, sim._antenna_off):
            a.set_shift(s, s < off)
    err = None
    try:
        (_restore_same if same else _restore_other)(sim, path, meta)
    except Exception as exc:  # noqa: BLE001 - reduced below
        err = exc
    if sim.comm.allreduce_scalar(float(err is not None), "max") > 0:
        if err is not None:
            raise RuntimeError(f"restoring the checkpoint {path} failed: {err!r}") from err
        raise RuntimeError(f"restoring the checkpoint {path} failed on another rank")
    if not same:  # ghosts by the exchange (collective)
        g = sim.grid
        g.apply_bc_E()
        g.apply_bc_B()
        g.fill_ghosts((g.Ex, g.Ey, g.Ez))
        g.fill_ghosts((g.Bx, g.By, g.Bz))
        # the CPML memory has its ghosts assembled from the files (_restore_other): only the ranks with
        # an absorbing layer hold it, so an exchange would not be made by both sides of a rank face
        if sim.implicit and sim.Eprev is not None:
            g.apply_bc_E(sim.Eprev)
            g.fill_ghosts(sim.Eprev)
    # scalars (the same on every rank)
    sim.step_count, sim.time = int(meta["step"]), float(meta["time"])
    sim.history = {k: list(v) for k, v in meta["history"].items()}
    sim.kinetic_species = {k: list(v) for k, v in meta["kinetic_species"].items()}
    sim.flux_out, sim.dissipated = float(meta["flux_out"]), float(meta["dissipated"])
    if meta["picard_its"] is not None:
        sim._picard_its = int(meta["picard_its"])
    if meta["nonconverged"] is not None:
        sim.nonconverged = int(meta["nonconverged"])
    out = sim.output
    om = meta["output"]
    out.last, out.acc_step = dict(om["last"]), dict(om["acc_step"])
    for k, a in om["acc"].items():
        dg = out.diags[int(k)]
        dg.n_acc, dg.t_acc = int(a["n"]), float(a["t"])
    out.history_seen = len(sim.history["step"])
    out.after_restore()
    # every particle kept exactly once (collective)
    S = len(sim.species)
    if S:
        tot = sim.comm.allreduce(cp.asarray([float(sp.n) for sp in sim.species]))
        if [int(round(v)) for v in cp.asnumpy(tot)] != meta["particles"]:
            raise RuntimeError(f"restoring the checkpoint {path}: the particle numbers {cp.asnumpy(tot).tolist()} "
                               f"differ from the stored {meta['particles']}")


def _set_counters(sim, flat):
    import cupy as cp

    off = 0
    for c in _counters(sim):
        c[...] = cp.asarray(flat[off: off + c.size]).reshape(c.shape)
        off += c.size


def _set_species(sim, sp, i, z, cap):
    """The store of species i from a file of the same layout, with the capacity of the run."""
    import cupy as cp

    n = int(z[f"s{i}.nc"][0])
    sp.cap = cap
    if sim.compact:
        sp.W = cp.zeros((2, cap), dtype=cp.float64)
        sp.C = cp.zeros(cap, dtype=cp.int32)
        sp.F = cp.zeros((3, cap), dtype=cp.uint32)
        sp.U = cp.zeros((3, cap), dtype=cp.float64)
        sp.Dsp = None
        sp.W[:, :n], sp.C[:n], sp.F[:, :n], sp.U[:, :n] = (cp.asarray(z[f"s{i}.{k}"]) for k in "WCFU")
    else:
        sp.P = cp.zeros((8, cap), dtype=cp.float64)
        sp.P[:, :n] = cp.asarray(z[f"s{i}.P"])
    if sp.Pn is not None:
        sp.Pn = cp.zeros((6, cap), dtype=cp.float64)
    sp._list = None  # the flag list has the size of the capacity, as in the run
    sp.n = n


def _set_accumulators(sim, values):
    """values: diagnostic index -> list of arrays (or None)."""
    import cupy as cp

    for k, dg in enumerate(sim.output.diags):
        if k in values:
            acc = [cp.asarray(a) for a in values[k][0]]
            dg.accum = acc if values[k][1] else acc[0]


def _restore_same(sim, path, meta):
    import cupy as cp

    with np.load(rank_file(path, sim.rank)) as z:
        if not np.array_equal(z["layout"], _rank_layout(sim)):
            raise ValueError("the layout of the rank differs from the stored one")
        for k, arr in _grid_arrays(sim).items():
            arr[...] = cp.asarray(z[k])
        if sim.implicit:
            sim.Eprev = list(cp.asarray(z["Eprev"])) if "Eprev" in z.files else None
        for i, sp in enumerate(sim.species):
            _set_species(sim, sp, i, z, int(z[f"s{i}.nc"][1]))
        _set_counters(sim, z["counters"])
        if z["bank"].shape != tuple(sim.bank.shape):
            raise ValueError("the injection bank differs in size from the stored one")
        sim.bank[...] = cp.asarray(z["bank"])
        acc = {}
        for k, a in meta["output"]["acc"].items():
            k = int(k)
            m = sum(1 for f in z.files if f.startswith(f"acc.{k}."))
            acc[k] = ([z[f"acc.{k}.{i}"] for i in range(m)], a["list"])
        _set_accumulators(sim, acc)
        _restore_tracks(sim, z)


def _restore_tracks(sim, z):
    if sim.rank != 0:
        return
    names = [sp.name for sp in sim.species]
    for f in z.files:
        if f.startswith("tracks."):
            _, k, si = f.split(".")
            sim.output.diags[int(k)].ids[names[int(si)]] = np.asarray(z[f])


def _restore_other(sim, path, meta):
    """Restore onto another layout: every rank reads every file and keeps its part."""
    import cupy as cp

    g, axes, cfg = sim.grid, sim.decomp.axes, sim.cfg
    mine = _layout_parts(_rank_layout(sim))
    targets = _grid_arrays(sim)
    eprev = None
    if sim.implicit:
        eprev = cp.zeros((3,) + g.shape)
    for arr in targets.values():
        arr.fill(0)
    # global cell of every local index along each axis, wrapped on a periodic axis (a rank without an
    # absorbing layer has no CPML memory: it stays zero, as in the run)
    gidx = []
    for a in range(3):
        x = np.arange(g.shape[a]) + mine["cell_lo"][a] - g.g
        gidx.append(x % cfg.n_cells[a] if cfg.boundaries[a][0] == "periodic" else x)
    S = len(sim.species)
    parts = [[] for _ in range(S)]
    # faces with a bank: the global plane of columns of each, (axes b, c of the face)
    planes = {}
    acc_meta = {int(k): a for k, a in meta["output"]["acc"].items()}
    acc_sum = {}
    has_eprev = False
    for q in range(meta["layout"]["world"]):
        with np.load(rank_file(path, q)) as z:
            lay = _layout_parts(z["layout"])
            gq = lay["ghost"]
            ov = [_overlap(lay["cell_lo"][a], lay["n_int"][a], mine["cell_lo"][a], mine["n_int"][a]) for a in range(3)]
            psi = [k for k in targets if k.startswith("psi") and k in z.files]
            if psi:  # with the ghosts: the value of the owner of each (periodic images wrapped)
                pos = [_ghost_overlap(gidx[a], lay["cell_lo"][a], lay["n_int"][a], gq) for a in range(3)]
                if all(p[0].size for p in pos):
                    dsel, ssel = np.ix_(*(p[0] for p in pos)), np.ix_(*(p[1] for p in pos))
                    for k in psi:
                        targets[k][dsel] = cp.asarray(z[k][ssel])
            if all(ov):
                src = tuple(slice(gq + o - lay["cell_lo"][a], gq + o - lay["cell_lo"][a] + m)
                            for a, (o, m) in enumerate(ov))
                dst = tuple(slice(g.g + o - mine["cell_lo"][a], g.g + o - mine["cell_lo"][a] + m)
                            for a, (o, m) in enumerate(ov))
                for k, arr in targets.items():
                    if not k.startswith("psi"):
                        arr[dst] = cp.asarray(z[k][src])
                if eprev is not None and "Eprev" in z.files:
                    has_eprev = True
                    eprev[(slice(None),) + dst] = cp.asarray(z["Eprev"][(slice(None),) + src])
            for i, sp in enumerate(sim.species):
                Q = _global_records(sim, z, i, lay)
                parts[i].append(Q[:, _mine(sim, sp, Q)])
            binfo = z["binfo"]
            for f in range(6):
                if binfo[f] < 0:
                    continue
                a = f >> 1
                b, c = (1, 2) if a == 0 else (0, 2) if a == 1 else (0, 1)
                pl = planes.setdefault(f, np.zeros((cfg.n_cells[b], cfg.n_cells[c])))
                nb, nc = binfo[9 + b], binfo[9 + c]
                pl[binfo[6 + b]: binfo[6 + b] + nb, binfo[6 + c]: binfo[6 + c] + nc] = \
                    z["bank"][binfo[f]: binfo[f] + nb * nc].reshape(nb, nc)
            for k, a in acc_meta.items():
                m = sum(1 for f in z.files if f.startswith(f"acc.{k}."))
                arrs = [z[f"acc.{k}.{i}"] for i in range(m)]
                if a["list"]:  # modes: x profiles of the rank's physical cells -> the global x axis
                    tot = acc_sum.setdefault(k, [np.zeros(cfg.n_cells[0], dtype=x.dtype) for x in arrs])
                    for t, x in zip(tot, arrs):
                        t[lay["phys_lo"][0]: lay["phys_lo"][0] + lay["n_phys"][0]] += x
                else:  # lineout: an array along the whole axis, the rank's part filled
                    tot = acc_sum.setdefault(k, [np.zeros_like(arrs[0])])
                    tot[0] += arrs[0]
            if q == 0:
                _restore_tracks(sim, z)
    for i, sp in enumerate(sim.species):
        Q = np.concatenate(parts[i], axis=1) if parts[i] else np.zeros((8, 0))
        sp.n = 0
        sp._set_records(cp.asarray(Q))
        if sp.Pn is not None and sp.Pn.shape[1] != sp.cap:
            sp.Pn = cp.zeros((6, sp.cap), dtype=cp.float64)
    if sim.implicit:
        sim.Eprev = list(eprev) if has_eprev else None
    bank = np.zeros(sim.bank.size)
    binfo = cp.asnumpy(sim.binfo)
    for f, pl in planes.items():
        if binfo[f] >= 0:
            a = f >> 1
            b, c = (1, 2) if a == 0 else (0, 2) if a == 1 else (0, 1)
            nb, nc = binfo[9 + b], binfo[9 + c]
            bank[binfo[f]: binfo[f] + nb * nc] = pl[binfo[6 + b]: binfo[6 + b] + nb,
                                                    binfo[6 + c]: binfo[6 + c] + nc].ravel()
    sim.bank[...] = cp.asarray(bank)
    _set_counters(sim, np.asarray(meta["reduced"]) if sim.rank == 0 else np.zeros(len(meta["reduced"])))
    acc = {}
    first_x = all(ax.coord == 0 for ax in axes[1:])  # one rank of each x range holds the modes sums
    for k, tot in acc_sum.items():
        if acc_meta[k]["list"]:
            x0, nx = axes[0].phys_lo, axes[0].n_phys
            acc[k] = ([t[x0: x0 + nx] if first_x else np.zeros(nx, dtype=t.dtype) for t in tot], True)
        else:
            acc[k] = ([tot[0] if sim.rank == 0 else np.zeros_like(tot[0])], False)
    _set_accumulators(sim, acc)
