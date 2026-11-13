"""HDF5 writers (CPU side, rank 0 only): appendable time series for reduced diagnostics and
openPMD 1.1 (ED-PIC) file-based series for field and particle dumps (OUT-D1, OUT-D3).

Arrays arrive as NumPy arrays already cast (and optionally mantissa-rounded) on the GPU.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os

import h5py
import numpy as np

try:  # zstd and other filters for HDF5 (registers them with the HDF5 library on import)
    import hdf5plugin
except ImportError:  # pragma: no cover - hdf5plugin is a hard dependency; keep deflate usable
    hdf5plugin = None

OPENPMD_VERSION = "1.1.0"


# ------------------------------------------------------------------ compression
def compression(codec: str, level: int) -> dict:
    """h5py dataset keywords for the lossless codec (byte shuffle + zstd/deflate)."""
    if codec == "none":
        return {}
    if codec == "deflate":
        return {"compression": "gzip", "compression_opts": int(level), "shuffle": True}
    if codec == "zstd":
        if hdf5plugin is None:
            raise RuntimeError("codec 'zstd' needs the hdf5plugin package (or set output.codec: deflate)")
        return {**hdf5plugin.Zstd(clevel=int(level)), "shuffle": True}
    raise ValueError(f"unknown codec {codec!r}")


def round_mantissa(x, keep_bits: int, xp=np):
    """Round a float32/float64 array to ``keep_bits`` explicit mantissa bits (round to nearest,
    ties to even). The relative error is at most 2^-(keep_bits+1); the zeroed low bits make the
    following lossless codec compress much better. Works with NumPy and CuPy (``xp``)."""
    x = xp.ascontiguousarray(x)
    if x.dtype == xp.float32:
        ut, nm = xp.uint32, 23
    elif x.dtype == xp.float64:
        ut, nm = xp.uint64, 52
    else:
        raise TypeError("round_mantissa needs float32 or float64 data")
    drop = nm - int(keep_bits)
    if drop <= 0:
        return x
    b = x.view(ut)
    one = ut(1)
    half = ut((1 << (drop - 1)) - 1)
    lsb = (b >> ut(drop)) & one
    mask = ~ut((1 << drop) - 1)
    r = ((b + half + lsb) & mask).view(x.dtype)
    # keep inf/nan untouched (rounding could turn a large finite value into inf only at the top binade)
    return xp.where(xp.isfinite(x), r, x)


# --------------------------------------------------------------- reduced series
class SeriesWriter:
    """One HDF5 file per reduced diagnostic: datasets ``t``, ``step`` and one dataset per
    quantity of shape (n_out, ...), extended along time.

    Rows are buffered and written one whole chunk at a time: rewriting a partly filled compressed
    chunk would leak its old file space (HDF5 does not reuse it), which for frequent small
    outputs made files several times larger than their data. ``flush()`` (end of a run, close)
    writes a partial chunk. The file is in SWMR mode, so it can be read while the run goes on
    (``h5py.File(path, 'r', swmr=True)``); readers see complete chunks."""

    def __init__(self, path: str, attrs: dict, codec: str, level: int, chunk_bytes: int = 256 << 10,
                 max_rows: int = 256, append: bool = False):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        self.path = path
        self.comp = compression(codec, level)
        self.chunk_bytes, self.max_rows = chunk_bytes, max_rows
        self.ds = None
        self.n = 0  # rows on disk
        self.buf = []  # pending rows: (step, t, data)
        if append and os.path.exists(path):  # a run continued from a checkpoint (OUT-D6)
            with h5py.File(path, "r") as f:
                committed = int(f["step"].shape[0]) if "step" in f else 0
            if committed == 0:
                # no committed row: the file may hold a schema that a killed writer left half made
                # (only ``t``, or ``t`` and ``step``); it is made anew (review of PR #50)
                os.remove(path)
                append = False
        if append and os.path.exists(path):
            self.f = h5py.File(path, "a", libver="latest")
            if "step" in self.f:
                self.ds = {k: v for k, v in self.f.items() if k not in self.FIXED}
                self.n = int(self.ds["step"].shape[0])
                rows = [v.chunks[0] for k, v in self.ds.items() if k not in ("t", "step")]
                self.rows = int(rows[0]) if rows else 1
                self.f.swmr_mode = True
            return
        self.f = h5py.File(path, "w", libver="latest")
        for k, v in attrs.items():
            self.f.attrs[k] = _attr(v)

    FIXED = ("id",)  # datasets written once, not along time

    def _create(self, data: dict):
        row = max([np.asarray(v).nbytes for v in data.values()] + [8])
        self.rows = int(max(1, min(self.max_rows, self.chunk_bytes // row)))
        ds = {}  # set as self.ds only once every dataset exists
        for k, dtype in (("t", np.float64), ("step", np.int64)):
            ds[k] = self.f.create_dataset(k, shape=(0,), maxshape=(None,), dtype=dtype,
                                          chunks=(max(self.rows, 64),), **self.comp)
        for k, v in data.items():
            v = np.asarray(v)
            comp = self.comp if v.dtype.kind in "fc" else {}
            # an empty dimension (e.g. no tracked particle) needs a positive chunk size, which HDF5
            # accepts only for an unlimited dimension
            ds[k] = self.f.create_dataset(k, shape=(0,) + v.shape,
                                          maxshape=(None,) + tuple(n if n else None for n in v.shape),
                                          dtype=v.dtype, chunks=(self.rows,) + tuple(max(n, 1) for n in v.shape),
                                          **comp)
        self._before_swmr()
        self.f.swmr_mode = True
        self.ds = ds

    def _before_swmr(self):
        pass

    def append(self, step: int, t: float, data: dict):
        if self.ds is None:
            self._create(data)
        self.buf.append((step, t, data))
        self.write_due()

    def buffer_rows(self, steps, ts, data: dict):
        """Take many rows (scalars history) into the buffer. The only file operation is creating
        the datasets on first use: if this raises, no row was taken. Once it returns, the rows are
        the writer's: write_due(), flush() or close() writes them, and a failed write keeps them
        buffered for the next attempt (a write is repeatable: it rewrites rows n.. of the file)."""
        rows = [(int(steps[i]), float(ts[i]), {k: np.asarray(v)[i] for k, v in data.items()})
                for i in range(len(steps))]
        if not rows:
            return
        if self.ds is None:
            self._create(rows[0][2])
        self.buf.extend(rows)

    def write_due(self):
        """Write the buffer once it holds at least a whole chunk."""
        if len(self.buf) >= self.rows:
            self._write()

    def _write(self):
        """Write the buffered rows at rows n.. of the file and flush. The rows leave the buffer
        only after the flush succeeded: a failure anywhere keeps them, and the next write repeats
        the same resize, assignment and flush (rewriting rows n.., so repeating is harmless)."""
        if not self.buf:
            return
        m = len(self.buf)
        n = self.n
        for ds in self.ds.values():
            ds.resize(n + m, axis=0)
        self.ds["t"][n:] = np.array([b[1] for b in self.buf], dtype=np.float64)
        self.ds["step"][n:] = np.array([b[0] for b in self.buf], dtype=np.int64)
        for k in self.buf[0][2]:
            self.ds[k][n:] = np.stack([np.asarray(b[2][k]) for b in self.buf])
        self.f.flush()  # SWMR readers see the rows only after this
        self.n = n + m
        # a partial chunk was written: the next write continues in a new chunk-aligned block
        self.buf = []

    def flush(self):
        """Write pending rows (a partial chunk) so that readers see everything recorded so far."""
        if self.f is not None:
            self._write()

    def set_attrs(self, name: str, attrs: dict):
        for k, v in attrs.items():
            self.ds[name].attrs[k] = _attr(v)

    def close(self):
        if self.f is not None:
            self._write()
            self.f.close()
            self.f = None


class TrackWriter(SeriesWriter):
    """Tracks: like a series but every quantity is (n_out, n_tracked) and the tracked ids are
    stored once in ``id`` (written before SWMR starts); a particle that has left the domain is
    NaN from then on."""

    def __init__(self, path: str, attrs: dict, codec: str, level: int, ids: np.ndarray, append: bool = False):
        self.ids = np.asarray(ids, dtype=np.uint64)
        if append and os.path.exists(path):
            # A continued run appends to tracks of the same particles. A file whose rows were all cut back
            # (the checkpoint came before its first output) belongs to no output of this run: it is made
            # anew for the ids of this run (review of PR #50). Rows of other particles are refused.
            with h5py.File(path, "r") as f:
                old = f["id"][...] if "id" in f else None
                n_rows = int(f["step"].shape[0]) if "step" in f else 0
            if old is None or not np.array_equal(old, self.ids):
                if n_rows:
                    raise ValueError(f"{path}: the tracked ids differ from those of the run being continued")
                os.remove(path)
                append = False
        super().__init__(path, attrs, codec, level, chunk_bytes=1 << 20, max_rows=64, append=append)

    def _before_swmr(self):
        self.f.create_dataset("id", data=self.ids)


_M32 = 0xFFFFFFFF


def _rot(x, k):
    return ((x << k) | (x >> (32 - k))) & _M32


def _lookup3(data: bytes) -> int:
    """Jenkins' lookup3 hash (hashlittle, initial value 0): the checksum of HDF5 metadata."""
    a = b = c = (0xDEADBEEF + len(data)) & _M32
    n, k = len(data), 0
    while n > 12:
        a = (a + int.from_bytes(data[k: k + 4], "little")) & _M32
        b = (b + int.from_bytes(data[k + 4: k + 8], "little")) & _M32
        c = (c + int.from_bytes(data[k + 8: k + 12], "little")) & _M32
        for x, y, z, r in ((0, 2, 1, 4), (1, 0, 2, 6), (2, 1, 0, 8), (0, 2, 1, 16), (1, 0, 2, 19), (2, 1, 0, 4)):
            v = [a, b, c]
            v[x] = (v[x] - v[y]) & _M32
            v[x] ^= _rot(v[y], r)
            v[y] = (v[y] + v[z]) & _M32
            a, b, c = v
        n, k = n - 12, k + 12
    if n == 0:
        return c
    tail = data[k:] + bytes(12 - n)
    a = (a + int.from_bytes(tail[0:4], "little")) & _M32
    b = (b + int.from_bytes(tail[4:8], "little")) & _M32
    c = (c + int.from_bytes(tail[8:12], "little")) & _M32
    for x, y, r in ((2, 1, 14), (0, 2, 11), (1, 0, 25), (2, 1, 16), (0, 2, 4), (1, 0, 14), (2, 1, 24)):
        v = [a, b, c]
        v[x] ^= v[y]
        v[x] = (v[x] - _rot(v[y], r)) & _M32
        a, b, c = v
    return c


_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def clear_open_flags(path: str) -> bool:
    """Clear the file consistency flags that a writer left in the superblock (``h5clear -s``).

    A file opened for (SWMR) writing carries these flags until it is closed; a writer that was
    killed leaves them set, and HDF5 then refuses to open the file ("file is already open for
    write"). A run continued from a checkpoint (OUT-D6) is the only writer of its directory, so the
    flags of its series files are those of the interrupted run. Only a version 2/3 superblock has
    them (the files of this module, ``libver="latest"``); its checksum is verified before the
    rewrite. Returns whether the flags were set."""
    with open(path, "r+b") as f:
        size = os.fstat(f.fileno()).st_size
        base = 0
        while True:  # the superblock is at 0 or after a user block of 512, 1024, 2048, ... bytes
            f.seek(base)
            if f.read(8) == _SIGNATURE:
                break
            base = 512 if base == 0 else 2 * base
            if base >= size:
                raise OSError(f"{path}: no HDF5 superblock")
        head = f.read(4)
        version, offsets, flags = head[0], head[1], head[3]
        if version < 2 or flags == 0:
            return False
        n = 12 + 4 * offsets  # signature .. root object header address, then the checksum
        f.seek(base)
        sb = bytearray(f.read(n + 4))
        if _lookup3(bytes(sb[:n])) != int.from_bytes(sb[n:], "little"):
            raise OSError(f"{path}: the superblock checksum does not match; not changing the file")
        sb[11] = 0
        sb[n:] = _lookup3(bytes(sb[:n])).to_bytes(4, "little")
        f.seek(base)
        f.write(sb)
    return True


def truncate_series(path: str, step: int) -> int:
    """Cut the rows of a series (or tracks) file at ``step`` and later (a run continued from the
    checkpoint of that step writes them again, OUT-D6). Returns the rows kept. The open-for-write
    flags of a writer that was killed are cleared first (:func:`clear_open_flags`)."""
    clear_open_flags(path)
    with h5py.File(path, "a", libver="latest") as f:
        if "step" not in f:
            return 0
        steps = f["step"][...]
        # The rows kept are the committed prefix: strictly increasing steps before ``step``. A writer
        # killed between extending the datasets and filling them leaves rows of fill zeros, so the
        # array need not be sorted (review of PR #50); the rows before the checkpoint step were flushed
        # at the checkpoint and are complete.
        keep = 0
        while keep < steps.size and steps[keep] < step and (keep == 0 or steps[keep] > steps[keep - 1]):
            keep += 1
        if keep < steps.size:
            for k, ds in f.items():
                if k not in SeriesWriter.FIXED:
                    ds.resize(keep, axis=0)
        return keep


def _attr(v):
    if isinstance(v, (dict, list, tuple)) and not (isinstance(v, (list, tuple)) and v and
                                                   all(isinstance(x, (int, float)) for x in v)):
        return json.dumps(v, default=str)
    if v is None:
        return "None"
    return v


# ---------------------------------------------------------------------- units
class Units:
    """Normalised units of pypic (c = eps0 = 1, n0 -> omega_pe = 1, lengths c/omega_pe,
    charge e, mass m_e) and their SI values when a reference density n0 [m^-3] is given."""

    C = 299792458.0
    QE = 1.602176634e-19
    ME = 9.1093837015e-31
    EPS0 = 8.8541878128e-12

    def __init__(self, n0_m3: float | None):
        self.n0 = n0_m3
        if n0_m3:
            self.wpe = math.sqrt(n0_m3 * self.QE ** 2 / (self.EPS0 * self.ME))
            self.L = self.C / self.wpe
            self.T = 1.0 / self.wpe
            self.E = self.ME * self.C * self.wpe / self.QE
            self.B = self.ME * self.wpe / self.QE
            self.rho = self.QE * n0_m3
            self.J = self.QE * n0_m3 * self.C
            self.n = n0_m3
            self.p = self.ME * self.C
            self.energy = self.ME * self.C ** 2
            self.weight = n0_m3 * self.L ** 3
            self.q = self.QE
            self.m = self.ME
        else:
            for k in ("L", "T", "E", "B", "rho", "J", "n", "p", "energy", "weight", "q", "m"):
                setattr(self, k, 1.0)
            self.wpe = None

    def as_dict(self):
        return {"n0_m3": self.n0, "normalised": self.n0 is None,
                "length": "c/omega_pe", "time": "1/omega_pe", "E": "m_e c omega_pe / e", "B": "m_e omega_pe / e",
                "density": "n0", "charge_density": "e n0", "current_density": "e n0 c", "momentum": "m_e c",
                "energy": "m_e c^2", "SI": None if self.n0 is None else {
                    "omega_pe": self.wpe, "length_m": self.L, "time_s": self.T, "E_V_per_m": self.E, "B_T": self.B,
                    "charge_density_C_per_m3": self.rho, "current_density_A_per_m2": self.J}}


# openPMD unitDimension: powers of (L, M, T, I, theta, N, J)
DIM = {
    "E": (1, 1, -3, -1, 0, 0, 0), "B": (0, 1, -2, -1, 0, 0, 0), "J": (-2, 0, 0, 1, 0, 0, 0),
    "rho": (-3, 0, 1, 1, 0, 0, 0), "n": (-3, 0, 0, 0, 0, 0, 0), "flux": (-2, 0, -1, 0, 0, 0, 0),
    "ekin": (-1, 1, -2, 0, 0, 0, 0), "position": (1, 0, 0, 0, 0, 0, 0), "momentum": (1, 1, -1, 0, 0, 0, 0),
    "weighting": (0,) * 7, "charge": (0, 0, 1, 1, 0, 0, 0), "mass": (0, 1, 0, 0, 0, 0, 0), "id": (0,) * 7,
    "E_at_particle": (1, 1, -3, -1, 0, 0, 0), "B_at_particle": (0, 1, -2, -1, 0, 0, 0),
}


def unit_si(u: Units, record: str) -> float:
    return {"E": u.E, "B": u.B, "J": u.J, "rho": u.rho, "n": u.n, "flux": u.n * (u.L / u.T),
            "ekin": u.energy * u.n, "position": u.L, "momentum": u.p, "weighting": 1.0,
            "charge": u.q, "mass": u.m, "id": 1.0}[record]


# ------------------------------------------------------------------- openPMD
def _root_attrs(f, name: str, meta: dict):
    f.attrs["openPMD"] = np.bytes_(OPENPMD_VERSION)
    f.attrs["openPMDextension"] = np.uint32(1)  # ED-PIC
    f.attrs["basePath"] = np.bytes_("/data/%T/")
    f.attrs["iterationEncoding"] = np.bytes_("fileBased")
    f.attrs["iterationFormat"] = np.bytes_(f"{name}_%T.h5")
    f.attrs["software"] = np.bytes_("pypic")
    f.attrs["softwareVersion"] = np.bytes_(meta.get("version", "?"))
    f.attrs["date"] = np.bytes_(_dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z"))
    f.attrs["author"] = np.bytes_(meta.get("author", "pypic"))


def _iteration(f, it: int, time: float, dt: float, u: Units):
    g = f.create_group(f"data/{it}")
    g.attrs["time"] = np.float64(time)
    g.attrs["dt"] = np.float64(dt)
    g.attrs["timeUnitSI"] = np.float64(u.T)
    return g


def _record_attrs(obj, record: str, time_offset: float, u: Units):
    obj.attrs["unitDimension"] = np.array(DIM[record], dtype=np.float64)
    obj.attrs["timeOffset"] = np.float32(time_offset)


def write_meshes(path: str, name: str, it: int, time: float, dt: float, meshes: dict, grid: dict, u: Units,
                 meta: dict, codec: str, level: int, edpic: dict):
    """meshes: {record: {"components": {comp or "": array}, "position": {comp: (px,py,pz)},
    "time_offset": float}}; grid: spacing, offset (global, in c/omega_pe)."""
    comp = compression(codec, level)
    with h5py.File(path, "w") as f:
        _root_attrs(f, name, meta)
        f.attrs["meshesPath"] = np.bytes_("meshes/")
        g = _iteration(f, it, time, dt, u)
        m = g.create_group("meshes")
        for k in ("fieldSolver", "fieldBoundary", "particleBoundary", "currentSmoothing", "chargeCorrection"):
            v = edpic[k]
            m.attrs[k] = np.array([np.bytes_(x) for x in v]) if isinstance(v, (list, tuple)) else np.bytes_(v)
        if edpic.get("currentSmoothing") != "none" and "currentSmoothingParameters" in edpic:
            m.attrs["currentSmoothingParameters"] = np.bytes_(edpic["currentSmoothingParameters"])
        for rec, spec in meshes.items():
            kind = rec.split("_", 1)[0]
            r = m.create_group(rec) if len(spec["components"]) > 1 or "" not in spec["components"] else None
            for c, arr in spec["components"].items():
                if r is None:
                    ds = m.create_dataset(rec, data=arr, chunks=_chunks(arr), **comp)
                    holder = ds
                else:
                    ds = r.create_dataset(c, data=arr, chunks=_chunks(arr), **comp)
                    holder = r
                ds.attrs["unitSI"] = np.float64(unit_si(u, kind))
                ds.attrs["position"] = np.array(spec["position"][c], dtype=np.float32)
            for obj in ([holder] if r is None else [r]):
                _record_attrs(obj, kind, spec["time_offset"], u)
                obj.attrs["geometry"] = np.bytes_("cartesian")
                obj.attrs["dataOrder"] = np.bytes_("C")
                obj.attrs["axisLabels"] = np.array([b"x", b"y", b"z"])
                obj.attrs["gridSpacing"] = np.array(grid["spacing"], dtype=np.float64)
                obj.attrs["gridGlobalOffset"] = np.array(grid["offset"], dtype=np.float64)
                obj.attrs["gridUnitSI"] = np.float64(u.L)
                obj.attrs["fieldSmoothing"] = np.bytes_("none")


def _chunks(arr):
    if arr.ndim == 0 or arr.size < 1024:
        return None
    # ~1-4 MB chunks: whole planes of the last two axes, a few along the first
    plane = max(1, arr[0].nbytes)
    n0 = int(max(1, min(arr.shape[0], (2 << 20) // plane)))
    return (n0,) + arr.shape[1:]


def _p_record(s, n, rname, comps, u, comp, time_offset=0.0, macro=1, wpow=1.0):
    """Particle record (scalar if comps has the single key '', else a group of components)."""
    chunks = None if n < 1024 else (min(n, 1 << 20),)
    if list(comps) == [""]:
        arr = comps[""]
        holder = s.create_dataset(rname, data=arr, chunks=chunks, **(comp if arr.dtype.kind == "f" else {}))
        holder.attrs["unitSI"] = np.float64(unit_si(u, rname))
    else:
        holder = s.create_group(rname)
        for c, arr in comps.items():
            ds = holder.create_dataset(c, data=arr, chunks=chunks, **comp)
            ds.attrs["unitSI"] = np.float64(unit_si(u, rname))
    _record_attrs(holder, rname, time_offset, u)
    holder.attrs["macroWeighted"] = np.uint32(macro)
    holder.attrs["weightingPower"] = np.float64(wpow)


def _p_constant(s, n, rname, value, u, comps=("",)):
    """Constant particle record (openPMD 'value' + 'shape')."""
    unit_rec = "position" if rname == "positionOffset" else rname
    r = s.create_group(rname)
    targets = [r] if comps == ("",) else [r.create_group(c) for c in comps]
    for t in targets:
        t.attrs["value"] = np.float64(value)
        t.attrs["shape"] = np.array([n], dtype=np.uint64)
        t.attrs["unitSI"] = np.float64(unit_si(u, unit_rec))
    _record_attrs(r, unit_rec, 0.0, u)
    r.attrs["macroWeighted"] = np.uint32(0)
    r.attrs["weightingPower"] = np.float64(1.0 if rname in ("charge", "mass") else 0.0)


def write_particles(path: str, name: str, it: int, time: float, dt: float, species: dict, u: Units, meta: dict,
                    codec: str, level: int, edpic: dict):
    """species: {name: {"x","y","z","ux","uy","uz","w": float arrays, "id": uint64, "charge", "mass",
    "momentum_offset": float}}; momentum is written as m u (m_e c units)."""
    comp = compression(codec, level)
    with h5py.File(path, "w") as f:
        _root_attrs(f, name, meta)
        f.attrs["particlesPath"] = np.bytes_("particles/")
        g = _iteration(f, it, time, dt, u)
        pg = g.create_group("particles")
        for sname, d in species.items():
            s = pg.create_group(sname)
            for k in ("particleShape", "currentDeposition", "particlePush", "particleInterpolation",
                      "particleSmoothing"):
                v = edpic[k]
                s.attrs[k] = np.float32(v) if k == "particleShape" else np.bytes_(v)
            n = len(d["id"])
            m = d["mass"]
            _p_record(s, n, "position", {"x": d["x"], "y": d["y"], "z": d["z"]}, u, comp, macro=0, wpow=0.0)
            _p_constant(s, n, "positionOffset", 0.0, u, ("x", "y", "z"))
            _p_record(s, n, "momentum", {c: d["u" + c] * d["u" + c].dtype.type(m) for c in "xyz"}, u, comp,
                      time_offset=d["momentum_offset"], macro=0, wpow=1.0)
            # openPMD: weighting = number of physical particles per macro-particle, unitSI 1
            # (normalised runs: in units of n0 (c/omega_pe)^3)
            _p_record(s, n, "weighting", {"": d["w"] * d["w"].dtype.type(u.weight)}, u, comp, macro=1, wpow=1.0)
            _p_record(s, n, "id", {"": d["id"]}, u, comp, macro=0, wpow=0.0)
            _p_constant(s, n, "charge", d["charge"], u)
            _p_constant(s, n, "mass", d["mass"], u)
