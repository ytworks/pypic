"""Reader for pypic output directories (CPU only: h5py + hdf5plugin, no CuPy).

    run = open_run("output/landau")
    run.scalars["kinetic"], run.scalars.t            energy history
    lo = run.series("lineout"); lo["Ex"], lo.t, lo.x("Ex")   reduced time series
    f = run.fields("fields"); f.get(f.iterations[-1], "Ez")  -> Field(data, coords, time)
    p = run.particles("sample").get(it, "electron")  -> dict x y z ux uy uz w id gamma
    tr = run.tracks("tracks", "electron"); tr["ux"]   (n_t, n_tracked), NaN once a particle is gone

Values are in the normalised units of the run (c = 1, lengths c/omega_pe, time 1/omega_pe,
fields m_e c omega_pe / e, densities n0); ``run.units`` converts to SI when the run had n0_m3.
"""
from __future__ import annotations

import glob
import json
import os
import re
from dataclasses import dataclass

import h5py
import numpy as np

try:
    import hdf5plugin  # noqa: F401  (registers the zstd filter for reading)
except ImportError:  # pragma: no cover
    pass

from ..output.h5 import Units
from ..output.names import expand_fields, record_of, stagger_of


def open_run(path: str) -> "Run":
    return Run(path)


@dataclass
class Field:
    data: np.ndarray
    coords: tuple  # 1-D coordinate arrays (x, y, z) of the samples, in c/omega_pe
    time: float
    time_offset: float  # time level relative to `time` (e.g. J at t - dt/2)

    @property
    def extent(self):
        """matplotlib imshow extent of a 2-D cut (first two coordinate axes)."""
        return _extent(self.coords[0], self.coords[1])


def _extent(a, b):
    da = a[1] - a[0] if len(a) > 1 else 1.0
    db = b[1] - b[0] if len(b) > 1 else 1.0
    return [a[0] - da / 2, a[-1] + da / 2, b[0] - db / 2, b[-1] + db / 2]


class Run:
    def __init__(self, path: str):
        self.path = path
        with open(os.path.join(path, "run.json")) as f:
            self.info = json.load(f)
        i = self.info
        self.dt = float(i["dt"])
        self.n_steps = int(i["n_steps"])
        self.n_cells = tuple(i["n_cells"])
        self.dx = tuple(i["dx"])
        self.L = tuple(i["L"])
        self.solver = i["solver"]
        self.species = [s["name"] for s in i["species"]]
        self.species_info = {s["name"]: s for s in i["species"]}
        self.units = Units(i["units"]["n0_m3"])
        self.diag = {d["name"]: d for d in i["diagnostics"]}
        self._scalars = None

    @property
    def meta(self) -> dict:
        return self.info.get("meta", {})

    @property
    def config(self) -> dict:
        return self.info["config"]

    @property
    def timescales(self) -> dict:
        return self.info.get("timescales", {})

    def reload(self):
        """Re-read run.json (it is rewritten when the experiment adds metadata)."""
        self.__init__(self.path)

    # ------------------------------------------------------------- access
    @property
    def scalars(self) -> "Series":
        if self._scalars is None:
            self._scalars = Series(self, os.path.join(self.path, "scalars.h5"), {"kind": "scalars"})
        return self._scalars

    def _layout(self, name, kinds):
        if name not in self.diag:
            raise KeyError(f"no diagnostic {name!r} in {self.path} (have {sorted(self.diag)})")
        d = self.diag[name]
        if d["kind"] not in kinds:
            raise TypeError(f"diagnostic {name!r} is a {d['kind']}, not one of {kinds}")
        return d

    def series(self, name: str) -> "Series":
        d = self._layout(name, ("lineout", "modes", "slice", "kspectrum", "histogram"))
        return Series(self, os.path.join(self.path, d["path"]), d)

    def fields(self, name: str) -> "FieldSeries":
        return FieldSeries(self, self._layout(name, ("fields",)))

    def particles(self, name: str) -> "ParticleSeries":
        return ParticleSeries(self, self._layout(name, ("particles",)))

    def tracks(self, name: str, species: str | None = None) -> "Tracks":
        d = self._layout(name, ("tracks",))
        sp = species or d["species"][0]
        return Tracks(os.path.join(self.path, d["path"] % sp), d, sp)

    def names(self, kind: str | None = None) -> list[str]:
        return [n for n, d in self.diag.items() if kind is None or d["kind"] == kind]

    def summary(self) -> str:
        """Diagnostics, cadence and the storage actually used (with the compression achieved)."""
        lines = [f"run {self.path}: solver={self.solver} cells={self.n_cells} dt={self.dt:.5g} "
                 f"steps={self.n_steps} GPUs={self.info['world_size']} commit={self.info.get('git_commit')}"]
        tot = 0
        rows = [("scalars", "scalars", self.info["config"]["diag_interval"], "scalars.h5", None, None)]
        rows += [(n, d["kind"], d["every_steps"], d["path"], d["bytes_per_output"], d["n_outputs"])
                 for n, d in self.diag.items()]
        lines.append(f"  {'diagnostic':14s} {'kind':10s} {'every':>6s} {'dt_out':>8s} {'files':>5s} "
                     f"{'on disk':>10s} {'raw est.':>10s} {'ratio':>6s}")
        for n, kind, every, path, bpo, nout in rows:
            pattern = os.path.join(self.path, path.replace("%08d", "*").replace("%s", "*"))
            files = sorted(glob.glob(pattern))
            size = sum(os.path.getsize(f) for f in files)
            tot += size
            raw = bpo * nout if bpo is not None else None
            ratio = f"{raw / size:.2f}" if raw and size else "-"
            lines.append(f"  {n:14s} {kind:10s} {every:6d} {every * self.dt:8.4g} {len(files):5d} "
                         f"{_fmt(size):>10s} {_fmt(raw) if raw else '-':>10s} {ratio:>6s}")
        lines.append(f"  total on disk {_fmt(tot)}")
        return "\n".join(lines)


def _fmt(b):
    if b is None:
        return "-"
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if b < 1000 or unit == "TB":
            return f"{b:.3g} {unit}"
        b /= 1000.0


class Series:
    """A reduced-diagnostic time series (lineout, modes, slice, kspectrum, histogram, scalars).
    Data are read on access, so a file that is still being written is seen up to date."""

    def __init__(self, run: Run, path: str, layout: dict):
        self.run, self.path, self.layout = run, path, layout
        self.f = h5py.File(path, "r", swmr=True)

    def close(self):
        self.f.close()

    def _ds(self, k):
        ds = self.f[k]
        ds.refresh()
        return ds

    @property
    def t(self) -> np.ndarray:
        return self._ds("t")[:]

    @property
    def step(self) -> np.ndarray:
        return self._ds("step")[:]

    def keys(self) -> list[str]:
        return [k for k in self.f.keys() if k not in ("t", "step")]

    def __getitem__(self, k) -> np.ndarray:
        return self._ds(k)[:]

    def __contains__(self, k) -> bool:
        return k in self.f

    def read(self, k, sel=slice(None)) -> np.ndarray:
        """Read part of a quantity, e.g. ``read('Ex', slice(0, 100))``."""
        return self._ds(k)[sel]

    @property
    def attrs(self) -> dict:
        return dict(self.f.attrs)

    @property
    def components(self) -> list[str]:
        return list(self.layout.get("components", self.keys()))

    # ---------------------------------------------------------- coordinates
    def x(self, comp: str) -> np.ndarray:
        """Lineout: sample positions along the profile axis (Yee staggering included)."""
        a = "xyz".index(self.layout["axis"])
        n = self.run.n_cells[a]
        return (np.arange(n) + stagger_of(comp)[a]) * self.run.dx[a]

    def plane_axes(self, comp: str) -> tuple:
        """Slice: coordinates of the two in-plane axes (block-averaged sample positions)."""
        a = "xyz".index(self.layout["axis"])
        c = self.run.config["output"]["diagnostics"]
        dcfg = next(d for d in c if d["name"] == self.layout["name"])
        co = dcfg.get("coarsen", (1, 1, 1))
        out = []
        for b in range(3):
            if b == a:
                continue
            n = -(-self.run.n_cells[b] // co[b])
            pos = ((co[b] - 1) / 2 + stagger_of(comp)[b]) / co[b]
            out.append((np.arange(n) + pos) * self.run.dx[b] * co[b])
        return tuple(out)

    def extent(self, comp: str):
        a, b = self.plane_axes(comp)
        return _extent(a, b)

    @property
    def edges(self) -> list[np.ndarray]:
        """Histogram: bin edges per quantity."""
        return [np.asarray(e) for e in json.loads(self.f.attrs["edges"])]

    @property
    def quantities(self) -> list[str]:
        return json.loads(self.f.attrs["quantities"])

    @property
    def k(self) -> np.ndarray:
        """kspectrum: centres of the |k| bins."""
        return np.asarray(self.f.attrs["k"])

    @property
    def modes(self) -> list:
        return self.layout.get("modes", [])

    def k_of_mode(self, j: int) -> np.ndarray:
        """Wave vector of mode j of a modes diagnostic."""
        return 2 * np.pi * np.asarray(self.layout["modes"][j], dtype=float) / np.asarray(self.run.L)


_IT_RE = re.compile(r"_(\d+)\.h5$")


class _DumpSeries:
    def __init__(self, run: Run, layout: dict):
        self.run, self.layout = run, layout
        pattern = os.path.join(run.path, layout["path"].replace("%08d", "*"))
        self.files = {int(_IT_RE.search(f).group(1)): f for f in glob.glob(pattern)}

    @property
    def iterations(self) -> list[int]:
        return sorted(self.files)

    @property
    def times(self) -> np.ndarray:
        return np.array([it * self.run.dt for it in self.iterations])

    def nearest(self, t: float) -> int:
        its = self.iterations
        return its[int(np.argmin(np.abs(np.asarray(its) * self.run.dt - t)))]


class FieldSeries(_DumpSeries):
    """openPMD mesh dumps of one fields diagnostic."""

    @property
    def components(self) -> list[str]:
        return list(self.layout["components"])

    def get(self, it: int, comp: str, si: bool = False) -> Field:
        rec, c = record_of(expand_fields([comp])[0])
        with h5py.File(self.files[it], "r") as f:
            g = f[f"data/{it}"]
            m = g["meshes"][rec]
            ds = m[c] if c else m
            data = ds[()]
            pos = ds.attrs["position"]
            spacing = m.attrs["gridSpacing"]
            off = m.attrs["gridGlobalOffset"]
            to = float(m.attrs["timeOffset"])
            time = float(g.attrs["time"])
            if si:
                data = data * ds.attrs["unitSI"]
        coords = tuple(off[a] + (np.arange(data.shape[a]) + pos[a]) * spacing[a] for a in range(3))
        return Field(data, coords, time, to)


class ParticleSeries(_DumpSeries):
    """openPMD particle dumps (decomposition-independent samples, sorted by id)."""

    def get(self, it: int, species: str) -> dict:
        with h5py.File(self.files[it], "r") as f:
            s = f[f"data/{it}/particles/{species}"]
            m = float(s["mass"].attrs["value"])
            out = {}
            for c in "xyz":
                out[c] = s["position"][c][()] + s["positionOffset"][c].attrs["value"]
                out["u" + c] = s["momentum"][c][()] / m
            out["w"] = s["weighting"][()] / self.run.units.weight
            out["id"] = s["id"][()]
            out["time"] = float(f[f"data/{it}"].attrs["time"])
            out["momentum_time_offset"] = float(s["momentum"].attrs["timeOffset"])
        out["gamma"] = np.sqrt(1.0 + out["ux"] ** 2 + out["uy"] ** 2 + out["uz"] ** 2)
        return out


class Tracks:
    """Tracked particles: quantities x y z ux uy uz (and Ex..Bz at the particle) as (n_t, n_tracked)."""

    def __init__(self, path: str, layout: dict, species: str):
        self.path, self.layout, self.species = path, layout, species
        self.f = h5py.File(path, "r", swmr=True)

    def _ds(self, k):
        ds = self.f[k]
        ds.refresh()
        return ds

    @property
    def t(self):
        return self._ds("t")[:]

    @property
    def ids(self):
        return self.f["id"][:]

    def keys(self):
        return [k for k in self.f.keys() if k not in ("t", "step", "id")]

    def __getitem__(self, k):
        if k == "gamma":
            return np.sqrt(1.0 + self["ux"] ** 2 + self["uy"] ** 2 + self["uz"] ** 2)
        return self._ds(k)[:]

    @property
    def time_offsets(self) -> dict:
        return json.loads(self.f.attrs["time_offsets"])

    def close(self):
        self.f.close()
