"""CPU tests of the output layer: cadence from physics, HDF5/openPMD writers + reader, lossy
mantissa rounding, and the analysis tools on synthetic data (no CuPy, no GPU)."""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.config import OutputConfig
from pypic.output import h5
from pypic.output.schedule import report, resolve


def _cfg(diags, ts=None, budget=None, n_steps=1000, dt=0.05):
    return SimConfig(n_cells=(64, 16, 16), dx=(0.1,) * 3, dt=dt, n_steps=n_steps, ppc=8,
                     species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 1e-3)])],
                     output=OutputConfig(dir="x", diagnostics=diags, timescales=ts or {}, budget_GB=budget))


# ------------------------------------------------------------------ schedule
def test_cadence_forms_and_auto_takes_the_strictest_requirement():
    ts = {"shortest_period": 2 * math.pi, "growth_time": 2.0}
    cfg = _cfg([{"name": "a", "kind": "lineout", "fields": ["Ex"], "every": 7},
                {"name": "b", "kind": "lineout", "fields": ["Ex"], "every": {"t": 0.5}},
                {"name": "c", "kind": "lineout", "fields": ["Ex"], "every": {"count": 10}},
                {"name": "d", "kind": "modes", "fields": ["Ex"], "modes": [[1, 0, 0]]},
                {"name": "e", "kind": "lineout", "fields": ["Ex"], "every": {"t": 5.0}}], ts)
    r = {x.name: x for x in resolve(cfg)}
    assert (r["a"].every, r["b"].every, r["c"].every) == (7, 10, 100)
    # auto: min(Nyquist of 2 pi oversampled 2x -> pi/2 = 31 steps, growth_time/10 = 0.2 -> 4 steps)
    assert r["d"].every == 4 and r["d"].average
    assert r["c"].n_outputs == 11 and r["c"].due(500) and not r["c"].due(550)
    assert any("aliasing" in n for n in r["e"].notes)  # 5.0 > pi without time averaging
    assert "lineout" in report(cfg, list(r.values()))


def test_budget_widens_only_auto_cadences():
    big = {"name": "f", "kind": "fields", "fields": ["E", "B"]}  # auto: 21 dumps of 393 kB (8.3 MB)
    r = resolve(_cfg([big], budget=0.002))[0]
    assert r.total_bytes <= 0.002e9 and "widened" in r.notes[0]
    assert r.n_outputs < resolve(_cfg([big]))[0].n_outputs
    with pytest.raises(ValueError, match="budget"):
        resolve(_cfg([{**big, "every": 1}], budget=0.002))


# ---------------------------------------------------------- writers / reader
def test_round_mantissa_bound():
    rng = np.random.default_rng(1)
    x = (rng.standard_normal(100000) * 10.0 ** rng.integers(-20, 20, 100000)).astype(np.float32)
    for bits in (3, 10, 18):
        r = h5.round_mantissa(x, bits)
        assert np.max(np.abs(r - x) / np.abs(x)) <= 2.0 ** -(bits + 1)


def _write_run(path, steps=(0, 10), n=(8, 4, 4), dx=0.5, dt=0.1):
    """A synthetic run directory as the engine lays it out (run.json + one dump of each kind)."""
    u = h5.Units(None)
    ed = dict(fieldSolver="Yee", fieldBoundary=["periodic"] * 6, particleBoundary=["periodic"] * 6,
              currentSmoothing="none", chargeCorrection="none", particleShape=1.0, currentDeposition="Esirkepov",
              particlePush="Boris", particleInterpolation="trilinear", particleSmoothing="none")
    x = (np.arange(n[0]) + 0.5) * dx
    for it in steps:
        os.makedirs(f"{path}/fields/f", exist_ok=True)
        ex = np.broadcast_to(np.sin(x + it * dt)[:, None, None], n).astype(np.float32)
        h5.write_meshes(f"{path}/fields/f/f_{it:08d}.h5", "f", it, it * dt, dt,
                        {"E": {"components": {"x": ex}, "position": {"x": (0.5, 0, 0)}, "time_offset": 0.0}},
                        {"spacing": (dx,) * 3, "offset": (0, 0, 0)}, u, {"version": "t"}, "zstd", 3, ed)
        os.makedirs(f"{path}/particles/p", exist_ok=True)
        m = 1500
        d = {k: np.linspace(0, 1, m, dtype=np.float32) for k in ("x", "y", "z", "ux", "uy", "uz", "w")}
        d.update(id=np.arange(m, dtype=np.uint64), charge=-1.0, mass=2.0, momentum_offset=-0.05)
        h5.write_particles(f"{path}/particles/p/p_{it:08d}.h5", "p", it, it * dt, dt, {"e": d}, u, {"version": "t"},
                           "zstd", 3, ed)
    w = h5.SeriesWriter(f"{path}/reduced/lo.h5", {"diagnostic": "lo"}, "zstd", 3)
    for k in range(300):
        w.append(k, k * dt, {"Ex": np.full(n[0], k, np.float32)})
    w.close()
    info = {"format": "pypic-output-1", "dt": dt, "n_steps": steps[-1], "n_cells": list(n), "dx": [dx] * 3,
            "L": [k * dx for k in n], "solver": "explicit", "world_size": 1, "species": [{"name": "e"}],
            "units": u.as_dict(), "config": {"diag_interval": 1, "output": {"diagnostics": []}},
            "diagnostics": [
                {"name": "f", "kind": "fields", "path": "fields/f/f_%08d.h5", "components": ["Ex"], "every_steps": 10,
                 "bytes_per_output": 512, "n_outputs": 2},
                {"name": "p", "kind": "particles", "path": "particles/p/p_%08d.h5", "species": ["e"],
                 "every_steps": 10, "bytes_per_output": 1, "n_outputs": 2},
                {"name": "lo", "kind": "lineout", "path": "reduced/lo.h5", "axis": "x", "components": ["Ex"],
                 "every_steps": 1, "bytes_per_output": 32, "n_outputs": 300}]}
    with open(f"{path}/run.json", "w") as f:
        json.dump(info, f)


def test_writers_reader_roundtrip_and_openpmd(tmp_path):
    from pypic.analysis import open_run

    _write_run(str(tmp_path))
    run = open_run(str(tmp_path))
    fs = run.fields("f")
    assert fs.iterations == [0, 10]
    fd = fs.get(10, "Ex")
    assert np.allclose(fd.coords[0], (np.arange(8) + 0.5) * 0.5) and abs(fd.time - 1.0) < 1e-12
    assert np.allclose(fd.data[:, 0, 0], np.sin(fd.coords[0] + 1.0), atol=1e-6)
    p = run.particles("p").get(0, "e")
    assert np.allclose(p["ux"], np.linspace(0, 1, 1500)) and p["id"][-1] == 1499  # momentum m*u -> u
    lo = run.series("lo")
    assert lo["Ex"].shape == (300, 8) and lo["Ex"][299, 0] == 299  # rows buffered into whole chunks
    assert "lo" in run.summary()
    checker = shutil.which("openPMD_check_h5", path=os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"])
    if checker is None:
        pytest.skip("openPMD-validator not installed")
    for f in (fs.files[10], run.particles("p").files[10]):
        out = subprocess.run([checker, "-i", f], capture_output=True, text=True).stdout
        assert "Result: 0 Errors" in out, out


# ------------------------------------------------------------------ analysis
def test_omega_k_ridge_and_damped_fit():
    from pypic.analysis.spectral import fit_damped, omega_k, ridge

    dt, dx, nt, nx = 0.1, 0.2, 2048, 64
    t = np.arange(nt) * dt
    x = np.arange(nx) * dx
    k = 2 * np.pi * np.fft.rfftfreq(nx, dx)[1:12]
    w = np.sqrt(1 + 3 * 0.1 * k ** 2)  # Bohm-Gross-like branch
    f = sum(np.cos(kk * x[None, :] - ww * t[:, None] + 0.3 * j) for j, (kk, ww) in enumerate(zip(k, w)))
    kk, ww, P = omega_k(f, dt, dx)
    rk, rw, _ = ridge(kk, ww, P, k_range=(k[0] * 0.99, k[-1] * 1.01))
    assert np.max(np.abs(rw - np.interp(rk, k, w))) < 2 * np.pi / (nt * dt) / 2
    z = 0.01 * np.exp(-0.07 * t[:600]) * np.cos(1.3 * t[:600]) + 1e-9
    r = fit_damped(t[:600], z.astype(complex), omega0=1.2, gamma0=-0.05)
    assert abs(r["gamma"] + 0.07) < 1e-3 and abs(r["omega"] - 1.3) < 1e-3


def test_power_law_fit_recovers_index():
    from pypic.analysis.particles import fit_power_law, spectrum

    rng = np.random.default_rng(7)
    p_true = 2.4
    E = (1.0 - rng.random(400000)) ** (-1.0 / (p_true - 1.0))  # Pareto: dN/dE ~ E^-2.4, E >= 1
    edges = np.geomspace(1, 1e3, 61)
    Ec, f = spectrum(edges, np.histogram(E, edges)[0])
    r = fit_power_law(Ec, f, lo=2.0, hi=200.0)
    assert abs(r["p"] - p_true) < 5 * r["p_err"] + 0.02


def test_movie_writes_an_animation(tmp_path):
    from pypic.analysis import Movie, open_run

    _write_run(str(tmp_path))
    run = open_run(str(tmp_path))
    out = Movie(run).lineout("lo", "Ex").xt("lo", "Ex").save(str(tmp_path / "m.gif"), fps=10, n_frames=6)
    from PIL import Image

    assert Image.open(out).n_frames == 6


def test_movie_shows_a_moving_window_at_its_laboratory_place(tmp_path):
    """Three outputs of a plane across y and of histograms over (x, u_x) and over (u_x, x); the window has
    gone 0, 4 and 10 cells of 0.5 at them: the pictures of the last frame lie 5 further in x than the first
    (OUT-D5), on whichever axis of the histogram x is (review of PR #42: only the first one moved)."""
    from types import SimpleNamespace

    from pypic.analysis import Movie
    from pypic.analysis.animate import PhasePanel, SlicePanel

    shift, t = np.array([0, 4, 10]), np.array([0.0, 1.0, 2.0])
    data = np.arange(32.0).reshape(8, 4) + 1.0

    def extent(comp, i=None):
        x0 = 0.0 if i is None else 0.5 * shift[i]
        return [x0, x0 + 4.0, 0.0, 2.0]

    plane = SimpleNamespace(t=t, layout={"axis": "y"}, window_shift=shift, extent=extent,
                            read=lambda c, i: data * (i + 1))
    hist = SimpleNamespace(t=t, window_shift=shift, quantities=["x", "ux"], read=lambda sp, i: data,
                           edges=[np.linspace(0.0, 4.0, 9), np.linspace(-1.0, 1.0, 5)])
    turned = SimpleNamespace(t=t, window_shift=shift, quantities=["ux", "x"], read=lambda sp, i: data.T,
                             edges=hist.edges[::-1])
    run = SimpleNamespace(series={"plane": plane, "phase": hist, "turned": turned}.get, dx=[0.5] * 3)
    a, b, c = SlicePanel(run, "plane", "Ex"), PhasePanel(run, "phase", "e"), PhasePanel(run, "turned", "e")
    Movie(run).add(a).add(b).add(c).save(str(tmp_path / "m.gif"), fps=5)
    assert list(a.im.get_extent()) == [5.0, 9.0, 0.0, 2.0] and a.ax.get_xlim() == (5.0, 9.0)
    assert list(b.im.get_extent()) == [5.0, 9.0, -1.0, 1.0] and b.ax.get_xlim() == (5.0, 9.0)
    assert list(c.im.get_extent()) == [-1.0, 1.0, 5.0, 9.0] and c.ax.get_ylim() == (5.0, 9.0)


class _Tracks:
    """Minimal stand-in for analysis.run.Tracks (float32 rows, as written by default)."""

    def __init__(self, t, **q):
        self.t = np.asarray(t, dtype=np.float64)
        self.q = {k: np.asarray(v, dtype=np.float32) for k, v in q.items()}

    def __getitem__(self, k):
        return self.q[k]


def test_parallel_perpendicular_work_keeps_the_work_where_B_vanishes():
    """Review of PR #2: float32 tracks with B = 0 gave NaN directions and the work vanished."""
    from pypic.analysis.particles import parallel_perpendicular_work, work_by_component

    t = np.array([0.0, 0.5, 1.0])
    one = np.ones((3, 2))
    zero = np.zeros((3, 2))
    # particle 0: no B; particle 1: B along x (E has parallel and perpendicular parts)
    tr = _Tracks(t, ux=np.stack([t, t], 1), uy=zero, uz=zero, Ex=one, Ey=0.5 * one, Ez=zero,
                 Bx=np.stack([zero[:, 0], one[:, 1]], 1), By=zero, Bz=zero)
    w = work_by_component(tr, 1.0, 1.0)
    total = w["W_x"] + w["W_y"] + w["W_z"]
    pp = parallel_perpendicular_work(tr, 1.0, 1.0)
    assert total[-1, 0] > 0.4
    assert np.allclose(pp["W_par"] + pp["W_perp"], total, rtol=1e-6, atol=0)
    assert pp["W_par"][-1, 0] == 0 and np.isclose(pp["W_par"][-1, 1], total[-1, 1], rtol=1e-6)


def test_two_wave_fit_returns_the_signal_and_handles_zero_samples():
    """Review of PR #2: a, b were returned for the window start instead of t = 0, and complex64
    modes with exact zeros got infinite weights (LinAlgError)."""
    from pypic.analysis.spectral import fit_damped, fit_two_wave

    t = np.linspace(3, 16, 261)
    g0, w0, a0, b0 = -0.04, 1.2, 2 + 0.3j, 0.4 - 0.1j
    z = np.exp(g0 * t) * (a0 * np.exp(-1j * w0 * t) + b0 * np.exp(1j * w0 * t))
    g, w, a, b = fit_two_wave(t, z, -0.03, 1.1)
    rec = np.exp(g * t) * (a * np.exp(-1j * w * t) + b * np.exp(1j * w * t))
    assert np.max(np.abs(rec - z)) < 1e-6 * np.max(np.abs(z))
    t = np.linspace(0, 12, 241)
    z = (np.exp(-0.08 * t) * np.sin(1.2 * t)).astype(np.complex64)  # standing wave, z(0) = 0 exactly
    r = fit_damped(t, z, omega0=1.2, gamma0=-0.08, t_max=8)
    assert abs(r["gamma"] + 0.08) < 1e-4 and abs(r["omega"] - 1.2) < 1e-4


def test_tracks_with_no_selected_particle(tmp_path):
    """Review of PR #2: a tracks diagnostic that selects no particle of a species (id stride, or all
    particles gone before t_start) failed with a zero chunk size; it records empty rows."""
    from pypic.analysis.run import Tracks

    path = str(tmp_path / "tr.h5")
    w = h5.TrackWriter(path, {"time_offsets": json.dumps({"x": 0.0})}, "zstd", 3, ids=np.zeros(0))
    for i in range(3):
        w.append(10 * i, 0.5 * i, {"x": np.zeros(0, np.float32), "ux": np.zeros(0, np.float32)})
    w.close()
    tr = Tracks(path, {}, "e")
    assert tr["x"].shape == (3, 0) and tr.ids.size == 0 and np.allclose(tr.t, [0, 0.5, 1])


def test_record_example_loads_only_repository_examples(tmp_path):
    """Review of PR #2: meta.example from run.json chose the directory that analyze.py was
    imported from, so a run directory could run its own code."""
    root = os.path.join(os.path.dirname(__file__), "..", "..")
    planted = tmp_path / "planted"
    planted.mkdir()
    marker = tmp_path / "marker"
    (planted / "analyze.py").write_text(f"open({str(marker)!r}, 'w').close()\nraise SystemExit(3)\n")
    run = tmp_path / "run"
    run.mkdir()
    (run / "run.json").write_text(json.dumps({"meta": {"example": str(planted)}}))
    r = subprocess.run([sys.executable, os.path.join(root, "scripts", "record_example.py"), str(run)],
                       capture_output=True, text=True)
    assert r.returncode == 2 and not marker.exists(), r.stderr


def test_series_rows_stay_buffered_until_the_file_is_flushed(tmp_path):
    """Review of PR #2: the rows left the buffer before the HDF5 flush, so after a failed flush the
    next flush() returned without flushing and a SWMR reader saw no row."""
    path = str(tmp_path / "s.h5")
    w = h5.SeriesWriter(path, {}, "zstd", 3)
    real, calls = w.f.flush, []

    def flush_failing_once():
        calls.append(1)
        if len(calls) == 1:
            raise OSError("injected flush failure")
        real()

    w.f.flush = flush_failing_once
    w.append(0, 0.0, {"x": np.zeros(2)})
    with pytest.raises(OSError):
        w.flush()
    w.flush()  # the retry
    reader = ("import sys, h5py, hdf5plugin\n"
              "with h5py.File(sys.argv[1], 'r', swmr=True) as f:\n"
              "    print(f['step'][:].tolist())")
    r = subprocess.run([sys.executable, "-c", reader, path], capture_output=True, text=True)
    w.close()
    assert r.stdout.strip() == "[0]", r.stderr
