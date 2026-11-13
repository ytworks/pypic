"""Restart checkpoints (PIC-D25) and the append mode of the output (OUT-D6).

A run of 2N steps against a run that stops after its checkpoint at N and a run that continues from
it: the histories, the fields, the particle numbers and ids, and every file of the output agree.
Two identical runs differ in the last bits (atomic deposits), so the criterion is relative 1e-10,
not bitwise. A checkpoint from 2 GPUs restored on 1 and on 4 GPUs agrees with the restore on 2."""
import json
import os
import shutil
import subprocess
import sys

import cupy as cp
import h5py
import hdf5plugin  # noqa: F401  (registers the zstd filter of the files)
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "mgpu_restart.py")
TOL = 1e-10
EXACT = ("step", "ids", "n_particles", "n_injected", "n_absorbed", "shift", "end")
NOISE = ("gauss", "divB", "rho", "fill")  # round-off quantities, and the fill of the arrays (depends on the layout)


def _run(n, dims, case, mode, d, *extra, check=True):
    r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", str(n), SCRIPT, dims, case, mode, d, *extra],
                       capture_output=True, text=True, timeout=900)
    if not check:
        return r
    assert r.returncode == 0, r.stderr[-3000:]
    return json.loads(r.stdout.strip().splitlines()[-1])


def _rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    assert a.shape == b.shape
    fin = np.isfinite(a)
    assert np.array_equal(fin, np.isfinite(b))
    if not fin.any():
        return 0.0
    scale = np.max(np.abs(a[fin]))
    return float(np.max(np.abs(a[fin] - b[fin])) / scale) if scale else float(np.max(np.abs(b[fin])))


def _files(d):
    names = (os.path.relpath(os.path.join(r, f), d) for r, _, fs in os.walk(d) for f in fs if f.endswith(".h5"))
    return sorted(p for p in names if not p.startswith("checkpoint" + os.sep))


def _compare(ref, o, dref, do):
    """Worst relative difference of the histories, the fields and the output files; exact items equal."""
    worst = 0.0
    for k, v in ref.items():
        if k == "start" or k in NOISE:
            continue
        if k in EXACT:
            assert o[k] == v, k
        else:
            worst = max(worst, _rel(v, o[k]))
    fa, fb = np.load(os.path.join(dref, "fields.npy")), np.load(os.path.join(do, "fields.npy"))
    worst = max([worst] + [_rel(fa[c], fb[c]) for c in range(6)])
    names = _files(dref)
    assert names and names == _files(do)
    for name in names:
        with h5py.File(os.path.join(dref, name), "r") as f, h5py.File(os.path.join(do, name), "r") as g:
            sets = []
            f.visit(sets.append)
            for k in (k for k in sets if isinstance(f[k], h5py.Dataset)):
                a, b = f[k][...], g[k][...]
                if k.split("/")[-1] in NOISE:
                    assert a.shape == b.shape, (name, k)
                elif a.dtype.kind in "fc":
                    worst = max(worst, _rel(np.abs(a) if a.dtype.kind == "c" else a,
                                            np.abs(b) if b.dtype.kind == "c" else b))
                else:
                    assert np.array_equal(a, b), (name, k)
    return worst


@pytest.mark.parametrize("case,n,dims", [("two_stream", 1, "1,1,1"), ("two_stream", 2, "2,1,1"),
                                         ("open_flow", 2, "2,1,1"), ("compact", 1, "1,1,1"), ("window", 1, "1,1,1"),
                                         ("implicit", 1, "1,1,1")])
def test_restart_continues_the_run(tmp_path, case, n, dims):
    """2N steps against N + checkpoint + 3 steps (cut) + restart + N: histories, fields, ids, output."""
    if cp.cuda.runtime.getDeviceCount() < n:
        pytest.skip(f"needs {n} GPUs")
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    ref = _run(n, dims, case, "full", a)
    part = _run(n, dims, case, "part", b)
    o = _run(n, dims, case, "resume", b)
    N = ref["end"] // 2
    assert part["end"] == N + 3 and o["start"] == N and o["end"] == 2 * N
    worst = _compare(ref, o, a, b)
    print(f"[restart {case} {n} GPU] largest relative difference {worst:.2e}")  # diagnostic
    assert worst <= TOL
    if case == "window":
        assert o["shift"] > 0
    if case == "open_flow":
        assert o["bank"] > 0.0 and o["n_injected"][-1] > 0


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 4, reason="needs >= 4 GPUs")
@pytest.mark.parametrize("case", ["two_stream", "open_flow", "compact", "window"])
def test_restart_onto_another_layout(tmp_path, case):
    """A checkpoint written on 2 GPUs, restored on 1 and on 4 GPUs: the same as the restore on 2."""
    base = str(tmp_path / "b")
    _run(2, "2,1,1", case, "part", base)
    dirs = {}
    for n in (2, 1, 4):
        dirs[n] = str(tmp_path / f"b{n}")
        shutil.copytree(base, dirs[n])
    ref = _run(2, "2,1,1", case, "resume", dirs[2])
    for n, dims in ((1, "1,1,1"), (4, "2,2,1" if case != "window" else "4,1,1")):
        o = _run(n, dims, case, "resume", dirs[n])
        assert o["start"] == ref["start"]
        worst = _compare(ref, o, dirs[2], dirs[n])
        print(f"[restart {case} 2 -> {n} GPUs] largest relative difference {worst:.2e}")  # diagnostic
        assert worst <= TOL, n


def test_partial_checkpoint_is_ignored_and_changed_configuration_refused(tmp_path):
    """A checkpoint whose writing did not end (step-<n>.tmp) is not read; a changed seed is refused."""
    a = str(tmp_path / "a")
    ref = _run(1, "1,1,1", "two_stream", "full", a)
    ref_dir = str(tmp_path / "ref")
    shutil.copytree(a, ref_dir)
    ck = os.path.join(a, "checkpoint")
    last = sorted(os.listdir(ck))[-1]
    os.rename(os.path.join(ck, last), os.path.join(ck, last + ".tmp"))  # crash before the rename
    r = _run(1, "1,1,1", "two_stream", "resume", a, "seed", check=False)
    assert r.returncode != 0 and "seed" in r.stderr and "PIC-D25" in r.stderr
    o = _run(1, "1,1,1", "two_stream", "resume", a)
    assert o["start"] == ref["end"] // 2
    worst = _compare(ref, o, ref_dir, a)
    print(f"[restart past a .tmp checkpoint] largest relative difference {worst:.2e}")  # diagnostic
    assert worst <= TOL


def test_restart_after_the_writer_was_killed(tmp_path):
    """A run killed (SIGKILL) after its checkpoint leaves its series files open for SWMR writing; the
    restart clears those flags (OUT-D6) and continues like the uninterrupted run."""
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    ref = _run(1, "1,1,1", "two_stream", "full", a)
    r = _run(1, "1,1,1", "two_stream", "part", b, "kill", check=False)
    assert r.returncode != 0
    with pytest.raises(OSError, match="already open for write"):
        h5py.File(os.path.join(b, "scalars.h5"), "r")
    o = _run(1, "1,1,1", "two_stream", "resume", b)
    assert o["start"] == ref["end"] // 2
    worst = _compare(ref, o, a, b)
    print(f"[restart after SIGKILL] largest relative difference {worst:.2e}")  # diagnostic
    assert worst <= TOL


def _checkpoints(d):
    return sorted(os.listdir(os.path.join(d, "checkpoint")))


def test_stale_checkpoints_are_removed(tmp_path):
    """The checkpoints of an earlier run in the directory are refused with its output and removed by
    output.overwrite; a restart from an earlier checkpoint removes the later ones (they belong to the
    continuation that it cuts). Left in place, they would count as newer when pruning, and a restart
    would take them (PIC-D25)."""
    b = str(tmp_path / "b")
    for step in (90, 100):  # the checkpoints of a longer earlier run
        p = os.path.join(b, "checkpoint", f"step-{step:08d}")
        os.makedirs(p)
        with open(os.path.join(p, "meta.json"), "w") as f:
            json.dump({"format": "pypic-checkpoint-1", "layout": {"world": 1}, "step": step}, f)
        open(os.path.join(p, "rank-0.npz"), "wb").close()
    r = _run(1, "1,1,1", "two_stream", "part", b, check=False)
    assert r.returncode != 0 and "checkpoint" in r.stderr and "overwrite" in r.stderr
    _run(1, "1,1,1", "two_stream", "part", b, "overwrite")
    assert _checkpoints(b) == ["step-00000020"]
    assert _run(1, "1,1,1", "two_stream", "resume", b)["start"] == 20
    a = str(tmp_path / "a")
    _run(1, "1,1,1", "two_stream", "full", a)
    assert _checkpoints(a) == ["step-00000020", "step-00000040"]
    o = _run(1, "1,1,1", "two_stream", "part", a, "from=" + os.path.join(a, "checkpoint", "step-00000020"))
    assert o["start"] == 20 and _checkpoints(a) == ["step-00000020"]


def _snapshot(d):
    out = {}
    for r, _, fs in os.walk(d):
        for f in fs:
            with open(os.path.join(r, f), "rb") as fh:
                out[os.path.relpath(os.path.join(r, f), d)] = fh.read()
    return out


def test_restart_refuses_the_output_of_another_run(tmp_path):
    """A continued run appends only to the run that wrote its checkpoint (the run id) or to an empty
    directory: the output of another run is refused before any file changes, and replaced with
    output.overwrite (OUT-D6)."""
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    _run(1, "1,1,1", "two_stream", "part", a)
    _run(1, "1,1,1", "two_stream", "part", b)
    before = _snapshot(b)
    src = "from=" + os.path.join(a, "checkpoint")
    r = _run(1, "1,1,1", "two_stream", "resume", b, src, check=False)
    assert r.returncode != 0 and "another run" in r.stderr
    assert _snapshot(b) == before
    o = _run(1, "1,1,1", "two_stream", "resume", b, src, "overwrite")
    assert o["start"] == 20
    with h5py.File(os.path.join(b, "reduced", "slice.h5"), "r") as f:
        assert f["step"][0] == 20  # only the continuation: the rows of the other run are gone
    assert _checkpoints(b) == ["step-00000040"]


def test_changed_average_schedule_drops_pending_samples(tmp_path):
    """The averaged lineout has the sample of step 19 pending at the checkpoint of step 20; the
    continued run starts the lineout at step 30. Its output at 30 averages the samples 25..30, as a
    run with that schedule from the start does, without step 19."""
    r, b = str(tmp_path / "r"), str(tmp_path / "b")
    _run(1, "1,1,1", "two_stream", "full", r, "tstart")
    _run(1, "1,1,1", "two_stream", "part", b)
    _run(1, "1,1,1", "two_stream", "resume", b, "tstart")
    name = os.path.join("reduced", "lineout.h5")
    with h5py.File(os.path.join(r, name), "r") as f, h5py.File(os.path.join(b, name), "r") as g:
        assert f["step"][...].tolist() == [30, 36]
        keep = g["step"][...] >= 30
        assert g["step"][keep].tolist() == [30, 36]
        worst = max(_rel(f[k][...], g[k][...][keep]) for k in f if k != "step")
    print(f"[restart with a later t_start] largest relative difference {worst:.2e}")  # diagnostic
    assert worst <= TOL


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_restart_without_particles_onto_another_layout(tmp_path):
    """No particle species, a laser antenna: the reduced counters (the antenna work) survive a
    checkpoint written on 2 GPUs and restored on 1 GPU (PIC-D25)."""
    base = str(tmp_path / "b")
    _run(2, "2,1,1", "vacuum", "part", base)
    dirs = {n: str(tmp_path / f"b{n}") for n in (2, 1)}
    for d in dirs.values():
        shutil.copytree(base, d)
    ref = _run(2, "2,1,1", "vacuum", "resume", dirs[2])
    o = _run(1, "1,1,1", "vacuum", "resume", dirs[1])
    assert o["start"] == ref["start"] == 20 and ref["antenna"][-1] > 0.0
    worst = _compare(ref, o, dirs[2], dirs[1])
    print(f"[restart without particles 2 -> 1 GPU] largest relative difference {worst:.2e}")  # diagnostic
    assert worst <= TOL
