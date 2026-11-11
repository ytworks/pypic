"""GPU tests of the file output (OUT-D1..D4): the files hold the simulation state, they do not
depend on the number of GPUs, and particle ids are unique and decomposition independent."""
from __future__ import annotations

import glob
import os
import signal
import subprocess
import sys

import cupy as cp
import h5py
import numpy as np
import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic.analysis import open_run
from pypic.config import Body, Dipole, ExternalFields, OutputConfig
from pypic.parallel.comm import SelfComm
from pypic.simulation import Simulation

HERE = os.path.dirname(__file__)


def launch(n: int, *args, timeout: float = 600) -> subprocess.CompletedProcess:
    """Run a script with ``pypic.launch -n n`` in its own process group; on a timeout the whole
    group (launcher and ranks) is killed, so a hung case cannot leave ranks spinning on the GPUs."""
    p = subprocess.Popen([sys.executable, "-m", "pypic.launch", "-n", str(n), *args], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        raise
    return subprocess.CompletedProcess(p.args, p.returncode, out, err)


def test_files_hold_the_state(tmp_path):
    """Fields, a coarsened field, a lineout, a histogram and a full particle dump written at the
    last step equal the in-memory state (float32)."""
    n = 8
    diags = [dict(name="f", kind="fields", fields=["Ex", "rho"], every=n),
             dict(name="fc", kind="fields", fields=["Ex"], coarsen=[2, 2, 4], every=n),
             dict(name="lo", kind="lineout", fields=["Ex"], axis="y", every=n),
             dict(name="h", kind="histogram", quantities=["x", "ux"], bins=[16, 8], ranges=[None, [-5, 5]], every=n),
             dict(name="p", kind="particles", every=n)]
    cfg = SimConfig(n_cells=(16, 8, 8), dx=(0.2,) * 3, dt=0.08, n_steps=n,
                    species=[SpeciesConfig("e", -1, 1, [Beam(0.5, (0.3, 0, 0), 0.01), Beam(0.5, (-0.3, 0, 0), 0.01)]),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                    ppc=4, diag_interval=2, output=OutputConfig(dir=str(tmp_path), diagnostics=diags))
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    sim.close()
    run = open_run(str(tmp_path))
    s = sim.grid.physical_slice()
    Ex = cp.asnumpy(sim.grid.Ex[s])
    f = run.fields("f").get(n, "Ex")
    assert np.array_equal(f.data, Ex.astype(np.float32))
    assert np.allclose(f.coords[0], (np.arange(16) + 0.5) * 0.2)  # Ex is staggered by half a cell in x
    fc = run.fields("fc").get(n, "Ex").data
    assert np.allclose(fc, Ex.reshape(8, 2, 4, 2, 2, 4).mean(axis=(1, 3, 5)), rtol=1e-6, atol=1e-7)
    assert np.allclose(run.series("lo")["Ex"][-1], Ex.mean(axis=(0, 2)), rtol=1e-6, atol=1e-8)
    h = run.series("h")
    for sp in sim.species:
        assert abs(h.read(sp.name, -1).sum() - float(cp.sum(sp.P[6, : sp.n]))) < 1e-6 * float(cp.sum(sp.P[6, : sp.n]))
        p = run.particles("p").get(n, sp.name)
        ids = cp.asnumpy(sp.P[7, : sp.n])
        assert len(np.unique(p["id"])) == sp.n and np.array_equal(np.sort(ids).astype(np.uint64), p["id"])
        order = np.argsort(ids)
        assert np.array_equal(p["x"], cp.asnumpy(sp.P[0, : sp.n])[order].astype(np.float32))
    assert len(run.scalars.t) == len(sim.history["t"])
    assert run.info["time_centering"]["momentum"] == -0.5  # explicit leapfrog: u at t - dt/2


def test_field_outputs_hold_the_total_field_with_external_fields(tmp_path):
    """With external fields (PIC-D12) E and B of a field dump and of the tracks are the fields the
    particles feel: the dipole at its place plus the uniform fields. The field energies of the
    history are those of the Yee fields, which stay zero with test particles only."""
    E, B0, m, c, core = (0.01, 0.02, -0.03), (0.1, 0.0, 0.05), np.array([0.0, 0.5, -2.0]), (8.1, 7.9, 8.2), 1.0
    n = 4
    pts = [[1.0, 2.0, 3.0, 0, 0, 0], [14.2, 8.0, 3.3, 0, 0, 0], [8.0, 14.6, 12.1, 0, 0, 0], [2.3, 13.1, 9.0, 0, 0, 0]]
    diags = [dict(name="f", kind="fields", fields=["E", "B"], every=n),
             dict(name="tr", kind="tracks", species=["t"], stride=1, every=n, with_fields=True, dtype="float64")]
    cfg = SimConfig(n_cells=(32, 32, 32), dx=(0.5,) * 3, dt=0.05, n_steps=n, diag_interval=1, gauss_interval=0,
                    species=[SpeciesConfig("t", -1, 1, [], test=True, particles=pts)],
                    external_fields=ExternalFields(E=E, B=B0, dipoles=[Dipole(tuple(m), c, core)]),
                    body=Body(center=c, radius=core),  # a dipole needs a body around its core (PIC-D14)
                    output=OutputConfig(dir=str(tmp_path), diagnostics=diags))
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    sim.run()
    sim.close()

    def total_B(x, y, z):  # regularised dipole + uniform field at the points
        r = np.stack([x - c[0], y - c[1], z - c[2]])
        r2 = np.sum(r * r, axis=0)
        mr = np.tensordot(m, r, axes=1)
        return [(3.0 * r[k] * mr - m[k] * (r2 - 2.0 * core**2)) / (r2 + core**2) ** 2.5 + B0[k] for k in range(3)], r2

    run = open_run(str(tmp_path))
    assert run.info["field_convention"] == {"E_B": "total", "field_energy": "yee", "external_fields": True}
    fs = run.fields("f")
    for k, name in enumerate("xyz"):
        assert np.array_equal(fs.get(n, "E" + name).data, np.full((32,) * 3, E[k], dtype=np.float32))
        f = fs.get(n, "B" + name)
        ref, r2 = total_B(*np.meshgrid(*f.coords, indexing="ij"))
        far = r2 > 64.0  # beyond 8 core radii: the face values are within 0.5 % of the formula
        dipole = np.sqrt(sum((ref[j] - B0[j]) ** 2 for j in range(3)))
        assert np.all(np.abs(f.data - ref[k])[far] < 0.01 * dipole[far])
    tr = run.tracks("tr", "t")
    ref, r2 = total_B(tr["x"][-1], tr["y"][-1], tr["z"][-1])
    dipole = np.sqrt(sum((ref[k] - B0[k]) ** 2 for k in range(3)))  # |B| of the dipole at every particle
    assert r2.min() > 36.0
    for k, name in enumerate("xyz"):
        assert np.allclose(tr["E" + name][-1], E[k], rtol=1e-14, atol=0)
        assert np.all(np.abs(tr["B" + name][-1] - ref[k]) < 0.01 * dipole)
    assert max(sim.history["magnetic"]) == 0.0 and max(sim.history["electric"]) == 0.0


def _datasets(path):
    out = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            out[name] = obj[()]
    with h5py.File(path, "r") as f:
        f.visititems(visit)
    return out


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_a_second_run_does_not_mix_with_the_first(tmp_path):
    """Review of PR #2: a re-run into the same directory kept the first run's other dumps, which
    the reader listed as part of the new run. The directory is refused unless output.overwrite,
    which removes the old run's files."""

    def run(steps, overwrite=False):
        cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.2,) * 3, dt=0.08, n_steps=steps,
                        species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.01)]),
                                 SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                        ppc=2, output=OutputConfig(dir=str(tmp_path), overwrite=overwrite, diagnostics=[
                            dict(name="f", kind="fields", fields=["Ex"], every=1)]))
        sim = Simulation(cfg, comm=SelfComm(), verbose=False)
        sim.run()
        sim.close()

    run(2)
    first = open_run(str(tmp_path)).info["run_id"]
    with pytest.raises(FileExistsError):
        run(1)
    run(1, overwrite=True)
    run_ = open_run(str(tmp_path))
    assert run_.fields("f").iterations == [0, 1]
    assert run_.info["run_id"] != first  # each run has its own id (one evidence record per run)


def test_a_failed_output_write_stops_every_rank(tmp_path):
    """Review of PR #2: when rank 0 could not prepare the output directory (e.g. old files it may
    not delete under output.overwrite) it raised alone, and the other ranks went on into the next
    collective. Two ranks run in threads, with a stand-in for the NCCL reduction."""
    import threading
    import types

    from pypic.output.manager import OutputError, OutputManager

    blocker = tmp_path / "a_file"
    blocker.write_text("")
    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.2,) * 3, dt=0.08, n_steps=1,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.01)]),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                    output=OutputConfig(dir=str(blocker / "run")))  # cannot be created
    barrier, vals, errors = threading.Barrier(2, timeout=30), [0.0, 0.0], [None, None]

    class Comm:
        size = 2

        def __init__(self, rank):
            self.rank = rank

        def allreduce_scalar(self, v, op="sum"):
            vals[self.rank] = v
            barrier.wait()
            out = max(vals) if op == "max" else sum(vals)
            barrier.wait()
            return float(out)

    def rank(r):
        try:
            OutputManager(types.SimpleNamespace(cfg=cfg, rank=r, comm=Comm(r)))
        except Exception as exc:  # noqa: BLE001 - checked below
            errors[r] = exc

    threads = [threading.Thread(target=rank, args=(r,)) for r in (0, 1)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(60)
    assert all(isinstance(e, OutputError) for e in errors), errors


def test_scalar_rows_of_a_failed_write_are_written_on_retry(tmp_path):
    """Review of PR #2: the scalar history rows were marked written before the write, so after a
    failed write (here: a directory where scalars.h5 goes) a retry wrote nothing."""
    from pypic.output.manager import OutputError

    cfg = SimConfig(n_cells=(8, 8, 8), dx=(0.2,) * 3, dt=0.08, n_steps=4, diag_interval=1,
                    species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0, 0, 0), 0.01)]),
                             SpeciesConfig("i", 1, 25, [Beam(1.0, (0, 0, 0), 0.001)])],
                    ppc=2, output=OutputConfig(dir=str(tmp_path)))
    sim = Simulation(cfg, comm=SelfComm(), verbose=False)
    (tmp_path / "scalars.h5").mkdir()  # the file cannot be created
    with pytest.raises(OutputError):
        sim.run()
    (tmp_path / "scalars.h5").rmdir()
    sim.close()
    with h5py.File(tmp_path / "scalars.h5", "r") as f:
        assert list(f["step"][:]) == list(sim.history["step"]) and len(sim.history["step"]) > 0


def test_output_is_independent_of_the_gpu_count(tmp_path):
    """Every output file of a mixed-boundary run (all diagnostic kinds, particles sorted by id)
    is the same on 1 and 2 GPUs to float32 round-off."""
    script = os.path.join(HERE, "mgpu_case.py")
    for n in (1, 2):
        r = launch(n, script, "explicit", str(tmp_path / f"g{n}"))
        assert r.returncode == 0, r.stderr[-2000:]
    import hdf5plugin  # noqa: F401

    files = sorted(glob.glob(str(tmp_path / "g1" / "**" / "*.h5"), recursive=True))
    assert len(files) > 10
    for fa in files:
        fb = fa.replace(os.sep + "g1" + os.sep, os.sep + "g2" + os.sep)
        A, B = _datasets(fa), _datasets(fb)
        assert A.keys() == B.keys(), fa
        for k in A:
            if os.path.basename(fa) == "scalars.h5" and k in ("gauss", "divB", "fill"):
                continue  # round-off residuals themselves; the fill of the fullest rank's arrays
            a, b = A[k], B[k]
            assert a.shape == b.shape, (fa, k)
            if a.dtype.kind in "fc":
                fin = np.isfinite(a)
                assert np.array_equal(fin, np.isfinite(b)), (fa, k)
                scale = np.max(np.abs(a[fin])) if fin.any() else 0.0
                assert np.max(np.abs(a[fin] - b[fin]), initial=0.0) <= 2e-6 * scale + 1e-12, (fa, k)
            else:
                assert np.array_equal(a, b), (fa, k)


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_particle_on_upper_reflecting_wall_is_kept():
    """Regression (PR #2): a particle exactly at x = L on a decomposed reflecting axis was dropped."""
    r = launch(2, os.path.join(HERE, "mgpu_regress.py"), "edge", timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "RESULT particles=1" in r.stdout, r.stdout


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs >= 2 GPUs")
def test_implicit_error_on_one_rank_is_raised_on_every_rank():
    """Regression (PR #2): a rank-local implicit-step error left the other ranks in a collective,
    so a caller that caught ImplicitStepError hung."""
    try:
        r = launch(2, os.path.join(HERE, "mgpu_regress.py"), "implicit_err", timeout=180)
    except subprocess.TimeoutExpired:
        pytest.fail("ranks hung after a rank-local implicit error")
    assert r.returncode == 0, r.stderr[-2000:]
    for rank in (0, 1):
        assert f"RESULT rank={rank} raised=True particles=4" in r.stdout, r.stdout
