"""Restart checkpoints (PIC-D25) and the append mode of the output (OUT-D6): the configuration, the
comparison of a continued run's configuration with the stored one, the directories of the
checkpoints, meta.json, and the cut and continuation of a series file. No GPU."""
import math
import os

import numpy as np
import pytest

from pypic import Beam, SimConfig, SpeciesConfig
from pypic import checkpoint as ck
from pypic.config import Checkpoint, Diagnostic, OutputConfig


def _cfg(**kw):
    base = dict(n_cells=(16, 8, 8), species=[SpeciesConfig("e", -1, 1, [Beam(1.0, (0.1, 0, 0), 0.01)])],
                output=OutputConfig(dir="/runs/a", diagnostics=[Diagnostic("line", "lineout", fields=("Ex",),
                                                                           every=5)]))
    base.update(kw)
    return SimConfig(**base)


@pytest.mark.parametrize("kw", [{"every": -1}, {"keep": 0}, {"every": True}, {"every": 2.5}])
def test_checkpoint_block_is_validated(kw):
    with pytest.raises(ValueError, match="checkpoint"):
        Checkpoint(**kw)


def test_checkpoint_needs_an_output_directory():
    cfg = SimConfig.from_dict({"checkpoint": {"every": 10}})
    assert cfg.checkpoint == Checkpoint(every=10, keep=2)
    with pytest.raises(ValueError, match="output.dir"):
        ck.check_config(cfg)
    ck.check_config(_cfg(checkpoint=Checkpoint(every=10)))


def test_allowed_changes_and_refused_ones():
    saved = ck.jsonable(_cfg(checkpoint=Checkpoint(every=10)).to_dict())
    same = _cfg(n_steps=5000, dims=(2, 1, 1), restart="/runs/a/checkpoint", diag_interval=3,
                checkpoint=Checkpoint(every=50, keep=4),
                output=OutputConfig(dir="/runs/b", overwrite=True,
                                    diagnostics=[Diagnostic("line", "lineout", fields=("Ex",), every=7)]))
    assert ck.config_differences(saved, same.to_dict()) == []
    seed = _cfg(seed=1)
    assert ck.config_differences(saved, seed.to_dict()) == ["seed"]
    dens = _cfg(species=[SpeciesConfig("e", -1, 1, [Beam(2.0, (0.1, 0, 0), 0.01)])])
    assert ck.config_differences(saved, dens.to_dict()) == ["species[0].beams[0].density"]
    fields = _cfg(output=OutputConfig(dir="/runs/a", diagnostics=[Diagnostic("line", "lineout", fields=("Ey",),
                                                                             every=5)]))
    assert ck.config_differences(saved, fields.to_dict()) == ["output.diagnostics[0].fields[0]"]
    seed.restart = "/runs/a/checkpoint/step-00000010"
    with pytest.raises(ValueError, match="seed"):
        ck.check_restart(seed, {"config": saved})


def _complete(root, step, world=1):
    p = os.path.join(root, ck.step_name(step))
    os.makedirs(p)
    ck.write_meta(p, {"format": ck.FORMAT, "layout": {"world": world}, "step": step})
    for r in range(world):
        open(ck.rank_file(p, r), "wb").close()
    return p


def test_partial_checkpoints_are_ignored(tmp_path):
    root = str(tmp_path / "checkpoint")
    p10 = _complete(root, 10)
    p20 = _complete(root, 20)
    os.rename(p20, p20 + ".tmp")  # written by every rank, but not renamed: a crash before the end
    os.makedirs(os.path.join(root, ck.step_name(30)))  # no meta.json
    assert ck.complete_steps(root) == [(10, p10)]
    assert ck.resolve(root) == p10
    assert ck.resolve(p10) == p10
    with pytest.raises(FileNotFoundError):
        ck.resolve(p20 + ".tmp")
    p40 = _complete(root, 40, world=2)
    os.remove(ck.rank_file(p40, 1))
    with pytest.raises(FileNotFoundError, match="ranks"):
        ck.read_meta(p40)


def test_prune_keeps_the_last_checkpoints(tmp_path):
    root = str(tmp_path)
    for s in (10, 20, 30):
        _complete(root, s)
    os.makedirs(os.path.join(root, ck.step_name(25) + ".tmp"))
    os.makedirs(os.path.join(root, ck.step_name(40) + ".tmp"))  # being written
    ck.prune(root, keep=2, current=30)
    assert sorted(os.listdir(root)) == [ck.step_name(20), ck.step_name(30), ck.step_name(40) + ".tmp"]


def test_meta_round_trip(tmp_path):
    meta = {"format": ck.FORMAT, "layout": {"world": 1, "dims": [1, 1, 1], "ghost": 2}, "step": 7,
            "time": 7 * 0.1, "history": {"t": [0.1 * 3, 1 / 3], "gauss": [float("nan"), 1e-300]},
            "config": ck.jsonable(_cfg().to_dict())}
    open(ck.rank_file(str(tmp_path), 0), "wb").close()
    ck.write_meta(str(tmp_path), meta)
    back = ck.read_meta(str(tmp_path))
    assert back["time"] == meta["time"] and back["history"]["t"] == meta["history"]["t"]
    assert math.isnan(back["history"]["gauss"][0]) and back["history"]["gauss"][1] == 1e-300
    assert back["config"] == meta["config"] and ck.config_differences(back["config"], _cfg().to_dict()) == []
    meta["format"] = "other"
    ck.write_meta(str(tmp_path), meta)
    with pytest.raises(ValueError, match="format"):
        ck.read_meta(str(tmp_path))


def test_series_is_cut_at_the_checkpoint_and_continued(tmp_path):
    """OUT-D6: rows at and after the checkpoint step go; the reopened writer appends after the rest."""
    from pypic.output.h5 import SeriesWriter, TrackWriter, truncate_series

    path = str(tmp_path / "s.h5")
    w = SeriesWriter(path, {"diagnostic": "s"}, "deflate", 1, max_rows=4)
    for s in range(10):
        w.append(s, 0.5 * s, {"v": np.full(3, float(s))})
    w.close()
    assert truncate_series(path, 6) == 6
    w = SeriesWriter(path, {"diagnostic": "s"}, "deflate", 1, max_rows=4, append=True)
    for s in range(6, 9):
        w.append(s, 0.5 * s, {"v": np.full(3, float(s))})
    w.close()
    import h5py

    with h5py.File(path, "r") as f:
        assert f["step"][...].tolist() == list(range(9))
        assert f["v"][:, 0].tolist() == [float(s) for s in range(9)] and f.attrs["diagnostic"] == "s"
    tp = str(tmp_path / "t.h5")
    t = TrackWriter(tp, {}, "deflate", 1, np.array([3, 5]))
    for s in range(4):
        t.append(s, s, {"x": np.array([s, -s], dtype=float)})
    t.close()
    assert truncate_series(tp, 2) == 2
    t = TrackWriter(tp, {}, "deflate", 1, np.array([3, 5]), append=True)
    t.append(2, 2.0, {"x": np.array([7.0, 8.0])})
    t.close()
    with h5py.File(tp, "r") as f:
        assert f["id"][...].tolist() == [3, 5] and f["step"][...].tolist() == [0, 1, 2]
        assert f["x"][2].tolist() == [7.0, 8.0]
