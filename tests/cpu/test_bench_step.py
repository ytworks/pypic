"""Phase accounting of the benchmark driver (scripts/bench_step.py --phases)."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "bench_step.py"


class Clock:
    t = 0.0

    def perf_counter(self):
        return self.t


def load(monkeypatch):
    spec = importlib.util.spec_from_file_location("bench_step", SCRIPT)
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    clock = Clock()
    monkeypatch.setattr(bench, "time", clock)
    return bench, clock


def test_phases_of_a_step_exclude_calls_made_outside_the_step(monkeypatch):
    bench, clock = load(monkeypatch)

    def write():
        clock.t += 1.0

    phases = bench.Phases(sync=lambda: None)
    output = phases.wrap("output", write)
    step = phases.wrap_step(output)
    phases.on = True
    step()    # Simulation.step: the output of time level n
    output()  # Simulation.run after its last step: the output of the final time level

    assert phases.step_time == 1.0
    assert phases.time["output"] == 1.0
    assert phases.calls["output"] == 1 and phases.series["output"] == [1.0]
    assert phases.outside["output"] == 1.0


def test_detail_of_a_phase_is_inclusive_and_exclusive_of_nested_details(monkeypatch):
    bench, clock = load(monkeypatch)
    phases = bench.Phases(sync=lambda: None)

    def gather():
        clock.t += 2.0

    gather = phases.wrap_detail("gather", gather)

    def compute():
        clock.t += 1.0
        gather()

    compute = phases.wrap_detail("compute", compute)
    step = phases.wrap_step(phases.wrap("output", compute))
    phases.on = True
    step()
    compute()  # outside a step: not part of the detail of a step

    assert phases.time["output"] == 3.0
    assert phases.detail["compute"] == [3.0, 1.0, 1]
    assert phases.detail["gather"] == [2.0, 2.0, 1]
