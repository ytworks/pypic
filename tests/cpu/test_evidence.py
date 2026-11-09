"""The evidence scripts (scripts/record_example.py, verify_bench_results.py) run on a throwaway
repository with a stand-in example, so their command-line behaviour is tested without touching
bench/results."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"


def _git(root: Path, *args: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", *args], cwd=root,
                          check=True, capture_output=True, text=True).stdout.strip()


def _repo(tmp_path: Path, analyze: str) -> Path:
    root = tmp_path / "repo"
    (root / "scripts").mkdir(parents=True)
    for name in ("record_example.py", "bench_records.py", "verify_bench_results.py"):
        shutil.copy(SCRIPTS / name, root / "scripts" / name)
    (root / "examples" / "demo").mkdir(parents=True)
    (root / "examples" / "demo" / "analyze.py").write_text(analyze)
    _git(root, "init", "-q")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "A")
    return root


def _run_dir(tmp_path: Path, name: str, run_id: str, value: float) -> Path:
    d = tmp_path / name
    d.mkdir()
    (d / "run.json").write_text(json.dumps({
        "git_commit": "0123456789abcdef", "run_id": run_id, "created": "2026-11-09T12:00:00+09:00",
        "world_size": 8, "gpu_name": "GPU", "meta": {"example": "demo", "config_file": "configs/big.yaml"}}))
    (d / "value.json").write_text(json.dumps(value))
    return d


def _record(root: Path, run_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(root / "scripts" / "record_example.py"), str(run_dir)], cwd=root,
                          capture_output=True, text=True)


def _verify(root: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(root / "scripts" / "verify_bench_results.py")], cwd=root,
                          capture_output=True, text=True)


def test_each_run_keeps_its_own_strict_json_record(tmp_path):
    """Review of PR #2: two runs of one configuration recorded on one day overwrote each other, and
    an unmeasured metric was written as NaN, which is not JSON."""
    root = _repo(tmp_path, "import json, os\n\n\ndef analyze(run_dir, show_plot=False):\n"
                           "    v = json.load(open(os.path.join(run_dir, 'value.json')))\n"
                           "    return {'passed': True, 'value': v, 'tail_index': float('nan')}\n")
    runs = [_run_dir(tmp_path, "r1", "00aa11bb", 1.0), _run_dir(tmp_path, "r2", "22cc33dd", 2.0)]
    for run in runs + runs[:1]:  # recording a run again replaces only its own record
        r = _record(root, run)
        assert r.returncode == 0, r.stderr
    res = root / "bench" / "results"
    files = sorted(res.glob("*-example-*.json"))
    assert len(files) == 2 and len(json.loads((res / "index.json").read_text())["entries"]) == 2
    recs = [json.loads(f.read_text(), parse_constant=_no_constant) for f in files]  # strict JSON
    assert sorted(r["metrics"]["value"] for r in recs) == [1.0, 2.0]
    for rec, f in zip(recs, files):
        assert rec["metrics"]["tail_index"] is None and f.name.endswith(f"-run-{rec['run_id']}.json")
    assert _verify(root).returncode == 0
    # the verifier rejects a record that is not strict JSON
    f = files[0]
    f.write_text(f.read_text().replace('"tail_index": null', '"tail_index": NaN'))
    r = _verify(root)
    assert r.returncode == 1 and "not strict JSON" in r.stderr


def test_the_analysis_code_must_not_change_while_it_runs(tmp_path):
    """Review of PR #2: the analysis commit was read after the analysis, so a checkout made while it
    ran was recorded as the code that produced the metrics."""
    root = _repo(tmp_path, "import subprocess\n\n\ndef analyze(run_dir, show_plot=False):\n"
                           "    subprocess.run(['git', 'checkout', '-q', 'other'], check=True)\n"
                           "    return {'passed': True, 'value': 1.0}\n")
    _git(root, "checkout", "-q", "-b", "other")
    (root / "examples" / "demo" / "analyze.py").write_text("def analyze(run_dir, show_plot=False):\n    return {}\n")
    _git(root, "commit", "-q", "-am", "B")
    _git(root, "checkout", "-q", "-")
    r = _record(root, _run_dir(tmp_path, "r1", "00aa11bb", 1.0))
    assert r.returncode == 1 and "changed while the analysis ran" in r.stderr
    assert not list((root / "bench" / "results").glob("*.json"))


def _no_constant(token):
    raise ValueError(token)
