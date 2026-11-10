"""GPU test of the launcher's CPU binding (needs the device's PCI address, no CUDA context)."""
from __future__ import annotations

import json
import os
import subprocess
import sys

import cupy as cp
import pytest

from pypic.launch import device_cpus


def test_launcher_binds_a_rank_to_the_numa_node_of_its_gpu(tmp_path):
    """A rank runs on the CPUs of its GPU's NUMA node, among those the launcher itself may use;
    --no-bind leaves the affinity as it is. The last device is used: on a multi-socket machine its
    node is not the one a process starts on."""
    dev = cp.cuda.runtime.getDeviceCount() - 1
    allowed = os.sched_getaffinity(0)
    node = device_cpus(dev)
    if node is None or not node & allowed:
        pytest.skip("no NUMA node known for the device, or none of its CPUs is allowed")
    script = tmp_path / "affinity.py"
    script.write_text("import json, os\nprint(json.dumps(sorted(os.sched_getaffinity(0))))\n")

    def affinity(*opts):
        r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", "1", "--devices", str(dev), *opts,
                            str(script)], capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr
        return set(json.loads(r.stdout))

    assert affinity() == node & allowed
    assert affinity("--no-bind") == allowed
