"""CPU test of the multi-GPU launcher's failure handling (no GPU or NCCL needed)."""
from __future__ import annotations

import subprocess
import sys


def test_launcher_stops_all_ranks_when_any_rank_fails(tmp_path):
    """Rank 1 fails at start-up while rank 0 blocks (as in an NCCL wait). The launcher must stop
    rank 0 and return rank 1's exit code; without that it would hang, which the subprocess timeout
    (hang detection only) turns into a failure."""
    script = tmp_path / "fail.py"
    script.write_text("import os, sys, time\n"
                      "if os.environ['PYPIC_RANK'] == '1':\n    sys.exit(17)\n"
                      "while True:\n    time.sleep(1)\n")
    r = subprocess.run([sys.executable, "-m", "pypic.launch", "-n", "2", str(script)], capture_output=True,
                       text=True, timeout=300)
    assert r.returncode == 17, r.stderr
