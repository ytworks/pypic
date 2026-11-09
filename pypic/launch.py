"""Multi-GPU launcher: one process per GPU, NCCL for communication.

Usage::

    python -m pypic.launch -n 8 examples/two_stream/run.py --config ...

Each child process gets PYPIC_RANK / PYPIC_WORLD_SIZE / PYPIC_DEVICE / PYPIC_NCCL_ID.
"""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time


def _terminate(procs, grace: float = 10.0):
    """SIGTERM the ranks that are still running, then SIGKILL whatever survives ``grace`` s."""
    for p in procs:
        if p.poll() is None:
            p.send_signal(signal.SIGTERM)
    t_end = time.monotonic() + grace
    for p in procs:
        try:
            p.wait(timeout=max(0.0, t_end - time.monotonic()))
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Launch a pypic script on several GPUs")
    ap.add_argument("-n", "--nproc", type=int, default=1, help="number of GPUs / processes")
    ap.add_argument("--devices", type=str, default=None, help="comma separated device ids (default 0..n-1)")
    ap.add_argument("script")
    ap.add_argument("args", nargs=argparse.REMAINDER)
    a = ap.parse_args(argv)

    devices = [int(x) for x in a.devices.split(",")] if a.devices else list(range(a.nproc))
    if len(devices) < a.nproc:
        ap.error("not enough devices")
    env_base = dict(os.environ)
    env_base["PYPIC_WORLD_SIZE"] = str(a.nproc)
    env_base.setdefault("PYTHONUNBUFFERED", "1")  # interleaved, immediate output from all ranks
    if a.nproc > 1:
        try:
            from cupy.cuda import nccl
        except ImportError:  # CPU-only environment: ranks that need NCCL will fail at init
            print("[pypic.launch] NCCL is not available; PYPIC_NCCL_ID is not set", file=sys.stderr)
        else:
            env_base["PYPIC_NCCL_ID"] = bytes(nccl.get_unique_id()).hex()
    procs = []
    for r in range(a.nproc):
        env = dict(env_base)
        env["PYPIC_RANK"] = str(r)
        env["PYPIC_DEVICE"] = str(devices[r])
        procs.append(subprocess.Popen([sys.executable, a.script, *a.args], env=env))
    code = 0
    try:
        # Poll every rank: a rank that dies (init error, OOM, ...) while the others block in NCCL
        # must be detected at once, whatever its rank number.
        while True:
            states = [p.poll() for p in procs]
            failed = [(r, rc) for r, rc in enumerate(states) if rc not in (None, 0)]
            if failed:
                r, rc = failed[0]
                code = rc if rc > 0 else 128 - rc  # killed by a signal -> 128 + signal number
                print(f"[pypic.launch] rank {r} exited with code {code}; terminating the other ranks",
                      file=sys.stderr, flush=True)
                _terminate(procs)
                break
            if all(rc == 0 for rc in states):
                break
            time.sleep(0.1)
    except KeyboardInterrupt:
        for p in procs:
            if p.poll() is None:
                p.send_signal(signal.SIGINT)
        _terminate(procs, grace=5.0)
        code = 130
    sys.exit(code)


if __name__ == "__main__":
    main()
