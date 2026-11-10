"""Multi-GPU launcher: one process per GPU, NCCL for communication.

Usage::

    python -m pypic.launch -n 8 examples/two_stream/run.py --config ...

Each child process gets PYPIC_RANK / PYPIC_WORLD_SIZE / PYPIC_DEVICE / PYPIC_NCCL_ID and runs on
the CPUs of the NUMA node its GPU is attached to (``--no-bind`` leaves the CPU affinity alone).
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


def parse_cpulist(text: str) -> set[int]:
    """CPU numbers of a Linux cpulist such as ``0-15,64-79``."""
    cpus: set[int] = set()
    for part in text.strip().split(","):
        if part:
            lo, _, hi = part.partition("-")
            cpus.update(range(int(lo), int(hi or lo) + 1))
    return cpus


def device_cpus(device: int) -> set[int] | None:
    """CPUs of the NUMA node a CUDA device is attached to, or None when that is not known (no
    CuPy, not Linux, no NUMA information for the device)."""
    try:
        from cupy.cuda import runtime

        bdf = runtime.deviceGetPCIBusId(device).lower()
        with open(f"/sys/bus/pci/devices/{bdf}/numa_node") as f:
            node = int(f.read())
        if node < 0:
            return None
        with open(f"/sys/devices/system/node/node{node}/cpulist") as f:
            return parse_cpulist(f.read())
    except Exception:  # binding is an optimisation: any failure means "do not bind"
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="Launch a pypic script on several GPUs")
    ap.add_argument("-n", "--nproc", type=int, default=1, help="number of GPUs / processes")
    ap.add_argument("--devices", type=str, default=None, help="comma separated device ids (default 0..n-1)")
    ap.add_argument("--no-bind", action="store_true",
                    help="do not bind each rank to the CPUs of its GPU's NUMA node")
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
    # A rank drives its GPU from one host thread; on the NUMA node of the GPU that thread and the
    # driver's memory stay next to the device. A child inherits the affinity set here.
    allowed = os.sched_getaffinity(0) if hasattr(os, "sched_setaffinity") and not a.no_bind else None
    procs = []
    for r in range(a.nproc):
        env = dict(env_base)
        env["PYPIC_RANK"] = str(r)
        env["PYPIC_DEVICE"] = str(devices[r])
        cpus = (device_cpus(devices[r]) or set()) & allowed if allowed else None
        if cpus:
            os.sched_setaffinity(0, cpus)
        try:
            procs.append(subprocess.Popen([sys.executable, a.script, *a.args], env=env))
        finally:
            if cpus:
                os.sched_setaffinity(0, allowed)
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
