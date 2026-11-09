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
        from cupy.cuda import nccl

        from pypic.parallel.comm import encode_uid

        env_base["PYPIC_NCCL_ID"] = encode_uid(nccl.get_unique_id())
    procs = []
    for r in range(a.nproc):
        env = dict(env_base)
        env["PYPIC_RANK"] = str(r)
        env["PYPIC_DEVICE"] = str(devices[r])
        procs.append(subprocess.Popen([sys.executable, a.script, *a.args], env=env))
    code = 0
    try:
        for p in procs:
            rc = p.wait()
            if rc != 0 and code == 0:
                code = rc
                for q in procs:
                    if q.poll() is None:
                        q.send_signal(signal.SIGTERM)
    except KeyboardInterrupt:
        for p in procs:
            p.send_signal(signal.SIGINT)
        code = 130
    sys.exit(code)


if __name__ == "__main__":
    main()
