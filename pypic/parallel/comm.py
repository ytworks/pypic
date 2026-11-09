"""Communication layer: single-process (SelfComm) and multi-GPU NCCL backend.

Point-to-point exchanges are issued as a list of (send_buffer, dest) and
(recv_buffer, src) pairs inside one NCCL group. Messages between the same pair
of ranks are matched in issue order, so callers must issue sends and receives in
a consistent order (see :mod:`pypic.parallel.decomp`).
"""
from __future__ import annotations

import os

import cupy as cp
import numpy as np


class SelfComm:
    rank = 0
    size = 1

    def exchange(self, sends, recvs):
        """sends: list of (array, dest); recvs: list of (array, src). Only self-messages."""
        assert len(sends) == len(recvs)
        for (s, d), (r, src) in zip(sends, recvs):
            assert d == 0 and src == 0
            r[...] = s

    def allreduce(self, arr: cp.ndarray, op: str = "sum") -> cp.ndarray:
        return arr

    def allreduce_scalar(self, v: float, op: str = "sum") -> float:
        return float(v)

    def bcast_obj(self, obj):
        return obj

    def barrier(self):
        cp.cuda.Device().synchronize()


def _nccl_dtype(dtype):
    from cupy.cuda import nccl

    return {
        np.dtype(np.float64): nccl.NCCL_FLOAT64,
        np.dtype(np.float32): nccl.NCCL_FLOAT32,
        np.dtype(np.int64): nccl.NCCL_INT64,
        np.dtype(np.int32): nccl.NCCL_INT32,
    }[np.dtype(dtype)]


class NcclComm:
    """NCCL communicator; one process per GPU. Created by :func:`init_from_env`."""

    def __init__(self, rank: int, size: int, uid):
        from cupy.cuda import nccl

        self.nccl = nccl
        self.rank = rank
        self.size = size
        self.comm = nccl.NcclCommunicator(size, uid, rank)

    @property
    def _stream(self):
        return cp.cuda.get_current_stream().ptr

    def exchange(self, sends, recvs):
        nccl = self.nccl
        local = []
        keep = []  # contiguous send copies must outlive groupEnd(), where NCCL actually launches
        nccl.groupStart()
        for s, d in sends:
            if d == self.rank:
                local.append(s)
                continue
            s = cp.ascontiguousarray(s)
            keep.append(s)
            if s.size:
                self.comm.send(s.data.ptr, s.size, _nccl_dtype(s.dtype), d, self._stream)
        tmp = []
        li = 0
        for r, src in recvs:
            if src == self.rank:
                r[...] = local[li]
                li += 1
                continue
            if r.size == 0:
                continue
            if r.flags.c_contiguous:
                self.comm.recv(r.data.ptr, r.size, _nccl_dtype(r.dtype), src, self._stream)
            else:
                buf = cp.empty(r.shape, r.dtype)
                self.comm.recv(buf.data.ptr, buf.size, _nccl_dtype(buf.dtype), src, self._stream)
                tmp.append((r, buf))
        nccl.groupEnd()
        for r, buf in tmp:
            r[...] = buf
        del keep

    def allreduce(self, arr: cp.ndarray, op: str = "sum") -> cp.ndarray:
        nccl = self.nccl
        o = {"sum": nccl.NCCL_SUM, "max": nccl.NCCL_MAX, "min": nccl.NCCL_MIN}[op]
        arr = cp.ascontiguousarray(arr)
        out = cp.empty_like(arr)
        self.comm.allReduce(arr.data.ptr, out.data.ptr, arr.size, _nccl_dtype(arr.dtype), o, self._stream)
        return out

    def allreduce_scalar(self, v: float, op: str = "sum") -> float:
        return float(self.allreduce(cp.asarray([v], dtype=cp.float64), op)[0])

    def bcast_obj(self, obj):
        # objects are constructed deterministically on every rank; nothing to do
        return obj

    def barrier(self):
        self.allreduce(cp.zeros(1))
        cp.cuda.Device().synchronize()


def init_from_env():
    """Create the communicator from environment variables set by ``pypic.launch``."""
    size = int(os.environ.get("PYPIC_WORLD_SIZE", "1"))
    rank = int(os.environ.get("PYPIC_RANK", "0"))
    if size == 1:
        dev = int(os.environ.get("PYPIC_DEVICE", "0"))
        cp.cuda.Device(dev).use()
        return SelfComm()
    dev = int(os.environ.get("PYPIC_DEVICE", str(rank)))
    cp.cuda.Device(dev).use()
    uid = bytes.fromhex(os.environ["PYPIC_NCCL_ID"])
    return NcclComm(rank, size, uid)


def encode_uid(uid) -> str:
    return bytes(uid).hex()
