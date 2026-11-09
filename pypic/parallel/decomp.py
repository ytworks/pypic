"""3-D Cartesian domain decomposition."""
from __future__ import annotations

from dataclasses import dataclass, field

BC_CODES = {"periodic": 0, "open": 1, "reflecting": 2}


def auto_dims(size: int, n_cells) -> tuple[int, int, int]:
    """Factor ``size`` into (px,py,pz) dividing the grid and minimising the halo surface."""
    best, best_cost = None, None
    for px in range(1, size + 1):
        if size % px:
            continue
        for py in range(1, size // px + 1):
            if (size // px) % py:
                continue
            pz = size // px // py
            dims = (px, py, pz)
            if any(n % p for n, p in zip(n_cells, dims)):
                continue
            loc = [n // p for n, p in zip(n_cells, dims)]
            cost = loc[0] * loc[1] * (pz > 1) + loc[0] * loc[2] * (py > 1) + loc[1] * loc[2] * (px > 1)
            if best_cost is None or cost < best_cost:
                best, best_cost = dims, cost
    if best is None:
        raise ValueError(f"cannot decompose grid {tuple(n_cells)} over {size} ranks")
    return best


@dataclass
class AxisLayout:
    n_phys: int  # physical cells owned by this rank
    phys_lo: int  # global index of first physical cell
    pml_lo: int  # absorbing cells appended at the low side (only on boundary ranks)
    pml_hi: int
    lo: int | None  # neighbour rank (None: global non-periodic face)
    hi: int | None
    coord: int
    nproc: int

    @property
    def n_int(self):
        return self.n_phys + self.pml_lo + self.pml_hi

    @property
    def cell_lo(self):
        return self.phys_lo - self.pml_lo


@dataclass
class Decomposition:
    n_cells: tuple
    dims: tuple
    rank: int
    bc: tuple  # ((lo,hi),(lo,hi),(lo,hi)) as strings
    npml: int
    axes: list = field(default_factory=list)

    @classmethod
    def build(cls, n_cells, dims, rank, bc, npml):
        n_cells = tuple(int(v) for v in n_cells)
        dims = tuple(int(v) for v in dims)
        px, py, pz = dims
        coords = (rank // (py * pz), (rank // pz) % py, rank % pz)

        def rank_of(c):
            return (c[0] * py + c[1]) * pz + c[2]

        axes = []
        for a in range(3):
            lo_bc, hi_bc = bc[a]
            if (lo_bc == "periodic") != (hi_bc == "periodic"):
                raise ValueError("periodic boundaries must be set on both faces of an axis")
            if n_cells[a] % dims[a]:
                raise ValueError(f"axis {a}: {n_cells[a]} cells not divisible by {dims[a]} ranks")
            n_phys = n_cells[a] // dims[a]
            c = coords[a]
            per = lo_bc == "periodic"

            def nb(cc, a=a, per=per):
                if 0 <= cc < dims[a]:
                    pass
                elif per:
                    cc %= dims[a]
                else:
                    return None
                cl = list(coords)
                cl[a] = cc
                return rank_of(cl)

            pml_lo = npml if (c == 0 and lo_bc == "open") else 0
            pml_hi = npml if (c == dims[a] - 1 and hi_bc == "open") else 0
            axes.append(AxisLayout(n_phys, c * n_phys, pml_lo, pml_hi, nb(c - 1), nb(c + 1), c, dims[a]))
        return cls(n_cells, dims, rank, tuple(tuple(b) for b in bc), npml, axes)
