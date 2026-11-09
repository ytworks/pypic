// Common definitions for pypic CUDA kernels.
// The following macros are injected at compile time by pypic.backend:
//   NX NY NZ        local array extents (including ghost layers)
//   GH              ghost width
//   OX OY OZ        global coordinate of local array index 0 (node position)
//   DX DY DZ DT     cell sizes and time step
//   LX LY LZ        physical (global) domain lengths, domain is [0, L)
//   SX0 SX1 ...     physical bounds of the local sub-domain
//   BCX0 BCX1 ...   particle/field boundary type at global faces (0 periodic, 1 open, 2 reflecting)
//   DECX DECY DECZ  1 if the axis is decomposed over more than one rank
#pragma once

#define BC_PERIODIC 0
#define BC_OPEN 1
#define BC_REFLECT 2

#define IDX(i, j, k) ((((long long)(i)) * NY + (long long)(j)) * NZ + (long long)(k))

__device__ __forceinline__ double cic(double d) {
    d = fabs(d);
    return d < 1.0 ? 1.0 - d : 0.0;
}

// Trilinear interpolation of a staggered Yee component.
// (X,Y,Z) are local index coordinates (node i at integer i); (sx,sy,sz) are the
// staggering offsets (0 or 0.5) of the component.
__device__ __forceinline__ double interp(const double* __restrict__ F, double X, double Y, double Z,
                                         double sx, double sy, double sz) {
    double fx = X - sx, fy = Y - sy, fz = Z - sz;
    int i = (int)floor(fx), j = (int)floor(fy), k = (int)floor(fz);
    double wx = fx - i, wy = fy - j, wz = fz - k;
    double c00 = F[IDX(i, j, k)] * (1 - wx) + F[IDX(i + 1, j, k)] * wx;
    double c10 = F[IDX(i, j + 1, k)] * (1 - wx) + F[IDX(i + 1, j + 1, k)] * wx;
    double c01 = F[IDX(i, j, k + 1)] * (1 - wx) + F[IDX(i + 1, j, k + 1)] * wx;
    double c11 = F[IDX(i, j + 1, k + 1)] * (1 - wx) + F[IDX(i + 1, j + 1, k + 1)] * wx;
    double c0 = c00 * (1 - wy) + c10 * wy;
    double c1 = c01 * (1 - wy) + c11 * wy;
    return c0 * (1 - wz) + c1 * wz;
}

// Warp-level reduction then one atomic per warp.
__device__ __forceinline__ void warp_atomic_add(double* dst, double v) {
    unsigned mask = 0xffffffffu;
    for (int off = 16; off > 0; off >>= 1) v += __shfl_down_sync(mask, v, off);
    if ((threadIdx.x & 31) == 0 && v != 0.0) atomicAdd(dst, v);
}

// ---- counter based RNG (splitmix64) for decomposition-independent loading ----
__device__ __forceinline__ unsigned long long splitmix64(unsigned long long& s) {
    unsigned long long z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
__device__ __forceinline__ double urand(unsigned long long& s) {
    // uniform in (0, 1)
    return ((splitmix64(s) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
__device__ __forceinline__ unsigned long long hash_key(unsigned long long a, unsigned long long b,
                                                       unsigned long long c) {
    unsigned long long s = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL);
    splitmix64(s);
    s ^= c * 0xD1B54A32D192ED03ULL;
    splitmix64(s);
    return s;
}
