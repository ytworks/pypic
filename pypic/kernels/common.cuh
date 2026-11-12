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
//   MOVING_WINDOW   defined for a run with a moving window (PIC-D22); its value numbers the
//                   simulation, so that two of them in one process never share window_state
#pragma once

#define BC_PERIODIC 0
#define BC_OPEN 1
#define BC_REFLECT 2

#define IDX(i, j, k) ((((long long)(i)) * NY + (long long)(j)) * NZ + (long long)(k))
// index steps of one node along x and y (along z: 1)
#define STRIDE_I ((long long)NY * NZ)
#define STRIDE_J ((long long)NZ)

// Moving window (PIC-D22). Positions stay laboratory positions; the arrays hold the cells
// s .. s + n along x after s shifts. XIDX is the index coordinate along x in the local arrays;
// FLO / FHI are the positions of the lower and upper face of the box along axis a; RX0 / RX1 are
// the bounds of the rank along x. Without a window they are the compile-time constants. With a
// window they read window_state, which the kernel window_set writes after every shift:
// [0] the shift s in cells, [1] [2] the faces x- and x+, [3] [4] the bounds of the rank along x
// (the doubles that the migration compares with).
//
// XIDX with a window does not depend on the rank or on s but for a whole number of cells, bit for
// bit: the deposit of a path starts exactly where the deposit of the step before has ended, also
// when the particle has changed rank or the window has shifted in between. It is built in three
// steps. The product x * (1 / dx) is rounded by itself (__dmul_rn: the compiler would otherwise
// fuse it with the subtraction). XLOW, the cells that the arrays of the first rank hold below the
// box, is added: the sum is the same on every rank and is not negative for a position inside the
// arrays. Then an integer is taken away, the global cell OXI of local index 0 plus XLOW plus s;
// it is not negative on any rank, so this subtraction is exact. (The form (x - OX) / dx - s of a
// box at rest rounds x - OX with the origin of the rank. Two ranks then give one position index
// coordinates that differ by a constant part of the rounding unit of x, which grows with the
// place of the window, and a particle that changes rank moves its charge by that much without a
// current.) Not window-aware (refused with a window): the injection kernels and the bank, the
// body, implicit.cu.
#ifdef MOVING_WINDOW
__device__ double window_state[5];
#define XIDX(x) ((__dmul_rn((x), 1.0 / DX) + XLOW) - ((double)(OXI + XLOW) + window_state[0]))
#define FLO(a) ((a) == 0 ? window_state[1] : 0.0)
#define FHI(a) ((a) == 0 ? window_state[2] : Lg[a])
#define RX0 window_state[3]
#define RX1 window_state[4]
#else
#define XIDX(x) (((x) - OX) * (1.0 / DX))
#define FLO(a) 0.0
#define FHI(a) Lg[a]
#define RX0 SX0
#define RX1 SX1
#endif

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

// The six Yee components at one particle share two sets of trilinear weights per axis: node
// (offset 0) and half-cell (offset 0.5). Computing them once gives the same values as interp()
// with fewer FP64 operations (issue #6).
struct YeeWeights {
    int in[3], ih[3];     // floor of X - 0 and X - 0.5 per axis
    double wn[3], wh[3];  // the matching fractions
};

__device__ __forceinline__ void yee_weights(double X, double Y, double Z, YeeWeights& w) {
    const double P[3] = {X, Y, Z};
#pragma unroll
    for (int a = 0; a < 3; a++) {
        const double fn = P[a], fh = P[a] - 0.5;
        const double gn = floor(fn), gh = floor(fh);
        w.in[a] = (int)gn;
        w.ih[a] = (int)gh;
        w.wn[a] = fn - gn;
        w.wh[a] = fh - gh;
    }
}

// interp(F, X, Y, Z, sx, sy, sz) with hx/hy/hz = 1 where the offset is 0.5.
__device__ __forceinline__ double interp_w(const double* __restrict__ F, const YeeWeights& w, int hx, int hy,
                                           int hz) {
    int i = hx ? w.ih[0] : w.in[0], j = hy ? w.ih[1] : w.in[1], k = hz ? w.ih[2] : w.in[2];
    double wx = hx ? w.wh[0] : w.wn[0], wy = hy ? w.wh[1] : w.wn[1], wz = hz ? w.wh[2] : w.wn[2];
    const double* __restrict__ f = F + IDX(i, j, k);  // the 8 corners at constant offsets from one index
    double c00 = f[0] * (1 - wx) + f[STRIDE_I] * wx;
    double c10 = f[STRIDE_J] * (1 - wx) + f[STRIDE_I + STRIDE_J] * wx;
    double c01 = f[1] * (1 - wx) + f[STRIDE_I + 1] * wx;
    double c11 = f[STRIDE_J + 1] * (1 - wx) + f[STRIDE_I + STRIDE_J + 1] * wx;
    double c0 = c00 * (1 - wy) + c10 * wy;
    double c1 = c01 * (1 - wy) + c11 * wy;
    return c0 * (1 - wz) + c1 * wz;
}

#ifdef MIXED_GATHER
// PIC-D11: E and B gathered from float32 copies. The weights come from global index coordinates
// (x / dx), computed in FP64 and rounded once to float, so they do not depend on the rank layout.
struct YeeWeightsF {
    int in[3], ih[3];
    float wn[3], wh[3];
};

__device__ __forceinline__ void yee_weights_global(const double* x, YeeWeightsF& w) {
    const double inv[3] = {1.0 / DX, 1.0 / DY, 1.0 / DZ};
#ifdef MOVING_WINDOW
    const int off[3] = {OXI + (int)window_state[0], OYI, OZI};
#else
    const int off[3] = {OXI, OYI, OZI};
#endif
#pragma unroll
    for (int a = 0; a < 3; a++) {
        double fn = x[a] * inv[a], fh = fn - 0.5;
        double gn = floor(fn), gh = floor(fh);
        w.in[a] = (int)gn - off[a];
        w.ih[a] = (int)gh - off[a];
        w.wn[a] = (float)(fn - gn);
        w.wh[a] = (float)(fh - gh);
    }
}

__device__ __forceinline__ float interp_f(const float* __restrict__ F, const YeeWeightsF& w, int hx, int hy,
                                          int hz) {
    int i = hx ? w.ih[0] : w.in[0], j = hy ? w.ih[1] : w.in[1], k = hz ? w.ih[2] : w.in[2];
    float wx = hx ? w.wh[0] : w.wn[0], wy = hy ? w.wh[1] : w.wn[1], wz = hz ? w.wh[2] : w.wn[2];
    const float* __restrict__ f = F + IDX(i, j, k);
    float c00 = f[0] * (1.0f - wx) + f[STRIDE_I] * wx;
    float c10 = f[STRIDE_J] * (1.0f - wx) + f[STRIDE_I + STRIDE_J] * wx;
    float c01 = f[1] * (1.0f - wx) + f[STRIDE_I + 1] * wx;
    float c11 = f[STRIDE_J + 1] * (1.0f - wx) + f[STRIDE_I + STRIDE_J + 1] * wx;
    float c0 = c00 * (1.0f - wy) + c10 * wy;
    float c1 = c01 * (1.0f - wy) + c11 * wy;
    return c0 * (1.0f - wz) + c1 * wz;
}
typedef float gather_t;
#else
typedef double gather_t;
#endif

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
