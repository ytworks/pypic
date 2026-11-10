// Particle migration between ranks. The particles flagged by the push kernels are taken out of
// P into a particle-major buffer Q[j * nrow + row] (a slice of particles is one contiguous
// message), the holes they leave are filled from the tail, and for each decomposed axis the
// buffer is split into the particles for the low / high neighbour and those that stay.
// nrow is the number of rows of P (8, or more while the initial displacement is applied).

__device__ __forceinline__ int lower_bound(const int* __restrict__ a, const int n, const int v) {
    int lo = 0, hi = n;
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (a[mid] < v) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

// Q[j] = particle idx[j] for the k flagged particles (full: every particle, idx is not read).
extern "C" __global__ void mig_take(const double* __restrict__ P, const long long cap, const int nrow,
                                    const int* __restrict__ idx, const int k, const int full,
                                    double* __restrict__ Q) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= k) return;
    const int p = full ? j : idx[j];
    for (int r = 0; r < nrow; r++) Q[(long long)j * nrow + r] = P[r * cap + p];
}

// idx (ascending) holds the k flagged slots of n = nk + k particles. nh of them are below nk:
// these holes take the nh particles of the tail [nk, n) that are not flagged, both in ascending
// order, which leaves the nk remaining particles in [0, nk).
extern "C" __global__ void mig_fill(double* P, const long long cap, const int nrow, const int* __restrict__ idx,
                                    const int k, const int nk) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= k) return;
    const int t = nk + j;
    const int pos = lower_bound(idx, k, t);
    if (pos < k && idx[pos] == t) return;  // a flagged slot of the tail
    const int nh = lower_bound(idx, k, nk);
    const int hole = idx[j - (pos - nh)];  // pos - nh flagged slots of the tail precede t
    for (int r = 0; r < nrow; r++) P[r * cap + hole] = P[r * cap + t];
}

// Class of the m buffered particles along one axis: 0 to the low neighbour (x < x_lo), 1 to the
// high neighbour (x >= x_hi), 2 stays, 3 dead (w == 0). At a global non-periodic face there is
// no neighbour (has_lo / has_hi = 0): a particle exactly on the face stays with the edge rank.
// axis < 0: dead or stays. key[j] = class * m + j: the sorted keys group the classes and keep the
// order inside each. cnt[class * stride] counts the classes.
extern "C" __global__ void mig_classify(const double* __restrict__ Q, const int nrow, const int m, const int axis,
                                        const int has_lo, const int has_hi, const double x_lo, const double x_hi,
                                        long long* __restrict__ key, int* cnt, const int stride) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= m) return;
    const double* q = Q + (long long)j * nrow;
    int c = 2;
    if (q[6] == 0.0) c = 3;
    else if (axis >= 0) {
        const double x = q[axis];
        if (has_lo && x < x_lo) c = 0;
        else if (has_hi && x >= x_hi) c = 1;
    }
    key[j] = (long long)c * m + j;
    atomicAdd(cnt + c * stride, 1);
}

// The sorted keys give the particles in the order low, high, stays, dead: the first n_send go to
// send, the next n_keep to keep, the dead ones are dropped. wrap_L > 0 (periodic axis): a sent
// particle is wrapped into [0, L), as it crosses the global boundary when its rank is at the edge.
extern "C" __global__ void mig_pack(const double* __restrict__ Q, const int nrow, const int m, const int axis,
                                    const long long* __restrict__ key, const int n_send, const int n_keep,
                                    const double wrap_L, double* __restrict__ send, double* __restrict__ keep) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n_send + n_keep) return;
    const double* q = Q + (key[j] % m) * nrow;
    double* d = j < n_send ? send + (long long)j * nrow : keep + (long long)(j - n_send) * nrow;
    for (int r = 0; r < nrow; r++) d[r] = q[r];
    if (j < n_send && wrap_L > 0.0) d[axis] -= wrap_L * floor(d[axis] / wrap_L);
}

// Append the m buffered particles behind the n particles of P.
extern "C" __global__ void mig_append(double* __restrict__ P, const long long cap, const int nrow, const int n,
                                      const double* __restrict__ Q, const int m) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= m) return;
    for (int r = 0; r < nrow; r++) P[r * cap + n + j] = Q[(long long)j * nrow + r];
}
