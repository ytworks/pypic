// Yee-lattice FDTD kernels (explicit leapfrog) with optional CPML.
// Layout (local index i,j,k):
//   Ex(i+1/2,j,k) Ey(i,j+1/2,k) Ez(i,j,k+1/2)
//   Bx(i,j+1/2,k+1/2) By(i+1/2,j,k+1/2) Bz(i+1/2,j+1/2,k)
//   Jx,Jy,Jz collocated with Ex,Ey,Ez; rho at nodes.
// Normalised units: c = eps0 = 1, so dE/dt = curl B - J, dB/dt = -curl E.
// Update range: I0X..I1X (exclusive) etc. are injected macros.

#define NIX (I1X - I0X)
#define NIY (I1Y - I0Y)
#define NIZ (I1Z - I0Z)

__device__ __forceinline__ bool range_index(long long t, int& i, int& j, int& k) {
    long long n = (long long)NIX * NIY * NIZ;
    if (t >= n) return false;
    k = (int)(t % NIZ) + I0Z;
    t /= NIZ;
    j = (int)(t % NIY) + I0Y;
    i = (int)(t / NIY) + I0X;
    return true;
}

extern "C" __global__ void update_B(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                    const double* __restrict__ Ez, double* Bx, double* By, double* Bz,
                                    const double dtb
#ifdef USE_PML
                                    , double* pbxy, double* pbxz, double* pbyx, double* pbyz,
                                    double* pbzx, double* pbzy,
                                    const double* __restrict__ bhx, const double* __restrict__ chx,
                                    const double* __restrict__ bhy, const double* __restrict__ chy,
                                    const double* __restrict__ bhz, const double* __restrict__ chz
#endif
) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    double dEz_dy = (Ez[IDX(i, j + 1, k)] - Ez[c]) * (1.0 / DY);
    double dEy_dz = (Ey[IDX(i, j, k + 1)] - Ey[c]) * (1.0 / DZ);
    double dEx_dz = (Ex[IDX(i, j, k + 1)] - Ex[c]) * (1.0 / DZ);
    double dEz_dx = (Ez[IDX(i + 1, j, k)] - Ez[c]) * (1.0 / DX);
    double dEy_dx = (Ey[IDX(i + 1, j, k)] - Ey[c]) * (1.0 / DX);
    double dEx_dy = (Ex[IDX(i, j + 1, k)] - Ex[c]) * (1.0 / DY);
#ifdef USE_PML
    // B-located derivatives sit at half-integer positions along the derivative axis.
    pbxy[c] = bhy[j] * pbxy[c] + chy[j] * dEz_dy;
    pbxz[c] = bhz[k] * pbxz[c] + chz[k] * dEy_dz;
    pbyz[c] = bhz[k] * pbyz[c] + chz[k] * dEx_dz;
    pbyx[c] = bhx[i] * pbyx[c] + chx[i] * dEz_dx;
    pbzx[c] = bhx[i] * pbzx[c] + chx[i] * dEy_dx;
    pbzy[c] = bhy[j] * pbzy[c] + chy[j] * dEx_dy;
    Bx[c] -= dtb * (dEz_dy - dEy_dz + pbxy[c] - pbxz[c]);
    By[c] -= dtb * (dEx_dz - dEz_dx + pbyz[c] - pbyx[c]);
    Bz[c] -= dtb * (dEy_dx - dEx_dy + pbzx[c] - pbzy[c]);
#else
    Bx[c] -= dtb * (dEz_dy - dEy_dz);
    By[c] -= dtb * (dEx_dz - dEz_dx);
    Bz[c] -= dtb * (dEy_dx - dEx_dy);
#endif
}

extern "C" __global__ void update_E(double* Ex, double* Ey, double* Ez, const double* __restrict__ Bx,
                                    const double* __restrict__ By, const double* __restrict__ Bz,
                                    const double* __restrict__ Jx, const double* __restrict__ Jy,
                                    const double* __restrict__ Jz, const double dte
#ifdef USE_PML
                                    , double* pexy, double* pexz, double* peyx, double* peyz,
                                    double* pezx, double* pezy,
                                    const double* __restrict__ bnx, const double* __restrict__ cnx,
                                    const double* __restrict__ bny, const double* __restrict__ cny,
                                    const double* __restrict__ bnz, const double* __restrict__ cnz
#endif
) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    double dBz_dy = (Bz[c] - Bz[IDX(i, j - 1, k)]) * (1.0 / DY);
    double dBy_dz = (By[c] - By[IDX(i, j, k - 1)]) * (1.0 / DZ);
    double dBx_dz = (Bx[c] - Bx[IDX(i, j, k - 1)]) * (1.0 / DZ);
    double dBz_dx = (Bz[c] - Bz[IDX(i - 1, j, k)]) * (1.0 / DX);
    double dBy_dx = (By[c] - By[IDX(i - 1, j, k)]) * (1.0 / DX);
    double dBx_dy = (Bx[c] - Bx[IDX(i, j - 1, k)]) * (1.0 / DY);
#ifdef USE_PML
    pexy[c] = bny[j] * pexy[c] + cny[j] * dBz_dy;
    pexz[c] = bnz[k] * pexz[c] + cnz[k] * dBy_dz;
    peyz[c] = bnz[k] * peyz[c] + cnz[k] * dBx_dz;
    peyx[c] = bnx[i] * peyx[c] + cnx[i] * dBz_dx;
    pezx[c] = bnx[i] * pezx[c] + cnx[i] * dBy_dx;
    pezy[c] = bny[j] * pezy[c] + cny[j] * dBx_dy;
    Ex[c] += dte * (dBz_dy - dBy_dz + pexy[c] - pexz[c] - Jx[c]);
    Ey[c] += dte * (dBx_dz - dBz_dx + peyz[c] - peyx[c] - Jy[c]);
    Ez[c] += dte * (dBy_dx - dBx_dy + pezx[c] - pezy[c] - Jz[c]);
#else
    Ex[c] += dte * (dBz_dy - dBy_dz - Jx[c]);
    Ey[c] += dte * (dBx_dz - dBz_dx - Jy[c]);
    Ez[c] += dte * (dBy_dx - dBx_dy - Jz[c]);
#endif
}

// out(B-located) = s * curl E   (used by the implicit Crank-Nicolson solver)
extern "C" __global__ void curl_E(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                  const double* __restrict__ Ez, double* Ox, double* Oy, double* Oz,
                                  const double s) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    Ox[c] = s * ((Ez[IDX(i, j + 1, k)] - Ez[c]) * (1.0 / DY) - (Ey[IDX(i, j, k + 1)] - Ey[c]) * (1.0 / DZ));
    Oy[c] = s * ((Ex[IDX(i, j, k + 1)] - Ex[c]) * (1.0 / DZ) - (Ez[IDX(i + 1, j, k)] - Ez[c]) * (1.0 / DX));
    Oz[c] = s * ((Ey[IDX(i + 1, j, k)] - Ey[c]) * (1.0 / DX) - (Ex[IDX(i, j + 1, k)] - Ex[c]) * (1.0 / DY));
}

// out(E-located) = s * curl B
extern "C" __global__ void curl_B(const double* __restrict__ Bx, const double* __restrict__ By,
                                  const double* __restrict__ Bz, double* Ox, double* Oy, double* Oz,
                                  const double s) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    Ox[c] = s * ((Bz[c] - Bz[IDX(i, j - 1, k)]) * (1.0 / DY) - (By[c] - By[IDX(i, j, k - 1)]) * (1.0 / DZ));
    Oy[c] = s * ((Bx[c] - Bx[IDX(i, j, k - 1)]) * (1.0 / DZ) - (Bz[c] - Bz[IDX(i - 1, j, k)]) * (1.0 / DX));
    Oz[c] = s * ((By[c] - By[IDX(i - 1, j, k)]) * (1.0 / DX) - (Bx[c] - Bx[IDX(i, j - 1, k)]) * (1.0 / DY));
}

// div E (at nodes) for Gauss-law diagnostics
extern "C" __global__ void div_E(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                 const double* __restrict__ Ez, double* out) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    out[c] = (Ex[c] - Ex[IDX(i - 1, j, k)]) * (1.0 / DX) + (Ey[c] - Ey[IDX(i, j - 1, k)]) * (1.0 / DY) +
             (Ez[c] - Ez[IDX(i, j, k - 1)]) * (1.0 / DZ);
}

// div B (at cell centres i+1/2,j+1/2,k+1/2)
extern "C" __global__ void div_B(const double* __restrict__ Bx, const double* __restrict__ By,
                                 const double* __restrict__ Bz, double* out) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    out[c] = (Bx[IDX(i + 1, j, k)] - Bx[c]) * (1.0 / DX) + (By[IDX(i, j + 1, k)] - By[c]) * (1.0 / DY) +
             (Bz[IDX(i, j, k + 1)] - Bz[c]) * (1.0 / DZ);
}

// ---- Binomial current filter (PIC-D10), one pass along one axis for up to 3 arrays.
// T_m = wl F_m[-1] + wc F_m[0] + wr F_m[+1] along `axis`, for local index lo <= i_axis < hi_m and
// every index of the other axes (ghosts included); copy_pass3 then writes T back to F. The
// explicit rounding intrinsics keep the order of the former array expression (no FMA
// contraction), so the result is bitwise identical to it.
__device__ __forceinline__ bool pass_index(long long t, int axis, int lo, int hi, long long& c, long long& s) {
    const int n[3] = {NX, NY, NZ};
    int m[3] = {NX, NY, NZ};
    m[axis] = hi - lo;
    if (m[axis] <= 0 || t >= (long long)m[0] * m[1] * m[2]) return false;
    int k = (int)(t % m[2]);
    t /= m[2];
    int j = (int)(t % m[1]);
    int i = (int)(t / m[1]);
    int ix[3] = {i, j, k};
    ix[axis] += lo;
    c = ((long long)ix[0] * n[1] + ix[1]) * n[2] + ix[2];
    s = axis == 0 ? (long long)n[1] * n[2] : (axis == 1 ? n[2] : 1);
    return true;
}

extern "C" __global__ void binomial_pass3(const double* __restrict__ F0, const double* __restrict__ F1,
                                          const double* __restrict__ F2, double* T0, double* T1, double* T2,
                                          const int narr, const int axis, const int lo, const int hi0,
                                          const int hi1, const int hi2, const double wl, const double wc,
                                          const double wr) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const double* F[3] = {F0, F1, F2};
    double* T[3] = {T0, T1, T2};
    const int hi[3] = {hi0, hi1, hi2};
    for (int m = 0; m < narr; m++) {
        long long c, s;
        if (!pass_index(t, axis, lo, hi[m], c, s)) continue;
        T[m][c] = __dadd_rn(__dadd_rn(__dmul_rn(wl, F[m][c - s]), __dmul_rn(wc, F[m][c])), __dmul_rn(wr, F[m][c + s]));
    }
}

extern "C" __global__ void copy_pass3(double* F0, double* F1, double* F2, const double* __restrict__ T0,
                                      const double* __restrict__ T1, const double* __restrict__ T2,
                                      const int narr, const int axis, const int lo, const int hi0, const int hi1,
                                      const int hi2) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    double* F[3] = {F0, F1, F2};
    const double* T[3] = {T0, T1, T2};
    const int hi[3] = {hi0, hi1, hi2};
    for (int m = 0; m < narr; m++) {
        long long c, s;
        if (pass_index(t, axis, lo, hi[m], c, s)) F[m][c] = T[m][c];
    }
}

// ---- Ghost layers (issue #8): up to 3 arrays per launch, `width` layers along `axis`, every
// index of the other axes. Thread t -> (array m, layer l, other indices); the layer layout is
// the C order of the slab F[..., lo:lo+width, ...], so a packed face buffer is (narr, slab).
__device__ __forceinline__ bool slab_index(long long t, int narr, int axis, int width, int& m, int& l,
                                           long long& off) {
    const int n[3] = {NX, NY, NZ};
    int s[3] = {NX, NY, NZ};
    s[axis] = width;
    long long per = (long long)s[0] * s[1] * s[2];
    if (t >= per * narr) return false;
    m = (int)(t / per);
    long long r = t - m * per;
    int k = (int)(r % s[2]);
    r /= s[2];
    int j = (int)(r % s[1]);
    int i = (int)(r / s[1]);
    int ix[3] = {i, j, k};
    l = ix[axis];
    ix[axis] = 0;
    off = ((long long)ix[0] * n[1] + ix[1]) * n[2] + ix[2];
    return true;
}

__device__ __forceinline__ long long axis_stride(int axis) {
    return axis == 0 ? (long long)NY * NZ : (axis == 1 ? (long long)NZ : 1LL);
}

// F[d0 + l] = F[s0 + l] (add = 0) or F[d0 + l] += F[s0 + l] (add = 1), l < width
extern "C" __global__ void slab_copy3(double* F0, double* F1, double* F2, const int narr, const int axis,
                                      const int s0, const int d0, const int width, const int add) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int m, l;
    long long off;
    if (!slab_index(t, narr, axis, width, m, l, off)) return;
    double* F = m == 0 ? F0 : (m == 1 ? F1 : F2);
    const long long st = axis_stride(axis);
    const double v = F[off + (s0 + l) * st];
    double* d = F + off + (d0 + l) * st;
    *d = add ? __dadd_rn(*d, v) : v;
}

// F[s0 + l] = 0, l < width
extern "C" __global__ void slab_zero3(double* F0, double* F1, double* F2, const int narr, const int axis,
                                      const int s0, const int width) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int m, l;
    long long off;
    if (!slab_index(t, narr, axis, width, m, l, off)) return;
    double* F = m == 0 ? F0 : (m == 1 ? F1 : F2);
    F[off + (s0 + l) * axis_stride(axis)] = 0.0;
}

// buf[m, slab] = F_m[s0 + l] (pack) or F_m[s0 + l] = / += buf[m, slab] (unpack, add)
extern "C" __global__ void slab_pack3(double* F0, double* F1, double* F2, double* buf, const int narr,
                                      const int axis, const int s0, const int width, const int unpack,
                                      const int add) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int m, l;
    long long off;
    if (!slab_index(t, narr, axis, width, m, l, off)) return;
    double* F = m == 0 ? F0 : (m == 1 ? F1 : F2);
    double* f = F + off + (s0 + l) * axis_stride(axis);
    if (!unpack) buf[t] = *f;
    else *f = add ? __dadd_rn(*f, buf[t]) : buf[t];
}
