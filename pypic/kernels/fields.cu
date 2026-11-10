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

// The update range and the ghost layer below it on every axis: the points a backward difference
// over the update range reads.
__device__ __forceinline__ bool range_index_low(long long t, int& i, int& j, int& k) {
    long long n = (long long)(NIX + 1) * (NIY + 1) * (NIZ + 1);
    if (t >= n) return false;
    k = (int)(t % (NIZ + 1)) + I0Z - 1;
    t /= NIZ + 1;
    j = (int)(t % (NIY + 1)) + I0Y - 1;
    i = (int)(t / (NIY + 1)) + I0X - 1;
    return true;
}

// The update range and one ghost layer around it on every axis: the points of B that the update of
// E, the gather and div B read. Computed there from the ghosts of E, B needs no exchange.
__device__ __forceinline__ bool range_index_wide(long long t, int& i, int& j, int& k) {
    long long n = (long long)(NIX + 2) * (NIY + 2) * (NIZ + 2);
    if (t >= n) return false;
    k = (int)(t % (NIZ + 2)) + I0Z - 1;
    t /= NIZ + 2;
    j = (int)(t % (NIY + 2)) + I0Y - 1;
    i = (int)(t / (NIY + 2)) + I0X - 1;
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
    if (!range_index_wide(t, i, j, k)) return;
    long long c = IDX(i, j, k);
    // Every product is rounded by itself (__dmul_rn, fma): a rank compiled with the absorbing layer
    // and one without give the same bits where the layer's terms are zero, which the ghost layer
    // computed here relies on.
    const double dEz_dy = __dmul_rn(Ez[IDX(i, j + 1, k)] - Ez[c], 1.0 / DY);
    const double dEy_dz = __dmul_rn(Ey[IDX(i, j, k + 1)] - Ey[c], 1.0 / DZ);
    const double dEx_dz = __dmul_rn(Ex[IDX(i, j, k + 1)] - Ex[c], 1.0 / DZ);
    const double dEz_dx = __dmul_rn(Ez[IDX(i + 1, j, k)] - Ez[c], 1.0 / DX);
    const double dEy_dx = __dmul_rn(Ey[IDX(i + 1, j, k)] - Ey[c], 1.0 / DX);
    const double dEx_dy = __dmul_rn(Ex[IDX(i, j + 1, k)] - Ex[c], 1.0 / DY);
#ifdef USE_PML
    // B-located derivatives sit at half-integer positions along the derivative axis.
    pbxy[c] = bhy[j] * pbxy[c] + chy[j] * dEz_dy;
    pbxz[c] = bhz[k] * pbxz[c] + chz[k] * dEy_dz;
    pbyz[c] = bhz[k] * pbyz[c] + chz[k] * dEx_dz;
    pbyx[c] = bhx[i] * pbyx[c] + chx[i] * dEz_dx;
    pbzx[c] = bhx[i] * pbzx[c] + chx[i] * dEy_dx;
    pbzy[c] = bhy[j] * pbzy[c] + chy[j] * dEx_dy;
    const double cx = dEz_dy - dEy_dz + pbxy[c] - pbxz[c];
    const double cy = dEx_dz - dEz_dx + pbyz[c] - pbyx[c];
    const double cz = dEy_dx - dEx_dy + pbzx[c] - pbzy[c];
#else
    const double cx = dEz_dy - dEy_dz, cy = dEx_dz - dEz_dx, cz = dEy_dx - dEx_dy;
#endif
    Bx[c] = fma(-dtb, cx, Bx[c]);
    By[c] = fma(-dtb, cy, By[c]);
    Bz[c] = fma(-dtb, cz, Bz[c]);
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
// low: also on the ghost layer below the update range, from the ghosts of E; curl_B of the result
// then needs no ghost exchange.
extern "C" __global__ void curl_E(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                  const double* __restrict__ Ez, double* Ox, double* Oy, double* Oz,
                                  const double s, const int low) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!(low ? range_index_low(t, i, j, k) : range_index(t, i, j, k))) return;
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

// ---- vector kernels of the implicit field solve: stacked (3, NX, NY, NZ) arrays
#define NVEC (3LL * NX * NY * NZ)

// Last kernel of the CG matrix product w = A u: w += Adiag u, and the inner products of the
// Chronopoulos-Gear recurrence with it. The sums of r u, w u and r r over the unknowns (W = 1) of
// each thread block go to part[3 * block]; the caller adds the blocks, a fixed order of summation.
// Launched with CG_BLOCK threads per block.
extern "C" __global__ void cg_matvec_end(double* __restrict__ w, const double* __restrict__ a,
                                         const double* __restrict__ u, const double* __restrict__ r,
                                         const double* __restrict__ W, double* __restrict__ part) {
    __shared__ double sum[3][CG_BLOCK];
    const int i = threadIdx.x;
    const long long t = (long long)blockIdx.x * blockDim.x + i;
    double ru = 0.0, wu = 0.0, rr = 0.0;
    if (t < NVEC) {
        const double wt = w[t] + a[t] * u[t];
        w[t] = wt;
        if (W[t] != 0.0) {  // ghost values are not defined: never multiplied, even by zero
            ru = r[t] * u[t];
            wu = wt * u[t];
            rr = r[t] * r[t];
        }
    }
    sum[0][i] = ru;
    sum[1][i] = wu;
    sum[2][i] = rr;
    __syncthreads();
    for (int h = CG_BLOCK / 2; h > 0; h >>= 1) {
        if (i < h)
            for (int q = 0; q < 3; q++) sum[q][i] += sum[q][i + h];
        __syncthreads();
    }
    if (i == 0)
        for (int q = 0; q < 3; q++) part[3LL * blockIdx.x + q] = sum[q][0];
}

// The CG matrix product in one launch, where no rank has a wall or an absorbing layer:
// w = a u + s curl curl u on the update range (curl_E on the range and the layer below it, curl_B and
// cg_matvec_end, with their expressions), and the inner products of cg_matvec_end. Outside the update
// range w is not written: those values are ghosts, which the CG never reads.
__device__ __forceinline__ double curl_curl(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                            const double* __restrict__ Ez, const int q, const int i, const int j,
                                            const int k) {
#define CURL_X(i, j, k) ((Ez[IDX(i, (j) + 1, k)] - Ez[IDX(i, j, k)]) * (1.0 / DY) - (Ey[IDX(i, j, (k) + 1)] - Ey[IDX(i, j, k)]) * (1.0 / DZ))
#define CURL_Y(i, j, k) ((Ex[IDX(i, j, (k) + 1)] - Ex[IDX(i, j, k)]) * (1.0 / DZ) - (Ez[IDX((i) + 1, j, k)] - Ez[IDX(i, j, k)]) * (1.0 / DX))
#define CURL_Z(i, j, k) ((Ey[IDX((i) + 1, j, k)] - Ey[IDX(i, j, k)]) * (1.0 / DX) - (Ex[IDX(i, (j) + 1, k)] - Ex[IDX(i, j, k)]) * (1.0 / DY))
    if (q == 0)
        return (CURL_Z(i, j, k) - CURL_Z(i, j - 1, k)) * (1.0 / DY) - (CURL_Y(i, j, k) - CURL_Y(i, j, k - 1)) * (1.0 / DZ);
    if (q == 1)
        return (CURL_X(i, j, k) - CURL_X(i, j, k - 1)) * (1.0 / DZ) - (CURL_Z(i, j, k) - CURL_Z(i - 1, j, k)) * (1.0 / DX);
    return (CURL_Y(i, j, k) - CURL_Y(i - 1, j, k)) * (1.0 / DX) - (CURL_X(i, j, k) - CURL_X(i, j - 1, k)) * (1.0 / DY);
#undef CURL_X
#undef CURL_Y
#undef CURL_Z
}

extern "C" __global__ void cg_matvec(const double* __restrict__ u, double* __restrict__ w,
                                     const double* __restrict__ a, const double* __restrict__ r,
                                     const double* __restrict__ W, double* __restrict__ part, const double s) {
    __shared__ double sum[3][CG_BLOCK];
    const int ti = threadIdx.x;
    const long long n = (long long)NX * NY * NZ;
    const long long t = (long long)blockIdx.x * blockDim.x + ti;
    double ru = 0.0, wu = 0.0, rr = 0.0;
    if (t < NVEC) {
        const int q = (int)(t / n);
        const long long c = t - q * n;
        const int k = (int)(c % NZ), j = (int)((c / NZ) % NY), i = (int)(c / ((long long)NY * NZ));
        if (i >= I0X && i < I1X && j >= I0Y && j < I1Y && k >= I0Z && k < I1Z) {
            const double wt = s * curl_curl(u, u + n, u + 2 * n, q, i, j, k) + a[t] * u[t];
            w[t] = wt;
            if (W[t] != 0.0) {
                ru = r[t] * u[t];
                wu = wt * u[t];
                rr = r[t] * r[t];
            }
        }
    }
    sum[0][ti] = ru;
    sum[1][ti] = wu;
    sum[2][ti] = rr;
    __syncthreads();
    for (int h = CG_BLOCK / 2; h > 0; h >>= 1) {
        if (ti < h)
            for (int v = 0; v < 3; v++) sum[v][ti] += sum[v][ti + h];
        __syncthreads();
    }
    if (ti == 0)
        for (int v = 0; v < 3; v++) part[3LL * blockIdx.x + v] = sum[v][0];
}

// Step lengths of one CG step from sc = ((r, u), (w, u), .) of all ranks, read and kept on the
// device: st = ((r, u) of the step before, alpha, beta), all zero before the first step.
//   beta = (r, u) / (r, u)_before,  alpha = (r, u) / ((w, u) - beta (r, u) / alpha_before)
// A zero denominator is a zero step: the residual is already exactly zero (the convergence test
// runs every few iterations only), or it is the first step (beta = 0).
extern "C" __global__ void cg_scalars(const double* __restrict__ sc, double* __restrict__ st) {
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    const double ru = sc[0];
    const double beta = st[0] != 0.0 ? ru / st[0] : 0.0;
    const double den = sc[1] - (beta != 0.0 && st[1] != 0.0 ? beta * ru / st[1] : 0.0);
    st[0] = ru;
    st[1] = den != 0.0 ? ru / den : 0.0;
    st[2] = beta;
}

// CG step with u = Minv r and w = A u given: p = u + beta p, s = w + beta s (= A p),
// x += alpha p, r -= alpha s, u = Minv r.
extern "C" __global__ void cg_step(double* __restrict__ x, double* __restrict__ r, double* __restrict__ u,
                                   double* __restrict__ p, double* __restrict__ s,
                                   const double* __restrict__ w, const double* __restrict__ Minv,
                                   const double* __restrict__ st) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= NVEC) return;
    const double alpha = st[1], beta = st[2];
    const double pt = u[t] + beta * p[t];
    const double sv = w[t] + beta * s[t];
    p[t] = pt;
    s[t] = sv;
    x[t] += alpha * pt;
    const double rt = r[t] - alpha * sv;
    r[t] = rt;
    u[t] = Minv[t] * rt;
}

// ---- half-step fields of the implicit scheme: Eh and Bh are stacked (3, NX, NY, NZ) arrays

// Eh = (E0 + E1) / 2 on every point, the ghosts included: the ghosts of E0 and of E1 already hold
// the values of their neighbours, so Eh needs no exchange of its own.
extern "C" __global__ void half_E(const double* __restrict__ Ex, const double* __restrict__ Ey,
                                  const double* __restrict__ Ez, const double* __restrict__ E1,
                                  double* __restrict__ Eh) {
    const long long n = (long long)NX * NY * NZ;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    Eh[t] = 0.5 * (Ex[t] + E1[t]);
    Eh[n + t] = 0.5 * (Ey[t] + E1[n + t]);
    Eh[2 * n + t] = 0.5 * (Ez[t] + E1[2 * n + t]);
}

// Bh = B0 + s curl Eh on the update range (s = -dt/2). The sum is not contracted with the product:
// the result is the one of curl_E followed by an addition.
extern "C" __global__ void half_B(const double* __restrict__ Eh, const double* __restrict__ Bx,
                                  const double* __restrict__ By, const double* __restrict__ Bz,
                                  double* __restrict__ Bh, const double s) {
    const long long n = (long long)NX * NY * NZ;
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    if (!range_index(t, i, j, k)) return;
    const double* Ex = Eh;
    const double* Ey = Eh + n;
    const double* Ez = Eh + 2 * n;
    long long c = IDX(i, j, k);
    Bh[c] = __dadd_rn(Bx[c], s * ((Ez[IDX(i, j + 1, k)] - Ez[c]) * (1.0 / DY) - (Ey[IDX(i, j, k + 1)] - Ey[c]) * (1.0 / DZ)));
    Bh[n + c] = __dadd_rn(By[c], s * ((Ex[IDX(i, j, k + 1)] - Ex[c]) * (1.0 / DZ) - (Ez[IDX(i + 1, j, k)] - Ez[c]) * (1.0 / DX)));
    Bh[2 * n + c] = __dadd_rn(Bz[c], s * ((Ey[IDX(i + 1, j, k)] - Ey[c]) * (1.0 / DX) - (Ex[IDX(i, j + 1, k)] - Ex[c]) * (1.0 / DY)));
}

// Change of the field over one Picard sweep: the largest |a - b| and the largest |a| over the
// unknowns (W = 1) of each thread block, to part[2 * block]; the caller takes the largest of the
// blocks. A NaN is kept, so that it can never pass the convergence test. Launched with CG_BLOCK
// threads per block.
extern "C" __global__ void picard_change(const double* __restrict__ a, const double* __restrict__ b,
                                         const double* __restrict__ W, double* __restrict__ part) {
    __shared__ double top[2][CG_BLOCK];
    const int i = threadIdx.x;
    const long long t = (long long)blockIdx.x * blockDim.x + i;
    double d = 0.0, m = 0.0;
    if (t < NVEC && W[t] != 0.0) {
        d = fabs(a[t] - b[t]);
        m = fabs(a[t]);
    }
    top[0][i] = d;
    top[1][i] = m;
    __syncthreads();
    for (int h = CG_BLOCK / 2; h > 0; h >>= 1) {
        if (i < h)
            for (int q = 0; q < 2; q++) {
                const double y = top[q][i + h];
                if (y > top[q][i] || y != y) top[q][i] = y;
            }
        __syncthreads();
    }
    if (i == 0)
        for (int q = 0; q < 2; q++) part[2LL * blockIdx.x + q] = top[q][0];
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

// buf[f, m, slab] = F_m[s_f + l] (pack) or F_m[s_f + l] = / += buf[f, m, slab] (unpack, add), for one
// face (s1 < 0) or two (s0, then s1) per launch. Two faces must not share a cell.
extern "C" __global__ void slab_pack3(double* F0, double* F1, double* F2, double* buf, const int narr,
                                      const int axis, const int s0, const int s1, const int width,
                                      const int unpack, const int add) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long face = (long long)NX * NY * NZ / (axis == 0 ? NX : (axis == 1 ? NY : NZ)) * width * narr;
    const bool second = t >= face;
    if (second && s1 < 0) return;
    int m, l;
    long long off;
    if (!slab_index(second ? t - face : t, narr, axis, width, m, l, off)) return;
    double* F = m == 0 ? F0 : (m == 1 ? F1 : F2);
    double* f = F + off + ((second ? s1 : s0) + l) * axis_stride(axis);
    if (!unpack) buf[t] = *f;
    else *f = add ? __dadd_rn(*f, buf[t]) : buf[t];
}

// Ghost layers of the axes that wrap onto this rank (bits of `mask`): every ghost cell takes the value of
// its periodic image, in one launch for all those axes. Thread t -> (region a, array, cell): region a has
// the 2 GH ghost layers of axis a, the owned range of the masked axes below a (their ghost layers are in
// those regions, so no cell is written twice) and every index of the other axes. An image is owned
// along every masked axis and is not written here.
extern "C" __global__ void slab_wrap3(double* F0, double* F1, double* F2, const int narr, const int mask) {
    long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const int n[3] = {NX, NY, NZ};
    for (int a = 0; a < 3; a++) {
        if (!(mask >> a & 1)) continue;
        int s[3], o[3];
        for (int b = 0; b < 3; b++) {
            const bool owned = b < a && (mask >> b & 1);
            s[b] = b == a ? 2 * GH : (owned ? n[b] - 2 * GH : n[b]);
            o[b] = owned ? GH : 0;
        }
        const long long per = (long long)s[0] * s[1] * s[2];
        if (r >= per * narr) {
            r -= per * narr;
            continue;
        }
        const int m = (int)(r / per);
        r -= m * per;
        int ix[3];
        ix[2] = (int)(r % s[2]);
        r /= s[2];
        ix[1] = (int)(r % s[1]);
        ix[0] = (int)(r / s[1]);
        long long dst = 0, src = 0;
        for (int b = 0; b < 3; b++) {
            const int ni = n[b] - 2 * GH;
            int d = ix[b] + o[b];
            if (b == a && d >= GH) d += ni;  // the upper GH layers of the region are the hi ghosts
            int q = d;
            if (mask >> b & 1) q = d < GH ? d + ni : (d >= GH + ni ? d - ni : d);
            dst = dst * n[b] + d;
            src = src * n[b] + q;
        }
        double* F = m == 0 ? F0 : (m == 1 ? F1 : F2);
        F[dst] = F[src];
        return;
    }
}
