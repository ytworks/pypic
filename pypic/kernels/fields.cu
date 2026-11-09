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
