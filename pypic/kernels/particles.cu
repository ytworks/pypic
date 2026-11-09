// Explicit particle kernels: relativistic Buneman-Boris push fused with
// Esirkepov charge-conserving current deposition (CIC = cubic super-particles).
// Particle storage: P[row * cap + p], rows = x, y, z, ux, uy, uz, w  (u = gamma v, c = 1).

#define PX(p) P[(p)]
#define PY(p) P[cap + (p)]
#define PZ(p) P[2 * cap + (p)]
#define PUX(p) P[3 * cap + (p)]
#define PUY(p) P[4 * cap + (p)]
#define PUZ(p) P[5 * cap + (p)]
#define PW(p) P[6 * cap + (p)]

// Esirkepov deposition of a straight segment (X0,Y0,Z0)->(X1,Y1,Z1) in local index
// coordinates. |X1-X0| < 1 is required. qw = charge * weight.
__device__ void esirkepov(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double qw,
                          double* Jx, double* Jy, double* Jz) {
    int i0 = (int)floor(X0) - 1, j0 = (int)floor(Y0) - 1, k0 = (int)floor(Z0) - 1;
    double S0x[4], S0y[4], S0z[4], dSx[4], dSy[4], dSz[4];
#pragma unroll
    for (int a = 0; a < 4; a++) {
        S0x[a] = cic(X0 - (i0 + a));
        dSx[a] = cic(X1 - (i0 + a)) - S0x[a];
        S0y[a] = cic(Y0 - (j0 + a));
        dSy[a] = cic(Y1 - (j0 + a)) - S0y[a];
        S0z[a] = cic(Z0 - (k0 + a));
        dSz[a] = cic(Z1 - (k0 + a)) - S0z[a];
    }
    const double cx = qw * (1.0 / (DY * DZ * DT));
    const double cy = qw * (1.0 / (DX * DZ * DT));
    const double cz = qw * (1.0 / (DX * DY * DT));
    const double third = 1.0 / 3.0;
    // Jx
    for (int b = 0; b < 4; b++) {
        if (S0y[b] == 0.0 && dSy[b] == 0.0) continue;
        for (int c = 0; c < 4; c++) {
            if (S0z[c] == 0.0 && dSz[c] == 0.0) continue;
            double wyz = S0y[b] * S0z[c] + 0.5 * dSy[b] * S0z[c] + 0.5 * S0y[b] * dSz[c] + third * dSy[b] * dSz[c];
            double acc = 0.0;
            for (int a = 0; a < 3; a++) {
                acc -= cx * dSx[a] * wyz;
                if (acc != 0.0) atomicAdd(&Jx[IDX(i0 + a, j0 + b, k0 + c)], acc);
            }
        }
    }
    // Jy
    for (int a = 0; a < 4; a++) {
        if (S0x[a] == 0.0 && dSx[a] == 0.0) continue;
        for (int c = 0; c < 4; c++) {
            if (S0z[c] == 0.0 && dSz[c] == 0.0) continue;
            double wxz = S0x[a] * S0z[c] + 0.5 * dSx[a] * S0z[c] + 0.5 * S0x[a] * dSz[c] + third * dSx[a] * dSz[c];
            double acc = 0.0;
            for (int b = 0; b < 3; b++) {
                acc -= cy * dSy[b] * wxz;
                if (acc != 0.0) atomicAdd(&Jy[IDX(i0 + a, j0 + b, k0 + c)], acc);
            }
        }
    }
    // Jz
    for (int a = 0; a < 4; a++) {
        if (S0x[a] == 0.0 && dSx[a] == 0.0) continue;
        for (int b = 0; b < 4; b++) {
            if (S0y[b] == 0.0 && dSy[b] == 0.0) continue;
            double wxy = S0x[a] * S0y[b] + 0.5 * dSx[a] * S0y[b] + 0.5 * S0x[a] * dSy[b] + third * dSx[a] * dSy[b];
            double acc = 0.0;
            for (int c = 0; c < 3; c++) {
                acc -= cz * dSz[c] * wxy;
                if (acc != 0.0) atomicAdd(&Jz[IDX(i0 + a, j0 + b, k0 + c)], acc);
            }
        }
    }
}

__device__ __forceinline__ void deposit_global(const double* p0, const double* p1, double qw, double* Jx,
                                               double* Jy, double* Jz) {
    esirkepov((p0[0] - OX) * (1.0 / DX), (p0[1] - OY) * (1.0 / DY), (p0[2] - OZ) * (1.0 / DZ),
              (p1[0] - OX) * (1.0 / DX), (p1[1] - OY) * (1.0 / DY), (p1[2] - OZ) * (1.0 / DZ), qw, Jx, Jy, Jz);
}

// diag: [0] kinetic energy (time centred), [1..3] momentum, [4] kinetic energy lost
// through open boundaries, [5] charge lost through open boundaries.
// flag[0]: number of particles that must be removed or migrated.
extern "C" __global__ void push_deposit(double* P, const int n, const int cap, const double* __restrict__ Ex,
                                        const double* __restrict__ Ey, const double* __restrict__ Ez,
                                        const double* __restrict__ Bx, const double* __restrict__ By,
                                        const double* __restrict__ Bz, double* Jx, double* Jy, double* Jz,
                                        const double q, const double m, double* diag, int* flag, int* list, const int listcap,
                                        const int do_diag) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0.0, mx = 0.0, my = 0.0, mz = 0.0, lke = 0.0, lq = 0.0;
    if (p < n) {
        const double w = PW(p);
        double pos0[3] = {PX(p), PY(p), PZ(p)};
        double u[3] = {PUX(p), PUY(p), PUZ(p)};
        const double X = (pos0[0] - OX) * (1.0 / DX), Y = (pos0[1] - OY) * (1.0 / DY),
                     Z = (pos0[2] - OZ) * (1.0 / DZ);
        const double h = q * DT / (2.0 * m);
        double ex = h * interp(Ex, X, Y, Z, 0.5, 0.0, 0.0);
        double ey = h * interp(Ey, X, Y, Z, 0.0, 0.5, 0.0);
        double ez = h * interp(Ez, X, Y, Z, 0.0, 0.0, 0.5);
        double bx = h * interp(Bx, X, Y, Z, 0.0, 0.5, 0.5);
        double by = h * interp(By, X, Y, Z, 0.5, 0.0, 0.5);
        double bz = h * interp(Bz, X, Y, Z, 0.5, 0.5, 0.0);
        double g0 = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        // --- relativistic Buneman-Boris ---
        double umx = u[0] + ex, umy = u[1] + ey, umz = u[2] + ez;  // half electric kick
        double gm = sqrt(1.0 + umx * umx + umy * umy + umz * umz);
        double tx = bx / gm, ty = by / gm, tz = bz / gm;
        double t2 = tx * tx + ty * ty + tz * tz;
        double sx = 2.0 * tx / (1.0 + t2), sy = 2.0 * ty / (1.0 + t2), sz = 2.0 * tz / (1.0 + t2);
        double upx = umx + (umy * tz - umz * ty);
        double upy = umy + (umz * tx - umx * tz);
        double upz = umz + (umx * ty - umy * tx);
        double urx = umx + (upy * sz - upz * sy);  // exact rotation: |u+| == |u-|
        double ury = umy + (upz * sx - upx * sz);
        double urz = umz + (upx * sy - upy * sx);
        u[0] = urx + ex;
        u[1] = ury + ey;
        u[2] = urz + ez;  // second half electric kick
        double g1 = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        double pos1[3] = {pos0[0] + DT * u[0] / g1, pos0[1] + DT * u[1] / g1, pos0[2] + DT * u[2] / g1};
        if (do_diag) {
            ke = w * m * (0.5 * (g0 + g1) - 1.0);
            mx = 0.5 * w * m * (PUX(p) + u[0]);
            my = 0.5 * w * m * (PUY(p) + u[1]);
            mz = 0.5 * w * m * (PUZ(p) + u[2]);
        }
        // --- reflecting walls: split the path at the wall so the deposited current
        //     follows the real (folded) trajectory and charge stays conserved ---
        const double Lg[3] = {LX, LY, LZ};
        const int b0[3] = {BCX0, BCY0, BCZ0};
        const int b1[3] = {BCX1, BCY1, BCZ1};
        const double qw = q * w;
        for (int it = 0; it < 3; it++) {
            double tmin = 2.0, wall = 0.0;
            int amin = -1;
            for (int a = 0; a < 3; a++) {
                double d = pos1[a] - pos0[a];
                if (b0[a] == BC_REFLECT && pos1[a] < 0.0 && d < 0.0) {
                    double tt = (0.0 - pos0[a]) / d;
                    if (tt < tmin) { tmin = tt; amin = a; wall = 0.0; }
                }
                if (b1[a] == BC_REFLECT && pos1[a] > Lg[a] && d > 0.0) {
                    double tt = (Lg[a] - pos0[a]) / d;
                    if (tt < tmin) { tmin = tt; amin = a; wall = Lg[a]; }
                }
            }
            if (amin < 0) break;
            tmin = fmin(fmax(tmin, 0.0), 1.0);
            double pw[3];
            for (int a = 0; a < 3; a++) pw[a] = pos0[a] + tmin * (pos1[a] - pos0[a]);
            pw[amin] = wall;
            deposit_global(pos0, pw, qw, Jx, Jy, Jz);
            pos1[amin] = 2.0 * wall - pos1[amin];
            u[amin] = -u[amin];
            for (int a = 0; a < 3; a++) pos0[a] = pw[a];
        }
        deposit_global(pos0, pos1, qw, Jx, Jy, Jz);

        // --- open boundaries: particle leaves and is removed ---
        bool dead = false;
        for (int a = 0; a < 3; a++) {
            if ((b0[a] == BC_OPEN && pos1[a] < 0.0) || (b1[a] == BC_OPEN && pos1[a] >= Lg[a])) dead = true;
        }
        if (dead) {
            if (do_diag) {
                ke = 0.0; mx = my = mz = 0.0;
            }
            lke = w * m * (g1 - 1.0);
            lq = qw;
            PW(p) = 0.0;
            { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        } else {
            // periodic wrap on non-decomposed axes; decomposed axes are wrapped after migration
            const int dec[3] = {DECX, DECY, DECZ};
            const double slo[3] = {SX0, SY0, SZ0};
            const double shi[3] = {SX1, SY1, SZ1};
            bool migrate = false;
            for (int a = 0; a < 3; a++) {
                if (!dec[a]) {
                    if (b0[a] == BC_PERIODIC) pos1[a] -= Lg[a] * floor(pos1[a] / Lg[a]);
                } else if (pos1[a] < slo[a] || pos1[a] >= shi[a]) {
                    migrate = true;
                }
            }
            if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        }
        PX(p) = pos1[0];
        PY(p) = pos1[1];
        PZ(p) = pos1[2];
        PUX(p) = u[0];
        PUY(p) = u[1];
        PUZ(p) = u[2];
    }
    if (do_diag) {
        warp_atomic_add(&diag[0], ke);
        warp_atomic_add(&diag[1], mx);
        warp_atomic_add(&diag[2], my);
        warp_atomic_add(&diag[3], mz);
    }
    warp_atomic_add(&diag[4], lke);
    warp_atomic_add(&diag[5], lq);
}

// CIC charge density at nodes (diagnostic / Gauss check).
extern "C" __global__ void deposit_rho(const double* __restrict__ P, const int n, const int cap, const double q,
                                       double* rho) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    double X = (PX(p) - OX) * (1.0 / DX), Y = (PY(p) - OY) * (1.0 / DY), Z = (PZ(p) - OZ) * (1.0 / DZ);
    int i = (int)floor(X), j = (int)floor(Y), k = (int)floor(Z);
    double wx = X - i, wy = Y - j, wz = Z - k;
    double c = q * PW(p) * (1.0 / (DX * DY * DZ));
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int d = 0; d < 2; d++) {
                double s = (a ? wx : 1 - wx) * (b ? wy : 1 - wy) * (d ? wz : 1 - wz);
                if (s != 0.0) atomicAdd(&rho[IDX(i + a, j + b, k + d)], c * s);
            }
}

// Kinetic energy and momentum of particles at the current time level (used by the implicit
// scheme where u is stored at integer time levels).
extern "C" __global__ void kinetic(const double* __restrict__ P, const int n, const int cap, const double m,
                                   double* diag) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0, mx = 0, my = 0, mz = 0;
    if (p < n) {
        double ux = PUX(p), uy = PUY(p), uz = PUZ(p), w = PW(p);
        ke = w * m * (sqrt(1.0 + ux * ux + uy * uy + uz * uz) - 1.0);
        mx = w * m * ux;
        my = w * m * uy;
        mz = w * m * uz;
    }
    warp_atomic_add(&diag[0], ke);
    warp_atomic_add(&diag[1], mx);
    warp_atomic_add(&diag[2], my);
    warp_atomic_add(&diag[3], mz);
}

// Cell index (for sorting) of each particle.
extern "C" __global__ void cell_index(const double* __restrict__ P, const int n, const int cap,
                                      long long* key) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    int i = (int)floor((PX(p) - OX) * (1.0 / DX)), j = (int)floor((PY(p) - OY) * (1.0 / DY)),
        k = (int)floor((PZ(p) - OZ) * (1.0 / DZ));
    key[p] = IDX(i, j, k);
}

// ---- Loading: deterministic (decomposition independent) uniform positions + drifting
//      Maxwell-Juttner momenta (Zenitani 2015). One thread per (cell, slot).
//      Positions depend only on (seed, global cell, slot) so that all species/beams loaded
//      with the same seed sit on identical positions -> exact charge neutrality.
__device__ __forceinline__ double gauss(unsigned long long& s) {
    double r1 = urand(s), r2 = urand(s);
    return sqrt(-2.0 * log(r1)) * cos(6.283185307179586 * r2);
}

// Kinetic energy eps = gamma - 1 of a Maxwell-Juttner distribution, exact for any theta.
//  * theta < 1: rejection with the envelope obtained from sqrt(1 + eps/2) <= 1 + eps/4,
//      f(eps) ~ sqrt(eps) sqrt(1+eps/2) (1+eps) e^{-eps/theta}
//            <= sqrt(eps) (1 + 5 eps/4 + eps^2/4) e^{-eps/theta},
//    i.e. a mixture of Gamma(3/2), Gamma(5/2), Gamma(7/2) (scale theta) with weights
//    1/2 : 15 theta/16 : 15 theta^2/32; Gamma(k/2, theta) = theta/2 * chi^2_k (sum of k squared
//    Gaussians). Acceptance sqrt(1+eps/2)/(1+eps/4) >= ~0.94 for theta <= 1.
//  * theta >= 1: Sobol's algorithm (efficient for relativistic temperatures).
// A sample is never returned unaccepted: if the theta < 1 branch ever hit its (practically
// unreachable) iteration limit it falls through to Sobol's algorithm, which is exact as well.
__device__ double sample_mj_energy(double theta, unsigned long long& s, double* ur_out) {
    if (theta < 1.0) {
        const double w1 = 0.5, w2 = 0.9375 * theta, w3 = 0.46875 * theta * theta;
        const double wt = w1 + w2 + w3;
        for (int it = 0; it < 1000000; it++) {
            double r = urand(s) * wt;
            int k = r < w1 ? 3 : (r < w1 + w2 ? 5 : 7);
            double chi2 = 0.0;
            for (int i = 0; i < k; i++) {
                double g = gauss(s);
                chi2 += g * g;
            }
            double eps = 0.5 * theta * chi2;
            double acc = sqrt(1.0 + 0.5 * eps) / (1.0 + 0.25 * eps);
            if (urand(s) <= acc) {
                *ur_out = sqrt(eps * (eps + 2.0));
                return eps;
            }
        }
    }
    // Sobol: u = -theta ln(x1 x2 x3) accepted when eta^2 - u^2 > 1, eta = -theta ln(x1 x2 x3 x4)
    for (;;) {
        double x1 = urand(s), x2 = urand(s), x3 = urand(s), x4 = urand(s);
        double ur = -theta * log(x1 * x2 * x3);
        double eta = -theta * log(x1 * x2 * x3 * x4);
        if (eta * eta - ur * ur > 1.0) {
            *ur_out = ur;
            return sqrt(1.0 + ur * ur) - 1.0;
        }
    }
}

__device__ void sample_mj(double theta, unsigned long long& s, double* u) {
    if (theta <= 0.0) {
        u[0] = u[1] = u[2] = 0.0;
        return;
    }
    double ur;
    sample_mj_energy(theta, s, &ur);
    double x5 = urand(s), x6 = urand(s);
    double ct = 2.0 * x5 - 1.0, st = sqrt(fmax(0.0, 1.0 - ct * ct)), ph = 6.283185307179586 * x6;
    u[0] = ur * st * cos(ph);
    u[1] = ur * st * sin(ph);
    u[2] = ur * ct;
}

extern "C" __global__ void load_particles(double* P, const int cap, const long long offset, const int ppc,
                                          const int nlx, const int nly, const int nlz, const int gx0,
                                          const int gy0, const int gz0, const int gnx, const int gny,
                                          const unsigned long long seed, const unsigned long long stream,
                                          const double w, const double theta, const double udx,
                                          const double udy, const double udz) {
    long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    long long ntot = (long long)nlx * nly * nlz * ppc;
    if (t >= ntot) return;
    int slot = (int)(t % ppc);
    long long c = t / ppc;
    int k = (int)(c % nlz);
    c /= nlz;
    int j = (int)(c % nly);
    int i = (int)(c / nly);
    long long gcell = ((long long)(gx0 + i) * gny + (gy0 + j)) * (long long)((long long)1 << 30) + (gz0 + k);
    unsigned long long sp = hash_key(seed, (unsigned long long)gcell, (unsigned long long)slot);
    double x = (gx0 + i + urand(sp)) * DX;
    double y = (gy0 + j + urand(sp)) * DY;
    double z = (gz0 + k + urand(sp)) * DZ;
    unsigned long long sm = hash_key(seed ^ (stream * 0x2545F4914F6CDD1DULL + 77ULL), (unsigned long long)gcell,
                                     (unsigned long long)slot);
    double u[3];
    sample_mj(theta, sm, u);
    // Lorentz boost with drift (flipping method of Zenitani 2015)
    double ud2 = udx * udx + udy * udy + udz * udz;
    if (ud2 > 0.0) {
        double Gd = sqrt(1.0 + ud2);
        double udn = sqrt(ud2);
        double nx = udx / udn, ny = udy / udn, nz = udz / udn;
        double beta = udn / Gd;
        double g = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        double upar = u[0] * nx + u[1] * ny + u[2] * nz;
        if (-beta * upar / g > urand(sm)) upar = -upar;  // flip
        // recompute perpendicular part (unchanged) and boost the parallel part
        double uperpx = u[0] - (u[0] * nx + u[1] * ny + u[2] * nz) * nx;
        double uperpy = u[1] - (u[0] * nx + u[1] * ny + u[2] * nz) * ny;
        double uperpz = u[2] - (u[0] * nx + u[1] * ny + u[2] * nz) * nz;
        double upar_b = Gd * (upar + beta * g);
        u[0] = uperpx + upar_b * nx;
        u[1] = uperpy + upar_b * ny;
        u[2] = uperpz + upar_b * nz;
    }
    long long p = offset + t;
    PX(p) = x;
    PY(p) = y;
    PZ(p) = z;
    PUX(p) = u[0];
    PUY(p) = u[1];
    PUZ(p) = u[2];
    PW(p) = w;
}
