// Energy- and charge-conserving implicit particle kernels (Chen-Chacon type).
// Particle push: relativistic implicit midpoint,
//   u1 - u0 = (q dt / m) (E_eff + vbar x B),  x1 - x0 = dt vbar,  vbar = (u0 + u1)/(g0 + g1)
// so that m (g1 - g0) = q dt vbar . E_eff exactly. The path x0 -> x1 is split at every cell
// face; the same path-integrated weights are used to gather E_eff and to deposit J, hence
// sum_g E_g . J_g V dt == sum_p q w dt vbar . E_eff (exact energy exchange) and the
// deposited J satisfies the discrete continuity equation for CIC charge (exact Gauss law).
// Pn rows: 0..2 displacement d = x1 - x0 (unfolded), 3..5 u1 (unfolded frame).

#ifndef MAXSEG
#define MAXSEG 32
#endif

__device__ __forceinline__ double fold1(double x, int a, double& s) {
    const double Lg[3] = {LX, LY, LZ};
    const int b0[3] = {BCX0, BCY0, BCZ0};
    const int b1[3] = {BCX1, BCY1, BCZ1};
    s = 1.0;
    if (b0[a] == BC_REFLECT && x < 0.0) { s = -1.0; return -x; }
    if (b1[a] == BC_REFLECT && x > Lg[a]) { s = -1.0; return 2.0 * Lg[a] - x; }
    return x;
}

// mode 0: gather E_eff (returned in Eo), mode 1: deposit J with charge qw.
__device__ void path_op(const double* p0, const double* d, int mode, const double* __restrict__ Ex,
                        const double* __restrict__ Ey, const double* __restrict__ Ez, double* Eo, double qw,
                        double* Jx, double* Jy, double* Jz, int* err) {
    const double inv[3] = {1.0 / DX, 1.0 / DY, 1.0 / DZ};
    const double org[3] = {OX, OY, OZ};
    double ts[MAXSEG + 2];
    int nt = 0;
    ts[nt++] = 0.0;
    for (int a = 0; a < 3; a++) {
        double Xa = p0[a] * inv[a], Xb = (p0[a] + d[a]) * inv[a];
        if (Xb > Xa) {
            for (double nn = floor(Xa) + 1.0; nn < Xb; nn += 1.0) {
                if (nt < MAXSEG + 1) ts[nt++] = (nn - Xa) / (Xb - Xa); else *err = 1;
            }
        } else if (Xb < Xa) {
            for (double nn = ceil(Xa) - 1.0; nn > Xb; nn -= 1.0) {
                if (nt < MAXSEG + 1) ts[nt++] = (nn - Xa) / (Xb - Xa); else *err = 1;
            }
        }
    }
    ts[nt++] = 1.0;
    // insertion sort
    for (int a = 1; a < nt; a++) {
        double v = ts[a];
        int b = a - 1;
        while (b >= 0 && ts[b] > v) { ts[b + 1] = ts[b]; b--; }
        ts[b + 1] = v;
    }
    if (mode == 0) Eo[0] = Eo[1] = Eo[2] = 0.0;
    const double third = 1.0 / 3.0;
    for (int sgi = 0; sgi + 1 < nt; sgi++) {
        double ta = ts[sgi], tb = ts[sgi + 1];
        double dt_ = tb - ta;
        if (dt_ <= 0.0) continue;
        double tm = 0.5 * (ta + tb);
        double A[3], B[3], s[3];
        int cell[3];
        for (int a = 0; a < 3; a++) {
            double sm;
            const double umid = p0[a] + tm * d[a];  // unfolded midpoint
            double mid = fold1(umid, a, sm);
            s[a] = sm;
            double xa = p0[a] + ta * d[a], xb = p0[a] + tb * d[a];
            // fold endpoints about the wall the segment lies beyond. The wall is chosen from the
            // UNFOLDED midpoint: the folded one can be past the box centre when |d| > L/2.
            if (sm < 0.0) {
                double Lg[3] = {LX, LY, LZ};
                double wall = (umid < 0.0) ? 0.0 : Lg[a];
                xa = 2.0 * wall - xa;
                xb = 2.0 * wall - xb;
            }
            A[a] = (xa - org[a]) * inv[a];
            B[a] = (xb - org[a]) * inv[a];
            cell[a] = (int)floor((mid - org[a]) * inv[a]);
            A[a] -= cell[a];
            B[a] -= cell[a];  // fractional position within the cell, in [0,1]
        }
        int i = cell[0], j = cell[1], k = cell[2];
        // linear weights of node (cell) and (cell+1) at the endpoints
        double S0x[2] = {1.0 - A[0], A[0]}, S0y[2] = {1.0 - A[1], A[1]}, S0z[2] = {1.0 - A[2], A[2]};
        double dSx[2] = {-(B[0] - A[0]), B[0] - A[0]}, dSy[2] = {-(B[1] - A[1]), B[1] - A[1]},
               dSz[2] = {-(B[2] - A[2]), B[2] - A[2]};
        double dX = B[0] - A[0], dY = B[1] - A[1], dZ = B[2] - A[2];  // folded displacement (cells)
        for (int b = 0; b < 2; b++)
            for (int c = 0; c < 2; c++) {
                double wyz = S0y[b] * S0z[c] + 0.5 * dSy[b] * S0z[c] + 0.5 * S0y[b] * dSz[c] + third * dSy[b] * dSz[c];
                double wxz = S0x[b] * S0z[c] + 0.5 * dSx[b] * S0z[c] + 0.5 * S0x[b] * dSz[c] + third * dSx[b] * dSz[c];
                double wxy = S0x[b] * S0y[c] + 0.5 * dSx[b] * S0y[c] + 0.5 * S0x[b] * dSy[c] + third * dSx[b] * dSy[c];
                long long ix = IDX(i, j + b, k + c), iy = IDX(i + b, j, k + c), iz = IDX(i + b, j + c, k);
                if (mode == 0) {
                    Eo[0] += dt_ * s[0] * wyz * Ex[ix];
                    Eo[1] += dt_ * s[1] * wxz * Ey[iy];
                    Eo[2] += dt_ * s[2] * wxy * Ez[iz];
                } else {
                    double jx = qw * dX * (1.0 / (DT * DY * DZ)) * wyz;
                    double jy = qw * dY * (1.0 / (DT * DX * DZ)) * wxz;
                    double jz = qw * dZ * (1.0 / (DT * DX * DY)) * wxy;
                    if (jx != 0.0) atomicAdd(&Jx[ix], jx);
                    if (jy != 0.0) atomicAdd(&Jy[iy], jy);
                    if (jz != 0.0) atomicAdd(&Jz[iz], jz);
                }
            }
    }
}

__device__ void wall_avg_B(const double* p0, const double* d, const double* __restrict__ Bx,
                           const double* __restrict__ By, const double* __restrict__ Bz, double* B) {
    const double Lg[3] = {LX, LY, LZ};
    const int b0[3] = {BCX0, BCY0, BCZ0};
    const int b1[3] = {BCX1, BCY1, BCZ1};
    double ts[5];
    int nt = 0;
    ts[nt++] = 0.0;
    for (int a = 0; a < 3; a++) {
        if (d[a] == 0.0) continue;
        double x1 = p0[a] + d[a];
        if (b0[a] == BC_REFLECT && (p0[a] < 0.0) != (x1 < 0.0)) ts[nt++] = (0.0 - p0[a]) / d[a];
        if (b1[a] == BC_REFLECT && (p0[a] > Lg[a]) != (x1 > Lg[a])) ts[nt++] = (Lg[a] - p0[a]) / d[a];
    }
    ts[nt++] = 1.0;
    for (int a = 1; a < nt; a++) {
        double v = ts[a];
        int b = a - 1;
        while (b >= 0 && ts[b] > v) { ts[b + 1] = ts[b]; b--; }
        ts[b + 1] = v;
    }
    B[0] = B[1] = B[2] = 0.0;
    for (int k = 0; k + 1 < nt; k++) {
        double f = ts[k + 1] - ts[k];
        if (f <= 0.0 && nt > 2) continue;
        if (nt == 2) f = 1.0;
        double tm = 0.5 * (ts[k] + ts[k + 1]), mid[3], s[3];
        for (int a = 0; a < 3; a++) mid[a] = fold1(p0[a] + tm * d[a], a, s[a]);
        double X = (mid[0] - OX) * (1.0 / DX), Y = (mid[1] - OY) * (1.0 / DY), Z = (mid[2] - OZ) * (1.0 / DZ);
        double det = s[0] * s[1] * s[2];
        B[0] += f * det * s[0] * interp(Bx, X, Y, Z, 0.0, 0.5, 0.5);
        B[1] += f * det * s[1] * interp(By, X, Y, Z, 0.5, 0.0, 0.5);
        B[2] += f * det * s[2] * interp(Bz, X, Y, Z, 0.5, 0.5, 0.0);
    }
}

// One evaluation of the particle map d -> dt * vbar(u1(E_eff(d), B(d))). Returns new displacement
// dn and u1 (both in the unfolded frame).
__device__ void particle_map(const double* p0, const double* u0, double g0, const double* d, double eps,
                             const double* __restrict__ Ex, const double* __restrict__ Ey,
                             const double* __restrict__ Ez, const double* __restrict__ Bx,
                             const double* __restrict__ By, const double* __restrict__ Bz, double* u1,
                             double* dn, int* err, double* E) {
    path_op(p0, d, 0, Ex, Ey, Ez, E, 0.0, nullptr, nullptr, nullptr, err);
    // B averaged over the pieces of the path on either side of reflecting walls, each piece
    // evaluated at its (folded) midpoint and mapped to the unfolded frame as a pseudovector
    // (B_unf = det(R) R B). Without the split the map would jump when the midpoint crosses a
    // wall with non-zero tangential B. B does no work, so energy conservation is unaffected.
    double B[3];
    wall_avg_B(p0, d, Bx, By, Bz, B);
    double av[3] = {u0[0] + eps * E[0], u0[1] + eps * E[1], u0[2] + eps * E[2]};
    double g1 = sqrt(1.0 + u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]);
    // Sweeps over g1, each u1 = 2 ub - u0 with ub = av + ub x t, t = b / m, b = 2 eps B, m = g0 + g1:
    //   ub = (m^2 av + m av x b + (av . b) b) / (m^2 + b^2),
    // one division per sweep (FP64 divisions are the most expensive operations here). They stop when a
    // sweep leaves g1 unchanged, or when the sweeps still to come provably change nothing above that
    // level: ub = M(t) av with |M| <= 1 gives |d u1 / d g1| <= 2 |b| |av| / m^2 = L, a bound that also
    // holds for the new g1 (|u1| < g1) and decreases with g1. If L <= 1/2 at the lowest g1 within the
    // last change dg of both ends of this sweep, all later sweeps stay in that range and u1 and g1 are
    // within 2 L dg of their limits. kb = (L m^2)^2: the test needs no division.
    const double b[3] = {2.0 * eps * B[0], 2.0 * eps * B[1], 2.0 * eps * B[2]};
    const double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    const double ab = av[0] * b[0] + av[1] * b[1] + av[2] * b[2];
    const double axb[3] = {av[1] * b[2] - av[2] * b[1], av[2] * b[0] - av[0] * b[2], av[0] * b[1] - av[1] * b[0]};
    const double kb = 4.0 * b2 * (av[0] * av[0] + av[1] * av[1] + av[2] * av[2]);
    for (int r = 0; r < 60; r++) {
        const double m = g0 + g1;
        const double inv = 2.0 / (m * m + b2);
        for (int a = 0; a < 3; a++) u1[a] = (m * (m * av[a] + axb[a]) + ab * b[a]) * inv - u0[a];
        const double usq = u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2];
        const double g1n = sqrt(1.0 + usq);
        const double dg = fabs(g1n - g1);
        bool done = dg <= 2e-16 * g1n;
        if (!done) {
            const double ml = g0 + fmax(1.0, fmin(g1, g1n) - dg);
            const double m4 = (ml * ml) * (ml * ml);
            done = 4.0 * kb <= m4 && kb * dg * dg <= 1e-32 * usq * m4;  // L <= 1/2 and L dg <= 1e-16 |u1|
        }
        g1 = g1n;
        if (done) break;
    }
    // vbar = (u0 + u1) / (g0 + g1): exact energy identity (u1 - u0).vbar = g1 - g0
    const double dtg = DT / (g0 + g1);
    for (int a = 0; a < 3; a++) dn[a] = dtg * (u0[a] + u1[a]);
}

extern "C" __global__ void implicit_push(const double* __restrict__ P, double* Pn, const int n, const long long cap,
                                         const double* __restrict__ Ex, const double* __restrict__ Ey,
                                         const double* __restrict__ Ez, const double* __restrict__ Bx,
                                         const double* __restrict__ By, const double* __restrict__ Bz,
                                         double* Jx, double* Jy, double* Jz, const double eps, const double qm,
                                         const double qdep, const int first, const double tol, const int maxit,
                                         int* err) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    const double p0[3] = {P[p], P[cap + p], P[2 * cap + p]};
    const double u0[3] = {P[3 * cap + p], P[4 * cap + p], P[5 * cap + p]};
    const double w = P[6 * cap + p];
    const double g0 = sqrt(1.0 + u0[0] * u0[0] + u0[1] * u0[1] + u0[2] * u0[2]);
    const double h[3] = {DX, DY, DZ};
    const double ih[3] = {1.0 / DX, 1.0 / DY, 1.0 / DZ};
    double d[3], u1[3], dn[3];
    if (first) {
        for (int a = 0; a < 3; a++) { d[a] = DT * u0[a] / g0; u1[a] = u0[a]; }
    } else {
        for (int a = 0; a < 3; a++) { d[a] = Pn[a * cap + p]; u1[a] = Pn[(3 + a) * cap + p]; }
    }
    double Eeff[3];
    bool conv = false;
    int it;
    // current iterate d with its evaluation (dn, u1, Eeff) and residual r = F(d) - d
    particle_map(p0, u0, g0, d, eps, Ex, Ey, Ez, Bx, By, Bz, u1, dn, err, Eeff);
    double r[3], nrm = 0.0;
    for (int a = 0; a < 3; a++) { r[a] = dn[a] - d[a]; nrm = fmax(nrm, fabs(r[a]) * ih[a]); }
    for (it = 0; it < maxit; it++) {
        if (nrm < tol) { conv = true; break; }  // d is the path used for both gather and deposit
        double dtr[3];
        bool accepted = false;
        if (it < 4) {
            for (int a = 0; a < 3; a++) dtr[a] = dn[a];  // plain Picard
            accepted = true;
        } else {
            // Newton step on G(d) = F(d) - d with a finite-difference Jacobian, trust region and
            // back-tracking. Needed for particles grazing a node plane or a wall, where the
            // 0th-order (charge-conserving) weight makes F stiff and only piecewise smooth.
            double JG[3][3];
            for (int b = 0; b < 3; b++) {
                double dp[3] = {d[0], d[1], d[2]}, up[3] = {u1[0], u1[1], u1[2]}, dnp[3], Et[3];
                double hb = 1e-7 * h[b];
                dp[b] += hb;
                particle_map(p0, u0, g0, dp, eps, Ex, Ey, Ez, Bx, By, Bz, up, dnp, err, Et);
                for (int a = 0; a < 3; a++) JG[a][b] = (a == b ? 1.0 : 0.0) - (dnp[a] - dn[a]) / hb;
            }
            double det = JG[0][0] * (JG[1][1] * JG[2][2] - JG[1][2] * JG[2][1]) -
                         JG[0][1] * (JG[1][0] * JG[2][2] - JG[1][2] * JG[2][0]) +
                         JG[0][2] * (JG[1][0] * JG[2][1] - JG[1][1] * JG[2][0]);
            double dl[3] = {r[0], r[1], r[2]};
            if (fabs(det) > 1e-14) {
                double inv = 1.0 / det;
                dl[0] = inv * (r[0] * (JG[1][1] * JG[2][2] - JG[1][2] * JG[2][1]) - JG[0][1] * (r[1] * JG[2][2] - JG[1][2] * r[2]) +
                               JG[0][2] * (r[1] * JG[2][1] - JG[1][1] * r[2]));
                dl[1] = inv * (JG[0][0] * (r[1] * JG[2][2] - JG[1][2] * r[2]) - r[0] * (JG[1][0] * JG[2][2] - JG[1][2] * JG[2][0]) +
                               JG[0][2] * (JG[1][0] * r[2] - r[1] * JG[2][0]));
                dl[2] = inv * (JG[0][0] * (JG[1][1] * r[2] - r[1] * JG[2][1]) - JG[0][1] * (JG[1][0] * r[2] - r[1] * JG[2][0]) +
                               r[0] * (JG[1][0] * JG[2][1] - JG[1][1] * JG[2][0]));
            }
            double mx = 0.0;
            for (int a = 0; a < 3; a++) mx = fmax(mx, fabs(dl[a]) * ih[a]);
            if (mx > 0.5) for (int a = 0; a < 3; a++) dl[a] *= 0.5 / mx;  // trust region: half a cell
            double lam = 1.0;
            for (int ls = 0; ls < 8 && !accepted; ls++, lam *= 0.5) {
                double dt2[3], u2[3] = {u1[0], u1[1], u1[2]}, dn2[3], E2[3], r2[3], n2 = 0.0;
                for (int a = 0; a < 3; a++) dt2[a] = d[a] + lam * dl[a];
                particle_map(p0, u0, g0, dt2, eps, Ex, Ey, Ez, Bx, By, Bz, u2, dn2, err, E2);
                for (int a = 0; a < 3; a++) { r2[a] = dn2[a] - dt2[a]; n2 = fmax(n2, fabs(r2[a]) * ih[a]); }
                if (n2 < nrm) {
                    for (int a = 0; a < 3; a++) { d[a] = dt2[a]; u1[a] = u2[a]; dn[a] = dn2[a]; Eeff[a] = E2[a]; r[a] = r2[a]; }
                    nrm = n2;
                    accepted = true;
                }
            }
            if (accepted) continue;
            // Robust fallback: 1-D bracketed root solve (Illinois) of G_a(d_a) = 0 on the axis with the
            // largest residual, other components fixed. Since |F_a| < c dt, [d_a-.., d_a+..] can always be
            // expanded up to [-DT, DT] to bracket a sign change, so this step always converges.
            int ax = 0;
            {
                double best = -1.0;
                for (int a = 0; a < 3; a++) if (fabs(r[a]) * ih[a] > best) { best = fabs(r[a]) * ih[a]; ax = a; }
            }
            {
                double dl0[3] = {d[0], d[1], d[2]}, ut[3] = {u1[0], u1[1], u1[2]}, dnt[3], Et[3];
                double ga = r[ax], xa = d[ax];
                double step = fmax(2.0 * fabs(ga), 1e-6 * h[ax]);
                double xb = xa, gb = ga;
                // expand in the direction indicated by the sign of the residual
                for (int e = 0; e < 60; e++) {
                    xb = xa + (ga > 0 ? step : -step);
                    if (xb > DT) xb = DT;
                    if (xb < -DT) xb = -DT;
                    dl0[ax] = xb;
                    particle_map(p0, u0, g0, dl0, eps, Ex, Ey, Ez, Bx, By, Bz, ut, dnt, err, Et);
                    gb = dnt[ax] - xb;
                    if ((gb > 0) != (ga > 0) || fabs(xb) >= DT) break;
                    xa = xb; ga = gb; step *= 2.0;
                }
                // Illinois iterations on [xa, xb]
                int side = 0;
                for (int e = 0; e < 200 && (gb > 0) != (ga > 0); e++) {
                    double xc = (xa * gb - xb * ga) / (gb - ga);
                    if (!(xc > fmin(xa, xb) && xc < fmax(xa, xb))) xc = 0.5 * (xa + xb);
                    dl0[ax] = xc;
                    particle_map(p0, u0, g0, dl0, eps, Ex, Ey, Ez, Bx, By, Bz, ut, dnt, err, Et);
                    double gc = dnt[ax] - xc;
                    if ((gc > 0) == (gb > 0)) {
                        xb = xc; gb = gc;
                        if (side == -1) ga *= 0.5;
                        side = -1;
                    } else {
                        xa = xc; ga = gc;
                        if (side == 1) gb *= 0.5;
                        side = 1;
                    }
                    if (fabs(xb - xa) < 0.1 * tol * h[ax] || gc == 0.0) break;
                }
                for (int a = 0; a < 3; a++) dtr[a] = d[a];
                dtr[ax] = (fabs(ga) < fabs(gb)) ? xa : xb;
            }
        }
        // |d| < c dt (speed of light) safeguard
        double dd = sqrt(dtr[0] * dtr[0] + dtr[1] * dtr[1] + dtr[2] * dtr[2]);
        if (dd >= DT) for (int a = 0; a < 3; a++) dtr[a] *= (1.0 - 1e-12) * DT / dd;
        for (int a = 0; a < 3; a++) d[a] = dtr[a];
        particle_map(p0, u0, g0, d, eps, Ex, Ey, Ez, Bx, By, Bz, u1, dn, err, Eeff);
        nrm = 0.0;
        for (int a = 0; a < 3; a++) { r[a] = dn[a] - d[a]; nrm = fmax(nrm, fabs(r[a]) * ih[a]); }
    }
    if (!conv) atomicAdd(err + 1, 1);  // non-converged particle counter
    // Exact energy closure: m (g1 - g0) = q d . E_eff along the deposited path d. For converged
    // particles this is a round-off level correction; for the rare non-converged ones it keeps
    // the total energy exact (charge conservation is unaffected since x1 = x0 + d).
    {
        double gt = g0 + qm * (d[0] * Eeff[0] + d[1] * Eeff[1] + d[2] * Eeff[2]);
        double un = sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]);
        if (gt >= 1.0 && un > 0.0) {
            double sc = sqrt(gt * gt - 1.0) / un;
            for (int a = 0; a < 3; a++) u1[a] *= sc;
        } else if (gt < 1.0) {
            err[0] = err[0] < 3 ? 3 : err[0];
        }
    }
    // ghost-width guard
    if (fabs(d[0]) > (GH - 1) * DX || fabs(d[1]) > (GH - 1) * DY || fabs(d[2]) > (GH - 1) * DZ) err[0] = 2;
    double dummy[3];
    if (qdep != 0.0) path_op(p0, d, 1, Ex, Ey, Ez, dummy, qdep * w, Jx, Jy, Jz, err);  // test species: no deposit
    for (int a = 0; a < 3; a++) { Pn[a * cap + p] = d[a]; Pn[(3 + a) * cap + p] = u1[a]; }
}

// Commit the converged step: x1 = x0 + d with reflections folded (u flipped), open-boundary
// removal, periodic wrap on non-decomposed axes and migration flagging.
extern "C" __global__ void implicit_finalize(double* P, const double* __restrict__ Pn, const int n, const long long cap,
                                             const double q, const double m, const double qdep, double* diag, int* flag,
                                             int* list, const int listcap) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    double lke = 0.0, lq = 0.0;
    if (p < n) {
        const double Lg[3] = {LX, LY, LZ};
        const int b0[3] = {BCX0, BCY0, BCZ0};
        const int b1[3] = {BCX1, BCY1, BCZ1};
        const int dec[3] = {DECX, DECY, DECZ};
        const double slo[3] = {SX0, SY0, SZ0};
        const double shi[3] = {SX1, SY1, SZ1};
        double x[3], u[3];
        bool dead = false, migrate = false;
        for (int a = 0; a < 3; a++) {
            double s;
            x[a] = fold1(P[a * cap + p] + Pn[a * cap + p], a, s);
            u[a] = s * Pn[(3 + a) * cap + p];
            if ((b0[a] == BC_OPEN && x[a] < 0.0) || (b1[a] == BC_OPEN && x[a] >= Lg[a])) dead = true;
        }
        double w = P[6 * cap + p];
        if (dead) {
            lke = w * m * (sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) - 1.0);
            lq = qdep * w;
            P[6 * cap + p] = 0.0;
            { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        } else {
            for (int a = 0; a < 3; a++) {
                if (!dec[a]) {
                    if (b0[a] == BC_PERIODIC) x[a] -= Lg[a] * floor(x[a] / Lg[a]);
                } else if (x[a] < slo[a] || x[a] >= shi[a]) {
                    migrate = true;
                }
            }
            if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        }
        for (int a = 0; a < 3; a++) { P[a * cap + p] = x[a]; P[(3 + a) * cap + p] = u[a]; }
    }
    warp_atomic_add(&diag[4], lke);
    warp_atomic_add(&diag[5], lq);
}
