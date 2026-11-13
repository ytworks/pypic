// Explicit particle kernels: relativistic Buneman-Boris push fused with
// Esirkepov charge-conserving current deposition (CIC = cubic super-particles).
//
// Particle storage. A kernel of this file reaches a particle only through the accessors below, so
// the layout behind them changes with the build (PIC-D19). cap is 64-bit (row * cap exceeds 2^31
// for large species); n and p stay int. Exceptions: the remaining displacement of displace_deposit
// (set-up only), and the row copies of migrate.cu.

#if defined(COMPACT_STORE) && defined(MOVING_WINDOW)
#error "the compact storage has no moving window (PIC-D19, PIC-D22)"
#endif
#ifdef COMPACT_STORE
// Compact layout, 56 bytes per particle in four arrays:
//   S_w[row * cap + p]  rows = w, id (double; the id is an integer-valued double)
//   S_c[p]              the cell of the particle in the arrays of this rank (IDX, int)
//   S_f[a * cap + p]    the position within that cell along axis a, in units of 2^-32 cell (unsigned int)
//   S_u[a * cap + p]    the momentum (double; u = gamma v, c = 1)
// A position is a point of a lattice of 2^-32 cell. Its index coordinate cell + f * 2^-32 is exact in
// a double while |cell| < 2^21, in the arrays of the rank and in the box (Simulation checks the grid).
#define TWO32 4294967296.0
#define ITWO32 (1.0 / 4294967296.0)
#define PSTORE_AS(S) double* S##_w, int* S##_c, unsigned int* S##_f, double* S##_u
#define PSTORE PSTORE_AS(P)
#define PSTORE_RO \
    const double* __restrict__ P_w, const int* __restrict__ P_c, const unsigned int* __restrict__ P_f, \
    const double* __restrict__ P_u
// index coordinate (cells from the first one of the arrays of this rank) of particle p of the storage S
#define S_XI(S, cap, p) ((double)((S##_c)[(p)] / (NY * NZ)) + (double)(S##_f)[(p)] * ITWO32)
#define S_YI(S, cap, p) ((double)((S##_c)[(p)] / NZ % NY) + (double)(S##_f)[(cap) + (p)] * ITWO32)
#define S_ZI(S, cap, p) ((double)((S##_c)[(p)] % NZ) + (double)(S##_f)[2 * (cap) + (p)] * ITWO32)
// values of particle p of the storage (S, cap): physical position, momentum, weight, id
#define S_X(S, cap, p) ((S_XI(S, cap, p) + (double)(OXI)) * DX)
#define S_Y(S, cap, p) ((S_YI(S, cap, p) + (double)(OYI)) * DY)
#define S_Z(S, cap, p) ((S_ZI(S, cap, p) + (double)(OZI)) * DZ)
#define S_UX(S, cap, p) (+(S##_u)[(p)])
#define S_UY(S, cap, p) (+(S##_u)[(cap) + (p)])
#define S_UZ(S, cap, p) (+(S##_u)[2 * (cap) + (p)])
#define S_W(S, cap, p) (+(S##_w)[(p)])
#define S_ID(S, cap, p) (+(S##_w)[(cap) + (p)])
#define S_SET_POS(S, cap, p, x, y, z) store_pos(S##_c, S##_f, cap, p, x, y, z)
#define S_SET_MOM(S, cap, p, ux, uy, uz) \
    { (S##_u)[(p)] = (ux); (S##_u)[(cap) + (p)] = (uy); (S##_u)[2 * (cap) + (p)] = (uz); }
#define S_SET_W(S, cap, p, w) { (S##_w)[(p)] = (w); }
#define S_SET_ID(S, cap, p, id) { (S##_w)[(cap) + (p)] = (id); }
#define XI(p) S_XI(P, cap, p)
#define YI(p) S_YI(P, cap, p)
#define ZI(p) S_ZI(P, cap, p)
// Working coordinate: what the kernels that move or place particles (push_deposit,
// displace_deposit, the loader, the injection) compute in, and what SET_POS takes. This layout:
// the index coordinate of the box (cells from its lower corner). Walls, the range of the rank and
// the period are integers there, so a reflection, a wrap and the change of rank keep a lattice
// point exact; every point handed to SET_POS is first put on the lattice (WSNAP, WRAND).
#define QX(p) (XI(p) + (double)(OXI))            // working coordinate of particle p
#define QY(p) (YI(p) + (double)(OYI))
#define QZ(p) (ZI(p) + (double)(OZI))
#define WBOX {(double)(GNX), (double)(GNY), (double)(GNZ)}     // the box
#define WCELL {1.0, 1.0, 1.0}                                  // one cell
#define WLO {(double)(SX0I), (double)(SY0I), (double)(SZ0I)}   // physical range of the rank
#define WHI {(double)(SX1I), (double)(SY1I), (double)(SZ1I)}
#define TO_XI(x) ((x) - (double)(OXI))           // working coordinate -> index coordinate
#define TO_YI(y) ((y) - (double)(OYI))
#define TO_ZI(z) ((z) - (double)(OZI))
#define WPHYS(x, d) ((x) * (d))                  // -> physical position (d: the cell of the axis)
#define WLEN(s, d) ((s) * (1.0 / (d)))           // physical length -> working coordinate
#define WDT(a) ((a) == 0 ? DT / DX : ((a) == 1 ? DT / DY : DT / DZ))  // time step per cell of axis a
#define WSNAP(x) (rint((x) * TWO32) * ITWO32)    // onto a point that the storage holds
#define WBELOW(wall) ((wall) - ITWO32)           // last point below an upper wall
// uniform point of cell c (advances the stream s): the upper 32 bits of the number urand takes
#define WRAND(c, s, d) ((double)(c) + (double)(splitmix64(s) >> 32) * ITWO32)
#define WPOINT(c, f, d) WSNAP((double)(c) + (f))
#define PDISP_ARG , double* Dsp                  // the remaining set-up displacement: its own array
#define DISP(a, p) Dsp[(a) * cap + p]

// Stores the lattice point (x, y, z) of the working coordinate: at most one cell outside the
// physical range of the rank, so its cell is one of the arrays of the rank.
__device__ __forceinline__ void store_pos(int* C, unsigned int* F, const long long cap, const long long p,
                                          const double x, const double y, const double z) {
    const double X[3] = {x, y, z};
    const int off[3] = {OXI, OYI, OZI};
    int c[3];
#pragma unroll
    for (int a = 0; a < 3; a++) {
        const double g = floor(X[a]);
        c[a] = (int)g - off[a];
        F[a * cap + p] = (unsigned int)((X[a] - g) * TWO32);
    }
    C[p] = (int)IDX(c[0], c[1], c[2]);
}

// The stored integers of particle p: c = its cell along the three axes (arrays of this rank), q = the
// position within that cell. The integer route of push_deposit computes in these.
__device__ __forceinline__ void load_cell(const int* C, const unsigned int* F, const long long cap, const long long p,
                                          int* c, unsigned int* q) {
    const int l = C[p];
    c[0] = l / (NY * NZ);
    c[1] = l / NZ % NY;
    c[2] = l % NZ;
#pragma unroll
    for (int a = 0; a < 3; a++) q[a] = F[a * cap + p];
}

__device__ __forceinline__ void store_cell(int* C, unsigned int* F, const long long cap, const long long p,
                                           const int* c, const unsigned int* q) {
    C[p] = (int)IDX(c[0], c[1], c[2]);
#pragma unroll
    for (int a = 0; a < 3; a++) F[a * cap + p] = q[a];
}
#define GET_CELL(p, c, q) load_cell(P_c, P_f, cap, p, c, q)
#define SET_CELL(p, c, q) store_cell(P_c, P_f, cap, p, c, q)
#else
// Layout of doubles: S[row * cap + p], rows = x, y, z, ux, uy, uz, w, id (u = gamma v, c = 1; the id
// is an integer-valued double); the remaining displacement of displace_deposit in rows 8..10.
#define PSTORE_AS(S) double* S                     // the storage of a species as kernel parameters
#define PSTORE PSTORE_AS(P)
#define PSTORE_RO const double* __restrict__ P
// values of particle p of the storage (S, cap): physical position, momentum, weight, id
#define S_X(S, cap, p) (+(S)[(p)])
#define S_Y(S, cap, p) (+(S)[(cap) + (p)])
#define S_Z(S, cap, p) (+(S)[2 * (cap) + (p)])
#define S_UX(S, cap, p) (+(S)[3 * (cap) + (p)])
#define S_UY(S, cap, p) (+(S)[4 * (cap) + (p)])
#define S_UZ(S, cap, p) (+(S)[5 * (cap) + (p)])
#define S_W(S, cap, p) (+(S)[6 * (cap) + (p)])
#define S_ID(S, cap, p) (+(S)[7 * (cap) + (p)])
#define S_SET_POS(S, cap, p, x, y, z) { (S)[(p)] = (x); (S)[(cap) + (p)] = (y); (S)[2 * (cap) + (p)] = (z); }
#define S_SET_MOM(S, cap, p, ux, uy, uz) \
    { (S)[3 * (cap) + (p)] = (ux); (S)[4 * (cap) + (p)] = (uy); (S)[5 * (cap) + (p)] = (uz); }
#define S_SET_W(S, cap, p, w) { (S)[6 * (cap) + (p)] = (w); }
#define S_SET_ID(S, cap, p, id) { (S)[7 * (cap) + (p)] = (id); }
// index coordinate: cells from the first one of the arrays of this rank
#define XI(p) XIDX(PX(p))
#define YI(p) ((PY(p) - OY) * (1.0 / DY))
#define ZI(p) ((PZ(p) - OZ) * (1.0 / DZ))
// Working coordinate: what the kernels that move or place particles (push_deposit,
// displace_deposit, the loader, the injection) compute in, and what SET_POS takes. This layout:
// the physical position.
#define QX(p) PX(p)                              // working coordinate of particle p
#define QY(p) PY(p)
#define QZ(p) PZ(p)
#define WBOX {LX, LY, LZ}                        // the box
#define WCELL {DX, DY, DZ}                       // one cell
#define WLO {RX0, SY0, SZ0}                      // physical range of the rank
#define WHI {RX1, SY1, SZ1}
#define TO_XI(x) XIDX(x)                         // working coordinate -> index coordinate
#define TO_YI(y) (((y) - OY) * (1.0 / DY))
#define TO_ZI(z) (((z) - OZ) * (1.0 / DZ))
#define WPHYS(x, d) (x)                          // -> physical position (d: the cell of the axis)
#define WLEN(s, d) (s)                           // physical length -> working coordinate
#define WDT(a) DT                                // time step, per cell of axis a where the unit is the cell
#define WSNAP(x) (x)                             // onto a point that the storage holds
#define WBELOW(wall) nextafter(wall, 0.0)        // last point below an upper wall
#define WRAND(c, s, d) (((c) + urand(s)) * (d))  // uniform point of cell c (advances the stream s)
#define WPOINT(c, f, d) (((c) + (f)) * (d))      // the point at the fraction f of cell c
#define PDISP_ARG                                // the remaining set-up displacement: rows 8..10
#define DISP(a, p) P[(8 + (a)) * cap + p]
#endif
// the same for the storage (P, cap) of the kernel
#define PX(p) S_X(P, cap, p)
#define PY(p) S_Y(P, cap, p)
#define PZ(p) S_Z(P, cap, p)
#define PUX(p) S_UX(P, cap, p)
#define PUY(p) S_UY(P, cap, p)
#define PUZ(p) S_UZ(P, cap, p)
#define PW(p) S_W(P, cap, p)
#define PID(p) S_ID(P, cap, p)
#define SET_POS(p, x, y, z) S_SET_POS(P, cap, p, x, y, z)
#define SET_MOM(p, ux, uy, uz) S_SET_MOM(P, cap, p, ux, uy, uz)
#define SET_W(p, w) S_SET_W(P, cap, p, w)
#define SET_ID(p, id) S_SET_ID(P, cap, p, id)

// CIC weights of the two ends of a move along one axis, on the 3 nodes from n = min(floor X0, floor X1)
// (|X1 - X0| < 1): (1 - w, w, 0) for an end in cell n, (0, 1 - w, w) for an end in cell n + 1, with
// w = X - floor X as in deposit_rho. Node 2 has weight only when the ends are in different cells (wide).
struct AxisShape {
    int n;
    bool wide;
    double S0[3], dS[3];
};

// The ends given as cell and position within the cell.
__device__ __forceinline__ void axis_shape_cells(const int c0, const double w0, const int c1, const double w1,
                                                 AxisShape& s) {
    const double v0 = 1.0 - w0, v1 = 1.0 - w1;
    const bool up0 = c0 > c1, up1 = c1 > c0;  // this end is in cell n + 1
    s.n = up0 ? c1 : c0;
    s.wide = c0 != c1;
    s.S0[0] = up0 ? 0.0 : v0;
    s.S0[1] = up0 ? v0 : w0;
    s.S0[2] = up0 ? w0 : 0.0;
    s.dS[0] = (up1 ? 0.0 : v1) - s.S0[0];
    s.dS[1] = (up1 ? v1 : w1) - s.S0[1];
    s.dS[2] = (up1 ? w1 : 0.0) - s.S0[2];
}

__device__ __forceinline__ void axis_shape(const double X0, const double X1, AxisShape& s) {
    const double f0 = floor(X0), f1 = floor(X1);
    axis_shape_cells((int)f0, X0 - f0, (int)f1, X1 - f1, s);
}

// Esirkepov deposition of a straight segment with the shapes x, y, z of its ends (axis_shape) in the
// arrays of the rank; mx, my, mz: the ends differ along that axis. qw = charge * weight.
// The current along an axis is the prefix sum of -c dS over its nodes: two faces carry current (the
// second only for a wide move; the third sum is sum dS = 0), each times the transverse weight
//   W(b, c) = S0 S0' + (dS S0' + S0 dS') / 2 + dS dS' / 3 = S0'[c] A[b] + dS'[c] B[b],
// A = S0 + dS / 2, B = S0 / 2 + dS / 3. An axis without displacement carries no current.
__device__ __forceinline__ void esirkepov_shapes(const AxisShape& x, const AxisShape& y, const AxisShape& z,
                                                 const bool mx, const bool my, const bool mz, const double qw,
                                                 double* Jx, double* Jy, double* Jz) {
    const double third = 1.0 / 3.0;
    double Ax[3], Bx[3], Ay[3], By[3];
#pragma unroll
    for (int a = 0; a < 3; a++) {
        Ax[a] = x.S0[a] + 0.5 * x.dS[a];
        Bx[a] = 0.5 * x.S0[a] + third * x.dS[a];
        Ay[a] = y.S0[a] + 0.5 * y.dS[a];
        By[a] = 0.5 * y.S0[a] + third * y.dS[a];
    }
    const long long base = IDX(x.n, y.n, z.n);
    if (mx) {
        const double cx = qw * (1.0 / (DY * DZ * DT));
        const double f0 = -cx * x.dS[0], f1 = f0 - cx * x.dS[1];
#pragma unroll
        for (int b = 0; b < 3; b++) {
            if (b == 2 && !y.wide) continue;
#pragma unroll
            for (int c = 0; c < 3; c++) {
                if (c == 2 && !z.wide) continue;
                const double W = z.S0[c] * Ay[b] + z.dS[c] * By[b];
                double* J = Jx + (base + b * STRIDE_J + c);
                atomicAdd(J, f0 * W);
                if (x.wide) atomicAdd(J + STRIDE_I, f1 * W);
            }
        }
    }
    if (my) {
        const double cy = qw * (1.0 / (DX * DZ * DT));
        const double f0 = -cy * y.dS[0], f1 = f0 - cy * y.dS[1];
#pragma unroll
        for (int a = 0; a < 3; a++) {
            if (a == 2 && !x.wide) continue;
#pragma unroll
            for (int c = 0; c < 3; c++) {
                if (c == 2 && !z.wide) continue;
                const double W = z.S0[c] * Ax[a] + z.dS[c] * Bx[a];
                double* J = Jy + (base + a * STRIDE_I + c);
                atomicAdd(J, f0 * W);
                if (y.wide) atomicAdd(J + STRIDE_J, f1 * W);
            }
        }
    }
    if (mz) {
        const double cz = qw * (1.0 / (DX * DY * DT));
        const double f0 = -cz * z.dS[0], f1 = f0 - cz * z.dS[1];
#pragma unroll
        for (int a = 0; a < 3; a++) {
            if (a == 2 && !x.wide) continue;
#pragma unroll
            for (int b = 0; b < 3; b++) {
                if (b == 2 && !y.wide) continue;
                const double W = y.S0[b] * Ax[a] + y.dS[b] * Bx[a];
                double* J = Jz + (base + a * STRIDE_I + b * STRIDE_J);
                atomicAdd(J, f0 * W);
                if (z.wide) atomicAdd(J + 1, f1 * W);
            }
        }
    }
}

// The same for the segment (X0,Y0,Z0)->(X1,Y1,Z1) in local index coordinates; |X1-X0| < 1 is required.
__device__ __forceinline__ void esirkepov(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double qw,
                                          double* Jx, double* Jy, double* Jz) {
    AxisShape x, y, z;
    axis_shape(X0, X1, x);
    axis_shape(Y0, Y1, y);
    axis_shape(Z0, Z1, z);
    esirkepov_shapes(x, y, z, X1 != X0, Y1 != Y0, Z1 != Z0, qw, Jx, Jy, Jz);
}

__device__ __forceinline__ void deposit_global(const double* p0, const double* p1, double qw, double* Jx,
                                               double* Jy, double* Jz) {
    esirkepov(TO_XI(p0[0]), TO_YI(p0[1]), TO_ZI(p0[2]), TO_XI(p1[0]), TO_YI(p1[1]), TO_ZI(p1[2]), qw, Jx, Jy, Jz);
}

// Deposit of the part of a path that lies inside the open faces of the box (PIC-D15): the current
// beyond an open face is not deposited. A particle that leaves through an open face leaves its
// charge on the nodes of the face, one that enters leaves the opposite charge there; nothing is
// deposited in the absorbing layer. Inside means 0 <= x < L, as for the removal. Without an open
// face this is deposit_global.
__device__ __forceinline__ void deposit_cut(const double* p0, const double* p1, double qw, double* Jx,
                                            double* Jy, double* Jz) {
    const double Lg[3] = WBOX;
    const bool o0[3] = {BCX0 == BC_OPEN, BCY0 == BC_OPEN, BCZ0 == BC_OPEN};
    const bool o1[3] = {BCX1 == BC_OPEN, BCY1 == BC_OPEN, BCZ1 == BC_OPEN};
    double t0 = 0.0, t1 = 1.0;
    for (int a = 0; a < 3; a++) {
        const double d = p1[a] - p0[a];
        if (o0[a]) {
            if (p0[a] < FLO(a) && p1[a] < FLO(a)) return;
            if (p0[a] < FLO(a)) t0 = fmax(t0, (FLO(a) - p0[a]) / d);
            else if (p1[a] < FLO(a)) t1 = fmin(t1, (FLO(a) - p0[a]) / d);
        }
        if (o1[a]) {
            if (p0[a] >= FHI(a) && p1[a] >= FHI(a)) return;
            if (p0[a] >= FHI(a)) t0 = fmax(t0, (FHI(a) - p0[a]) / d);
            else if (p1[a] >= FHI(a)) t1 = fmin(t1, (FHI(a) - p0[a]) / d);
        }
    }
    if (t0 == 0.0 && t1 == 1.0) {
        deposit_global(p0, p1, qw, Jx, Jy, Jz);
        return;
    }
    if (!(t0 < t1)) return;  // the path passes outside an edge of the box
    double q0[3], q1[3];
    for (int a = 0; a < 3; a++) {
        const double d = p1[a] - p0[a];
        q0[a] = t0 > 0.0 ? p0[a] + t0 * d : p0[a];
        q1[a] = t1 < 1.0 ? p0[a] + t1 * d : p1[a];
        if (o0[a]) { q0[a] = fmax(q0[a], FLO(a)); q1[a] = fmax(q1[a], FLO(a)); }
        if (o1[a]) { q0[a] = fmin(q0[a], FHI(a)); q1[a] = fmin(q1[a], FHI(a)); }
    }
    deposit_global(q0, q1, qw, Jx, Jy, Jz);
}

// Injection at open faces (PIC-D13); defined with inject_face below.
__device__ void sample_inflow(const double* rb, const int nb, const int axis, const double sgn,
                              unsigned long long& s, double* u, const double* dd);
// Fluctuation of the inflow (PIC-D26): du = sum over the nfl modes of fl (11 doubles each:
// kx, ky, kz, omega, phase, cos part ux, uy, uz, sin part ux, uy, uz) of
// c cos(th) + s sin(th), th = k . r - omega t + phase; added to the drift of an entering particle.
__device__ __forceinline__ void fluct_du(const double* fl, const int nfl, const double* r, const double t,
                                         double* du) {
    du[0] = du[1] = du[2] = 0.0;
    for (int m = 0; m < nfl; m++) {
        const double* q = fl + 11 * m;
        const double th = q[0] * r[0] + q[1] * r[1] + q[2] * r[2] - q[3] * t + q[4];
        const double c = cos(th), s = sin(th);
        for (int a = 0; a < 3; a++) du[a] += q[5 + a] * c + q[8 + a] * s;
    }
}
__device__ __forceinline__ int bank_slot(const int* binfo, const int face, const int* gc);
// salt of the random keys of the particles returned at a paired face: (seed, step, id, wall hit)
static const unsigned long long RET_SALT = 0x14057B7EF767814FULL;

// diag: [0] kinetic energy (time centred), [1..3] momentum, [4] kinetic energy lost
// through open boundaries (net of the energy of the particles returned at a paired face),
// [5] charge lost through open boundaries.
// flag[0]: number of particles that must be removed or migrated.
// h = q DT / (2 m): the half kick factor of the species.
// fmask (PIC-D13; 0 without injection): bit f = this species is returned at the paired face f
// instead of leaving (face f = 2 * axis + side); bit 8 + f = a particle of this species that
// leaves through the paired face f is added to the bank of its face column (bank, binfo).
// Only the faces of the compile-time mask PAIRED_FACES count: without a paired face the
// compiler removes these paths (and, without a reflecting wall, the whole wall loop).
// rb: inflow tables of this species, nrb rows of 6 per face (see sample_inflow).
extern "C" __global__ void push_deposit(PSTORE, const int n, const long long cap, const gather_t* __restrict__ Ex,
                                        const gather_t* __restrict__ Ey, const gather_t* __restrict__ Ez,
                                        const gather_t* __restrict__ Bx, const gather_t* __restrict__ By,
                                        const gather_t* __restrict__ Bz, double* Jx, double* Jy, double* Jz,
                                        const double h, const double m, const double qdep, double* diag, int* flag,
                                        int* list, const int listcap, const int do_diag,
                                        const unsigned long long seed, const long long step, const int fmask,
                                        const double* rb, const int nrb, double* bank, const int* binfo,
                                        const double* flu, const int nflu, const int flface, const double tnow) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0.0, mx = 0.0, my = 0.0, mz = 0.0, lke = 0.0, lq = 0.0;
#ifdef BODY_SPHERE
    double ake = 0.0, aq = 0.0, an = 0.0;  // absorbed by the body: kinetic energy, charge, particles
#endif
    if (p < n) {
        const double w = PW(p);
#ifdef COMPACT_STORE
        // the stored integers, and the working coordinate of the general path (exact)
        int c0[3];
        unsigned int q0[3];
        GET_CELL(p, c0, q0);
        const int g0[3] = {c0[0] + OXI, c0[1] + OYI, c0[2] + OZI};  // cell in the box
        const double w0[3] = {(double)q0[0] * ITWO32, (double)q0[1] * ITWO32, (double)q0[2] * ITWO32};
        double pos0[3] = {(double)g0[0] + w0[0], (double)g0[1] + w0[1], (double)g0[2] + w0[2]};
#else
        double pos0[3] = {QX(p), QY(p), QZ(p)};
#endif
        double u[3] = {PUX(p), PUY(p), PUZ(p)};
#ifdef MIXED_GATHER
        YeeWeightsF yw;
#ifdef COMPACT_STORE
        yee_weights_cell(c0, q0, yw);
#else
        yee_weights_global(pos0, yw);
#endif
        double ex = h * (double)interp_f(Ex, yw, 1, 0, 0);
        double ey = h * (double)interp_f(Ey, yw, 0, 1, 0);
        double ez = h * (double)interp_f(Ez, yw, 0, 0, 1);
        double bx = h * (double)interp_f(Bx, yw, 0, 1, 1);
        double by = h * (double)interp_f(By, yw, 1, 0, 1);
        double bz = h * (double)interp_f(Bz, yw, 1, 1, 0);
#else
        const double X = XIDX(pos0[0]), Y = (pos0[1] - OY) * (1.0 / DY),
                     Z = (pos0[2] - OZ) * (1.0 / DZ);
        YeeWeights yw;
        yee_weights(X, Y, Z, yw);
        double ex = h * interp_w(Ex, yw, 1, 0, 0);
        double ey = h * interp_w(Ey, yw, 0, 1, 0);
        double ez = h * interp_w(Ez, yw, 0, 0, 1);
        double bx = h * interp_w(Bx, yw, 0, 1, 1);
        double by = h * interp_w(By, yw, 1, 0, 1);
        double bz = h * interp_w(Bz, yw, 1, 1, 0);
#endif
        // --- relativistic Buneman-Boris (one division per reciprocal: FP64 divisions are the
        //     most expensive operations of the kernel) ---
        double umx = u[0] + ex, umy = u[1] + ey, umz = u[2] + ez;  // half electric kick
        const double igm = 1.0 / sqrt(1.0 + umx * umx + umy * umy + umz * umz);
        double tx = bx * igm, ty = by * igm, tz = bz * igm;
        double t2 = tx * tx + ty * ty + tz * tz;
        const double sf = 2.0 / (1.0 + t2);
        double sx = sf * tx, sy = sf * ty, sz = sf * tz;
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
#ifdef COMPACT_STORE
        // The end point, put on the lattice once: s1 = its distance from the lower face of the cell of
        // the start in units of 2^-32 cell. It depends on the stored fraction and on u only, not on
        // the cell, so not on the rank layout. A step is shorter than a cell: fl = -1, 0 or 1.
        const double rg = 1.0 / g1;
        const double dtc[3] = {rg * (DT / DX), rg * (DT / DY), rg * (DT / DZ)};
        long long s1[3];
        int fl[3];
        double pos1[3];
#pragma unroll
        for (int a = 0; a < 3; a++) {
            s1[a] = __double2ll_rn((w0[a] + dtc[a] * u[a]) * TWO32);
            fl[a] = (int)(s1[a] >> 32);
            pos1[a] = (double)g0[a] + (double)s1[a] * ITWO32;
        }
#else
        const double dtg = DT / g1;
        double pos1[3] = {pos0[0] + dtg * u[0], pos0[1] + dtg * u[1], pos0[2] + dtg * u[2]};
#endif
        if (do_diag) {
            const double g0 = sqrt(1.0 + PUX(p) * PUX(p) + PUY(p) * PUY(p) + PUZ(p) * PUZ(p));
            ke = w * m * (0.5 * (g0 + g1) - 1.0);
            mx = 0.5 * w * m * (PUX(p) + u[0]);
            my = 0.5 * w * m * (PUY(p) + u[1]);
            mz = 0.5 * w * m * (PUZ(p) + u[2]);
        }
        // --- reflecting walls: split the path at the wall so the deposited current
        //     follows the real (folded) trajectory and charge stays conserved. A paired face
        //     returns its neutraliser the same way, with a velocity of the inflow distribution.
        //     Three hits are enough: a reflecting wall is hit once per axis and step, and a returned
        //     neutraliser meets no second paired face (they lie on one axis, injection.check_config) ---
        const double Lg[3] = WBOX;
        const int b0[3] = {BCX0, BCY0, BCZ0};
        const int b1[3] = {BCX1, BCY1, BCZ1};
        // deposit charge: q for plasma species, 0 for test species (they feel the fields only)
        const double qw = qdep * w;
#ifdef COMPACT_STORE
        // --- integer route (PIC-D19): a step that ends inside the faces of the box that are not
        //     periodic, and outside the body, meets no wall, no cut and no removal. Its deposit, wrap
        //     and store are done on the stored integers; they give what the general path below gives
        //     for the same step (the shapes of the two ends are the same doubles). The choice depends
        //     on the cell in the box only, not on the rank layout ---
        const int gn[3] = {GNX, GNY, GNZ};
        bool inner = true;
#pragma unroll
        for (int a = 0; a < 3; a++) {
            if (b0[a] != BC_PERIODIC && (g0[a] + fl[a] < 0 || g0[a] + fl[a] >= gn[a])) inner = false;
        }
#ifdef BODY_SPHERE
        {
            const double rx = WPHYS(pos1[0], DX) - BODY_CX, ry = WPHYS(pos1[1], DY) - BODY_CY,
                         rz = WPHYS(pos1[2], DZ) - BODY_CZ;
            if (rx * rx + ry * ry + rz * rz < BODY_R2) inner = false;
        }
#endif
        if (inner) {
            const unsigned int q1[3] = {(unsigned int)s1[0], (unsigned int)s1[1], (unsigned int)s1[2]};
            if (qw != 0.0) {
                AxisShape x, y, z;
                axis_shape_cells(c0[0], w0[0], c0[0] + fl[0], (double)q1[0] * ITWO32, x);
                axis_shape_cells(c0[1], w0[1], c0[1] + fl[1], (double)q1[1] * ITWO32, y);
                axis_shape_cells(c0[2], w0[2], c0[2] + fl[2], (double)q1[2] * ITWO32, z);
                esirkepov_shapes(x, y, z, s1[0] != (long long)q0[0], s1[1] != (long long)q0[1],
                                 s1[2] != (long long)q0[2], qw, Jx, Jy, Jz);
            }
            // periodic wrap on non-decomposed axes, change of rank on decomposed ones (as below)
            const int dec[3] = {DECX, DECY, DECZ};
            const int slo[3] = {SX0I, SY0I, SZ0I};
            const int shi[3] = {SX1I, SY1I, SZ1I};
            int c1[3];
            bool migrate = false;
#pragma unroll
            for (int a = 0; a < 3; a++) {
                const int g = g0[a] + fl[a];
                c1[a] = c0[a] + fl[a];
                if (!dec[a]) {
                    if (b0[a] == BC_PERIODIC) {
                        if (g < 0) c1[a] += gn[a];
                        else if (g >= gn[a]) c1[a] -= gn[a];
                    }
                } else if (g < slo[a] || g >= shi[a]) {
                    migrate = true;
                }
            }
            if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
            SET_CELL(p, c1, q1);
        } else {
#endif
        double trem = 1.0;  // fraction of the step left after the wall hits so far
        const int fm = fmask & (PAIRED_FACES | (PAIRED_FACES << 8));
        for (int it = 0; it < 3; it++) {
            double tmin = 2.0, wall = 0.0;
            int amin = -1;
            for (int a = 0; a < 3; a++) {
                double d = pos1[a] - pos0[a];
                if ((b0[a] == BC_REFLECT || ((fm >> (2 * a)) & 1)) && pos1[a] < FLO(a) && d < 0.0) {
                    double tt = (FLO(a) - pos0[a]) / d;
                    if (tt < tmin) { tmin = tt; amin = a; wall = FLO(a); }
                }
                if (((b1[a] == BC_REFLECT && pos1[a] > FHI(a)) || (((fm >> (2 * a + 1)) & 1) && pos1[a] >= FHI(a))) &&
                    d > 0.0) {
                    double tt = (FHI(a) - pos0[a]) / d;
                    if (tt < tmin) { tmin = tt; amin = a; wall = FHI(a); }
                }
            }
            if (amin < 0) break;
            tmin = fmin(fmax(tmin, 0.0), 1.0);
            double pw[3];
            for (int a = 0; a < 3; a++) pw[a] = pos0[a] + tmin * (pos1[a] - pos0[a]);
            pw[amin] = wall;
            if (qw != 0.0) deposit_cut(pos0, pw, qw, Jx, Jy, Jz);
            trem *= 1.0 - tmin;
            const int face = 2 * amin + (wall > 0.0 ? 1 : 0);
            if ((fm >> face) & 1) {
                unsigned long long sr = hash_key(hash_key(seed ^ RET_SALT, (unsigned long long)step, 0ULL),
                                                 (unsigned long long)PID(p), (unsigned long long)it);
                lke += w * m * (g1 - 1.0);
                if (nflu > 0 && face == flface) {  // returned with the fluctuation at the wall point (PIC-D26)
                    double dd[3];
                    fluct_du(flu, nflu, pw, tnow, dd);
                    sample_inflow(rb + (long long)face * nrb * 6, nrb, amin, wall > 0.0 ? -1.0 : 1.0, sr, u, dd);
                } else {
                    sample_inflow(rb + (long long)face * nrb * 6, nrb, amin, wall > 0.0 ? -1.0 : 1.0, sr, u, nullptr);
                }
                g1 = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
                lke -= w * m * (g1 - 1.0);
                for (int a = 0; a < 3; a++) pos1[a] = WSNAP(pw[a] + (trem * WDT(a) / g1) * u[a]);
                if (wall > 0.0 && pos1[amin] >= wall) pos1[amin] = WBELOW(wall);
            } else {
                pos1[amin] = 2.0 * wall - pos1[amin];
                u[amin] = -u[amin];
            }
            for (int a = 0; a < 3; a++) pos0[a] = pw[a];
        }
        if (qw != 0.0) deposit_cut(pos0, pos1, qw, Jx, Jy, Jz);

        // --- open boundaries: particle leaves and is removed ---
        bool dead = false;
        for (int a = 0; a < 3; a++) {
            if ((b0[a] == BC_OPEN && pos1[a] < FLO(a)) || (b1[a] == BC_OPEN && pos1[a] >= FHI(a))) dead = true;
        }
#ifdef BODY_SPHERE
        // --- absorbing body (PIC-D14): a particle that ends its step inside is removed; its whole
        //     path was deposited, so its charge stays at the nodes of the body ---
        const double rbx = WPHYS(pos1[0], DX) - BODY_CX, rby = WPHYS(pos1[1], DY) - BODY_CY,
                     rbz = WPHYS(pos1[2], DZ) - BODY_CZ;
        if (!dead && rbx * rbx + rby * rby + rbz * rbz < BODY_R2) {
            if (do_diag) {
                ke = 0.0; mx = my = mz = 0.0;
            }
            ake = w * m * (g1 - 1.0);
            aq = qw;  // zero for test species
            an = 1.0;
            SET_W(p, 0.0);
            { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        } else
#endif
        if (dead) {
            if (do_diag) {
                ke = 0.0; mx = my = mz = 0.0;
            }
            lke += w * m * (g1 - 1.0);
            lq = qw;  // charge lost: zero for test species
            SET_W(p, 0.0);
            { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
            if (fm & 0x3f00) {  // a carrier that leaves through a paired face: the column owes one
                const double D[3] = WCELL;
                for (int a = 0; a < 3; a++) {
                    const int face = pos1[a] < 0.0 ? 2 * a : (pos1[a] >= Lg[a] ? 2 * a + 1 : -1);
                    if (face < 0 || !((fm >> (8 + face)) & 1)) continue;
                    int gc[3];
                    for (int b = 0; b < 3; b++) {
                        gc[b] = (int)floor(pos0[b] / D[b]);
                        gc[b] = min(max(gc[b], binfo[6 + b]), binfo[6 + b] + binfo[9 + b] - 1);
                    }
                    const int col = bank_slot(binfo, face, gc);
                    if (col >= 0) atomicAdd(&bank[col], fabs(qw));
                    break;
                }
            }
        } else {
            // periodic wrap on non-decomposed axes (a step moves a particle by less than a cell,
            // so by at most one period); decomposed axes are wrapped after migration
            const int dec[3] = {DECX, DECY, DECZ};
            const double slo[3] = WLO;
            const double shi[3] = WHI;
            bool migrate = false;
            for (int a = 0; a < 3; a++) {
                if (!dec[a]) {
                    if (b0[a] == BC_PERIODIC) {
                        if (pos1[a] < 0.0) pos1[a] += Lg[a];
                        else if (pos1[a] >= Lg[a]) pos1[a] -= Lg[a];
                    }
                } else if (pos1[a] < slo[a] || pos1[a] >= shi[a]) {
                    migrate = true;
                }
            }
            if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
        }
        SET_POS(p, pos1[0], pos1[1], pos1[2]);
#ifdef COMPACT_STORE
        }  // general path
#endif
        SET_MOM(p, u[0], u[1], u[2]);
    }
    if (do_diag) {
        warp_atomic_add(&diag[0], ke);
        warp_atomic_add(&diag[1], mx);
        warp_atomic_add(&diag[2], my);
        warp_atomic_add(&diag[3], mz);
    }
    warp_atomic_add(&diag[4], lke);
    warp_atomic_add(&diag[5], lq);
#ifdef BODY_SPHERE
    warp_atomic_add(&diag[9], ake);
    warp_atomic_add(&diag[10], aq);
    warp_atomic_add(&diag[11], an);
#endif
}

// Set-up of a body (PIC-D14): the particles loaded inside it get weight zero. The beams of all
// species share their positions, so what is left is as neutral as what was loaded.
extern "C" __global__ void mark_body(PSTORE, const int n, const long long cap) {
#ifdef BODY_SPHERE
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    const double rx = PX(p) - BODY_CX, ry = PY(p) - BODY_CY, rz = PZ(p) - BODY_CZ;
    if (rx * rx + ry * ry + rz * rz < BODY_R2) SET_W(p, 0.0);
#endif
}

// The particles of weight zero, into the list that migrate() takes out.
extern "C" __global__ void flag_dead(PSTORE_RO, const int n, const long long cap, int* flag,
                                     int* list, const int listcap) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    if (PW(p) == 0.0) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
}

#ifdef MOVING_WINDOW
// ---- Moving window (PIC-D22). window_state: the shift in cells, the two faces of the window along
//      x and the two bounds of this rank along x, in laboratory coordinates (the doubles that the
//      migration gets).
extern "C" __global__ void window_set(const double s, const double lo, const double hi, const double r0,
                                      const double r1) {
    if (blockIdx.x || threadIdx.x) return;
    window_state[0] = s;
    window_state[1] = lo;
    window_state[2] = hi;
    window_state[3] = r0;
    window_state[4] = r1;
}

// After a shift: the list that migrate() takes out, from the particles alone (the lists of the
// step are void, the bounds have moved). Dead slots; the particles below the new x- face, which
// leave with their plane: kinetic energy, charge and number into wd[0..2], weight zero, nothing
// deposited (their charge stays on the nodes that the Gauss check leaves out at an open face,
// PIC-D15); the particles outside this rank. wd[3] counts the charged particles at x >= xnear,
// within the reach of the x+ wall, where the rule for the field of the new wall cell does not hold.
extern "C" __global__ void window_drop(PSTORE, const int n, const long long cap, const double m,
                                       const double qdep, const double xnear, double* wd, int* flag, int* list,
                                       const int listcap) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    double ke = 0.0, q = 0.0, gone = 0.0, near = 0.0;
    if (p < n) {
        const double w = PW(p);
        bool take = w == 0.0;
        if (!take) {
            const double pos[3] = {PX(p), PY(p), PZ(p)};
            if (pos[0] < FLO(0)) {
                const double ux = PUX(p), uy = PUY(p), uz = PUZ(p);
                ke = w * m * (sqrt(1.0 + ux * ux + uy * uy + uz * uz) - 1.0);
                q = qdep * w;
                gone = 1.0;
                SET_W(p, 0.0);
                take = true;
            } else {
                if (qdep != 0.0 && pos[0] >= xnear) near = 1.0;
                const int dec[3] = {DECX, DECY, DECZ};
                const double slo[3] = WLO;
                const double shi[3] = WHI;
                for (int a = 0; a < 3; a++) {
                    if (dec[a] && (pos[a] < slo[a] || pos[a] >= shi[a])) take = true;
                }
            }
        }
        if (take) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
    }
    warp_atomic_add(&wd[0], ke);
    warp_atomic_add(&wd[1], q);
    warp_atomic_add(&wd[2], gone);
    warp_atomic_add(&wd[3], near);
}

// A plane that enters the window is loaded as it was at the shift before. One round of its move
// to the present: the particles [n0, n) go on for the time tau (less than a cell on every axis)
// and are reflected at the upper face xw of the plane, as the plane ahead, which is loaded later
// in the same state, holds them back. The current of the move is deposited like that of a step
// (the caller sets E -= DT J after the last round); the walls and open faces across the window
// act as in the push, with the lost energy and charge in diag[4..5]; a particle that crosses a
// rank bound is listed for migrate(), which the caller runs after every round. first: the kinetic
// energy of the plane into wd[4] (a particle lost in the same round is in both sums), and the
// slots that a density profile left without weight are listed, so that migrate() drops them.
// tau = 0: nothing moves (a cold plane is loaded where it stands).
extern "C" __global__ void window_move(PSTORE, const int n0, const int n, const long long cap, const double tau,
                                       const double xw, double* Jx, double* Jy, double* Jz, const double m,
                                       const double qdep, const int first, double* wd, double* diag, int* flag,
                                       int* list, const int listcap) {
    int p = n0 + blockIdx.x * blockDim.x + threadIdx.x;
    double kin = 0.0, lke = 0.0, lq = 0.0;
    if (first && p < n && PW(p) == 0.0) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
    if (p < n && PW(p) != 0.0) {
        const double w = PW(p), qw = qdep * w;
        double u[3] = {PUX(p), PUY(p), PUZ(p)};
        const double g = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (first) kin = w * m * (g - 1.0);
        if (tau != 0.0) {
            const double tg = tau / g;
            double pos0[3] = {PX(p), PY(p), PZ(p)};
            double pos1[3] = {pos0[0] + tg * u[0], pos0[1] + tg * u[1], pos0[2] + tg * u[2]};
            const double Lg[3] = WBOX;
            const int b0[3] = {BCX0, BCY0, BCZ0};
            const int b1[3] = {BCX1, BCY1, BCZ1};
            // the path is split at a wall as in the push; each of the three walls is hit once at most
            for (int it = 0; it < 3; it++) {
                double tmin = 2.0, wall = 0.0;
                int amin = -1;
                for (int a = 0; a < 3; a++) {
                    const double d = pos1[a] - pos0[a];
                    const double hi = a == 0 ? xw : Lg[a];
                    if (a > 0 && b0[a] == BC_REFLECT && pos1[a] < 0.0 && d < 0.0) {
                        const double tt = -pos0[a] / d;
                        if (tt < tmin) { tmin = tt; amin = a; wall = 0.0; }
                    }
                    if ((a == 0 || b1[a] == BC_REFLECT) && pos1[a] > hi && d > 0.0) {
                        const double tt = (hi - pos0[a]) / d;
                        if (tt < tmin) { tmin = tt; amin = a; wall = hi; }
                    }
                }
                if (amin < 0) break;
                tmin = fmin(fmax(tmin, 0.0), 1.0);
                double pw[3];
                for (int a = 0; a < 3; a++) pw[a] = pos0[a] + tmin * (pos1[a] - pos0[a]);
                pw[amin] = wall;
                if (qw != 0.0) deposit_cut(pos0, pw, qw, Jx, Jy, Jz);
                pos1[amin] = 2.0 * wall - pos1[amin];
                u[amin] = -u[amin];
                for (int a = 0; a < 3; a++) pos0[a] = pw[a];
            }
            if (qw != 0.0) deposit_cut(pos0, pos1, qw, Jx, Jy, Jz);
            bool dead = false;
            for (int a = 1; a < 3; a++) {
                if ((b0[a] == BC_OPEN && pos1[a] < 0.0) || (b1[a] == BC_OPEN && pos1[a] >= Lg[a])) dead = true;
            }
            if (dead) {
                lke = w * m * (g - 1.0);
                lq = qw;
                SET_W(p, 0.0);
                { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
            } else {
                const int dec[3] = {DECX, DECY, DECZ};
                const double slo[3] = WLO;
                const double shi[3] = WHI;
                bool migrate = false;
                for (int a = 0; a < 3; a++) {
                    if (!dec[a]) {
                        if (b0[a] == BC_PERIODIC) {
                            if (pos1[a] < 0.0) pos1[a] += Lg[a];
                            else if (pos1[a] >= Lg[a]) pos1[a] -= Lg[a];
                        }
                    } else if (pos1[a] < slo[a] || pos1[a] >= shi[a]) {
                        migrate = true;
                    }
                }
                if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
            }
            SET_POS(p, pos1[0], pos1[1], pos1[2]);
            SET_MOM(p, u[0], u[1], u[2]);
        }
    }
    warp_atomic_add(&wd[4], kin);
    warp_atomic_add(&diag[4], lke);
    warp_atomic_add(&diag[5], lq);
}
#endif

// ---- Initial perturbation (PIC-D9): virtual, charge-conserving displacement. The remaining
//      displacement of each particle is stored in rows 8..10 (present only during set-up). One
//      call moves every particle by at most h_a (< 1 cell) per axis and deposits the current of
//      that move with Esirkepov; the caller sets E -= DT * J after all rounds, so div E = rho.
extern "C" __global__ void displace_deposit(PSTORE PDISP_ARG, const int n, const long long cap, const double hx,
                                            const double hy, const double hz, double* Jx, double* Jy, double* Jz,
                                            const double qdep, int* flag, int* list, const int listcap) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    const double h[3] = {hx, hy, hz};
    double rem[3], s[3];
    bool moved = false;
    for (int a = 0; a < 3; a++) {
        rem[a] = DISP(a, p);
        s[a] = fmax(-h[a], fmin(h[a], rem[a]));
        moved |= (s[a] != 0.0);
    }
    if (!moved) return;
    double p0[3] = {QX(p), QY(p), QZ(p)};
    double p1[3] = {WSNAP(p0[0] + WLEN(s[0], DX)), WSNAP(p0[1] + WLEN(s[1], DY)), WSNAP(p0[2] + WLEN(s[2], DZ))};
    if (qdep != 0.0) deposit_global(p0, p1, qdep * PW(p), Jx, Jy, Jz);
    for (int a = 0; a < 3; a++) DISP(a, p) = rem[a] - s[a];
    // periodic wrap on non-decomposed axes; decomposed axes are wrapped after migration
    const double Lg[3] = WBOX;
    const int dec[3] = {DECX, DECY, DECZ};
    const double slo[3] = WLO;
    const double shi[3] = WHI;
    bool migrate = false;
    for (int a = 0; a < 3; a++) {
        if (!dec[a]) {
            p1[a] -= Lg[a] * floor(p1[a] / Lg[a]);
        } else if (p1[a] < slo[a] || p1[a] >= shi[a]) {
            migrate = true;
        }
    }
    if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = p; }
    SET_POS(p, p1[0], p1[1], p1[2]);
}

// Boris velocity step of length tau without moving the particle: with tau = -DT/2 it turns
// momenta given at t = 0 into the u^{-1/2} of the explicit leapfrog when E or B is non-zero.
extern "C" __global__ void boris_kick(PSTORE, const int n, const long long cap, const double* __restrict__ Ex,
                                      const double* __restrict__ Ey, const double* __restrict__ Ez,
                                      const double* __restrict__ Bx, const double* __restrict__ By,
                                      const double* __restrict__ Bz, const double q, const double m,
                                      const double tau) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    const double X = XI(p), Y = YI(p), Z = ZI(p);
    const double h = q * tau / (2.0 * m);
    double ex = h * interp(Ex, X, Y, Z, 0.5, 0.0, 0.0);
    double ey = h * interp(Ey, X, Y, Z, 0.0, 0.5, 0.0);
    double ez = h * interp(Ez, X, Y, Z, 0.0, 0.0, 0.5);
    double bx = h * interp(Bx, X, Y, Z, 0.0, 0.5, 0.5);
    double by = h * interp(By, X, Y, Z, 0.5, 0.0, 0.5);
    double bz = h * interp(Bz, X, Y, Z, 0.5, 0.5, 0.0);
    double umx = PUX(p) + ex, umy = PUY(p) + ey, umz = PUZ(p) + ez;
    double gm = sqrt(1.0 + umx * umx + umy * umy + umz * umz);
    double tx = bx / gm, ty = by / gm, tz = bz / gm;
    double t2 = tx * tx + ty * ty + tz * tz;
    double sx = 2.0 * tx / (1.0 + t2), sy = 2.0 * ty / (1.0 + t2), sz = 2.0 * tz / (1.0 + t2);
    double upx = umx + (umy * tz - umz * ty);
    double upy = umy + (umz * tx - umx * tz);
    double upz = umz + (umx * ty - umy * tx);
    SET_MOM(p, umx + (upy * sz - upz * sy) + ex, umy + (upz * sx - upx * sz) + ey,
            umz + (upx * sy - upy * sx) + ez);
}

// CIC charge density at nodes (diagnostic / Gauss check). shell > 0: only the particles within `shell`
// cells of an open face, which is all that the Marder zone needs (PIC-D16).
extern "C" __global__ void deposit_rho(PSTORE_RO, const int n, const long long cap, const double q,
                                       const double shell, double* rho) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    if (shell > 0.0) {
        const bool o0[3] = {BCX0 == BC_OPEN, BCY0 == BC_OPEN, BCZ0 == BC_OPEN};
        const bool o1[3] = {BCX1 == BC_OPEN, BCY1 == BC_OPEN, BCZ1 == BC_OPEN};
        bool in = false;
#ifdef COMPACT_STORE
        // a position is its cell plus a fraction below one: the cell decides (the shell in whole cells)
        const int l = P_c[p], sh = (int)ceil(shell), gn[3] = {GNX, GNY, GNZ};
        const int c[3] = {l / (NY * NZ) + OXI, l / NZ % NY + OYI, l % NZ + OZI};  // cell of the box
        for (int a = 0; a < 3; a++) in = in || (o0[a] && c[a] < sh) || (o1[a] && c[a] >= gn[a] - sh);
#else
        const double x[3] = {PX(p), PY(p), PZ(p)}, Lg[3] = {LX, LY, LZ}, d[3] = {DX, DY, DZ};
        for (int a = 0; a < 3; a++)
            in = in || (o0[a] && x[a] - FLO(a) < shell * d[a]) || (o1[a] && x[a] >= FHI(a) - shell * d[a]);
#endif
        if (!in) return;
    }
    double X = XI(p), Y = YI(p), Z = ZI(p);
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

// CIC moment of a species at nodes (output diagnostics): sum of w * f(u) / cell volume with
// kind 0: f = 1 (number density), 1..3: f = v_x, v_y, v_z (particle flux), 4: f = m (gamma - 1).
extern "C" __global__ void deposit_moment(PSTORE_RO, const int n, const long long cap, const int kind,
                                          const double m, double* out) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    double X = XI(p), Y = YI(p), Z = ZI(p);
    int i = (int)floor(X), j = (int)floor(Y), k = (int)floor(Z);
    double wx = X - i, wy = Y - j, wz = Z - k;
    double f = 1.0;
    if (kind > 0) {
        double ux = PUX(p), uy = PUY(p), uz = PUZ(p);
        double g = sqrt(1.0 + ux * ux + uy * uy + uz * uz);
        f = kind == 1 ? ux / g : kind == 2 ? uy / g : kind == 3 ? uz / g : m * (g - 1.0);
    }
    double c = f * PW(p) * (1.0 / (DX * DY * DZ));
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int d = 0; d < 2; d++) {
                double s = (a ? wx : 1 - wx) * (b ? wy : 1 - wy) * (d ? wz : 1 - wz);
                if (s != 0.0) atomicAdd(&out[IDX(i + a, j + b, k + d)], c * s);
            }
}

// Field the particles gather with external fields (PIC-D12): Yee (difference) field + external
// array (if has_ext) + uniform value, in the type of the gather (float32 with MIXED_GATHER).
extern "C" __global__ void ext_sum(const double* d, const double* ext, const int has_ext, const double uni,
                                   gather_t* __restrict__ out, const long long n) {
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double v = d[i] + uni;
    if (has_ext) v += ext[i];
    out[i] = (gather_t)v;
}

// Fields at the positions of selected particles (track output): out[c * m + t] for
// c = Ex Ey Ez Bx By Bz, interpolated exactly as in the explicit push.
extern "C" __global__ void gather_at_particles(PSTORE_RO, const long long cap,
                                               const long long* __restrict__ sel, const int m,
                                               const double* __restrict__ Ex, const double* __restrict__ Ey,
                                               const double* __restrict__ Ez, const double* __restrict__ Bx,
                                               const double* __restrict__ By, const double* __restrict__ Bz,
                                               double* out) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= m) return;
    long long p = sel[t];
    const double X = XI(p), Y = YI(p), Z = ZI(p);
    out[t] = interp(Ex, X, Y, Z, 0.5, 0.0, 0.0);
    out[m + t] = interp(Ey, X, Y, Z, 0.0, 0.5, 0.0);
    out[2 * m + t] = interp(Ez, X, Y, Z, 0.0, 0.0, 0.5);
    out[3 * m + t] = interp(Bx, X, Y, Z, 0.0, 0.5, 0.5);
    out[4 * m + t] = interp(By, X, Y, Z, 0.5, 0.0, 0.5);
    out[5 * m + t] = interp(Bz, X, Y, Z, 0.5, 0.5, 0.0);
}

// Kinetic energy and momentum of particles at the current time level (used by the implicit
// scheme where u is stored at integer time levels).
extern "C" __global__ void kinetic(PSTORE_RO, const int n, const long long cap, const double m,
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
extern "C" __global__ void cell_index(PSTORE_RO, const int n, const long long cap,
                                      long long* key) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= n) return;
    int i = (int)floor(XI(p)), j = (int)floor(YI(p)), k = (int)floor(ZI(p));
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

// Lorentz boost of a rest-frame sample with the drift ud (flipping method of Zenitani 2015).
__device__ __forceinline__ void boost_drift(double* u, const double udx, const double udy, const double udz,
                                            unsigned long long& sm) {
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
}

extern "C" __global__ void load_particles(PSTORE, const long long cap, const long long offset, const int ppc,
                                          const int nlx, const int nly, const int nlz, const int gx0,
                                          const int gy0, const int gz0, const int gnx, const int gny,
                                          const int gnz, const unsigned long long seed,
                                          const unsigned long long stream, const double w, const double theta,
                                          const double udx, const double udy, const double udz,
                                          const double ur, const double ox, const double oy, const double oz,
                                          const double id_base) {
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
#ifdef LATTICE_X
    // PIC-D23: the slots of a cell sit on LATTICE_X x LATTICE_Y x LATTICE_Z places (their product is ppc)
    double x = WPOINT(gx0 + i, (slot % LATTICE_X + 0.5) / LATTICE_X, DX);
    double y = WPOINT(gy0 + j, (slot / LATTICE_X % LATTICE_Y + 0.5) / LATTICE_Y, DY);
    double z = WPOINT(gz0 + k, (slot / (LATTICE_X * LATTICE_Y) + 0.5) / LATTICE_Z, DZ);
#else
    unsigned long long sp = hash_key(seed, (unsigned long long)gcell, (unsigned long long)slot);
    double x = WRAND(gx0 + i, sp, DX);
    double y = WRAND(gy0 + j, sp, DY);
    double z = WRAND(gz0 + k, sp, DZ);
#endif
    unsigned long long sm = hash_key(seed ^ (stream * 0x2545F4914F6CDD1DULL + 77ULL), (unsigned long long)gcell,
                                     (unsigned long long)slot);
    double u[3];
    sample_mj(theta, sm, u);
    double vx = udx, vy = udy, vz = udz;
    if (ur != 0.0) {  // radial beam (PIC-D17): the drift ur along the unit vector from the origin
        const double rx = WPHYS(x, DX) - ox, ry = WPHYS(y, DY) - oy, rz = WPHYS(z, DZ) - oz;
        const double r = sqrt(rx * rx + ry * ry + rz * rz);
        const double f = r > 0.0 ? ur / r : 0.0;
        vx = f * rx;
        vy = f * ry;
        vz = f * rz;
    }
    boost_drift(u, vx, vy, vz, sm);
    long long p = offset + t;
    SET_POS(p, x, y, z);
    SET_MOM(p, u[0], u[1], u[2]);
    SET_W(p, w);
    // persistent id (PIC-D8): id_base + linear global cell * ppc + slot (exact below 2^53)
    const long long lin = ((long long)(gx0 + i) * gny + (gy0 + j)) * (long long)gnz + (gz0 + k);
    SET_ID(p, id_base + (double)(lin * ppc + slot));
}

// ---- Injection at open faces (PIC-D13) ----
// Velocity of a particle that crosses a face inwards, of a species with nb beams: distributed as
// v_n f(u), v_n = sgn * v[axis] > 0. rb: per beam [cumulative density * vcap, theta, udx, udy, udz,
// vcap], vcap an upper bound of v_n for the beam (pypic.injection.inflow_table). A beam is picked
// with probability density * vcap and its sample accepted with probability v_n / vcap.
__device__ void sample_inflow(const double* rb, const int nb, const int axis, const double sgn,
                              unsigned long long& s, double* u, const double* dd) {
    const double tot = rb[6 * (nb - 1)];
    for (int it = 0; it < 4096; it++) {
        const double r = urand(s) * tot;
        int b = 0;
        while (b < nb - 1 && r >= rb[6 * b]) b++;
        const double* q = rb + 6 * b;
        sample_mj(q[1], s, u);
        if (dd) boost_drift(u, q[2] + dd[0], q[3] + dd[1], q[4] + dd[2], s);
        else boost_drift(u, q[2], q[3], q[4], s);
        const double vn = sgn * u[axis] / sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (vn > 0.0 && urand(s) * q[5] <= vn) return;
    }
    u[axis] = sgn * fabs(u[axis]);  // practically unreachable: an inward velocity all the same
}

// Bank of a paired face: one entry per face column of the rank (a column = a physical cell of the
// two axes of the face). binfo: [0..5] offset of the bank of face f (-1: none), [6..8] first
// physical cell of the rank per axis, [9..11] number of physical cells. gc: global cell.
// Returns -1 for a cell outside the columns of the rank (an edge or corner cell of the layer).
__device__ __forceinline__ int bank_slot(const int* binfo, const int face, const int* gc) {
    const int a = face >> 1, b = a == 0 ? 1 : 0, c = a == 2 ? 1 : 2;
    const int jb = gc[b] - binfo[6 + b], jc = gc[c] - binfo[6 + c];
    if (binfo[face] < 0 || jb < 0 || jb >= binfo[9 + b] || jc < 0 || jc >= binfo[9 + c]) return -1;
    return binfo[face] + jb * binfo[9 + c] + jc;
}

// salt of the random keys of the injected particles: (seed, step, global cell, slot); the keys of
// the loader are (seed, global cell, slot)
static const unsigned long long INJ_SALT = 0x5851F42D4C957F2DULL;

// One sample of the layer outside the box (PIC-D13): the position r0 of (cell gc, slot), shared by
// all species and beams (key: seed, step, cell, slot), the momentum u of the beam row tb, drawn per
// stream as in load_particles, and the end point r1 of its first step. True if the particle ends
// that step inside the box. sm: the random stream of the momentum, left after the sample.
__device__ bool layer_sample(const unsigned long long ks, const unsigned long long gcell, const int* gc, const int slot,
                             const double* tb, double* r0, double* u, double* r1, double& g, unsigned long long& sm,
                             const double* flu, const int nflu, const double tnow) {
    unsigned long long sp = hash_key(ks, gcell, (unsigned long long)slot);
    r0[0] = WRAND(gc[0], sp, DX);
    r0[1] = WRAND(gc[1], sp, DY);
    r0[2] = WRAND(gc[2], sp, DZ);
    sm = hash_key(ks ^ ((unsigned long long)tb[5] * 0x2545F4914F6CDD1DULL + 77ULL), gcell, (unsigned long long)slot);
    sample_mj(tb[1], sm, u);
    if (nflu > 0) {  // the inflow fluctuation at the birth point (PIC-D26)
        double dd[3];
        fluct_du(flu, nflu, r0, tnow, dd);
        boost_drift(u, tb[2] + dd[0], tb[3] + dd[1], tb[4] + dd[2], sm);
    } else {
        boost_drift(u, tb[2], tb[3], tb[4], sm);
    }
    g = sqrt(1.0 + u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    const double Lg[3] = WBOX;
    const int b0[3] = {BCX0, BCY0, BCZ0};
    bool inside = true;
    for (int a = 0; a < 3; a++) {
        r1[a] = WSNAP(r0[a] + (WDT(a) / g) * u[a]);
        if (b0[a] != BC_PERIODIC && (r1[a] < 0.0 || r1[a] >= Lg[a])) inside = false;
    }
    return inside;
}

// Fill a block of cells outside the box (nlx x nly x nlz cells from the global cell gx0, gy0, gz0;
// one thread per cell, slot and beam) with the beams of one species and keep the particles that
// end their first step inside the box (layer_sample). The particle is appended at
// n0 + ninj[0]++ (the caller made room) at the end of its step, and flagged if it belongs to
// another rank. bt: per beam [weight, theta, udx, udy, udz, stream, id offset, partner id offset];
// id = id_step + id offset + (cell of the extended grid (gn + 2)^3) * ppc + slot.
//   paired = 0 (flux face): the current of the whole path, from the cell outside, is deposited.
//   paired = 1: the particle enters together with a particle of the neutraliser species (Pe, qe,
//     ...) born at the same point with the opposite charge: the pair appears inside the box,
//     nothing is deposited. The partner gets its momentum from inject_partners, which must run
//     before anything else reads Pe: pstate / pface (one entry per particle appended to Pe in
//     this step) keep the random stream and the face + 1 of the partner until then. What
//     its face column owed before this launch (bank0, a copy of bank; see push_deposit) is paid
//     by the first particles that enter from the cell, in the order (slot, beam): such a particle
//     enters alone, deposits its path and pays the bank. A thread of a column that owes replays
//     the samples before its own to find its turn, so the result does not depend on the order
//     of the threads.
// diag (per species): [6] kinetic energy, [7] charge, [8] number of the injected particles; with
// do_diag their kinetic energy and momentum are added to [0..3], as they are part of the plasma
// from this step on.
extern "C" __global__ void inject_face(PSTORE, const long long cap, const int n0, int* ninj, double* diag, int* flag,
                                       int* list, const int listcap, double* Jx, double* Jy, double* Jz,
                                       const double q, const double m, const int gx0, const int gy0, const int gz0,
                                       const int nlx, const int nly, const int nlz, const int gnx, const int gny,
                                       const int gnz, const int ppc, const unsigned long long seed,
                                       const long long step, const double id_step, const double* bt,
                                       const int nbeam, const int face, const int paired, PSTORE_AS(Pe),
                                       const long long cape, const int ne0, int* ninje, double* diage, int* flage,
                                       int* liste, const int listcape, const double qe,
                                       unsigned long long* pstate, int* pface, double* bank,
                                       const double* bank0, const int* binfo, const int do_diag,
                                       const double* flu, const int nflu, const double tnow) {
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long per = (long long)ppc * nbeam;
    if (t >= (long long)nlx * nly * nlz * per) return;
    const long long cell = t / per;
    const int turn = (int)(t % per), slot = turn / nbeam;
    const int k = (int)(cell % nlz), j = (int)((cell / nlz) % nly), i = (int)(cell / nlz / nly);
    const int gc[3] = {gx0 + i, gy0 + j, gz0 + k};
    const unsigned long long ks = hash_key(seed ^ INJ_SALT, (unsigned long long)step, 0ULL);
    const unsigned long long gcell =
        ((unsigned long long)(gc[0] + 1) * (gny + 2) + (gc[1] + 1)) * (1ULL << 30) + (gc[2] + 1);
    const long long lin = ((long long)(gc[0] + 1) * (gny + 2) + (gc[1] + 1)) * (long long)(gnz + 2) + (gc[2] + 1);
    const double* tb = bt + 8 * (turn % nbeam);
    double r0[3], u[3], r1[3], g;
    unsigned long long sm;
    if (!layer_sample(ks, gcell, gc, slot, tb, r0, u, r1, g, sm, flu, nflu, tnow)) return;
    const double w = tb[0], qw = q * w;
    bool pair = paired != 0;
    const int column = pair ? bank_slot(binfo, face, gc) : -1;
    if (column >= 0 && bank0[column] > 0.0) {
        double rem = bank0[column];  // what is left for this particle after those before it
        for (int o = 0; o < turn && rem > 0.0; o++) {
            const double* to = bt + 8 * (o % nbeam);
            const double co = fabs(q * to[0]);
            double a0[3], au[3], a1[3], ag;
            unsigned long long as;
            if (rem >= co * (1.0 - 1e-9) && layer_sample(ks, gcell, gc, o / nbeam, to, a0, au, a1, ag, as, flu, nflu,
                                                         tnow))
                rem = fmax(rem - co, 0.0);
        }
        if (rem >= fabs(qw) * (1.0 - 1e-9)) {
            atomicAdd(&bank[column], -fmin(fabs(qw), rem));
            pair = false;
        }
    }
    if (!pair) deposit_cut(r0, r1, qw, Jx, Jy, Jz);
    const double Lg[3] = WBOX;
    const int b0[3] = {BCX0, BCY0, BCZ0};
    const int dec[3] = {DECX, DECY, DECZ};
    const double slo[3] = WLO;
    const double shi[3] = WHI;
    bool migrate = false;
    for (int a = 0; a < 3; a++) {
        if (!dec[a]) {
            if (b0[a] == BC_PERIODIC) {
                if (r1[a] < 0.0) r1[a] += Lg[a];
                else if (r1[a] >= Lg[a]) r1[a] -= Lg[a];
            }
        } else if (r1[a] < slo[a] || r1[a] >= shi[a]) {
            migrate = true;
        }
    }
    const double cellslot = (double)(lin * ppc + slot);
    {
        const long long p = n0 + atomicAdd(ninj, 1);
        SET_POS(p, r1[0], r1[1], r1[2]);
        SET_MOM(p, u[0], u[1], u[2]);
        SET_W(p, w);
        SET_ID(p, id_step + tb[6] + cellslot);
        if (migrate) { int sl = atomicAdd(flag, 1); if (sl < listcap) list[sl] = (int)p; }
    }
    const double ke = w * m * (g - 1.0);
    atomicAdd(&diag[6], ke);
    atomicAdd(&diag[7], qw);
    atomicAdd(&diag[8], 1.0);
    if (do_diag) {
        atomicAdd(&diag[0], ke);
        for (int a = 0; a < 3; a++) atomicAdd(&diag[1 + a], w * m * u[a]);
    }
    if (pair) {
        const int ie = atomicAdd(ninje, 1);
        const long long p = ne0 + ie;
        S_SET_POS(Pe, cape, p, r1[0], r1[1], r1[2]);
        S_SET_MOM(Pe, cape, p, 0.0, 0.0, 0.0);
        S_SET_W(Pe, cape, p, -qw / qe);
        S_SET_ID(Pe, cape, p, id_step + tb[7] + cellslot);
        if (migrate) { int sl = atomicAdd(flage, 1); if (sl < listcape) liste[sl] = (int)p; }
        pstate[ie] = sm;
        pface[ie] = face + 1;
    }
}

// Momenta of the partners that inject_face appended to the neutraliser species at paired faces:
// one thread per particle appended to Pe in this step (ninje[0] of them, from ne0), drawn from the
// neutraliser's inflow distribution at the face (rb, see sample_inflow) with the stream the
// carrier left. A kernel of its own because the sampler is a rejection loop: here every thread of
// a warp runs it. diage: as diag of inject_face.
extern "C" __global__ void inject_partners(PSTORE_AS(Pe), const long long cape, const int ne0, const int* ninje,
                                           double* diage, const double qe, const double me, const double* rb,
                                           const int nrb, const unsigned long long* pstate, int* pface,
                                           const int do_diag, const double* fle, const int nfle, const int flface,
                                           const double tnow) {
    const int ie = blockIdx.x * blockDim.x + threadIdx.x;
    if (ie >= ninje[0] || pface[ie] == 0) return;
    const int face = pface[ie] - 1;
    pface[ie] = 0;
    unsigned long long sm = pstate[ie];
    double ue[3];
    const long long p = ne0 + ie;
    if (nfle > 0 && face == flface) {  // the partner takes the neutraliser's fluctuation at its point (PIC-D26)
        const double r[3] = {S_X(Pe, cape, p), S_Y(Pe, cape, p), S_Z(Pe, cape, p)};
        double dd[3];
        fluct_du(fle, nfle, r, tnow, dd);
        sample_inflow(rb + (long long)face * nrb * 6, nrb, face >> 1, (face & 1) ? -1.0 : 1.0, sm, ue, dd);
    } else {
        sample_inflow(rb + (long long)face * nrb * 6, nrb, face >> 1, (face & 1) ? -1.0 : 1.0, sm, ue, nullptr);
    }
    S_SET_MOM(Pe, cape, p, ue[0], ue[1], ue[2]);
    const double we = S_W(Pe, cape, p);
    const double eke = we * me * (sqrt(1.0 + ue[0] * ue[0] + ue[1] * ue[1] + ue[2] * ue[2]) - 1.0);
    atomicAdd(&diage[6], eke);
    atomicAdd(&diage[7], qe * we);
    atomicAdd(&diage[8], 1.0);
    if (do_diag) {
        atomicAdd(&diage[0], eke);
        for (int a = 0; a < 3; a++) atomicAdd(&diage[1 + a], we * me * ue[a]);
    }
}
