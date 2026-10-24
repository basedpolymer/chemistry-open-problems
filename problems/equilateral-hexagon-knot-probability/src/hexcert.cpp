// hexcert: certified inner/outer counts of knotted equilateral hexagons on (d-box) x (theta grid).
//
// Coordinates (see result.md):  d1 in [0,2], s in [0,1], t in [0,1] parametrise the sorted
// moment-polytope chamber d1 >= d2 >= d3 (d2 = d1(1+s)/2, d3 = d1((1-s)/2 + t s)), and
// theta = (th1, th2, th3) in (0,pi)^3 are the T135 dihedral angles (theta = 0: flap folded onto the
// central triangle, theta = pi: planar convex hexagon, z(v_{2k}) = h_k sin(theta_k) > 0).
//
// For a d-box D and theta level L (n = 2^L cells per angle) we compute, for the three angle pairs
// (th1,th2), (th2,th3), (th1,th3), bitmaps of grid cells on which the six "pierce" conditions
//   A12: e4=[v4,v5] pierces T2=(v1,v2,v3)    B12: e1=[v1,v2] pierces T4=(v3,v4,v5)
//   A23: e6=[v6,v1] pierces T4=(v3,v4,v5)    B23: e3=[v3,v4] pierces T6=(v5,v6,v1)
//   A13: e2=[v2,v3] pierces T6=(v5,v6,v1)    B13: e5=[v5,v6] pierces T2=(v1,v2,v3)
// are certified TRUE on all of D x cell (inner) or not certified FALSE (outer), by interval arithmetic.
// The pattern R = A12 & A23 & A13 and L = B12 & B23 & B13 contain all knotted hexagons of the
// region (upper bound); certified R- (L-) cells are embedded hexagons, and every connected component
// of certified cells is certified knotted by an exact Fox 3-colouring of the xy-diagram at one
// sample point (lower bound).
//
// Input (stdin), one box per line:  id L1 a1 L2 a2 L3 a3 Lth
//   d1 in [a1, a1+1] * 2 / 2^L1,  s in [a2, a2+1] / 2^L2,  t in [a3, a3+1] / 2^L3,  n = 2^Lth.
// Output, one line per box:
//   id Lth outR inR certR outL inL certL ncompR ncompL failR failL unknotR unknotL
// (counts of fine 3D theta cells; cert = cells in components certified knotted by Fox colouring;
//  fail = components where the sample-point certification failed; unknot = components whose
//  sample point was certified UNKNOTTED, which would contradict Calvo's theorem; must be 0).
//
// Build: g++ -O2 -std=c++20 -ffp-contract=off -fno-fast-math -o hexcert hexcert.cpp
#include <bits/stdc++.h>
#include "interval.hpp"
using namespace std;

static int LT = 0;                // trig table level
static vector<double> CLO, CHI, SLO, SHI;
static I PI_I;
static bool USE_G = true;

static void load_table(const char* fn) {
    FILE* f = fopen(fn, "rb");
    if (!f) { fprintf(stderr, "cannot open table %s\n", fn); exit(2); }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    long N1 = sz / 32; long N = N1 - 1;
    LT = 0; while ((1L << LT) < N) LT++;
    if ((1L << LT) != N) { fprintf(stderr, "bad table size\n"); exit(2); }
    CLO.resize(N1); CHI.resize(N1); SLO.resize(N1); SHI.resize(N1);
    for (long k = 0; k < N1; k++) {
        double v[4];
        if (fread(v, 8, 4, f) != 4) { fprintf(stderr, "read error\n"); exit(2); }
        CLO[k] = v[0]; CHI[k] = v[1]; SLO[k] = v[2]; SHI[k] = v[3];
        if (!(CLO[k] <= CHI[k] && SLO[k] <= SHI[k])) { fprintf(stderr, "bad table entry\n"); exit(2); }
    }
    fclose(f);
    PI_I = I(ndn(M_PI), nup(M_PI));
}

// ---------------------------------------------------------------------------------------------
// geometry over a d-box: rows R[v][k][0..3] (homogeneous (1,x,y,z)); for moving vertex v (1,3,5):
// k=0 constant part, k=1 coefficient of cos(theta), k=2 coefficient of sin(theta).
struct Geo { A R[6][3][4]; };

static void build_geo(const A& d1, const A& s, const A& t, Geo& G) {
    A one(1.0), half(0.5), two(2.0), zero(0.0);
    A r2 = (one + s) * half;
    A r3 = (one - s) * half + t * s;
    A d2 = d1 * r2, d3 = d1 * r3;
    A u = s * (one - t);
    A Q2 = t * s * (one - u * u) * (two + t * s);     // = 4 r3^2 - (1 - r2^2 + r3^2)^2 >= 0
    A sq = asqrt(Q2);
    A pp = d1 * (one - r2 * r2 + r3 * r3) * half;      // v5 = (pp, qq, 0)
    A qq = d1 * sq * half;
    A quarter(0.25);
    A h1 = asqrt(one - d1 * d1 * quarter);
    A h2 = asqrt(one - d2 * d2 * quarter);
    A h3 = asqrt(one - d3 * d3 * quarter);
    // in-plane inward unit normals
    A n1x = zero, n1y = one;
    A n2x = -(sq / (two * r2)), n2y = (r3 * r3 - r2 * r2 - one) / (two * r2);
    A n3x = sq / (two * r3), n3y = -((one - r2 * r2 + r3 * r3) / (two * r3));
    // midpoints
    A m1x = d1 * half, m1y = zero;
    A m2x = (d1 + pp) * half, m2y = qq * half;
    A m3x = pp * half, m3y = qq * half;
    for (int v = 0; v < 6; v++) for (int k = 0; k < 3; k++) for (int c = 0; c < 4; c++) G.R[v][k][c] = zero;
    // fixed vertices v1=0 (idx0), v3=(d1,0,0) (idx2), v5=(pp,qq,0) (idx4)
    G.R[0][0][0] = one;
    G.R[2][0][0] = one; G.R[2][0][1] = d1;
    G.R[4][0][0] = one; G.R[4][0][1] = pp; G.R[4][0][2] = qq;
    // moving vertices v2 (idx1), v4 (idx3), v6 (idx5)
    A mx[3] = {m1x, m2x, m3x}, my[3] = {m1y, m2y, m3y}, nx[3] = {n1x, n2x, n3x}, ny[3] = {n1y, n2y, n3y}, hh[3] = {h1, h2, h3};
    for (int k = 0; k < 3; k++) {
        int v = 2 * k + 1;
        G.R[v][0][0] = one; G.R[v][0][1] = mx[k]; G.R[v][0][2] = my[k];
        G.R[v][1][1] = hh[k] * nx[k]; G.R[v][1][2] = hh[k] * ny[k];
        G.R[v][2][3] = hh[k];
    }
}

static A det4(const A* a, const A* b, const A* c, const A* d) {
    A s0 = a[0] * b[1] - a[1] * b[0], s1 = a[0] * b[2] - a[2] * b[0], s2 = a[0] * b[3] - a[3] * b[0];
    A s3 = a[1] * b[2] - a[2] * b[1], s4 = a[1] * b[3] - a[3] * b[1], s5 = a[2] * b[3] - a[3] * b[2];
    A c5 = c[2] * d[3] - c[3] * d[2], c4 = c[1] * d[3] - c[3] * d[1], c3 = c[1] * d[2] - c[2] * d[1];
    A c2 = c[0] * d[3] - c[3] * d[0], c1 = c[0] * d[2] - c[2] * d[0], c0 = c[0] * d[1] - c[1] * d[0];
    return s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
}

// ---------------------------------------------------------------------------------------------
// pair conditions
struct PierceSpec { int p, q, a, b, c; };
// pairs: 0 = (th1,th2) moving (v2,v4) ; 1 = (th2,th3) moving (v4,v6) ; 2 = (th1,th3) moving (v2,v6)
static const int PAIR_MV[3][2] = {{1, 3}, {3, 5}, {1, 5}};
static const PierceSpec SPEC_A[3] = {{3, 4, 0, 1, 2}, {5, 0, 2, 3, 4}, {1, 2, 4, 5, 0}};
static const PierceSpec SPEC_B[3] = {{0, 1, 2, 3, 4}, {2, 3, 4, 5, 0}, {4, 5, 0, 1, 2}};

struct Coef { I M[3][3]; I Mc[3][3]; I G[3][3][3]; I dD[3]; bool useG; };       // F(ti,tj) = sum_kl M[k][l] u_k(ti) u_l(tj), u = (1, cos, sin)

static void det_coef(const Geo& Gb, const Geo& Gc, const double* lo, const double* hi, const double* cen,
                     int pair, int w0, int w1, int w2, int w3, Coef& C) {
    int mi = PAIR_MV[pair][0], mj = PAIR_MV[pair][1];
    int w[4] = {w0, w1, w2, w3};
    for (int m = 0; m < 3; m++) C.dD[m] = I(ndn(lo[m] - cen[m]), nup(hi[m] - cen[m]));
    C.useG = USE_G;
    bool hasI = false, hasJ = false;
    for (int r = 0; r < 4; r++) {
        if (w[r] == mi) hasI = true;
        if (w[r] == mj) hasJ = true;
        if ((w[r] & 1) && w[r] != mi && w[r] != mj) { fprintf(stderr, "determinant involves a third moving vertex\n"); exit(3); }
    }
    for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
        if ((k > 0 && !hasI) || (l > 0 && !hasJ)) { C.M[k][l] = I(0.0); C.Mc[k][l] = I(0.0); for (int m = 0; m < 3; m++) C.G[m][k][l] = I(0.0); continue; }
        const A* rb[4]; const A* rc[4];
        for (int r = 0; r < 4; r++) {
            int kk = (w[r] == mi) ? k : (w[r] == mj) ? l : 0;
            rb[r] = Gb.R[w[r]][kk]; rc[r] = Gc.R[w[r]][kk];
        }
        A vb = det4(rb[0], rb[1], rb[2], rb[3]);
        A vc = det4(rc[0], rc[1], rc[2], rc[3]);
        I cf = vc.v;
        for (int m = 0; m < 3; m++) cf = cf + vb.g[m] * I(ndn(lo[m] - cen[m]), nup(hi[m] - cen[m]));
        I r = meet(vb.v, cf);
        if (r.lo > r.hi) r = vb.v;   // cannot happen for correct enclosures; keep the naive one
        C.M[k][l] = r;
        C.Mc[k][l] = vc.v;
        for (int m = 0; m < 3; m++) C.G[m][k][l] = vb.g[m];
    }
}

struct Cell { I C, S, c, s, H; };   // cos/sin over cell, at centre, half-width
static inline Cell make_cell(int L, int k) {
    Cell z;
    int sh = LT - L;
    long k0 = (long)k << sh, k1 = (long)(k + 1) << sh, kc = ((long)(2 * k + 1)) << (sh - 1);
    long half = 1L << (LT - 1);
    z.C = I(CLO[k1], CHI[k0]);
    if (k1 <= half) z.S = I(SLO[k0], SHI[k1]);
    else if (k0 >= half) z.S = I(SLO[k1], SHI[k0]);
    else z.S = I(min(SLO[k0], SLO[k1]), 1.0);
    z.c = I(CLO[kc], CHI[kc]); z.s = I(SLO[kc], SHI[kc]);
    z.H = PI_I * I(ldexp(1.0, -(L + 1)));
    z.H = I(-z.H.hi, z.H.hi);
    return z;
}

// Per-level cell data for the column angle (theta_j): U = (1,C,S) box, uc = centre, up = derivative box.
struct ColData { I U1, U2, c1, c2, p1, p2, H; };
// Per-level, per-cell data of one determinant for the row angle (theta_i):
//   ac = Mc^T u(centre_i), ap = Mc^T u'(cell_i), bm[m] = G[m]^T U(cell_i), ab = M^T U(cell_i)
struct RowData { I ac[3], ap[3], bm[3][3], ab[3]; };

static inline I dot3(const I* a, I u0, I u1, I u2) { return a[0] * u0 + a[1] * u1 + a[2] * u2; }

struct DetTab {
    const Coef* C; int Lmax;
    vector<vector<RowData>> rows;    // rows[L][k]
    void build(const Coef* C_, int Lmax_, const vector<vector<Cell>>& cells) {
        C = C_; Lmax = Lmax_;
        rows.assign(Lmax + 1, {});
        for (int L = 0; L <= Lmax; L++) {
            int n = 1 << L; rows[L].resize(n);
            for (int k = 0; k < n; k++) {
                const Cell& z = cells[L][k];
                RowData& r = rows[L][k];
                for (int l = 0; l < 3; l++) {
                    r.ac[l] = C->Mc[0][l] + C->Mc[1][l] * z.c + C->Mc[2][l] * z.s;
                    r.ap[l] = C->Mc[1][l] * (-z.S) + C->Mc[2][l] * z.C;
                    r.ab[l] = C->M[0][l] + C->M[1][l] * z.C + C->M[2][l] * z.S;
                    for (int m = 0; m < 3; m++) r.bm[m][l] = C->G[m][0][l] + C->G[m][1][l] * z.C + C->G[m][2][l] * z.S;
                }
            }
        }
    }
    // enclosure of F over D x cell_i x cell_j
    inline I eval(int L, int ki, const ColData& cj, const Cell& ci) const {
        const RowData& r = rows[L][ki];
        I F2 = r.ab[0] + r.ab[1] * cj.U1 + r.ab[2] * cj.U2;
        if (!C->useG) return F2;
        I Fc = r.ac[0] + r.ac[1] * cj.c1 + r.ac[2] * cj.c2;
        I dj = r.ac[1] * cj.p1 + r.ac[2] * cj.p2;
        I di = r.ap[0] + r.ap[1] * cj.U1 + r.ap[2] * cj.U2;
        I F1 = Fc + ci.H * di + cj.H * dj;
        for (int m = 0; m < 3; m++) F1 = F1 + (r.bm[m][0] + r.bm[m][1] * cj.U1 + r.bm[m][2] * cj.U2) * C->dD[m];
        I F = meet(F1, F2);
        if (F.lo > F.hi) return F2;   // cannot happen for correct enclosures
        return F;
    }
};

// pierce status: +1 certified true, -1 certified false, 0 unknown; sg[5] = known signs (0 = unknown), updated
static inline int sgn(I a) { return a.lo > 0 ? 1 : (a.hi < 0 ? -1 : 0); }
static inline int pierce_status(const DetTab* T5, int L, int ki, int kj, const Cell& ci, const ColData& cj, signed char* sg) {
    auto get = [&](int q) { if (!sg[q]) sg[q] = (signed char)sgn(T5[q].eval(L, ki, cj, ci)); return (int)sg[q]; };
    int sP = get(0), sQ = get(1);
    if (sP != 0 && sQ != 0 && sP == sQ) return -1;
    int s1 = get(2), s2 = get(3);
    if (s1 != 0 && s2 != 0 && s1 != s2) return -1;
    int s3 = get(4);
    if (s3 != 0 && ((s1 != 0 && s1 != s3) || (s2 != 0 && s2 != s3))) return -1;
    if (sP != 0 && sQ != 0 && sP != sQ && s1 != 0 && s1 == s2 && s2 == s3) return 1;
    return 0;
}

// ---------------------------------------------------------------------------------------------
struct Bitmap {
    int n, W; vector<uint64_t> w;
    void init(int n_) { n = n_; W = (n + 63) / 64; w.assign((size_t)n * W, 0); }
    inline void set_block(int r0, int r1, int c0, int c1) {   // [r0,r1) x [c0,c1)
        for (int r = r0; r < r1; r++) {
            uint64_t* row = &w[(size_t)r * W];
            for (int c = c0; c < c1;) {
                int wi = c >> 6, b = c & 63;
                int len = min(64 - b, c1 - c);
                uint64_t m = (len == 64) ? ~0ULL : (((1ULL << len) - 1) << b);
                row[wi] |= m; c += len;
            }
        }
    }
    inline const uint64_t* row(int r) const { return &w[(size_t)r * W]; }
    inline bool get(int r, int c) const { return (w[(size_t)r * W + (c >> 6)] >> (c & 63)) & 1; }
};

struct PairMaps { Bitmap Ain, Aout, Bin, Bout; };

struct QuadCtx { const DetTab* TA; const DetTab* TB; const vector<vector<Cell>>* cells; const vector<vector<ColData>>* cols; int Lmax; PairMaps* PM; };
static void quad(const QuadCtx& Q, int L, int ki, int kj, bool needA, bool needB, const signed char* sgA0, const signed char* sgB0) {
    const Cell& ci = (*Q.cells)[L][ki];
    const ColData& cj = (*Q.cols)[L][kj];
    signed char sgA[5], sgB[5];
    memcpy(sgA, sgA0, 5); memcpy(sgB, sgB0, 5);
    int stA = 0, stB = 0;
    if (needA) stA = pierce_status(Q.TA, L, ki, kj, ci, cj, sgA);
    if (needB) stB = pierce_status(Q.TB, L, ki, kj, ci, cj, sgB);
    int sh = Q.Lmax - L;
    int r0 = ki << sh, r1 = (ki + 1) << sh, c0 = kj << sh, c1 = (kj + 1) << sh;
    PairMaps& PM = *Q.PM;
    if (needA && stA == 1) { PM.Ain.set_block(r0, r1, c0, c1); PM.Aout.set_block(r0, r1, c0, c1); }
    if (needB && stB == 1) { PM.Bin.set_block(r0, r1, c0, c1); PM.Bout.set_block(r0, r1, c0, c1); }
    bool uA = needA && stA == 0, uB = needB && stB == 0;
    if (!uA && !uB) return;
    if (L == Q.Lmax) {
        if (uA) PM.Aout.set_block(r0, r1, c0, c1);
        if (uB) PM.Bout.set_block(r0, r1, c0, c1);
        return;
    }
    for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++)
        quad(Q, L + 1, 2 * ki + a, 2 * kj + b, uA, uB, sgA, sgB);
}

// count sum_{a,b} P12[a][b] popcount(P23[b] & P13[a])
static uint64_t trace3(const Bitmap& P12, const Bitmap& P23, const Bitmap& P13) {
    int n = P12.n, W = P12.W; uint64_t tot = 0;
    for (int a = 0; a < n; a++) {
        const uint64_t* r12 = P12.row(a); const uint64_t* r13 = P13.row(a);
        bool any13 = false; for (int k = 0; k < W; k++) if (r13[k]) { any13 = true; break; }
        if (!any13) continue;
        for (int wb = 0; wb < W; wb++) {
            uint64_t m = r12[wb];
            while (m) {
                int b = wb * 64 + __builtin_ctzll(m); m &= m - 1;
                const uint64_t* r23 = P23.row(b);
                for (int k = 0; k < W; k++) tot += __builtin_popcountll(r23[k] & r13[k]);
            }
        }
    }
    return tot;
}

// ---------------------------------------------------------------------------------------------
// Point evaluation of the hexagon (thin intervals) and exact Fox 3-colouring count of the xy-diagram.
struct P3 { I x, y, z; };
static void point_hexagon(const Geo& Gc, int Lc, const int* kc, P3* V) {
    // theta_k = (kc[k] + 1/2) pi / 2^Lc  (cell centres at level Lc)
    for (int v = 0; v < 6; v++) {
        I c(1.0), s(0.0);
        bool mv = (v % 2 == 1);
        if (mv) {
            int k = v / 2;
            long idx = ((long)(2 * kc[k] + 1)) << (LT - Lc - 1);
            c = I(CLO[idx], CHI[idx]); s = I(SLO[idx], SHI[idx]);
        }
        I comp[4];
        for (int q = 0; q < 4; q++) {
            comp[q] = Gc.R[v][0][q].v;
            if (mv) comp[q] = comp[q] + c * Gc.R[v][1][q].v + s * Gc.R[v][2][q].v;
        }
        V[v].x = comp[1]; V[v].y = comp[2]; V[v].z = comp[3];
    }
}
static inline I orient2(const P3& p, const P3& q, const P3& r) {
    return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
}
// returns number of Fox 3-colourings (3 or 9 or ...), or -1 if some sign could not be certified
static int fox3_certified(const P3* V) {
    struct Cr { int over, under; I tov, tun; };
    vector<Cr> cr;
    for (int i = 0; i < 6; i++) for (int j = i + 2; j < 6; j++) {
        if (i == 0 && j == 5) continue;
        const P3 &p = V[i], &q = V[(i + 1) % 6], &r = V[j], &s = V[(j + 1) % 6];
        I o1 = orient2(p, q, r), o2 = orient2(p, q, s), o3 = orient2(r, s, p), o4 = orient2(r, s, q);
        int a1 = sgn(o1), a2 = sgn(o2), a3 = sgn(o3), a4 = sgn(o4);
        if (!a1 || !a2 || !a3 || !a4) return -1;
        if (a1 == a2 || a3 == a4) continue;
        I tt = o3 / (o3 - o4), uu = o1 / (o1 - o2);
        I zi = p.z + tt * (q.z - p.z), zj = r.z + uu * (s.z - r.z);
        int dz = sgn(zi - zj);
        if (!dz) return -1;
        if (dz > 0) cr.push_back({i, j, I(i) + tt, I(j) + uu});
        else cr.push_back({j, i, I(j) + uu, I(i) + tt});
    }
    int m = cr.size();
    if (m == 0) return 3;
    // all crossing parameters (over and under) must be certified distinct
    vector<I> allp;
    for (auto& c : cr) { allp.push_back(c.tov); allp.push_back(c.tun); }
    for (size_t a = 0; a < allp.size(); a++) for (size_t b = a + 1; b < allp.size(); b++)
        if (!(allp[a].hi < allp[b].lo || allp[b].hi < allp[a].lo)) return -1;
    // undercrossing parameters sorted
    vector<double> us; for (auto& c : cr) us.push_back(c.tun.mid());
    sort(us.begin(), us.end());
    auto arc_of = [&](double x) { int idx = 0; for (double u : us) if (u < x) idx++; return ((idx - 1) % m + m) % m; };
    vector<vector<int>> M(m, vector<int>(m, 0));
    for (int k = 0; k < m; k++) {
        double x = cr[k].tun.mid();
        int idx = 0; for (double u : us) if (u < x) idx++;
        int aout = idx % m, ain = ((idx - 1) % m + m) % m;
        int aov = arc_of(cr[k].tov.mid());
        M[k][aov] = (M[k][aov] + 2) % 3;
        M[k][ain] = (M[k][ain] + 2) % 3;
        M[k][aout] = (M[k][aout] + 2) % 3;
    }
    int rank = 0;
    for (int col = 0; col < m && rank < m; col++) {
        int piv = -1;
        for (int r = rank; r < m; r++) if (M[r][col] % 3) { piv = r; break; }
        if (piv < 0) continue;
        swap(M[piv], M[rank]);
        int inv = (M[rank][col] == 1) ? 1 : 2;
        for (int c = 0; c < m; c++) M[rank][c] = (M[rank][c] * inv) % 3;
        for (int r = 0; r < m; r++) if (r != rank && M[r][col]) {
            int f = M[r][col];
            for (int c = 0; c < m; c++) M[r][c] = ((M[r][c] - f * M[rank][c]) % 3 + 3) % 3;
        }
        rank++;
    }
    int cnt = 1; for (int k = 0; k < m - rank; k++) cnt *= 3;
    return cnt;
}

// connected components (face adjacency) of the 3D cell set {(a,b,c): P12[a][b], P23[b][c], P13[a][c]}
// using runs along c; for each component try up to `tries` cells for a certified Fox count.
struct CompResult { uint64_t cert = 0; int ncomp = 0, fail = 0, unknot = 0; };
struct DSU { vector<int> p; int f(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; } void u(int a, int b) { a = f(a); b = f(b); if (a != b) p[a] = b; } };

static CompResult components(const Bitmap& P12, const Bitmap& P23, const Bitmap& P13, const Geo& Gc, int Lth) {
    int n = P12.n, W = P12.W;
    struct Run { int a, b, c0, c1; };
    vector<Run> runs;
    vector<int> rowstart((size_t)n * n + 1, -1);   // index of first run of row (a,b)
    vector<int> rowcount((size_t)n * n, 0);
    for (int a = 0; a < n; a++) {
        const uint64_t* r12 = P12.row(a); const uint64_t* r13 = P13.row(a);
        for (int b = 0; b < n; b++) {
            if (!((r12[b >> 6] >> (b & 63)) & 1)) continue;
            const uint64_t* r23 = P23.row(b);
            int start = runs.size();
            int c = 0;
            while (c < n) {
                // find next set bit
                int wi = c >> 6; uint64_t m = (r23[wi] & r13[wi]) & (~0ULL << (c & 63));
                while (!m) { wi++; if (wi >= W) break; m = r23[wi] & r13[wi]; }
                if (wi >= W) break;
                int c0 = wi * 64 + __builtin_ctzll(m);
                // find next zero bit after c0
                int c1 = c0;
                while (c1 + 1 < n) {
                    int w2 = (c1 + 1) >> 6, b2 = (c1 + 1) & 63;
                    if ((((r23[w2] & r13[w2]) >> b2) & 1)) c1++; else break;
                }
                runs.push_back({a, b, c0, c1});
                c = c1 + 1;
            }
            if ((int)runs.size() > start) { rowstart[(size_t)a * n + b] = start; rowcount[(size_t)a * n + b] = runs.size() - start; }
        }
    }
    CompResult res;
    int R = runs.size();
    if (R == 0) return res;
    DSU d; d.p.resize(R); for (int i = 0; i < R; i++) d.p[i] = i;
    auto link = [&](int i, size_t row) {
        if (rowcount[row] == 0) return;
        int s = rowstart[row], e = s + rowcount[row];
        for (int j = s; j < e; j++) if (runs[j].c0 <= runs[i].c1 && runs[i].c0 <= runs[j].c1) d.u(i, j);
    };
    for (int i = 0; i < R; i++) {
        int a = runs[i].a, b = runs[i].b;
        if (b > 0) link(i, (size_t)a * n + (b - 1));
        if (a > 0) link(i, (size_t)(a - 1) * n + b);
    }
    unordered_map<int, vector<int>> comp;
    for (int i = 0; i < R; i++) comp[d.f(i)].push_back(i);
    res.ncomp = comp.size();
    for (auto& kv : comp) {
        uint64_t cells = 0;
        for (int i : kv.second) cells += runs[i].c1 - runs[i].c0 + 1;
        int verdict = -1;
        int tries = 0;
        for (size_t q = 0; q < kv.second.size() && tries < 8; q += max<size_t>(1, kv.second.size() / 8), tries++) {
            const Run& r = runs[kv.second[q]];
            int kc[3] = {r.a, r.b, (r.c0 + r.c1) / 2};
            P3 V[6]; point_hexagon(Gc, Lth, kc, V);
            int f = fox3_certified(V);
            if (f >= 0) { verdict = f; break; }
        }
        if (verdict == 9) res.cert += cells;
        else if (verdict == 3) res.unknot++;
        else res.fail++;
    }
    return res;
}

// ---------------------------------------------------------------------------------------------
int main(int argc, char** argv) {
    const char* table = argc > 1 ? argv[1] : "trig_table.bin";
    if (argc > 2) USE_G = atoi(argv[2]) != 0;
    load_table(table);
    char line[512];
    while (fgets(line, sizeof line, stdin)) {
        long long id; int L1, L2, L3, Lth; long long a1, a2, a3;
        if (sscanf(line, "%lld %d %lld %d %lld %d %lld %d", &id, &L1, &a1, &L2, &a2, &L3, &a3, &Lth) != 8) continue;
        if (Lth + 1 > LT) { fprintf(stderr, "theta level too large for table\n"); return 2; }
        double lo[3] = {ldexp(2.0 * a1, -L1), ldexp((double)a2, -L2), ldexp((double)a3, -L3)};
        double hi[3] = {ldexp(2.0 * (a1 + 1), -L1), ldexp((double)(a2 + 1), -L2), ldexp((double)(a3 + 1), -L3)};
        double cen[3]; for (int m = 0; m < 3; m++) cen[m] = 0.5 * (lo[m] + hi[m]);   // exact (dyadic)
        A db[3], dc[3];
        for (int m = 0; m < 3; m++) {
            db[m] = A(I(lo[m], hi[m])); db[m].g[m] = I(1.0);
            dc[m] = A(I(cen[m])); dc[m].g[m] = I(1.0);
        }
        Geo Gb, Gc;
        build_geo(db[0], db[1], db[2], Gb);
        build_geo(dc[0], dc[1], dc[2], Gc);
        int n = 1 << Lth;
        vector<vector<Cell>> cells(Lth + 1);
        vector<vector<ColData>> cols(Lth + 1);
        for (int L = 0; L <= Lth; L++) {
            cells[L].resize(1 << L); cols[L].resize(1 << L);
            for (int k = 0; k < (1 << L); k++) {
                Cell z = make_cell(L, k); cells[L][k] = z;
                cols[L][k] = ColData{z.C, z.S, z.c, z.s, -z.S, z.C, z.H};
            }
        }
        PairMaps PM[3];
        for (int pr = 0; pr < 3; pr++) {
            Coef DA[5], DB[5];
            const PierceSpec* sp[2] = {&SPEC_A[pr], &SPEC_B[pr]};
            for (int ab = 0; ab < 2; ab++) {
                const PierceSpec& S = *sp[ab];
                Coef* D = ab ? DB : DA;
                det_coef(Gb, Gc, lo, hi, cen, pr, S.a, S.b, S.c, S.p, D[0]);
                det_coef(Gb, Gc, lo, hi, cen, pr, S.a, S.b, S.c, S.q, D[1]);
                det_coef(Gb, Gc, lo, hi, cen, pr, S.p, S.q, S.a, S.b, D[2]);
                det_coef(Gb, Gc, lo, hi, cen, pr, S.p, S.q, S.b, S.c, D[3]);
                det_coef(Gb, Gc, lo, hi, cen, pr, S.p, S.q, S.c, S.a, D[4]);
            }
            DetTab TA[5], TB[5];
            for (int q = 0; q < 5; q++) { TA[q].build(&DA[q], Lth, cells); TB[q].build(&DB[q], Lth, cells); }
            PM[pr].Ain.init(n); PM[pr].Aout.init(n); PM[pr].Bin.init(n); PM[pr].Bout.init(n);
            QuadCtx Q{TA, TB, &cells, &cols, Lth, &PM[pr]};
            signed char z5[5] = {0, 0, 0, 0, 0};
            quad(Q, 0, 0, 0, true, true, z5, z5);
        }
        uint64_t outR = trace3(PM[0].Aout, PM[1].Aout, PM[2].Aout);
        uint64_t inR = trace3(PM[0].Ain, PM[1].Ain, PM[2].Ain);
        uint64_t outL = trace3(PM[0].Bout, PM[1].Bout, PM[2].Bout);
        uint64_t inL = trace3(PM[0].Bin, PM[1].Bin, PM[2].Bin);
        CompResult cR, cL;
        if (inR) cR = components(PM[0].Ain, PM[1].Ain, PM[2].Ain, Gc, Lth);
        if (inL) cL = components(PM[0].Bin, PM[1].Bin, PM[2].Bin, Gc, Lth);
        printf("%lld %d %llu %llu %llu %llu %llu %llu %d %d %d %d %d %d\n", id, Lth,
               (unsigned long long)outR, (unsigned long long)inR, (unsigned long long)cR.cert,
               (unsigned long long)outL, (unsigned long long)inL, (unsigned long long)cL.cert,
               cR.ncomp, cL.ncomp, cR.fail, cL.fail, cR.unknot, cL.unknot);
        fflush(stdout);
    }
    return 0;
}
