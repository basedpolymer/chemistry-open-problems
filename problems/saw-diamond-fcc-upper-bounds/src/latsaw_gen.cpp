// latsaw_gen.cpp -- forbidden-set automata dominating self-avoiding walks on the DIAMOND and FACE-CENTRED CUBIC
// lattices; optimisation of the child choice; integer Collatz-Wielandt certificate with an explicit child table.
// Adapted from the Z^d generators of problems/cubic-saw-upper-bound (saw_prod.cpp / saw_big.cpp): same states
// ("loop size <= K" walks, Ponitz-Tittmann), same candidate children (self-avoiding paths through the remembered
// set), new lattice geometry and symmetry groups.
//
// Build:  g++ -O3 -march=native -fopenmp -o latsaw_gen latsaw_gen.cpp
// Usage:  latsaw_gen diamond|fcc K outprefix [mode] [maxit]
//   mode = count : enumerate the states only (counts by walk length)
//   mode = pt    : Ponitz-Tittmann child only (longest admissible suffix of the walk)
//   mode = opt   : all candidate children, min-sum value iteration (default)
//   mode = pi    : all candidate children, policy iteration with candidates cached on disk (large K);
//                  extra argument: a certificate of level K-1 or K-2 used as proxy for the initial policy
//   environment LATSAW_MONO = 1 | 2 | 3 restricts the candidate paths (1: only to older remembered vertices,
//                  2: also to the next younger one -- used for fcc K = 11, 3: at most one move to a younger vertex)
//
// GEOMETRY (normalised frame: the current vertex is the origin o).
//   fcc     : vertices {x in Z^3 : x1+x2+x3 even}; the 12 neighbour vectors u_t = (+-1,+-1,0) and permutations.
//   diamond : vertices A = {x in (2Z)^3 : x1+x2+x3 = 0 mod 4} and B = A + (1,1,1); a_t = the four vectors
//             (+-1,+-1,+-1) with product +1; neighbours of an A vertex x: x + a_t; of a B vertex y: y - a_t.
//             In the normalised frame the current vertex is always of type A.
// A state is a self-avoiding walk w = (v_l, ..., v_1, v_0 = o) stored newest-first by step codes t_0..t_{l-1}:
// the step v_{i+1} -> v_i has displacement u_{t_i} (fcc) or (-1)^{i+1} a_{t_i} (diamond); positions
// p_0 = 0, p_{i+1} = p_i - displacement_i.  Forbidden set F(w) = {p_1, ..., p_l}.
// Symmetry: the point group fixing o and the normalised frame (fcc: the 48 signed permutation matrices;
// diamond: the 24 signed permutation matrices preserving {a_t}, = all permutations of the 4 codes).
// Canonical form = lexicographically smallest image of (t_0, t_1, ...) under the group.
// Admissible: l + delta(w) <= K, delta = graph distance from v_0 to v_l avoiding v_1..v_{l-1} (delta = l if l <= 1).
// CHILDREN of (w, step to e = neighbour of o not in F(w)): the visited set seen from the new vertex, re-normalised by
// T(x) = x - e (fcc) or T(x) = e - x (diamond), i.e. Q = {T(p_i)} u {0}; every self-avoiding path from 0 through
// Q, read as a walk ending at the new vertex, is a valid child (its forbidden set maps into F(w) u {o}).
//
// OUTPUT: <outprefix>.lcert (format LATCERT1, see the checker), <outprefix>.hdr; path counts a_n on stdout.
// Validity of the output never depends on how the automaton was found: everything is re-checked by
// latsaw_check.cpp (written separately).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
#include <unistd.h>
using u64 = uint64_t;
using u32 = uint32_t;
using u128 = unsigned __int128;

static int LAT, K, NDIR, BITS, MAXSTEPS;
static const int MAXL = 32;
static int DV[12][3];          // neighbour vectors of the origin in the normalised frame
static int NG;                 // group order
static int GP[48][12];         // action of the group on step codes
static const u32 NONE = 0xffffffffu;

struct P3 { int8_t c[3]; };
static inline bool eqp(const P3 &a, const P3 &b) { return a.c[0] == b.c[0] && a.c[1] == b.c[1] && a.c[2] == b.c[2]; }
// displacement (from v_{i+1} to v_i) of step i with code t
static inline void disp(int i, int t, int *d) {
    int s = (LAT == 0) ? ((i & 1) ? 1 : -1) : 1;
    for (int a = 0; a < 3; a++) d[a] = s * DV[t][a];
}
static inline int code_of_disp(int i, const int *d) {
    int s = (LAT == 0) ? ((i & 1) ? 1 : -1) : 1;
    for (int t = 0; t < NDIR; t++) if (s * DV[t][0] == d[0] && s * DV[t][1] == d[1] && s * DV[t][2] == d[2]) return t;
    return -1;
}
static inline u64 encode(int l, const int *t) {
    u64 c = 0;
    for (int i = 0; i < l; i++) c |= (u64)t[i] << (BITS * i);
    return c | ((u64)l << 57) | (1ull << 63);
}
static inline int decode(u64 c, int *t) {
    int l = (int)((c >> 57) & 63);
    u64 m = (1ull << BITS) - 1;
    for (int i = 0; i < l; i++) t[i] = (int)((c >> (BITS * i)) & m);
    return l;
}
static inline void positions(int l, const int *t, P3 *p) {
    p[0].c[0] = p[0].c[1] = p[0].c[2] = 0;
    for (int i = 0; i < l; i++) { int d[3]; disp(i, t[i], d); for (int a = 0; a < 3; a++) p[i + 1].c[a] = (int8_t)(p[i].c[a] - d[a]); }
}
// lexicographically smallest image under the group
static inline void canon(int l, int *t) {
    int cand[48], nc = NG;
    for (int g = 0; g < NG; g++) cand[g] = g;
    for (int i = 0; i < l; i++) {
        int m = 99;
        for (int k = 0; k < nc; k++) m = std::min(m, GP[cand[k]][t[i]]);
        int nn = 0;
        for (int k = 0; k < nc; k++) if (GP[cand[k]][t[i]] == m) cand[nn++] = cand[k];
        nc = nn; t[i] = m;
        if (nc == 1) { for (int j = i + 1; j < l; j++) t[j] = GP[cand[0]][t[j]]; return; }
    }
}

// ---------------- lattice set-up ------------------------------------------------------------------------
static void setup_lattice() {
    if (LAT == 1) {
        NDIR = 12; BITS = 4;
        int n = 0;
        for (int a = 0; a < 3; a++) for (int b = a + 1; b < 3; b++)
            for (int sa = -1; sa <= 1; sa += 2) for (int sb = -1; sb <= 1; sb += 2) {
                DV[n][0] = DV[n][1] = DV[n][2] = 0; DV[n][a] = sa; DV[n][b] = sb; n++;
            }
    } else {
        NDIR = 4; BITS = 2;
        int A[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
        memcpy(DV, A, sizeof(A));
    }
    MAXSTEPS = 57 / BITS;
    // group: signed permutation matrices g with g(DV) = DV (as a set)
    NG = 0;
    int perm[3] = {0, 1, 2};
    std::vector<std::vector<int>> perms;
    do { perms.push_back({perm[0], perm[1], perm[2]}); } while (std::next_permutation(perm, perm + 3));
    for (auto &pm : perms) for (int sm = 0; sm < 8; sm++) {
        int sg[3]; for (int a = 0; a < 3; a++) sg[a] = (sm >> a & 1) ? -1 : 1;
        int act[12]; bool ok = true;
        for (int t = 0; t < NDIR && ok; t++) {
            int y[3] = {0, 0, 0};
            for (int a = 0; a < 3; a++) y[pm[a]] = sg[a] * DV[t][a];   // (g x)_{pm[a]} = sg[a] x_a
            int f = -1;
            for (int u = 0; u < NDIR; u++) if (DV[u][0] == y[0] && DV[u][1] == y[1] && DV[u][2] == y[2]) f = u;
            if (f < 0) ok = false; else act[t] = f;
        }
        if (ok) { memcpy(GP[NG], act, sizeof(act)); NG++; }
    }
}
static inline bool is_vertex(const int *x) {
    if (LAT == 1) return ((x[0] + x[1] + x[2]) & 1) == 0;
    bool ev = !(x[0] & 1) && !(x[1] & 1) && !(x[2] & 1), od = (x[0] & 1) && (x[1] & 1) && (x[2] & 1);
    int s = ((x[0] + x[1] + x[2]) % 4 + 4) % 4;
    return (ev && s == 0) || (od && s == 3);
}
// neighbour j of vertex x (diamond: depends on the type of x)
static inline void nbr(const int8_t *x, int j, int8_t *y) {
    int s = 1;
    if (LAT == 0 && (x[0] & 1)) s = -1;   // B vertex
    for (int a = 0; a < 3; a++) y[a] = (int8_t)(x[a] + s * DV[j][a]);
}

// ---------------- free graph distance from the origin (table) ------------------------------------------
static int RB, WB;   // box radius / width for distance tables and BFS grids
static std::vector<uint8_t> D0;
static inline int bidx(const int8_t *x) { return ((x[0] + RB) * WB + (x[1] + RB)) * WB + (x[2] + RB); }
static inline bool inbox(const int8_t *x, int margin) { for (int a = 0; a < 3; a++) if (x[a] < -RB + margin || x[a] > RB - margin) return false; return true; }
static void build_D0() {
    RB = 2 * K + 6; WB = 2 * RB + 1;   // large margin: D0 is exact for |x|_inf <= 2K (shortest paths stay close to the box of 0, x)
    D0.assign((size_t)WB * WB * WB, 255);
    std::vector<P3> q; P3 o; o.c[0] = o.c[1] = o.c[2] = 0; q.push_back(o); D0[bidx(o.c)] = 0;
    for (size_t h = 0; h < q.size(); h++) {
        P3 x = q[h]; int dx = D0[bidx(x.c)];
        if (dx >= 250) continue;
        for (int j = 0; j < NDIR; j++) {
            P3 y; nbr(x.c, j, y.c);
            if (!inbox(y.c, 0)) continue;
            int id = bidx(y.c);
            if (D0[id] == 255) { D0[id] = (uint8_t)(dx + 1); q.push_back(y); }
        }
    }
}
// free graph distance between two vertices
static inline int dfree(const P3 &f, const P3 &t) {
    int8_t d[3];
    bool flip = (LAT == 0) && (f.c[0] & 1);   // diamond, f of type B: use the automorphism x -> f - x
    for (int a = 0; a < 3; a++) d[a] = (int8_t)(flip ? f.c[a] - t.c[a] : t.c[a] - f.c[a]);
    if (!inbox(d, 0)) return 255;
    return D0[bidx(d)];
}
static inline int dorig(const P3 &t) { if (!inbox(t.c, 0)) return 255; return D0[bidx(t.c)]; }

// ---------------- exact avoiding distance --------------------------------------------------------------
struct Grid { std::vector<u32> st; u32 stamp = 0; std::vector<P3> q; };
static thread_local Grid G;
// P[0] = A, P[m] = B: true iff there is a lattice path A -> B of length <= lim avoiding P[1..m-1]
static bool avoid_ok(int m, const P3 *P, int lim) {
    if (m <= 1) return true;
    if (dfree(P[0], P[m]) > lim) return false;
    if (G.st.empty()) { G.st.assign(D0.size(), 0); G.q.resize(D0.size()); }
    G.stamp += 2;
    if (G.stamp > 0xfffffff0u) { std::fill(G.st.begin(), G.st.end(), 0); G.stamp = 2; }
    u32 BLK = G.stamp, SEEN = G.stamp + 1;
    for (int i = 1; i < m; i++) G.st[bidx(P[i].c)] = BLK;
    const P3 &B = P[m];
    size_t head = 0, tail = 0, layerEnd;
    G.q[tail++] = P[0]; G.st[bidx(P[0].c)] = SEEN; layerEnd = tail;
    int dist = 0;
    while (head < tail) {
        if (head == layerEnd) { dist++; layerEnd = tail; if (dist >= lim) return false; }
        P3 x = G.q[head++];
        for (int j = 0; j < NDIR; j++) {
            P3 y; nbr(x.c, j, y.c);
            if (eqp(y, B)) return true;
            if (dist + 1 + dfree(y, B) > lim) continue;
            if (!inbox(y.c, 1)) continue;
            int id = bidx(y.c);
            if (G.st[id] == BLK || G.st[id] == SEEN) continue;
            G.st[id] = SEEN; G.q[tail++] = y;
        }
    }
    return false;
}
static inline bool admissible(int m, const P3 *P) {
    if (m == 0) return true;
    if (m + dfree(P[0], P[m]) > K) return false;
    return avoid_ok(m, P, K - m);
}

// ---------------- enumeration of all admissible canonical walks ----------------------------------------
struct EnumCtx { int t[MAXL]; P3 p[MAXL + 1]; std::vector<u64> out; };
static void enum_dfs(EnumCtx &C, int l, u64 stab) {
    if (admissible(l, C.p)) C.out.push_back(encode(l, C.t));
    if (l + 1 > MAXSTEPS || l + 1 > K - 1) return;
    for (int t = 0; t < NDIR; t++) {
        // canonical constraint: for every g in the stabiliser of the prefix, g(t) >= t
        bool ok = true; u64 ns = 0;
        for (u64 s = stab; s; s &= s - 1) { int g = __builtin_ctzll(s); if (GP[g][t] < t) { ok = false; break; } if (GP[g][t] == t) ns |= 1ull << g; }
        if (!ok) continue;
        int d[3]; disp(l, t, d);
        P3 q; for (int a = 0; a < 3; a++) q.c[a] = (int8_t)(C.p[l].c[a] - d[a]);
        bool hit = false;
        for (int i = 0; i <= l; i++) if (eqp(q, C.p[i])) { hit = true; break; }
        if (hit) continue;
        if (l + 1 + dorig(q) > K) continue;   // l + d(o, p_l) is non-decreasing along extensions
        C.t[l] = t; C.p[l + 1] = q;
        enum_dfs(C, l + 1, ns);
    }
}

// ---------------- lookup table ------------------------------------------------------------------------
struct Lookup {
    const std::vector<u64> *codes; std::vector<u32> tab; u64 mask;
    static inline u64 h(u64 x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ull; x ^= x >> 27; x *= 0x81dadef4bc2dd44dull; x ^= x >> 33; return x; }
    void build(const std::vector<u64> &c) {
        codes = &c; u64 n = 1; while (n < c.size() * 3 / 2 + 16) n <<= 1;
        tab.assign(n, NONE); mask = n - 1;
        for (size_t i = 0; i < c.size(); i++) { u64 j = h(c[i]) & mask; while (tab[j] != NONE) j = (j + 1) & mask; tab[j] = (u32)i; }
    }
    inline long find(u64 k) const { u64 j = h(k) & mask; while (tab[j] != NONE) { if ((*codes)[tab[j]] == k) return tab[j]; j = (j + 1) & mask; } return -1; }
};
static const Lookup *GLK = nullptr;
static int MODE;   // 0 count, 1 pt, 2 opt, 3 pi
static thread_local u64 STAT_NODES = 0, STAT_LEAVES = 0;
static u64 TOT_NODES = 0, TOT_LEAVES = 0;
static int MONO = 0;   // env LATSAW_MONO (experiments): 1 = candidate paths visit remembered vertices in decreasing age only,
                       // 2 = also allow a move to the next younger vertex, 3 = at most one move to a younger vertex

// ---------------- children ----------------------------------------------------------------------------
// fills out with the candidate child indices for state (l,t) and direction j; returns false if j is forbidden
static bool children(int l, const int *t, int j, std::vector<u32> &out) {
    out.clear();
    P3 p[MAXL + 2]; positions(l, t, p);
    P3 e; for (int a = 0; a < 3; a++) e.c[a] = (int8_t)DV[j][a];
    for (int i = 1; i <= l; i++) if (eqp(p[i], e)) return false;
    int n = l + 2;
    P3 q[MAXL + 3];
    q[0].c[0] = q[0].c[1] = q[0].c[2] = 0;
    for (int i = 0; i <= l; i++) for (int a = 0; a < 3; a++)
        q[i + 1].c[a] = (int8_t)(LAT == 1 ? p[i].c[a] - e.c[a] : e.c[a] - p[i].c[a]);
    static thread_local std::vector<std::pair<u64, u32>> mc; mc.clear();   // (vertex mask, child index)
    // codes of a path idx[0..m] (idx[0] = 0): steps newest-first, canonicalised; returns deepest listed prefix
    auto best_prefix = [&](int m, const int *idx, u64 &mask) -> long {
        int tt[MAXL + 2];
        for (int i = 0; i < m; i++) {
            int d[3]; for (int a = 0; a < 3; a++) d[a] = q[idx[i]].c[a] - q[idx[i + 1]].c[a];
            tt[i] = code_of_disp(i, d);
            if (tt[i] < 0) { fprintf(stderr, "internal error: not a lattice step\n"); exit(7); }
        }
        canon(m, tt);
        for (int len = std::min(m, MAXSTEPS); len >= 0; len--) {
            long id = GLK->find(encode(len, tt));
            if (id >= 0) { mask = 0; for (int i = 0; i <= len; i++) mask |= 1ull << idx[i]; return id; }
        }
        return -1;
    };
    if (MODE == 1) {
        int idx[MAXL + 3]; for (int i = 0; i < n; i++) idx[i] = i;
        u64 mask; long id = best_prefix(n - 1, idx, mask);
        if (id < 0) { fprintf(stderr, "no child found\n"); exit(7); }
        out.push_back((u32)id);
        return true;
    }
    // adjacency inside Q (two vertices are adjacent iff their difference is a neighbour vector, +-a_t for diamond),
    // through a thread-local position grid: O(|Q| * #neighbours) instead of O(|Q|^2)
    static thread_local std::vector<int8_t> qgrid;
    if (qgrid.empty()) qgrid.assign(D0.size(), -1);
    for (int i = 0; i < n; i++) qgrid[bidx(q[i].c)] = (int8_t)i;
    u64 adj[MAXL + 3];
    for (int i = 0; i < n; i++) {
        adj[i] = 0;
        for (int u = 0; u < NDIR; u++) for (int sg = 1; sg >= (LAT == 0 ? -1 : 1); sg -= 2) {
            int8_t y[3]; for (int a = 0; a < 3; a++) y[a] = (int8_t)(q[i].c[a] + sg * DV[u][a]);
            if (!inbox(y, 0)) continue;
            int k = qgrid[bidx(y)];
            if (k >= 0) adj[i] |= 1ull << k;
        }
    }
    for (int i = 0; i < n; i++) qgrid[bidx(q[i].c)] = -1;
    int dq[MAXL + 3]; for (int i = 0; i < n; i++) dq[i] = dorig(q[i]);
    int path[MAXL + 3]; u64 rem[MAXL + 3]; int nback[MAXL + 3]; nback[0] = 0;
    auto ext = [&](int depth, u64 used) {
        u64 a = adj[path[depth]] & ~used, r = 0;
        if (MONO == 1) a &= ~((2ull << path[depth]) - 1);   // only vertices older than the current one
        else if (MONO == 2) a &= ~((2ull << path[depth]) - 1) | (path[depth] ? (1ull << (path[depth] - 1)) : 0);   // or the next younger one
        else if (MONO == 3 && nback[depth] >= 1) a &= ~((2ull << path[depth]) - 1);   // at most one move to a younger vertex
        while (a) { int k = __builtin_ctzll(a); a &= a - 1; if (depth + 1 + dq[k] <= K && depth + 1 <= MAXSTEPS) r |= 1ull << k; }
        return r;
    };
    // Incremental DFS: along the current path we keep the canonical code of every prefix (lexicographically minimal
    // image; the canonical form of a prefix is the prefix of the canonical form, so it can be extended step by
    // step with the set of group elements still achieving the minimum), its index in the state table (or -1),
    // and the deepest listed prefix so far.  At a leaf (no admissible extension) the deepest listed prefix of
    // the branch is a candidate, exactly as in the non-incremental version (best_prefix).
    static thread_local int rawc[2][MAXL + 3][MAXL + 3];
    for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) if (adj[i] >> k & 1)
        for (int par = 0; par < 2; par++) {
            int d[3]; for (int a = 0; a < 3; a++) d[a] = q[i].c[a] - q[k].c[a];
            rawc[par][i][k] = code_of_disp(par, d);   // -1 when the parity does not fit (diamond); never used then
        }
    u64 gmask[MAXL + 3], ccode[MAXL + 3]; int bestd[MAXL + 3]; long besti[MAXL + 3];
    const u64 ALLG = (NG == 64) ? ~0ull : ((1ull << NG) - 1);
    gmask[0] = ALLG; ccode[0] = 0; bestd[0] = 0; besti[0] = GLK->find(1ull << 63);
    auto leaf = [&](int depth) {
        STAT_LEAVES++;
        if (besti[depth] < 0) { fprintf(stderr, "no child found (leaf)\n"); exit(7); }
        u64 mask = 0; for (int i = 0; i <= bestd[depth]; i++) mask |= 1ull << path[i];
        mc.push_back({mask, (u32)besti[depth]});
    };
    auto push = [&](int d1, int nx) {   // extend the path to depth d1 with vertex nx
        int traw = rawc[(d1 - 1) & 1][path[d1 - 1]][nx];
        if (traw < 0) { fprintf(stderr, "internal error: not a lattice step (push)\n"); exit(7); }
        STAT_NODES++;
        int m = 99; u64 nm = 0;
        for (u64 s = gmask[d1 - 1]; s; s &= s - 1) { int g = __builtin_ctzll(s); int v = GP[g][traw]; if (v < m) { m = v; nm = 0; } if (v == m) nm |= 1ull << g; }
        gmask[d1] = nm; ccode[d1] = ccode[d1 - 1] | ((u64)m << (BITS * (d1 - 1)));
        long id = GLK->find(ccode[d1] | ((u64)d1 << 57) | (1ull << 63));
        if (id >= 0) { bestd[d1] = d1; besti[d1] = id; } else { bestd[d1] = bestd[d1 - 1]; besti[d1] = besti[d1 - 1]; }
    };
    int depth = 0; path[0] = 0; u64 used = 1;
    rem[0] = ext(0, used);
    if (!rem[0]) leaf(0);
    while (depth >= 0) {
        if (!rem[depth]) { used &= ~(1ull << path[depth]); depth--; continue; }
        int nx = __builtin_ctzll(rem[depth]);
        rem[depth] &= rem[depth] - 1;
        depth++; path[depth] = nx; used |= 1ull << nx;
        nback[depth] = nback[depth - 1] + (nx < path[depth - 1] ? 1 : 0);
        push(depth, nx);
        rem[depth] = ext(depth, used);
        if (!rem[depth]) leaf(depth);
    }
    // dominance pruning: drop candidates whose vertex set is strictly contained in another candidate's set
    std::sort(mc.begin(), mc.end());
    mc.erase(std::unique(mc.begin(), mc.end()), mc.end());
    std::stable_sort(mc.begin(), mc.end(), [](const std::pair<u64, u32> &a, const std::pair<u64, u32> &b) {
        return __builtin_popcountll(a.first) > __builtin_popcountll(b.first); });
    static thread_local std::vector<std::pair<u64, u32>> keep; keep.clear();
    for (size_t i = 0; i < mc.size(); i++) {
        bool dom = false;
        for (size_t k = 0; k < keep.size() && !dom; k++)
            if ((mc[i].first & keep[k].first) == mc[i].first && mc[i].first != keep[k].first) dom = true;
        if (!dom) keep.push_back(mc[i]);
    }
    for (auto &x : keep) out.push_back(x.second);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return true;
}

static double wall() { return omp_get_wtime(); }

// ---------------- certificate output ----------------------------------------------------------------------
static void write_cert(const std::string &outp, const std::string &lat, const std::string &mode, const std::vector<u64> &codes,
                       const std::vector<u64> &V, const std::vector<u32> &child, u64 P, u64 Q, u64 bn, u64 bd, double lb, double ub) {
    size_t N = codes.size();
    std::string cp = outp + ".lcert";
    FILE *f = fopen(cp.c_str(), "wb");
    u32 ll = (u32)LAT, kk = (u32)K; u64 NN = N;
    fwrite("LATCERT1", 1, 8, f); fwrite(&ll, 4, 1, f); fwrite(&kk, 4, 1, f); fwrite(&NN, 8, 1, f);
    fwrite(&P, 8, 1, f); fwrite(&Q, 8, 1, f);
    fwrite(codes.data(), 8, N, f); fwrite(V.data(), 8, N, f); fwrite(child.data(), 4, child.size(), f);
    fclose(f);
    std::string hp = outp + ".hdr";
    f = fopen(hp.c_str(), "w");
    fprintf(f, "lattice %s\nK %d\nmode %s\nstates %zu\nlambda_bar %llu/%llu\nmax_ratio %llu/%llu\nfloat_bracket %.10f %.10f\n", lat.c_str(), K,
            mode.c_str(), N, (unsigned long long)P, (unsigned long long)Q, (unsigned long long)bn, (unsigned long long)bd, lb, ub);
    fclose(f);
}

// ---------------- policy iteration with a disk cache of the candidates (large K) ---------------------------
// Pass 0 computes the candidate lists once and appends them, chunk by chunk, to one file per thread; the initial
// policy takes, for every (state, step), the candidate whose longest prefix listed in the proxy certificate (level
// K-1 or K-2, optional) has the smallest proxy weight.  Each later pass: evaluation of the current policy (power
// iteration), then improvement (argmin of the current vector over the cached candidates; ties keep the old child).
static int run_pi(const std::vector<u64> &codes, Lookup &LK, long startIdx, const std::string &outp, const std::string &lat, const char *proxy,
                  int maxpass, double T0) {
    size_t N = codes.size();
    // proxy weights
    std::vector<u64> pcodes, pV; Lookup PL; bool haveProxy = false;
    if (proxy) {
        FILE *f = fopen(proxy, "rb"); char mg[8]; u32 hl[2]; u64 hq[3];
        if (f && fread(mg, 1, 8, f) == 8 && memcmp(mg, "LATCERT1", 8) == 0 && fread(hl, 4, 2, f) == 2 && fread(hq, 8, 3, f) == 3 && (int)hl[0] == LAT) {
            pcodes.resize(hq[0]); pV.resize(hq[0]);
            if (fread(pcodes.data(), 8, hq[0], f) == hq[0] && fread(pV.data(), 8, hq[0], f) == hq[0]) haveProxy = true;
        }
        if (f) fclose(f);
        if (haveProxy) PL.build(pcodes);
        printf("proxy %s: %s (%zu states)\n", proxy, haveProxy ? "loaded" : "NOT loaded", pcodes.size());
    }
    auto proxyw = [&](u64 c) -> u64 {
        int l = (int)((c >> 57) & 63);
        for (int m = l; m >= 0; m--) {
            u64 cm = (c & ((m * BITS >= 57) ? ~0ull : ((1ull << (m * BITS)) - 1))) & ((1ull << 57) - 1);
            cm |= ((u64)m << 57) | (1ull << 63);
            long id = PL.find(cm);
            if (id >= 0) return pV[id];
        }
        return ~0ull;
    };
    std::vector<u32> pol((size_t)N * NDIR, NONE);
    const size_t CH = 1 << 13;
    size_t nch = (N + CH - 1) / CH;
    int nth = omp_get_max_threads();
    std::vector<FILE *> cf(nth);
    std::vector<int> cfd(nth);
    for (int t = 0; t < nth; t++) {
        std::string fn = outp + ".cand" + std::to_string(t);
        cf[t] = fopen(fn.c_str(), "w+b");
        if (!cf[t]) { fprintf(stderr, "cannot open candidate cache\n"); return 3; }
        cfd[t] = fileno(cf[t]);
    }
    std::vector<int> chTh(nch); std::vector<u64> chOff(nch), chLen(nch);
    std::vector<u64> off(nth, 0);
    double t0 = wall();
    long totCand = 0, totGrp = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : totCand, totGrp)
    for (size_t ch = 0; ch < nch; ch++) {
        int tid = omp_get_thread_num();
        std::vector<u32> kids; std::vector<u32> blob; int t[MAXL];
        for (size_t s = ch * CH; s < std::min(N, (ch + 1) * CH); s++) {
            int l = decode(codes[s], t);
            for (int j = 0; j < NDIR; j++) {
                if (!children(l, t, j, kids)) continue;
                blob.push_back(((u32)j << 24) | (u32)kids.size());
                blob.insert(blob.end(), kids.begin(), kids.end());
                u32 best = kids[0]; u64 bw = ~0ull;
                if (haveProxy) for (u32 k : kids) { u64 w = proxyw(codes[k]); if (w < bw) { bw = w; best = k; } }
                pol[s * NDIR + j] = best;
                totCand += kids.size(); totGrp++;
            }
        }
        size_t bytes = blob.size() * 4;
        if (fwrite(blob.data(), 1, bytes, cf[tid]) != bytes) { fprintf(stderr, "cache write error\n"); exit(6); }
        chTh[ch] = tid; chOff[ch] = off[tid]; chLen[ch] = bytes; off[tid] += bytes;
        if (ch % 256 == 0 && tid == 0) { printf("    pass 0: chunk %zu/%zu  %.1fs\n", ch, nch, wall() - t0); fflush(stdout); }
    }
    for (int t = 0; t < nth; t++) fflush(cf[t]);
    std::vector<u32>().swap(LK.tab); GLK = nullptr;   // the state table is not needed any more
    { std::vector<u32>().swap(PL.tab); std::vector<u64>().swap(pcodes); std::vector<u64>().swap(pV); }
    u64 cacheBytes = 0; for (int t = 0; t < nth; t++) cacheBytes += off[t];
    printf("pass 0: candidates %ld in %ld groups (%.2f per group), cache %.2f GB, %.1fs\n", totCand, totGrp, (double)totCand / std::max(1L, totGrp),
           cacheBytes / 1e9, wall() - t0);
    fflush(stdout);
    std::vector<double> v(N, 1.0), nv(N);
    double lb = 0, ub = 0;
    auto evaluate = [&](int maxit) {
        int it;
        for (it = 1; it <= maxit; it++) {
            double mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
            for (size_t s = 0; s < N; s++) {
                double sum = 0;
                for (int j = 0; j < NDIR; j++) { u32 c = pol[s * NDIR + j]; if (c != NONE) sum += v[c]; }
                nv[s] = sum; mx = std::max(mx, sum);
            }
            double hi = 0, lo = 1e300;
#pragma omp parallel for schedule(static) reduction(max : hi) reduction(min : lo)
            for (size_t s = 0; s < N; s++) if (v[s] > 0) { double r = nv[s] / v[s]; hi = std::max(hi, r); if (nv[s] > 0) lo = std::min(lo, r); }
            lb = lo; ub = hi;
#pragma omp parallel for schedule(static)
            for (size_t s = 0; s < N; s++) v[s] = nv[s] / mx;
            if (ub - lb < 2e-11) break;
        }
        return it;
    };
    // improvement step over the cached candidates; w = weight array (double or integer via lambda)
    auto improve = [&](auto weight, double tol) {
        long changed = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : changed)
        for (size_t ch = 0; ch < nch; ch++) {
            std::vector<u32> blob(chLen[ch] / 4);
            if (chLen[ch] && pread(cfd[chTh[ch]], blob.data(), chLen[ch], chOff[ch]) != (ssize_t)chLen[ch]) { fprintf(stderr, "cache read error\n"); exit(6); }
            size_t p = 0;
            for (size_t s = ch * CH; s < std::min(N, (ch + 1) * CH); s++) {
                // groups of state s are consecutive; the number of groups = allowed steps (pol != NONE)
                for (int j = 0; j < NDIR; j++) {
                    if (pol[s * NDIR + j] == NONE) continue;
                    u32 hdr = blob[p++]; int jj = (int)(hdr >> 24); u32 cnt = hdr & 0xffffff;
                    if (jj != j) { fprintf(stderr, "cache desync\n"); exit(6); }
                    u32 cur = pol[s * NDIR + j]; auto bw = weight(cur); u32 best = cur; bool moved = false;
                    for (u32 k = 0; k < cnt; k++) {
                        u32 c = blob[p + k]; auto w = weight(c);
                        // leave the current child unless the improvement exceeds the relative tolerance
                        if (moved ? (w < bw) : ((double)w < (double)bw * (1.0 - tol))) { bw = w; best = c; moved = true; }
                    }
                    p += cnt;
                    if (best != cur) { pol[s * NDIR + j] = best; changed++; }
                }
            }
        }
        return changed;
    };
    for (int pass = 1; pass <= maxpass; pass++) {
        double t1 = wall();
        int its = evaluate(1500);
        printf("pass %d: evaluation %d its %.1fs, lambda(policy) in [%.10f, %.10f]\n", pass, its, wall() - t1, lb, ub); fflush(stdout);
        t1 = wall();
        long changed = improve([&](u32 c) { return v[c]; }, 1e-12);
        printf("pass %d: improvement changed %ld  %.1fs\n", pass, changed, wall() - t1); fflush(stdout);
        if (changed == 0) break;
    }
    int its = evaluate(1500);
    printf("final evaluation %d its, lambda in [%.10f, %.10f]  %.1fs\n", its, lb, ub, wall() - T0); fflush(stdout);
    double vmax = 0; for (size_t s = 0; s < N; s++) vmax = std::max(vmax, v[s]);
    std::vector<double>().swap(nv);   // release before allocating V (memory peak for the largest runs)
    std::vector<u64> V(N);
    for (size_t s = 0; s < N; s++) { double x = std::ceil(v[s] / vmax * 1125899906842624.0); V[s] = (u64)std::max(1.0, x); }
    std::vector<double>().swap(v);
    long ch2 = improve([&](u32 c) { return V[c]; }, 0.0);   // argmin of the integer weights (V <= 2^50 + 1: exact in double)
    printf("integer argmin changed %ld children\n", ch2);
    for (int t = 0; t < nth; t++) { fclose(cf[t]); std::string fn = outp + ".cand" + std::to_string(t); unlink(fn.c_str()); }
    u64 bestNum = 0, bestDen = 1;
#pragma omp parallel
    {
        u64 lNum = 0, lDen = 1;
#pragma omp for schedule(static)
        for (size_t s = 0; s < N; s++) {
            u64 sum = 0;
            for (int j = 0; j < NDIR; j++) { u32 c = pol[s * NDIR + j]; if (c != NONE) sum += V[c]; }
            if ((u128)sum * lDen > (u128)lNum * V[s]) { lNum = sum; lDen = V[s]; }
        }
#pragma omp critical
        { if ((u128)lNum * bestDen > (u128)bestNum * lDen) { bestNum = lNum; bestDen = lDen; } }
    }
    u64 Q = 10000000ull;
    u64 P = (u64)(((u128)bestNum * Q + bestDen - 1) / bestDen);
    printf("integer CW: max ratio = %llu/%llu = %.12f  => certified lambda_bar = %llu/%llu\n", (unsigned long long)bestNum,
           (unsigned long long)bestDen, (double)bestNum / bestDen, (unsigned long long)P, (unsigned long long)Q);
    write_cert(outp, lat, "pi", codes, V, pol, P, Q, bestNum, bestDen, lb, ub);
    (void)startIdx;
    printf("RESULT lattice=%s K=%d mode=pi states=%zu lambda_float=[%.10f,%.10f] certified<=%llu/%llu total %.1fs\n", lat.c_str(), K, N, lb, ub,
           (unsigned long long)P, (unsigned long long)Q, wall() - T0);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s diamond|fcc K outprefix [count|pt|opt|pi] [maxit]\n", argv[0]); return 1; }
    std::string lat = argv[1];
    if (lat == "diamond") LAT = 0; else if (lat == "fcc") LAT = 1; else { fprintf(stderr, "lattice?\n"); return 1; }
    K = atoi(argv[2]); std::string outp = argv[3];
    std::string mode = argc > 4 ? argv[4] : "opt";
    MONO = getenv("LATSAW_MONO") ? atoi(getenv("LATSAW_MONO")) : 0;
    MODE = mode == "count" ? 0 : mode == "pt" ? 1 : mode == "opt" ? 2 : mode == "pi" ? 3 : -1;
    if (MODE < 0) { fprintf(stderr, "mode?\n"); return 1; }
    int MAXIT = argc > 5 ? atoi(argv[5]) : 20000;
    setup_lattice();
    if (K - 1 > MAXSTEPS || K > MAXL - 2) { fprintf(stderr, "K too large for the encoding\n"); return 1; }
    build_D0();
    double T0 = wall();
    printf("lattice=%s K=%d mode=%s ndir=%d group=%d threads=%d\n", lat.c_str(), K, mode.c_str(), NDIR, NG, omp_get_max_threads());
    // sanity: every DV is a vertex, origin's neighbours have distance 1
    for (int j = 0; j < NDIR; j++) if (!is_vertex(DV[j])) { fprintf(stderr, "geometry error\n"); return 3; }
    // ---- enumerate states ----
    std::vector<u64> codes;
    {
        const int S = std::min(5, K - 1);
        struct Item { EnumCtx c; int l; u64 stab; };
        std::vector<Item> layer, seeds;
        Item it0; it0.c.p[0].c[0] = it0.c.p[0].c[1] = it0.c.p[0].c[2] = 0; it0.l = 0; it0.stab = (NG == 64) ? ~0ull : ((1ull << NG) - 1);
        layer.push_back(it0);
        std::vector<u64> shortcodes;
        for (int l = 0; l < S; l++) {
            std::vector<Item> nxt;
            for (auto &it : layer) {
                if (admissible(it.l, it.c.p)) shortcodes.push_back(encode(it.l, it.c.t));
                for (int t = 0; t < NDIR; t++) {
                    bool ok = true; u64 ns = 0;
                    for (u64 s = it.stab; s; s &= s - 1) { int g = __builtin_ctzll(s); if (GP[g][t] < t) { ok = false; break; } if (GP[g][t] == t) ns |= 1ull << g; }
                    if (!ok) continue;
                    int d[3]; disp(it.l, t, d);
                    P3 q; for (int a = 0; a < 3; a++) q.c[a] = (int8_t)(it.c.p[it.l].c[a] - d[a]);
                    bool hit = false;
                    for (int i = 0; i <= it.l; i++) if (eqp(q, it.c.p[i])) { hit = true; break; }
                    if (hit || it.l + 1 + dorig(q) > K) continue;
                    Item n2 = it; n2.c.t[it.l] = t; n2.c.p[it.l + 1] = q; n2.l = it.l + 1; n2.stab = ns;
                    nxt.push_back(n2);
                }
            }
            layer.swap(nxt);
        }
        int nth = omp_get_max_threads();
        std::vector<std::vector<u64>> outs(nth);
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = 0; i < layer.size(); i++) {
            EnumCtx c = layer[i].c; c.out.clear();
            enum_dfs(c, layer[i].l, layer[i].stab);
            auto &o = outs[omp_get_thread_num()]; o.insert(o.end(), c.out.begin(), c.out.end());
        }
        codes = shortcodes;
        for (auto &o : outs) { codes.insert(codes.end(), o.begin(), o.end()); std::vector<u64>().swap(o); }
        std::sort(codes.begin(), codes.end());
        size_t before = codes.size();
        codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
        if (codes.size() != before) { fprintf(stderr, "duplicate states?!\n"); return 3; }
    }
    size_t N = codes.size();
    {
        std::vector<u64> byl(MAXL + 1, 0); for (u64 c : codes) byl[(c >> 57) & 63]++;
        printf("states=%zu  enum %.1fs  by length:", N, wall() - T0);
        for (int l = 0; l <= MAXL; l++) if (byl[l]) printf(" %d:%llu", l, (unsigned long long)byl[l]);
        printf("\n"); fflush(stdout);
    }
    if (MODE == 0) return 0;
    if (N >= NONE) { fprintf(stderr, "too many states\n"); return 3; }
    Lookup LK; LK.build(codes); GLK = &LK;
    long startIdx;
    { int tt[1]; startIdx = LK.find(encode(0, tt)); }
    if (startIdx < 0) { fprintf(stderr, "no start state\n"); return 3; }
    if (MODE == 3) return run_pi(codes, LK, startIdx, outp, lat, argc > 6 ? argv[6] : nullptr, MAXIT, T0);
    // ---- candidate children (CSR) ----
    double T1 = wall();
    std::vector<u64> rowg(N + 1);       // groups (allowed directions) per state
    std::vector<uint8_t> gdir;          // direction of each group
    std::vector<u64> gstart;            // candidate range of each group
    std::vector<u32> cand;
    {
        const size_t CH = 1 << 13;
        size_t nch = (N + CH - 1) / CH;
        std::vector<std::vector<u32>> cg(nch), cc(nch); std::vector<std::vector<uint8_t>> cd(nch);
        std::vector<u32> ng(N);
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t ch = 0; ch < nch; ch++) {
            std::vector<u32> kids; int t[MAXL];
            for (size_t s = ch * CH; s < std::min(N, (ch + 1) * CH); s++) {
                int l = decode(codes[s], t); u32 cnt = 0;
                for (int j = 0; j < NDIR; j++) {
                    if (!children(l, t, j, kids)) continue;
                    cg[ch].push_back((u32)kids.size()); cd[ch].push_back((uint8_t)j);
                    cc[ch].insert(cc[ch].end(), kids.begin(), kids.end()); cnt++;
                }
                ng[s] = cnt;
            }
#pragma omp atomic
            TOT_NODES += STAT_NODES;
#pragma omp atomic
            TOT_LEAVES += STAT_LEAVES;
            STAT_NODES = STAT_LEAVES = 0;
        }
        rowg[0] = 0; for (size_t s = 0; s < N; s++) rowg[s + 1] = rowg[s] + ng[s];
        size_t totg = rowg[N], totc = 0; for (auto &v : cc) totc += v.size();
        gstart.resize(totg + 1); gdir.reserve(totg); cand.reserve(totc);
        size_t g = 0; gstart[0] = 0;
        for (size_t ch = 0; ch < nch; ch++) {
            for (size_t k = 0; k < cg[ch].size(); k++) { gstart[g + 1] = gstart[g] + cg[ch][k]; g++; }
            gdir.insert(gdir.end(), cd[ch].begin(), cd[ch].end());
            cand.insert(cand.end(), cc[ch].begin(), cc[ch].end());
            std::vector<u32>().swap(cg[ch]); std::vector<u32>().swap(cc[ch]); std::vector<uint8_t>().swap(cd[ch]);
        }
    }
    printf("DFS nodes %llu, leaves %llu\n", (unsigned long long)TOT_NODES, (unsigned long long)TOT_LEAVES);
    printf("groups=%llu candidates=%zu (%.2f per group)  children %.1fs\n", (unsigned long long)rowg[N], cand.size(),
           (double)cand.size() / std::max<u64>(1, rowg[N]), wall() - T1);
    fflush(stdout);
    // ---- min-sum value iteration ----
    std::vector<double> v(N, 1.0), nv(N);
    double ub = 0, lb = 0; int it;
    for (it = 1; it <= MAXIT; it++) {
        double mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
        for (size_t s = 0; s < N; s++) {
            double sum = 0;
            for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                double m = 1e300;
                for (u64 c = gstart[g]; c < gstart[g + 1]; c++) m = std::min(m, v[cand[c]]);
                sum += m;
            }
            nv[s] = sum; mx = std::max(mx, sum);
        }
        double hi = 0, lo = 1e300;
#pragma omp parallel for schedule(static) reduction(max : hi) reduction(min : lo)
        for (size_t s = 0; s < N; s++) if (v[s] > 0) { double r = nv[s] / v[s]; hi = std::max(hi, r); if (nv[s] > 0) lo = std::min(lo, r); }
        ub = hi; lb = lo;
#pragma omp parallel for schedule(static)
        for (size_t s = 0; s < N; s++) v[s] = nv[s] / mx;
        if (it % 100 == 0) { printf("  it %d  [%.10f, %.10f]\n", it, lb, ub); fflush(stdout); }
        if (ub - lb < 1e-10) break;
    }
    printf("value iteration: it=%d  lambda in [%.10f, %.10f]  %.1fs\n", it, lb, ub, wall() - T0);
    fflush(stdout);
    // ---- integer vector, explicit policy, exact ratio ----
    double vmax = 0, vmin = 1e300; for (size_t s = 0; s < N; s++) { vmax = std::max(vmax, v[s]); vmin = std::min(vmin, v[s]); }
    printf("vmin/vmax = %.3e\n", vmin / vmax);
    std::vector<u64> V(N);
    for (size_t s = 0; s < N; s++) { double x = std::ceil(v[s] / vmax * 1125899906842624.0); V[s] = (u64)std::max(1.0, x); }
    std::vector<double>().swap(nv); std::vector<double>().swap(v);
    std::vector<u32> child((size_t)N * NDIR, NONE);
    u64 bestNum = 0, bestDen = 1;
#pragma omp parallel
    {
        u64 lNum = 0, lDen = 1;
#pragma omp for schedule(static)
        for (size_t s = 0; s < N; s++) {
            u64 sum = 0;
            for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                u64 m = ~0ull; u32 arg = NONE;
                for (u64 c = gstart[g]; c < gstart[g + 1]; c++) if (V[cand[c]] < m) { m = V[cand[c]]; arg = cand[c]; }
                sum += m; child[(size_t)s * NDIR + gdir[g]] = arg;
            }
            if ((u128)sum * lDen > (u128)lNum * V[s]) { lNum = sum; lDen = V[s]; }
        }
#pragma omp critical
        { if ((u128)lNum * bestDen > (u128)bestNum * lDen) { bestNum = lNum; bestDen = lDen; } }
    }
    u64 Q = 10000000ull;
    u64 P = (u64)(((u128)bestNum * Q + bestDen - 1) / bestDen);
    printf("integer CW: max ratio = %llu/%llu = %.12f  => certified lambda_bar = %llu/%llu\n", (unsigned long long)bestNum,
           (unsigned long long)bestDen, (double)bestNum / bestDen, (unsigned long long)P, (unsigned long long)Q);
    // ---- path counts of the explicit automaton from the start state ----
    {
        std::vector<u128> x(N, 1), y(N);
        printf("path counts a_n (explicit automaton):\n");
        for (int n = 1; n <= 32; n++) {
#pragma omp parallel for schedule(static)
            for (size_t s = 0; s < N; s++) { u128 acc = 0; for (int j = 0; j < NDIR; j++) { u32 c = child[(size_t)s * NDIR + j]; if (c != NONE) acc += x[c]; } y[s] = acc; }
            x.swap(y);
            u128 a = x[startIdx]; char buf[64]; int k = 63; buf[k] = 0; if (!a) buf[--k] = '0'; while (a) { buf[--k] = char('0' + (int)(a % 10)); a /= 10; }
            printf("  n=%2d a_n=%s\n", n, buf + k);
            if ((double)x[startIdx] > 1e36) break;
        }
    }
    // ---- write certificate ----
    {
        std::string cp = outp + ".lcert";
        FILE *f = fopen(cp.c_str(), "wb");
        u32 ll = (u32)LAT, kk = (u32)K; u64 NN = N;
        fwrite("LATCERT1", 1, 8, f); fwrite(&ll, 4, 1, f); fwrite(&kk, 4, 1, f); fwrite(&NN, 8, 1, f);
        fwrite(&P, 8, 1, f); fwrite(&Q, 8, 1, f);
        fwrite(codes.data(), 8, N, f); fwrite(V.data(), 8, N, f); fwrite(child.data(), 4, child.size(), f);
        fclose(f);
        std::string hp = outp + ".hdr";
        f = fopen(hp.c_str(), "w");
        fprintf(f, "lattice %s\nK %d\nmode %s\nstates %zu\nlambda_bar %llu/%llu\nmax_ratio %llu/%llu\nfloat_bracket %.10f %.10f\n", lat.c_str(), K,
                mode.c_str(), N, (unsigned long long)P, (unsigned long long)Q, (unsigned long long)bestNum, (unsigned long long)bestDen, lb, ub);
        fclose(f);
    }
    printf("RESULT lattice=%s K=%d mode=%s states=%zu lambda_float=[%.10f,%.10f] certified<=%llu/%llu total %.1fs\n", lat.c_str(), K,
           mode.c_str(), N, lb, ub, (unsigned long long)P, (unsigned long long)Q, wall() - T0);
    return 0;
}
