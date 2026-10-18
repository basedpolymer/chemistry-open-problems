// saw_automaton.cpp -- exploration version of the finite-memory automaton for SAWs on Z^d.
//
// States: self-avoiding walks w = (B = v_0, ..., v_l = A) taken up to translation and up to the
// hyperoctahedral group (canonical form: scanning steps from the newest one, axes are relabelled in order
// of first appearance and the first step along each axis is made positive). The forbidden set of a state is
// {v_0, ..., v_{l-1}} - A.  A state is admissible when loopsize(w) = l + dist(A, B) <= K.
//
// Modes:
//   mode 0  Ponitz-Tittmann: unique child = longest admissible suffix of w+e.
//   mode 1  every admissible suffix of w+e is a candidate child.
//   mode 2  mode 1 + "bridge" shortcuts: if v_i ~ v_j (adjacent, j >= i+3) the vertices strictly between are
//           removed; every admissible suffix of every shortcut walk is a candidate.
// dist: 0 = Manhattan, 1 = shortest path from A to B avoiding the other vertices of the walk.
//
// Value iteration: v <- T v, T(v)_s = sum_{allowed e} min_{c in C(s,e)} v_c ; reports Collatz-Wielandt
// upper bound max_s (Tv)_s / v_s (for the argmin policy) -- floating point, exploration only.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <array>

#include <omp.h>
using u64 = uint64_t;
using u32 = uint32_t;

static int D, K, MODE, DISTMODE, BITS;
static const int MAXL = 40;

struct Walk {
    int l;
    unsigned char s[MAXL];  // s[0] oldest step (B -> v1), s[l-1] newest step (v_{l-1} -> A); dir = 2*axis + neg
};

static inline int dsign(int dir) { return (dir & 1) ? -1 : 1; }
static inline int daxis(int dir) { return dir >> 1; }

static u64 encode(const Walk &w) {
    u64 c = 0;
    for (int i = 0; i < w.l; i++) c |= (u64)w.s[i] << (BITS * i);
    c |= (u64)w.l << 57;
    return c | (1ull << 63);
}
static void decode(u64 c, Walk &w) {
    w.l = (int)((c >> 57) & 63);
    u64 m = (1ull << BITS) - 1;
    for (int i = 0; i < w.l; i++) w.s[i] = (unsigned char)((c >> (BITS * i)) & m);
}
static void canon(Walk &w) {
    int map[8], flip[8], nxt = 0;
    for (int a = 0; a < D; a++) map[a] = -1;
    for (int i = w.l - 1; i >= 0; i--) {
        int a = daxis(w.s[i]), sg = w.s[i] & 1;
        if (map[a] < 0) { map[a] = nxt++; flip[a] = sg; }
        w.s[i] = (unsigned char)((map[a] << 1) | (sg ^ flip[a]));
    }
}
// positions relative to A: pos[i][*] for vertex v_i, i=0..l (v_l = A = 0)
static void positions(const Walk &w, int pos[][8]) {
    for (int a = 0; a < D; a++) pos[w.l][a] = 0;
    for (int i = w.l - 1; i >= 0; i--) {
        for (int a = 0; a < D; a++) pos[i][a] = pos[i + 1][a];
        pos[i][daxis(w.s[i])] -= dsign(w.s[i]);
    }
}
static inline int l1(const int *p) { int s = 0; for (int a = 0; a < D; a++) s += abs(p[a]); return s; }

// ---- exact avoiding distance via BFS on a small box --------------------------------------------------
static const int BR = 24;           // box radius
static thread_local std::vector<u32> stampv;
static thread_local u32 stampc = 0;
static inline int boxidx(const int *p) {
    int idx = 0;
    for (int a = 0; a < D; a++) idx = idx * (2 * BR + 1) + (p[a] + BR);
    return idx;
}
// shortest path length from A (=0) to v_0 avoiding v_1..v_{l-1}; returns > limit if larger than limit
static int avoid_dist(int l, int pos[][8], int limit) {
    int dd = l1(pos[0]);
    if (dd > limit) return dd;
    size_t boxn = 1; for (int a = 0; a < D; a++) boxn *= (2 * BR + 1);
    if (stampv.size() != boxn) { stampv.assign(boxn, 0); stampc = 0; }
    stampc += 2;
    if (stampc > 0xfffffff0u) { std::fill(stampv.begin(), stampv.end(), 0); stampc = 2; }
    u32 BLOCK = stampc, SEEN = stampc + 1;
    for (int i = 1; i < l; i++) stampv[boxidx(pos[i])] = BLOCK;
    int target = boxidx(pos[0]);
    std::vector<std::array<int, 8>> cur, nx;
    std::array<int, 8> z{}; cur.push_back(z);
    stampv[boxidx(z.data())] = SEEN;
    for (int dist = 0; dist < limit; dist++) {
        nx.clear();
        for (auto &p : cur) for (int dir = 0; dir < 2 * D; dir++) {
            std::array<int, 8> q = p; q[daxis(dir)] += dsign(dir);
            if (abs(q[daxis(dir)]) >= BR) continue;
            int id = boxidx(q.data());
            if (id == target) return dist + 1;
            if (stampv[id] == BLOCK || stampv[id] == SEEN) continue;
            // prune: remaining budget
            if (dist + 1 + l1(q.data()) > limit + 0 && false) continue;
            stampv[id] = SEEN; nx.push_back(q);
        }
        std::swap(cur, nx);
        if (cur.empty()) break;
    }
    return limit + 1;
}


static int EXTRA = 0, EXTRA_L = 0, PRUNE = 0;
static bool has_bridge(const Walk &w, int pos[][8]) {
    for (int i = 0; i <= w.l; i++) for (int j = i + 3; j <= w.l; j++) {
        int s = 0; for (int c = 0; c < D; c++) s += abs(pos[i][c] - pos[j][c]);
        if (s == 1) return true;
    }
    return false;
}
static bool line_start(const Walk &w) {
    if (w.l < EXTRA_L) return false;
    for (int i = 1; i < EXTRA_L; i++) if (w.s[i] != w.s[0]) return false;
    return true;
}
static int budget(const Walk &w, int pos[][8]) {
    int b = K;
    if (EXTRA == 1 && !has_bridge(w, pos)) b = K + 2;
    if (EXTRA == 2 && line_start(w)) b = K + 2;
    if (EXTRA == 3 && (line_start(w) || !has_bridge(w, pos))) b = K + 2;
    if (EXTRA == 4 && !has_bridge(w, pos)) b = K + 2 + (line_start(w) ? 2 : 0);
    return b;
}
static bool admissible(const Walk &w, int pos[][8]) {
    if (w.l == 0) return true;
    int dd = l1(pos[0]);
    if (w.l + dd > K + 4) return false;
    int B = budget(w, pos);
    if (w.l + dd > B) return false;
    if (DISTMODE == 0) return true;
    int lim = B - w.l;
    return avoid_dist(w.l, pos, lim) <= lim;
}

// ---- hash table ----------------------------------------------------------------------------------
struct HT {
    std::vector<u64> keys; std::vector<u32> vals; u64 mask;
    void init(u64 cap) { u64 n = 1; while (n < cap) n <<= 1; keys.assign(n, 0); vals.assign(n, 0); mask = n - 1; }
    static inline u64 h(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return x; }
    // returns index or -1
    long find(u64 k) const { u64 i = h(k) & mask; while (keys[i]) { if (keys[i] == k) return vals[i]; i = (i + 1) & mask; } return -1; }
    bool insert(u64 k, u32 v) { u64 i = h(k) & mask; while (keys[i]) { if (keys[i] == k) return false; i = (i + 1) & mask; } keys[i] = k; vals[i] = v; return true; }
};

// make the walk consisting of the vertices idx[0..m-1] (indices into pos) -- consecutive must be adjacent
static void walk_from_vertices(int m, const int *idx, int pos[][8], Walk &out) {
    out.l = m - 1;
    for (int t = 0; t + 1 < m; t++) {
        const int *p = pos[idx[t]], *q = pos[idx[t + 1]];
        int dir = -1;
        for (int a = 0; a < D; a++) if (q[a] != p[a]) dir = 2 * a + (q[a] - p[a] < 0 ? 1 : 0);
        out.s[t] = (unsigned char)dir;
    }
}

// children candidates of state w in direction dir.  Returns false if the step is forbidden.
static bool children(const Walk &w, int dir, std::vector<u64> &out) {
    out.clear();
    static thread_local int pos[MAXL + 2][8];
    positions(w, pos);
    // forbidden?
    int e[8] = {0}; e[daxis(dir)] = dsign(dir);
    for (int i = 0; i < w.l; i++) {
        bool eq = true;
        for (int a = 0; a < D; a++) if (pos[i][a] != e[a]) { eq = false; break; }
        if (eq) return false;
    }
    // augmented walk a = w + e, positions relative to new A
    Walk a = w; a.s[a.l] = (unsigned char)dir; a.l++;
    int apos[MAXL + 2][8];
    for (int i = 0; i <= w.l; i++) for (int c = 0; c < D; c++) apos[i][c] = pos[i][c] - e[c];
    for (int c = 0; c < D; c++) apos[a.l][c] = 0;
    // list of base walks (as vertex index lists) : the walk itself and bridge shortcuts
    std::vector<std::vector<int>> bases;
    { std::vector<int> v(a.l + 1); for (int i = 0; i <= a.l; i++) v[i] = i; bases.push_back(v); }
    if (MODE >= 2) {
        // bridges: v_i adjacent to v_j (positions i<j in the current vertex list, j >= i+3); remove the vertices
        // strictly between.  MODE 2: single bridge applied to the full walk; MODE 3: closure under repeated bridges.
        for (size_t b = 0; b < bases.size() && bases.size() < 256; b++) {
            if (MODE == 2 && b > 0) break;
            std::vector<int> cur = bases[b];
            int m = (int)cur.size();
            for (int i = 0; i < m; i++) for (int j = i + 3; j < m; j++) {
                int dsum = 0; for (int c = 0; c < D; c++) dsum += abs(apos[cur[i]][c] - apos[cur[j]][c]);
                if (dsum != 1) continue;
                std::vector<int> v; for (int t = 0; t <= i; t++) v.push_back(cur[t]); for (int t = j; t < m; t++) v.push_back(cur[t]);
                bool dup = false;
                for (auto &q : bases) if (q == v) { dup = true; break; }
                if (!dup) bases.push_back(v);
            }
        }
    }
    if (MODE == 4) {
        // all maximal self-avoiding paths starting at A' inside the induced subgraph on the visited vertices
        bases.clear();
        int n = a.l + 1;
        std::vector<unsigned> adj(n, 0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) {
            int s = 0; for (int c = 0; c < D; c++) s += abs(apos[i][c] - apos[j][c]);
            if (s == 1) adj[i] |= 1u << j;
        }
        std::vector<int> path; path.push_back(a.l);
        unsigned used = 1u << a.l;
        // iterative DFS
        std::vector<unsigned> rem; rem.push_back(adj[a.l]);
        while (!rem.empty() && bases.size() < 4096) {
            unsigned &r = rem.back();
            unsigned avail = r & ~used;
            if (!avail) {
                // maximal if no extension at all from this node
                if ((adj[path.back()] & ~used) == 0) {
                    std::vector<int> v(path.rbegin(), path.rend());
                    bases.push_back(v);
                }
                used &= ~(1u << path.back()); path.pop_back(); rem.pop_back();
                continue;
            }
            int nxt = __builtin_ctz(avail);
            r &= ~(1u << nxt);
            path.push_back(nxt); used |= 1u << nxt; rem.push_back(adj[nxt]);
        }
    }
    int tpos[MAXL + 2][8];
    std::vector<std::pair<unsigned, u64>> mc;  // (vertex mask, code)
    for (auto &bv : bases) {
        int m = (int)bv.size();
        // suffixes: start index st from 0 .. m-1 ; first admissible = longest
        for (int st = 0; st < m; st++) {
            Walk sw; walk_from_vertices(m - st, bv.data() + st, apos, sw);
            for (int t = 0; t < m - st; t++) for (int c = 0; c < D; c++) tpos[t][c] = apos[bv[st + t]][c];
            if (!admissible(sw, tpos)) continue;
            canon(sw);
            unsigned mask = 0; for (int t = st; t < m; t++) mask |= 1u << bv[t];
            mc.push_back({mask, encode(sw)});
            if (MODE != 1) break;   // modes 0,2,3: longest admissible suffix of each base only
        }
        if (MODE == 0) break;
    }
    if (PRUNE) {
        // drop candidates whose vertex set is a strict subset of another candidate's vertex set
        std::vector<std::pair<unsigned, u64>> keep;
        for (size_t i = 0; i < mc.size(); i++) {
            bool dom = false;
            for (size_t j = 0; j < mc.size() && !dom; j++)
                if (j != i && (mc[i].first & mc[j].first) == mc[i].first && mc[i].first != mc[j].first) dom = true;
            if (!dom) keep.push_back(mc[i]);
        }
        mc.swap(keep);
    }
    for (auto &p : mc) out.push_back(p.second);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return true;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s d k mode distmode [iters]\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]); MODE = atoi(argv[3]); DISTMODE = atoi(argv[4]);
    int ITERS = argc > 5 ? atoi(argv[5]) : 2000;
    EXTRA = argc > 6 ? atoi(argv[6]) : 0; EXTRA_L = argc > 7 ? atoi(argv[7]) : 0;
    if (MODE >= 5) { PRUNE = 1; MODE -= 1; }
    BITS = (2 * D <= 8) ? 3 : 4;
    double t0 = omp_get_wtime();
    HT ht; ht.init(1 << 20);
    std::vector<u64> states;
    Walk w0; w0.l = 0; u64 k0 = encode(w0);
    ht.insert(k0, 0); states.push_back(k0);
    // rows: for each state, list of groups (one per allowed direction), each group = list of candidate indices
    std::vector<u64> rowptr; std::vector<u32> grpptr; std::vector<u64> cand_keys;
    std::vector<u64> rowg; // group start per state
    std::vector<u32> gstart; // candidate start per group (into cand)
    std::vector<u32> cand;
    rowg.push_back(0); gstart.push_back(0);
    std::vector<u64> ch;
    for (size_t si = 0; si < states.size(); si++) {
        Walk w; decode(states[si], w);
        for (int dir = 0; dir < 2 * D; dir++) {
            if (!children(w, dir, ch)) continue;
            for (u64 c : ch) {
                long id = ht.find(c);
                if (id < 0) {
                    id = (long)states.size(); states.push_back(c);
                    if (states.size() * 2 > ht.mask) {
                        HT n2; n2.init((ht.mask + 1) * 4);
                        for (size_t q = 0; q < states.size(); q++) n2.insert(states[q], (u32)q);
                        std::swap(ht, n2);
                    } else ht.insert(c, (u32)id);
                }
                cand.push_back((u32)id);
            }
            gstart.push_back((u32)cand.size());
        }
        rowg.push_back(gstart.size() - 1);
    }
    size_t N = states.size();
    double t1 = omp_get_wtime();
    printf("d=%d k=%d mode=%d dist=%d : states=%zu groups=%zu cands=%zu  gen %.1fs\n", D, K, MODE, DISTMODE, N,
           gstart.size() - 1, cand.size(), t1 - t0);
    fflush(stdout);
    // value iteration
    std::vector<double> v(N, 1.0), nv(N);
    double ub = 0, lb = 0;
    for (int it = 1; it <= ITERS; it++) {
        double mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
        for (size_t s = 0; s < N; s++) {
            double sum = 0;
            for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                double m = 1e300;
                for (u32 c = gstart[g]; c < gstart[g + 1]; c++) m = std::min(m, v[cand[c]]);
                sum += m;
            }
            nv[s] = sum; mx = std::max(mx, sum);
        }
        // CW bounds over states with v>0 (skip start state 0, transient)
        double hi = 0, lo = 1e300;
#pragma omp parallel for schedule(static) reduction(max : hi) reduction(min : lo)
        for (size_t s = 1; s < N; s++) {
            if (v[s] > 0) { double r = nv[s] / v[s]; hi = std::max(hi, r); if (nv[s] > 0) lo = std::min(lo, r); }
        }
        ub = hi; lb = lo;
        for (size_t s = 0; s < N; s++) v[s] = nv[s] / mx;
        if (it % 100 == 0 || it == ITERS) { printf("  it %d  CW lower %.9f upper %.9f\n", it, lb, ub); fflush(stdout); }
        if (ub - lb < 1e-9) { printf("  converged it %d  CW lower %.10f upper %.10f\n", it, lb, ub); break; }
    }
    printf("RESULT extra=%d L=%d d=%d k=%d mode=%d dist=%d states=%zu lambda in [%.9f, %.9f]  total %.1fs\n",
           EXTRA, EXTRA_L, D, K, MODE, DISTMODE, N, lb, ub,
           omp_get_wtime() - t0);
    return 0;
}
