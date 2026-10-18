// saw_prod.cpp -- production code: builds the "loop-size <= K" automaton for self-avoiding walks on Z^d with
// path-shortcut children, optimises the child choice (min-sum value iteration), and writes an integer
// Collatz-Wielandt certificate vector.
//
// Build:  g++ -O3 -march=native -fopenmp -o saw_prod saw_prod.cpp          (2 <= d <= 6)
// Usage:  saw_prod d K outprefix [maxiter] [policy]
//           policy = 0 : Ponitz-Tittmann child (longest admissible suffix of the walk)  -> reproduces EJC 7 R21
//           policy = 1 : all maximal self-avoiding paths from the new position inside the visited set
//           policy = 2 : (experiment) single heuristic child "most sites" -- worse than policy 0, not used
//           policy = 3 : same candidate set as 1, faster search (admissibility by BFS)
//           policy = 4 : same candidate set as 1, admissibility by table lookup -- used for all results
// Environment options used only for experiments reported in PROGRESS.md: SAW_PI_TEST (policy-iteration
// convergence test), SAW_IMPORTANCE_OUT / SAW_REFINE / SAW_REFINE_FRAC (hybrid state sets), SAW_NOBFS (timing).
//
// Validity of the output never depends on these choices: the certificate is checked by verify_cert.cpp.
//
// STATES.  A state is a self-avoiding walk w = (v_l = B, ..., v_1, v_0 = A) stored newest-first: t[i] is the
// direction (2*axis + neg) of the step v_{i+1} -> v_i.  Positions are taken relative to A = 0, so
// p_0 = 0, p_{i+1} = p_i - unit(t[i]).  Forbidden set F(w) = {p_1, ..., p_l}.
// Canonical form: scanning t[0], t[1], ..., axes are relabelled 0,1,2,... in order of first appearance and
// the first step along each axis is made positive.  Admissible: l + |p_l|_1 <= K and the shortest lattice path
// from A to B avoiding p_1..p_{l-1} has length <= K - l.  The state set is the set of ALL admissible
// canonical walks (enumerated directly), sorted by code.
//
// CHILDREN.  From state w and direction e with e not in F(w): the visited set relative to A' = e is
// {p_i - e} u {-e}... i.e. q_0 = 0 (A'), q_{i+1} = p_i - e.  Any self-avoiding path in the induced subgraph on
// {q_0..q_{l+1}} starting at q_0 is a sub-walk whose forbidden set is contained in (F(w) u {0}) - e, hence a
// valid child.  We use the longest admissible prefix of every maximal such path and drop candidates whose
// vertex set is strictly contained in another candidate's vertex set.
//
// OUTPUT.  <outprefix>.cert (SAWCERT1: header, sorted codes, uint64 weights), <outprefix>.hdr (text summary),
// <outprefix>.v64 (weights only).  Also prints the automaton path counts a_n (n <= 24) for comparison with c_n.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
using u64 = uint64_t;
using u32 = uint32_t;
using u128 = unsigned __int128;

static int D, K, BITS, POLICY, NOBFS = 0;
static const int MAXL = 24;
static const int MAXD = 6;

struct Pos { int8_t c[MAXD]; };
static inline int l1(const Pos &p) { int s = 0; for (int a = 0; a < D; a++) s += abs(p.c[a]); return s; }
static inline int l1d(const Pos &p, const Pos &q) { int s = 0; for (int a = 0; a < D; a++) s += abs(p.c[a] - q.c[a]); return s; }
static inline bool eqp(const Pos &p, const Pos &q) { for (int a = 0; a < D; a++) if (p.c[a] != q.c[a]) return false; return true; }

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
static inline void canon(int l, int *t) {
    int map[MAXD], flip[MAXD], nxt = 0;
    for (int a = 0; a < D; a++) map[a] = -1;
    for (int i = 0; i < l; i++) {
        int a = t[i] >> 1, sg = t[i] & 1;
        if (map[a] < 0) { map[a] = nxt++; flip[a] = sg; }
        t[i] = (map[a] << 1) | (sg ^ flip[a]);
    }
}
static inline void positions(int l, const int *t, Pos *p) {
    for (int a = 0; a < D; a++) p[0].c[a] = 0;
    for (int i = 0; i < l; i++) { p[i + 1] = p[i]; p[i + 1].c[t[i] >> 1] -= (t[i] & 1) ? -1 : 1; }
}

// ---------------- exact avoiding distance (BFS with A* pruning) ------------------------------------------
struct Grid {
    int R, W; std::vector<u32> st; u32 stamp = 0; std::vector<int> q1, q2;
    void init(int r) { R = r; W = 2 * r + 1; size_t n = 1; for (int a = 0; a < D; a++) n *= W; st.assign(n, 0); stamp = 0; }
    inline int idx(const Pos &p) const { int x = 0; for (int a = 0; a < D; a++) x = x * W + (p.c[a] + R); return x; }
};
static thread_local Grid G;
// path P[0..m] (P[0] = A, P[m] = B); returns true iff there is a lattice path A -> B of length <= lim avoiding P[1..m-1]
static bool avoid_ok(int m, const Pos *P, int lim) {
    if (m <= 1) return true;
    int dAB = l1d(P[0], P[m]);
    if (dAB > lim) return false;
    if (G.st.empty()) G.init(K / 2 + 3);
    G.stamp += 2;
    if (G.stamp > 0xfffffff0u) { std::fill(G.st.begin(), G.st.end(), 0); G.stamp = 2; }
    u32 BLK = G.stamp, SEEN = G.stamp + 1;
    for (int i = 1; i < m; i++) G.st[G.idx(P[i])] = BLK;
    {
        // fast exact test: is there a monotone (shortest) lattice path A -> B avoiding the obstacles?
        // If yes the distance is dAB <= lim.  If no, every avoiding path has length >= dAB + 2.
        int n[MAXD], sg[MAXD], tot = 1;
        for (int a = 0; a < D; a++) { int df = P[m].c[a] - P[0].c[a]; sg[a] = df >= 0 ? 1 : -1; n[a] = abs(df) + 1; tot *= n[a]; }
        static thread_local std::vector<unsigned char> R;
        if ((int)R.size() < tot) R.resize(tot);
        // iterate cells in lexicographic order of offsets (predecessors along each axis come first)
        int off[MAXD] = {0, 0, 0, 0, 0, 0};
        for (int id = 0; id < tot; id++) {
            // decode id -> off (axis 0 slowest)
            int x = id; for (int a = D - 1; a >= 0; a--) { off[a] = x % n[a]; x /= n[a]; }
            Pos c; for (int a = 0; a < D; a++) c.c[a] = (int8_t)(P[0].c[a] + sg[a] * off[a]);
            unsigned char r = 0;
            if (id == 0) r = 1;
            else if (G.st[G.idx(c)] != BLK) {
                int stride = 1;
                for (int a = D - 1; a >= 0 && !r; a--) { if (off[a] > 0 && R[id - stride]) r = 1; stride *= n[a]; }
            }
            R[id] = r;
        }
        if (R[tot - 1]) return true;
        if (lim < dAB + 2) return false;
    }
    // layered BFS from A with A* pruning (dist + |q - B|_1 <= lim); fixed thread-local queue, no allocation
    static thread_local std::vector<Pos> Qb;
    if (Qb.size() < G.st.size() + 8) Qb.resize(G.st.size() + 8);
    const Pos &B = P[m];
    size_t head = 0, tail = 0, layerEnd;
    Qb[tail++] = P[0]; G.st[G.idx(P[0])] = SEEN; layerEnd = tail;
    int dist = 0;
    while (head < tail) {
        if (head == layerEnd) { dist++; layerEnd = tail; if (dist >= lim) return false; }
        Pos p = Qb[head++];
        for (int dir = 0; dir < 2 * D; dir++) {
            Pos q = p; q.c[dir >> 1] += (dir & 1) ? -1 : 1;
            if (eqp(q, B)) return true;
            if (dist + 1 + l1d(q, B) > lim) continue;
            if (abs(q.c[dir >> 1]) > G.R - 1) continue;
            int id = G.idx(q);
            if (G.st[id] == BLK || G.st[id] == SEEN) continue;
            G.st[id] = SEEN; Qb[tail++] = q;
        }
    }
    return false;
}
static inline bool admissible(int m, const Pos *P) {  // P[0] = current position (A), P[m] = B
    if (m == 0) return true;
    int d = l1d(P[0], P[m]);
    if (m + d > K) return false;
    if (NOBFS) return true;
    return avoid_ok(m, P, K - m);
}

// ---------------- enumeration of all admissible canonical walks ------------------------------------------
// DFS from A into the past; canonical constraint is checked incrementally (newest-first scanning).
struct EnumCtx {
    int t[MAXL]; Pos p[MAXL + 1]; int map_next; int axis_used[MAXD];
    std::vector<u64> out;
};
static void enum_dfs(EnumCtx &C, int l, int nused) {
    // current walk: steps t[0..l-1], positions p[0..l]
    if (!admissible(l, C.p)) {
        // admissibility can fail through the exact avoid-distance while longer walks are admissible again;
        // we only prune on the monotone Manhattan part.
        if (l + l1(C.p[l]) > K) return;
    } else {
        C.out.push_back(encode(l, C.t));
    }
    if (l + 1 + std::max(0, l1(C.p[l]) - 1) > K) return;
    for (int dir = 0; dir < 2 * D; dir++) {
        int a = dir >> 1;
        if (a > nused) continue;                       // new axes must appear in order
        if (a == nused && (dir & 1)) continue;          // first step along a new axis is positive
        Pos q = C.p[l]; q.c[a] -= (dir & 1) ? -1 : 1;   // older vertex = current oldest - unit(dir)
        bool hit = false;
        for (int i = 0; i <= l; i++) if (eqp(q, C.p[i])) { hit = true; break; }
        if (hit) continue;
        if (l + 1 + l1(q) > K) continue;
        C.t[l] = dir; C.p[l + 1] = q;
        enum_dfs(C, l + 1, a == nused ? nused + 1 : nused);
    }
}

// ---------------- lookup table -----------------------------------------------------------------------------
struct Lookup {
    const std::vector<u64> *codes; std::vector<u32> tab; u64 mask;
    static inline u64 h(u64 x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ull; x ^= x >> 27; x *= 0x81dadef4bc2dd44dull; x ^= x >> 33; return x; }
    void build(const std::vector<u64> &c) {
        codes = &c; u64 n = 1; while (n < c.size() * 3 / 2 + 16) n <<= 1;
        tab.assign(n, 0xffffffffu); mask = n - 1;
        for (size_t i = 0; i < c.size(); i++) { u64 j = h(c[i]) & mask; while (tab[j] != 0xffffffffu) j = (j + 1) & mask; tab[j] = (u32)i; }
    }
    inline long find(u64 k) const { u64 j = h(k) & mask; while (tab[j] != 0xffffffffu) { if ((*codes)[tab[j]] == k) return tab[j]; j = (j + 1) & mask; } return -1; }
};

// ---------------- children ---------------------------------------------------------------------------------
static const Lookup *GLK = nullptr;
// admissibility of the walk P[d2] (oldest) ... P[0] (newest) by table lookup: the state table contains exactly the
// admissible canonical walks, and admissibility is invariant under the symmetry group.
static bool adm_lookup(int m, const Pos *P) {
    if (m == 0) return true;
    if (m + l1d(P[0], P[m]) > K) return false;
    int tt[MAXL + 2];
    for (int i = 0; i < m; i++) {
        int d2 = -1;
        for (int a = 0; a < D; a++) if (P[i].c[a] != P[i + 1].c[a]) d2 = 2 * a + (P[i].c[a] - P[i + 1].c[a] < 0 ? 1 : 0);
        tt[i] = d2;
    }
    canon(m, tt);
    return GLK->find(encode(m, tt)) >= 0;
}
// fills out with candidate child codes for state (l,t) and direction dir; returns false if forbidden.
static bool children(int l, const int *t, int dir, std::vector<u64> &out) {
    out.clear();
    Pos p[MAXL + 2]; positions(l, t, p);
    Pos e; for (int a = 0; a < D; a++) e.c[a] = 0; e.c[dir >> 1] = (dir & 1) ? -1 : 1;
    for (int i = 1; i <= l; i++) if (eqp(p[i], e)) return false;
    int n = l + 2;
    Pos q[MAXL + 3];
    for (int a = 0; a < D; a++) q[0].c[a] = 0;
    for (int i = 0; i <= l; i++) for (int a = 0; a < D; a++) q[i + 1].c[a] = p[i].c[a] - e.c[a];
    // vertex-index paths from q[0]
    static thread_local std::vector<std::pair<u32, u64>> mc; mc.clear();
    auto add_path = [&](int m, const int *idx) {
        // longest admissible prefix of path idx[0..m-1] (idx[0] = 0 = A')
        Pos P[MAXL + 3];
        for (int i = 0; i < m; i++) P[i] = q[idx[i]];
        int best = -1;
        for (int len = m - 1; len >= 0; len--) if (admissible(len, P)) { best = len; break; }
        int tt[MAXL + 2]; u32 mask = 1;
        for (int i = 0; i < best; i++) {
            const Pos &newer = P[i], &older = P[i + 1];
            int d2 = -1;
            for (int a = 0; a < D; a++) if (newer.c[a] != older.c[a]) d2 = 2 * a + (newer.c[a] - older.c[a] < 0 ? 1 : 0);
            tt[i] = d2; mask |= 1u << idx[i + 1];
        }
        canon(best, tt);
        mc.push_back({mask, encode(best, tt)});
    };
    if (POLICY == 0) {
        int idx[MAXL + 3]; for (int i = 0; i < n; i++) idx[i] = i;
        add_path(n, idx);
    } else if (POLICY == 3 || POLICY == 4) {
        // same candidate set as POLICY 1 (after dominance pruning), computed faster: DFS over paths from A'
        // inside the visited set, cut by the monotone Manhattan bound m + |P_m|_1 <= K, admissibility evaluated
        // once per DFS node; at each leaf the deepest admissible node of the branch is a candidate.
        u32 adj[MAXL + 3];
        for (int i = 0; i < n; i++) { adj[i] = 0; for (int j = 0; j < n; j++) if (i != j && l1d(q[i], q[j]) == 1) adj[i] |= 1u << j; }
        int path[MAXL + 3], deep[MAXL + 3]; u32 rem[MAXL + 3]; Pos P[MAXL + 3];
        auto ext = [&](int depth, u32 used) {
            u32 a = adj[path[depth]] & ~used, r = 0;
            while (a) { int j = __builtin_ctz(a); a &= a - 1; if (depth + 1 + l1(q[j]) <= K) r |= 1u << j; }
            return r;
        };
        auto emit = [&](int best) {
            int tt[MAXL + 2]; u32 mask = 1;
            for (int i = 0; i < best; i++) {
                const Pos &newer = P[i], &older = P[i + 1];
                int d2 = -1;
                for (int a = 0; a < D; a++) if (newer.c[a] != older.c[a]) d2 = 2 * a + (newer.c[a] - older.c[a] < 0 ? 1 : 0);
                tt[i] = d2; mask |= 1u << path[i + 1];
            }
            canon(best, tt);
            mc.push_back({mask, encode(best, tt)});
        };
        // deep[i] : -1 = admissibility of prefix i not yet known, 0/1 = known
        auto leaf = [&](int depth) {
            for (int d2 = depth; d2 >= 0; d2--) {
                if (deep[d2] < 0) deep[d2] = (POLICY == 4 ? adm_lookup(d2, P) : admissible(d2, P)) ? 1 : 0;
                if (deep[d2]) { emit(d2); return; }
            }
        };
        int depth = 0; path[0] = 0; P[0] = q[0]; deep[0] = 1; u32 used = 1;
        rem[0] = ext(0, used);
        if (!rem[0]) leaf(0);
        while (depth >= 0) {
            if (!rem[depth]) { used &= ~(1u << path[depth]); depth--; continue; }
            int nx = __builtin_ctz(rem[depth]);
            rem[depth] &= rem[depth] - 1;
            depth++; path[depth] = nx; P[depth] = q[nx]; used |= 1u << nx;
            deep[depth] = -1;
            rem[depth] = ext(depth, used);
            if (!rem[depth]) leaf(depth);
        }
        static thread_local std::vector<std::pair<u32, u64>> keep;
        // dominance pruning: drop (mask, index) pairs whose mask is strictly contained in another candidate's mask.
        // Exact duplicates are removed first; then candidates are scanned by decreasing cardinality and compared with
        // the maximal ones kept so far (strict inclusion is transitive, so this gives the same result as all pairs).
        std::sort(mc.begin(), mc.end());
        mc.erase(std::unique(mc.begin(), mc.end()), mc.end());
        std::stable_sort(mc.begin(), mc.end(), [](const std::pair<u32, u64> &a, const std::pair<u32, u64> &b) {
            return __builtin_popcount(a.first) > __builtin_popcount(b.first); });
        keep.clear();
        for (size_t i = 0; i < mc.size(); i++) {
            bool dom = false;
            for (size_t j = 0; j < keep.size() && !dom; j++)
                if ((mc[i].first & keep[j].first) == mc[i].first && mc[i].first != keep[j].first) dom = true;
            if (!dom) keep.push_back(mc[i]);
        }
        mc.swap(keep);
    } else {
        u32 adj[MAXL + 3];
        for (int i = 0; i < n; i++) { adj[i] = 0; for (int j = 0; j < n; j++) if (i != j && l1d(q[i], q[j]) == 1) adj[i] |= 1u << j; }
        int path[MAXL + 3]; u32 rem[MAXL + 3]; int depth = 0; path[0] = 0; rem[0] = adj[0]; u32 used = 1;
        int npaths = 0;
        while (depth >= 0 && npaths < 100000) {
            u32 avail = rem[depth] & ~used;
            if (!avail) {
                if ((adj[path[depth]] & ~used) == 0) { add_path(depth + 1, path); npaths++; }
                used &= ~(1u << path[depth]); depth--;
                continue;
            }
            int nx = __builtin_ctz(avail);
            rem[depth] &= ~(1u << nx);
            depth++; path[depth] = nx; rem[depth] = adj[nx]; used |= 1u << nx;
        }
        if (npaths >= 100000) { fprintf(stderr, "too many paths\n"); exit(2); }
        // prune dominated
        static thread_local std::vector<std::pair<u32, u64>> keep2; std::vector<std::pair<u32, u64>> &keep = keep2; keep.clear();
        for (size_t i = 0; i < mc.size(); i++) {
            bool dom = false;
            for (size_t j = 0; j < mc.size() && !dom; j++)
                if ((mc[i].first & mc[j].first) == mc[i].first && mc[i].first != mc[j].first) dom = true;
            if (!dom) keep.push_back(mc[i]);
        }
        mc.swap(keep);
        if (POLICY == 2 && mc.size() > 1) {
            // heuristic single child: most vertices; tie -> most vertices among the 8 most recent; tie -> smallest code
            auto score = [&](const std::pair<u32, u64> &x) { return (u64)__builtin_popcount(x.first) * 1024 + __builtin_popcount(x.first & 0xffu); };
            size_t b = 0;
            for (size_t i = 1; i < mc.size(); i++) {
                u64 si = score(mc[i]), sb = score(mc[b]);
                if (si > sb || (si == sb && mc[i].second < mc[b].second)) b = i;
            }
            std::pair<u32, u64> x = mc[b]; mc.clear(); mc.push_back(x);
        }
    }
    for (auto &x : mc) out.push_back(x.second);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return true;
}

static void sha_like(const std::vector<u64> &codes, u64 &h1, u64 &h2) {
    // simple order-dependent 128-bit fingerprint (FNV-style + mixing); not cryptographic
    h1 = 1469598103934665603ull; h2 = 0x9e3779b97f4a7c15ull;
    for (u64 c : codes) { h1 ^= c; h1 *= 1099511628211ull; h2 += c * 0xbf58476d1ce4e5b9ull; h2 ^= h2 >> 29; }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s d K outprefix [maxiter] [policy]\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]); std::string outp = argv[3];
    int MAXIT = argc > 4 ? atoi(argv[4]) : 5000;
    POLICY = argc > 5 ? atoi(argv[5]) : 1;
    NOBFS = getenv("SAW_NOBFS") ? 1 : 0;
    BITS = (2 * D <= 8) ? 3 : 4;
    if ((K - 1) * BITS > 57 || K > MAXL) { fprintf(stderr, "K too large for encoding\n"); return 1; }
    double T0 = omp_get_wtime();
    // ---- enumerate states (parallel over the first few backward steps) ----
    std::vector<u64> codes;
    {
        // seeds: canonical walks of length <= S enumerated serially, then DFS in parallel from seeds of length S
        const int S = std::min(4, K - 1);
        std::vector<EnumCtx> seeds;
        std::vector<u64> shortcodes;
        // serial enumeration up to depth S
        EnumCtx c0; for (int a = 0; a < D; a++) c0.p[0].c[a] = 0;
        // BFS by depth with explicit nused
        struct Item { EnumCtx c; int l, nused; };
        std::vector<Item> layer; layer.push_back({c0, 0, 0});
        for (int l = 0; l < S; l++) {
            std::vector<Item> nxt;
            for (auto &it : layer) {
                if (admissible(it.l, it.c.p)) shortcodes.push_back(encode(it.l, it.c.t));
                for (int dir = 0; dir < 2 * D; dir++) {
                    int a = dir >> 1;
                    if (a > it.nused) continue;
                    if (a == it.nused && (dir & 1)) continue;
                    Pos q = it.c.p[it.l]; q.c[a] -= (dir & 1) ? -1 : 1;
                    bool hit = false;
                    for (int i = 0; i <= it.l; i++) if (eqp(q, it.c.p[i])) { hit = true; break; }
                    if (hit || it.l + 1 + l1(q) > K) continue;
                    Item n2 = it; n2.c.t[it.l] = dir; n2.c.p[it.l + 1] = q; n2.l = it.l + 1; n2.nused = a == it.nused ? it.nused + 1 : it.nused;
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
            enum_dfs(c, layer[i].l, layer[i].nused);
            auto &o = outs[omp_get_thread_num()]; o.insert(o.end(), c.out.begin(), c.out.end());
        }
        codes = shortcodes;
        for (auto &o : outs) { codes.insert(codes.end(), o.begin(), o.end()); std::vector<u64>().swap(o); }
        std::sort(codes.begin(), codes.end());
        size_t before = codes.size();
        codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
        if (codes.size() != before) { fprintf(stderr, "duplicate states?!\n"); return 3; }
    }
    if (getenv("SAW_REFINE")) {
        // hybrid state set: keep the level-(K-2) states (listed in the .imp file) and the level-K walks whose
        // level-(K-2) truncation is among the top fraction (by importance pi) of level-(K-2) states.
        double frac = getenv("SAW_REFINE_FRAC") ? atof(getenv("SAW_REFINE_FRAC")) : 0.1;
        FILE *f = fopen(getenv("SAW_REFINE"), "rb");
        u64 n2 = 0; if (!f || fread(&n2, 8, 1, f) != 1) { fprintf(stderr, "bad imp file\n"); return 3; }
        std::vector<u64> c2(n2); std::vector<double> pi2(n2);
        if (fread(c2.data(), 8, n2, f) != n2 || fread(pi2.data(), 8, n2, f) != n2) { fprintf(stderr, "bad imp file\n"); return 3; }
        fclose(f);
        std::vector<double> srt(pi2); std::sort(srt.begin(), srt.end(), std::greater<double>());
        size_t kk = std::min(n2 - 1, (size_t)(frac * n2));
        double thr = srt[kk];
        Lookup L2; L2.build(c2);
        std::vector<unsigned char> keep(codes.size(), 0);
#pragma omp parallel for schedule(dynamic, 4096)
        for (size_t i = 0; i < codes.size(); i++) {
            int t[MAXL]; int l = decode(codes[i], t);
            if (L2.find(codes[i]) >= 0) { keep[i] = 1; continue; }
            Pos p[MAXL + 1]; positions(l, t, p);
            for (int m = l - 1; m >= 0; m--) {
                int tt[MAXL]; for (int j = 0; j < m; j++) tt[j] = t[j];
                canon(m, tt);   // prefix of a canonical walk is canonical, kept for safety
                long id = L2.find(encode(m, tt));
                if (id >= 0) { keep[i] = pi2[id] >= thr ? 1 : 0; break; }
            }
        }
        std::vector<u64> nc; nc.reserve(codes.size());
        for (size_t i = 0; i < codes.size(); i++) if (keep[i]) nc.push_back(codes[i]);
        printf("refine: level-%d states %llu, top fraction %.4f (pi >= %.3e): kept %zu of %zu level-%d walks\n", K - 2,
               (unsigned long long)n2, frac, thr, nc.size(), codes.size(), K);
        codes.swap(nc);
    }
    size_t N = codes.size();
    u64 h1, h2; sha_like(codes, h1, h2);
    printf("d=%d K=%d policy=%d states=%zu  fingerprint=%016llx%016llx  enum %.1fs\n", D, K, POLICY, N,
           (unsigned long long)h1, (unsigned long long)h2, omp_get_wtime() - T0);
    fflush(stdout);
    Lookup LK; LK.build(codes);
    GLK = &LK;
    int startCode = -1;
    { int tt[1]; long s = LK.find(encode(0, tt)); startCode = (int)s; }
    if (startCode < 0) { fprintf(stderr, "no start state\n"); return 3; }
    // ---- children (parallel, chunked) ----
    double T1 = omp_get_wtime();
    std::vector<u32> ngroups(N);           // number of allowed directions per state
    std::vector<u64> rowg(N + 1);
    std::vector<u32> gsize;                // candidates per group
    std::vector<u32> cand;
    {
        const size_t CH = 1 << 14;
        size_t nch = (N + CH - 1) / CH;
        std::vector<std::vector<u32>> cg(nch), cc(nch);
        long missing = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : missing)
        for (size_t ch = 0; ch < nch; ch++) {
            std::vector<u64> kids;
            int t[MAXL];
            for (size_t s = ch * CH; s < std::min(N, (ch + 1) * CH); s++) {
                int l = decode(codes[s], t);
                u32 ng = 0;
                for (int dir = 0; dir < 2 * D; dir++) {
                    if (!children(l, t, dir, kids)) continue;
                    u32 cnt = 0;
                    for (u64 k : kids) { long id = LK.find(k); if (id < 0) { missing++; continue; } cc[ch].push_back((u32)id); cnt++; }
                    cg[ch].push_back(cnt); ng++;
                }
                ngroups[s] = ng;
            }
        }
        if (missing) { fprintf(stderr, "missing children: %ld\n", missing); return 3; }
        rowg[0] = 0; for (size_t s = 0; s < N; s++) rowg[s + 1] = rowg[s] + ngroups[s];
        size_t totg = rowg[N], totc = 0; for (auto &v : cc) totc += v.size();
        gsize.reserve(totg); cand.reserve(totc);
        for (size_t ch = 0; ch < nch; ch++) {
            gsize.insert(gsize.end(), cg[ch].begin(), cg[ch].end()); std::vector<u32>().swap(cg[ch]);
            cand.insert(cand.end(), cc[ch].begin(), cc[ch].end()); std::vector<u32>().swap(cc[ch]);
        }
    }
    std::vector<u64> gstart(rowg[N] + 1); gstart[0] = 0;
    for (size_t g = 0; g < rowg[N]; g++) gstart[g + 1] = gstart[g] + gsize[g];
    std::vector<u32>().swap(gsize);
    printf("groups=%llu candidates=%zu  children %.1fs\n", (unsigned long long)rowg[N], cand.size(), omp_get_wtime() - T1);
    fflush(stdout);
    if (getenv("SAW_PI_TEST")) {
        // experiment: policy iteration from the "first candidate" policy (or from v = proxy), full evaluation per pass
        std::vector<double> w(N, 1.0), nw(N);
        std::vector<u32> pol(rowg[N]);
        for (int pass = 0; pass < 8; pass++) {
#pragma omp parallel for schedule(static)
            for (size_t s = 0; s < N; s++)
                for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                    double m = 1e300; u32 arg = cand[gstart[g]];
                    for (u64 c = gstart[g]; c < gstart[g + 1]; c++) if (w[cand[c]] < m) { m = w[cand[c]]; arg = cand[c]; }
                    pol[g] = arg;
                }
            double lo = 0, hi = 0; int it;
            for (it = 0; it < 3000; it++) {
                double mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
                for (size_t s = 0; s < N; s++) { double sum = 0; for (u64 g = rowg[s]; g < rowg[s + 1]; g++) sum += w[pol[g]]; nw[s] = sum; mx = std::max(mx, sum); }
                hi = 0; lo = 1e300;
                for (size_t s = 0; s < N; s++) if (w[s] > 0) { double r = nw[s] / w[s]; hi = std::max(hi, r); if (nw[s] > 0) lo = std::min(lo, r); }
                for (size_t s = 0; s < N; s++) w[s] = nw[s] / mx;
                if (hi - lo < 1e-10) break;
            }
            printf("  PI pass %d: lambda(policy) in [%.10f, %.10f] after %d power its\n", pass, lo, hi, it);
            fflush(stdout);
        }
    }
    // ---- min-sum value iteration ----
    std::vector<double> v(N, 1.0), nv(N);
    double ub = 0, lb = 0;
    int it;
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
        for (size_t s = 0; s < N; s++) {
            if (v[s] > 0) { double r = nv[s] / v[s]; hi = std::max(hi, r); if (nv[s] > 0) lo = std::min(lo, r); }
        }
        ub = hi; lb = lo;
#pragma omp parallel for schedule(static)
        for (size_t s = 0; s < N; s++) v[s] = nv[s] / mx;
        if (it % 50 == 0) { printf("  it %d  [%.10f, %.10f]\n", it, lb, ub); fflush(stdout); }
        if (ub - lb < 2e-10) break;
    }
    printf("value iteration: it=%d  lambda in [%.10f, %.10f]  %.1fs\n", it, lb, ub, omp_get_wtime() - T0);
    fflush(stdout);
    // ---- integer vector and exact check ----
    // V_s = max(1, ceil(v_s / vmax * 2^50))
    double vmax = 0; for (size_t s = 0; s < N; s++) vmax = std::max(vmax, v[s]);
    std::vector<u64> V(N);
    for (size_t s = 0; s < N; s++) { double x = std::ceil(v[s] / vmax * 1125899906842624.0); V[s] = (u64)std::max(1.0, x); }
    std::vector<double>().swap(nv);
    // exact max ratio: find state maximising sum/V as rational; compare fractions with 128-bit
    u64 bestNum = 0, bestDen = 1;
#pragma omp parallel
    {
        u64 lNum = 0, lDen = 1;
#pragma omp for schedule(static)
        for (size_t s = 0; s < N; s++) {
            u64 sum = 0;
            for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                u64 m = ~0ull;
                for (u64 c = gstart[g]; c < gstart[g + 1]; c++) m = std::min(m, V[cand[c]]);
                sum += m;
            }
            if ((u128)sum * lDen > (u128)lNum * V[s]) { lNum = sum; lDen = V[s]; }
        }
#pragma omp critical
        { if ((u128)lNum * bestDen > (u128)bestNum * lDen) { bestNum = lNum; bestDen = lDen; } }
    }
    double ratio = (double)bestNum / (double)bestDen;
    // certified bound with 6 decimals: smallest P/10^6 >= bestNum/bestDen
    u64 Q = 1000000ull;
    u64 P = (u64)(((u128)bestNum * Q + bestDen - 1) / bestDen);
    printf("integer CW: max ratio = %llu/%llu = %.12f  => certified lambda_bar = %llu/%llu\n", (unsigned long long)bestNum,
           (unsigned long long)bestDen, ratio, (unsigned long long)P, (unsigned long long)Q);
    // ---- path counts from start (argmin policy w.r.t. V) vs c_n ----
    {
        std::vector<u32> pol(rowg[N]);
#pragma omp parallel for schedule(static)
        for (size_t s = 0; s < N; s++)
            for (u64 g = rowg[s]; g < rowg[s + 1]; g++) {
                u64 m = ~0ull; u32 arg = 0;
                for (u64 c = gstart[g]; c < gstart[g + 1]; c++) if (V[cand[c]] < m) { m = V[cand[c]]; arg = cand[c]; }
                pol[g] = arg;
            }
        std::vector<u128> cnt(N, 0), nc(N, 0);
        // counts of paths of length n from start = sum over states of number of length-n paths ending there
        cnt[startCode] = 1;
        printf("path counts a_n (automaton) :\n");
        for (int n = 1; n <= 24; n++) {
            std::fill(nc.begin(), nc.end(), (u128)0);
            for (size_t s = 0; s < N; s++) if (cnt[s]) for (u64 g = rowg[s]; g < rowg[s + 1]; g++) nc[pol[g]] += cnt[s];
            cnt.swap(nc);
            u128 tot = 0; for (size_t s = 0; s < N; s++) tot += cnt[s];
            char buf[64]; int bp = 63; buf[bp] = 0; u128 x = tot; do { buf[--bp] = '0' + (int)(x % 10); x /= 10; } while (x);
            printf("  n=%2d a_n=%s\n", n, buf + bp);
        }
        if (getenv("SAW_IMPORTANCE_OUT")) {
            // left Perron vector u of the argmin policy (power iteration on the transpose); pi_s = u_s V_s
            std::vector<double> u(N, 1.0), nu(N);
            for (int it2 = 0; it2 < 400; it2++) {
                std::fill(nu.begin(), nu.end(), 0.0);
                for (size_t s = 0; s < N; s++) for (u64 g = rowg[s]; g < rowg[s + 1]; g++) nu[pol[g]] += u[s];
                double mx = 0; for (size_t s = 0; s < N; s++) mx = std::max(mx, nu[s]);
                for (size_t s = 0; s < N; s++) u[s] = nu[s] / mx;
            }
            double tot = 0; for (size_t s = 0; s < N; s++) { nu[s] = u[s] * (double)V[s]; tot += nu[s]; }
            for (size_t s = 0; s < N; s++) nu[s] /= tot;
            FILE *f = fopen(getenv("SAW_IMPORTANCE_OUT"), "wb");
            u64 NN = N; fwrite(&NN, 8, 1, f); fwrite(codes.data(), 8, N, f); fwrite(nu.data(), 8, N, f); fclose(f);
            std::vector<double> sorted(nu); std::sort(sorted.begin(), sorted.end(), std::greater<double>());
            double acc = 0; size_t k1 = 0;
            for (double fr : {0.01, 0.05, 0.1, 0.2, 0.5}) { while (k1 < (size_t)(fr * N)) acc += sorted[k1++]; printf("  importance: top %.0f%% of states carry %.4f of pi\n", fr * 100, acc); }
        }
    }
    // ---- write certificate ----
    {
        std::string hp = outp + ".hdr", vp = outp + ".v64";
        FILE *f = fopen(hp.c_str(), "w");
        fprintf(f, "d %d\nK %d\npolicy %d\nstates %zu\nfingerprint %016llx%016llx\nlambda_num %llu\nlambda_den %llu\nmax_ratio %llu/%llu\n",
                D, K, POLICY, N, (unsigned long long)h1, (unsigned long long)h2, (unsigned long long)P, (unsigned long long)Q,
                (unsigned long long)bestNum, (unsigned long long)bestDen);
        fclose(f);
        f = fopen(vp.c_str(), "wb"); fwrite(V.data(), 8, N, f); fclose(f);
        // self-contained certificate: magic, d, K, N, P, Q, codes, weights
        std::string cp = outp + ".cert";
        f = fopen(cp.c_str(), "wb");
        u32 dd = (u32)D, kk = (u32)K; u64 NN = N;
        fwrite("SAWCERT1", 1, 8, f); fwrite(&dd, 4, 1, f); fwrite(&kk, 4, 1, f); fwrite(&NN, 8, 1, f);
        fwrite(&P, 8, 1, f); fwrite(&Q, 8, 1, f);
        fwrite(codes.data(), 8, N, f); fwrite(V.data(), 8, N, f);
        fclose(f);
    }
    printf("RESULT d=%d K=%d policy=%d states=%zu lambda_float=[%.10f,%.10f] certified<=%llu/%llu total %.1fs\n", D, K, POLICY, N, lb,
           ub, (unsigned long long)P, (unsigned long long)Q, omp_get_wtime() - T0);
    return 0;
}
