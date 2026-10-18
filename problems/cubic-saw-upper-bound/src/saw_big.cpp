// saw_big.cpp -- memory-lean version of saw_prod.cpp for large K (policy iteration, one stored child per
// (state, direction)).  Same states, same candidate children as saw_prod.cpp policy 4, same certificate format.
//
// Build:  g++ -O3 -march=native -fopenmp -o saw_big saw_big.cpp
// Usage:  saw_big d K outprefix passes [proxy.cert]
//   proxy.cert: certificate of level K-2 (same d); the initial policy picks, for every (state, direction), the
//   candidate whose level-(K-2) truncation has the smallest proxy weight.  Without proxy: first candidate.
//
// Pass structure: improvement (choose, for every (state, direction), the candidate of smallest current weight;
// keep the current child on ties) followed by evaluation (float power iteration for the fixed policy).  Pass 0
// computes the candidate lists and caches them on disk (<outprefix>.cand<thread>, deleted at the end); later
// passes read the cache.  The final evaluation is done in double precision, then V = ceil(v / vmax * 2^50) and
// the Collatz-Wielandt ratio is computed exactly.
// Output: <outprefix>.cert (SAWCERT1), <outprefix>.cert2 (SAWCERT2 = SAWCERT1 + explicit child table),
// <outprefix>.hdr; path counts a_n (n <= 22) on stdout.  Memory for d = 3, K = 18: about 10.5 GB.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>
#include <string>
#include <omp.h>
#include <sys/mman.h>
#include <unistd.h>
using u64 = uint64_t;
using u32 = uint32_t;
using u128 = unsigned __int128;

static int D, K, BITS;
static const int MAXL = 24;
static const int MAXD = 6;
static const u32 NONE = 0xffffffffu;

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

// ---------------- exact admissibility (used only for the enumeration of the state set) -------------------
struct Grid {
    int R, W; std::vector<u32> st; u32 stamp = 0;
    void init(int r) { R = r; W = 2 * r + 1; size_t n = 1; for (int a = 0; a < D; a++) n *= W; st.assign(n, 0); stamp = 0; }
    inline int idx(const Pos &p) const { int x = 0; for (int a = 0; a < D; a++) x = x * W + (p.c[a] + R); return x; }
};
static thread_local Grid G;
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
        int n[MAXD], sg[MAXD], tot = 1;
        for (int a = 0; a < D; a++) { int df = P[m].c[a] - P[0].c[a]; sg[a] = df >= 0 ? 1 : -1; n[a] = abs(df) + 1; tot *= n[a]; }
        static thread_local std::vector<unsigned char> R;
        if ((int)R.size() < tot) R.resize(tot);
        int off[MAXD] = {0, 0, 0, 0, 0};
        for (int id = 0; id < tot; id++) {
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
static inline bool admissible(int m, const Pos *P) {
    if (m == 0) return true;
    if (m + l1d(P[0], P[m]) > K) return false;
    return avoid_ok(m, P, K - m);
}
struct EnumCtx { int t[MAXL]; Pos p[MAXL + 1]; std::vector<u64> out; };
static void enum_dfs(EnumCtx &C, int l, int nused) {
    if (admissible(l, C.p)) C.out.push_back(encode(l, C.t));
    else if (l + l1(C.p[l]) > K) return;
    if (l + 1 + std::max(0, l1(C.p[l]) - 1) > K) return;
    for (int dir = 0; dir < 2 * D; dir++) {
        int a = dir >> 1;
        if (a > nused) continue;
        if (a == nused && (dir & 1)) continue;
        Pos q = C.p[l]; q.c[a] -= (dir & 1) ? -1 : 1;
        bool hit = false;
        for (int i = 0; i <= l; i++) if (eqp(q, C.p[i])) { hit = true; break; }
        if (hit) continue;
        if (l + 1 + l1(q) > K) continue;
        C.t[l] = dir; C.p[l + 1] = q;
        enum_dfs(C, l + 1, a == nused ? nused + 1 : nused);
    }
}

// ---------------- lookup ---------------------------------------------------------------------------------
// Open-addressing table of 64-bit entries: high (64 - sbits) bits = tag (high bits of a bijective hash of the
// code), low sbits bits = state index.  A probe reads the code array only when the tag matches, so a lookup costs
// about one cache miss (two for a hit).  Allocated with transparent huge pages.
static void *big_alloc(size_t bytes) {
    size_t al = (size_t)1 << 21, b = (bytes + al - 1) / al * al;
    void *p = aligned_alloc(al, b);
    if (!p) { fprintf(stderr, "out of memory (%zu bytes)\n", bytes); exit(5); }
    madvise(p, b, MADV_HUGEPAGE);
    return p;
}
struct Lookup {
    const u64 *codes = nullptr; size_t n = 0; u64 *tab = nullptr; u64 mask = 0, size = 0; int sbits = 0;
    static const u64 EMPTY = ~0ull;
    static inline u64 h(u64 x) { x ^= x >> 31; x *= 0x7fb5d329728ea185ull; x ^= x >> 27; x *= 0x81dadef4bc2dd44dull; x ^= x >> 33; return x; }
    void build(const u64 *c, size_t nn, double load) {
        codes = c; n = nn; size = 1; sbits = 0;
        while ((double)size * load < (double)nn + 16) { size <<= 1; sbits++; }
        mask = size - 1;
        tab = (u64 *)big_alloc(size * 8);
#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < size; j++) tab[j] = EMPTY;
        for (size_t i = 0; i < nn; i++) {
            u64 x = h(c[i]), j = x & mask;
            while (tab[j] != EMPTY) j = (j + 1) & mask;
            tab[j] = (x & ~mask) | (u64)i;
        }
    }
    void release() { if (tab) free(tab); tab = nullptr; }
    inline long find(u64 k) const {
        u64 x = h(k), j = x & mask, tag = x & ~mask;
        for (;;) {
            u64 e = tab[j];
            if (e == EMPTY) return -1;
            if ((e & ~mask) == tag && codes[e & mask] == k) return (long)(e & mask);
            j = (j + 1) & mask;
        }
    }
};
static Lookup LK;          // level K
static Lookup PLK;         // proxy level K-2 (optional)
static std::vector<u64> pcodes, pV;

static long lookup_prefix(const Lookup &L, int m, const Pos *P, int Klev) {
    if (m + l1d(P[0], P[m]) > Klev) return -1;
    int tt[MAXL + 2];
    for (int i = 0; i < m; i++) {
        int d2 = -1;
        for (int a = 0; a < D; a++) if (P[i].c[a] != P[i + 1].c[a]) d2 = 2 * a + (P[i].c[a] - P[i + 1].c[a] < 0 ? 1 : 0);
        tt[i] = d2;
    }
    canon(m, tt);
    return L.find(encode(m, tt));
}

// candidate children (indices into the level-K state table) for state (l,t), direction dir.
// returns false if the step is forbidden.  Candidate set = saw_prod.cpp policy 3/4.
// If proxyval != nullptr, also returns, per candidate, the proxy weight of its level-(K-2) truncation.
static bool children_idx(int l, const int *t, int dir, std::vector<u32> &out, std::vector<u64> *proxyval) {
    out.clear(); if (proxyval) proxyval->clear();
    Pos p[MAXL + 2]; positions(l, t, p);
    Pos e; for (int a = 0; a < D; a++) e.c[a] = 0; e.c[dir >> 1] = (dir & 1) ? -1 : 1;
    for (int i = 1; i <= l; i++) if (eqp(p[i], e)) return false;
    int n = l + 2;
    Pos q[MAXL + 3];
    for (int a = 0; a < D; a++) q[0].c[a] = 0;
    for (int i = 0; i <= l; i++) for (int a = 0; a < D; a++) q[i + 1].c[a] = p[i].c[a] - e.c[a];
    static thread_local std::vector<std::pair<u32, u32>> mc;   // (vertex mask, state index)
    static thread_local std::vector<std::pair<u32, u32>> keep;
    mc.clear();
    u32 adj[MAXL + 3];
    for (int i = 0; i < n; i++) { adj[i] = 0; for (int j = 0; j < n; j++) if (i != j && l1d(q[i], q[j]) == 1) adj[i] |= 1u << j; }
    int path[MAXL + 3]; long deep[MAXL + 3]; u32 rem[MAXL + 3]; Pos P[MAXL + 3];
    auto ext = [&](int depth, u32 used) {
        u32 a = adj[path[depth]] & ~used, r = 0;
        while (a) { int j = __builtin_ctz(a); a &= a - 1; if (depth + 1 + l1(q[j]) <= K) r |= 1u << j; }
        return r;
    };
    // deep[i]: -2 unknown, -1 not a state, >= 0 state index of prefix i
    auto leaf = [&](int depth) {
        for (int d2 = depth; d2 >= 0; d2--) {
            if (deep[d2] == -2) deep[d2] = lookup_prefix(LK, d2, P, K);
            if (deep[d2] >= 0) {
                u32 mask = 1; for (int i = 1; i <= d2; i++) mask |= 1u << path[i];
                mc.push_back({mask, (u32)deep[d2]});
                return;
            }
        }
    };
    int depth = 0; path[0] = 0; P[0] = q[0]; deep[0] = -2; u32 used = 1;
    rem[0] = ext(0, used);
    if (!rem[0]) leaf(0);
    while (depth >= 0) {
        if (!rem[depth]) { used &= ~(1u << path[depth]); depth--; continue; }
        int nx = __builtin_ctz(rem[depth]);
        rem[depth] &= rem[depth] - 1;
        depth++; path[depth] = nx; P[depth] = q[nx]; used |= 1u << nx;
        deep[depth] = -2;
        rem[depth] = ext(depth, used);
        if (!rem[depth]) leaf(depth);
    }
    // dominance pruning: drop (mask, index) pairs whose mask is strictly contained in another candidate's mask.
    // Exact duplicates are removed first; then candidates are scanned by decreasing cardinality and compared with
    // the maximal ones kept so far (strict inclusion is transitive, so this gives the same result as all pairs).
    std::sort(mc.begin(), mc.end());
    mc.erase(std::unique(mc.begin(), mc.end()), mc.end());
    std::stable_sort(mc.begin(), mc.end(), [](const std::pair<u32, u32> &a, const std::pair<u32, u32> &b) {
        return __builtin_popcount(a.first) > __builtin_popcount(b.first); });
    keep.clear();
    for (size_t i = 0; i < mc.size(); i++) {
        bool dom = false;
        for (size_t j = 0; j < keep.size() && !dom; j++)
            if ((mc[i].first & keep[j].first) == mc[i].first && mc[i].first != keep[j].first) dom = true;
        if (!dom) keep.push_back(mc[i]);
    }
    // dedupe by index
    std::sort(keep.begin(), keep.end(), [](const std::pair<u32, u32> &a, const std::pair<u32, u32> &b) { return a.second < b.second; });
    for (size_t i = 0; i < keep.size(); i++) {
        if (i && keep[i].second == keep[i - 1].second) continue;
        out.push_back(keep[i].second);
        if (proxyval) {
            // proxy: weight of the longest prefix (from the newest vertex) of the candidate that is a level-(K-2) state
            u64 code = LK.codes[keep[i].second]; int tt[MAXL]; int m = decode(code, tt);
            Pos R[MAXL + 1]; positions(m, tt, R);
            u64 pv = ~0ull;
            for (int mm = m; mm >= 0; mm--) { long id = lookup_prefix(PLK, mm, R, K - 2); if (id >= 0) { pv = pV[id]; break; } }
            proxyval->push_back(pv);
        }
    }
    return true;
}

static void u128str(u128 x, char *buf) { char tmp[64]; int p = 0; do { tmp[p++] = '0' + (int)(x % 10); x /= 10; } while (x); for (int i = 0; i < p; i++) buf[i] = tmp[p - 1 - i]; buf[p] = 0; }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s d K outprefix passes [proxy.cert]\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]); std::string outp = argv[3]; int PASSES = atoi(argv[4]);
    const char *proxyfile = argc > 5 ? argv[5] : nullptr;
    BITS = (2 * D <= 8) ? 3 : 4;
    const int ND = 2 * D - 1;   // stored directions for states with l >= 1: all except dir 1 (reversal)
    if ((K - 1) * BITS > 57 || K > MAXL) { fprintf(stderr, "K too large\n"); return 1; }
    double T0 = omp_get_wtime();
    auto stamp = [&]() { return omp_get_wtime() - T0; };
    // ---------- enumeration ----------
    std::vector<u64> codes;
    {
        const int S = std::min(5, K - 1);
        struct Item { EnumCtx c; int l, nused; };
        std::vector<Item> layer; Item it0; for (int a = 0; a < D; a++) it0.c.p[0].c[a] = 0; it0.l = 0; it0.nused = 0;
        for (int i = 0; i < MAXL; i++) it0.c.t[i] = 0;
        layer.push_back(it0);
        std::vector<u64> shortcodes;
        for (int l = 0; l < S; l++) {
            std::vector<Item> nxt;
            for (auto &it : layer) {
                if (admissible(it.l, it.c.p)) shortcodes.push_back(encode(it.l, it.c.t));
                for (int dir = 0; dir < 2 * D; dir++) {
                    int a = dir >> 1;
                    if (a > it.nused || (a == it.nused && (dir & 1))) continue;
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
        size_t tot = shortcodes.size(); for (auto &o : outs) tot += o.size();
        codes.reserve(tot);
        codes = shortcodes;
        for (auto &o : outs) { codes.insert(codes.end(), o.begin(), o.end()); std::vector<u64>().swap(o); }
        __gnu_parallel::sort(codes.begin(), codes.end());
        for (size_t i = 1; i < codes.size(); i++) if (codes[i] == codes[i - 1]) { fprintf(stderr, "duplicate\n"); return 3; }
    }
    const size_t N = codes.size();
    printf("d=%d K=%d states=%zu  enum %.1fs\n", D, K, N, stamp()); fflush(stdout);
    if (N >= NONE) { fprintf(stderr, "too many states\n"); return 3; }
    if (codes[0] != encode(0, nullptr)) { fprintf(stderr, "empty walk not first\n"); return 3; }
    LK.build(codes.data(), N, 0.8);
    printf("lookup built (%llu slots) %.1fs\n", (unsigned long long)LK.size, stamp()); fflush(stdout);
    // ---------- proxy ----------
    bool haveProxy = false;
    if (proxyfile) {
        FILE *f = fopen(proxyfile, "rb");
        char mg[8]; u32 dd, kk; u64 nn, pp, qq;
        if (f && fread(mg, 1, 8, f) == 8 && fread(&dd, 4, 1, f) == 1 && fread(&kk, 4, 1, f) == 1 && fread(&nn, 8, 1, f) == 1 &&
            fread(&pp, 8, 1, f) == 1 && fread(&qq, 8, 1, f) == 1 && (int)dd == D && (int)kk == K - 2) {
            pcodes.resize(nn); pV.resize(nn);
            if (fread(pcodes.data(), 8, nn, f) == nn && fread(pV.data(), 8, nn, f) == nn) haveProxy = true;
        }
        if (f) fclose(f);
        if (!haveProxy) { fprintf(stderr, "bad proxy file\n"); return 3; }
        PLK.build(pcodes.data(), pcodes.size(), 0.5);
        printf("proxy level %d loaded: %zu states\n", K - 2, pcodes.size()); fflush(stdout);
    }
    // ---------- policy storage ----------
    // pol[s*ND + j], j indexes directions {0, 2, 3, ..., 2D-1}; NONE = forbidden.  Start state (index 0, empty walk)
    // has all 2D directions; its children are all the one-step walk, handled separately.
    u32 *pol = (u32 *)big_alloc((size_t)N * ND * 4);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < (size_t)N * ND; i++) pol[i] = NONE;
    std::vector<float> vf(N, 1.0f), nvf(N);
    long oneStep = LK.find(([&]() { int tt[1] = {0}; return encode(1, tt); })());
    if (oneStep < 0) { fprintf(stderr, "no one-step state\n"); return 3; }
    auto dirOf = [&](int j) { return j == 0 ? 0 : j + 1; };
    double lamLo = 0, lamHi = 0;
    // Candidate cache: pass 0 computes the candidate lists once and appends them, chunk by chunk, to one file per
    // thread; later passes read them back (pread) instead of recomputing.  Record per (state, direction):
    // u16 count (0xffff = forbidden step) followed by count u32 state indices.
    const size_t CH = 4096;
    const size_t nch = (N + CH - 1) / CH;
    const int nth = omp_get_max_threads();
    std::vector<FILE *> cf(nth);
    for (int th = 0; th < nth; th++) {
        cf[th] = fopen((outp + ".cand" + std::to_string(th)).c_str(), "w+b");
        if (!cf[th]) { fprintf(stderr, "cannot open candidate cache\n"); return 3; }
    }
    struct ChunkLoc { int thr; u64 off, len; };
    std::vector<ChunkLoc> loc(nch);
    for (int pass = 0; pass <= PASSES; pass++) {
        // ---------- improvement ----------
        double t1 = stamp();
        long changed = 0, totalCand = 0;
        static long done = 0; done = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : changed, totalCand)
        for (size_t ch = 0; ch < nch; ch++) {
            std::vector<u32> kids; std::vector<u64> pvals;
            std::vector<unsigned char> buf;
            int t[MAXL];
            int tid = omp_get_thread_num();
            size_t rp = 0;
            if (pass > 0) {
                buf.resize(loc[ch].len);
                if (loc[ch].len && pread(fileno(cf[loc[ch].thr]), buf.data(), loc[ch].len, (off_t)loc[ch].off) != (ssize_t)loc[ch].len) {
                    fprintf(stderr, "cache read error\n"); exit(6);
                }
            }
            for (size_t s = std::max<size_t>(1, ch * CH); s < std::min(N, (ch + 1) * CH); s++) {
                int l = (pass == 0) ? decode(codes[s], t) : 0;
                for (int j = 0; j < ND; j++) {
                    int dir = dirOf(j);
                    bool useProxy = (pass == 0 && haveProxy);
                    if (pass == 0) {
                        bool ok = children_idx(l, t, dir, kids, useProxy ? &pvals : nullptr);
                        uint16_t cnt = ok ? (uint16_t)kids.size() : (uint16_t)0xffff;
                        if (ok && kids.empty()) { fprintf(stderr, "no child for state %zu dir %d\n", s, dir); exit(4); }
                        if (ok && kids.size() >= 0xffff) { fprintf(stderr, "too many candidates\n"); exit(4); }
                        size_t o = buf.size(); buf.resize(o + 2 + (ok ? 4 * kids.size() : 0));
                        memcpy(&buf[o], &cnt, 2); if (ok) memcpy(&buf[o + 2], kids.data(), 4 * kids.size());
                        if (!ok) { pol[s * ND + j] = NONE; continue; }
                    } else {
                        uint16_t cnt; memcpy(&cnt, &buf[rp], 2); rp += 2;
                        if (cnt == 0xffff) { pol[s * ND + j] = NONE; continue; }
                        kids.resize(cnt); memcpy(kids.data(), &buf[rp], 4 * (size_t)cnt); rp += 4 * (size_t)cnt;
                    }
                    totalCand += kids.size();
                    u32 cur = pol[s * ND + j], best = kids[0];
                    if (useProxy) {
                        u64 bv = pvals[0];
                        for (size_t i = 1; i < kids.size(); i++) if (pvals[i] < bv) { bv = pvals[i]; best = kids[i]; }
                    } else if (pass > 0) {
                        float bv = vf[kids[0]];
                        for (size_t i = 1; i < kids.size(); i++) if (vf[kids[i]] < bv) { bv = vf[kids[i]]; best = kids[i]; }
                        if (cur != NONE && vf[cur] <= bv) best = cur;   // keep current child on ties
                    }
                    if (best != cur) changed++;
                    pol[s * ND + j] = best;
                }
            }
            if (pass == 0) {
                loc[ch].thr = tid; loc[ch].off = (u64)ftello(cf[tid]); loc[ch].len = buf.size();
                if (buf.size() && fwrite(buf.data(), 1, buf.size(), cf[tid]) != buf.size()) { fprintf(stderr, "cache write error\n"); exit(6); }
            }
            { long dn = __atomic_add_fetch(&done, 1, __ATOMIC_RELAXED); long step = (long)std::max<size_t>(1, nch / 20);
              if (dn % step == 0) { printf("    pass %d: %ld/%zu chunks  %.1fs\n", pass, dn, nch, stamp()); fflush(stdout); } }
        }
        if (pass == 0) {
            for (int th = 0; th < nth; th++) fflush(cf[th]);
            LK.release(); PLK.release(); std::vector<u64>().swap(pcodes); std::vector<u64>().swap(pV);   // not needed any more
        }
        printf("pass %d: improvement %.1fs, changed %ld, candidates %ld\n", pass, stamp() - t1, changed, totalCand); fflush(stdout);
        // ---------- evaluation (float) ----------
        double t2 = stamp();
        int it;
        for (it = 1; it <= 3000; it++) {
            float mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
            for (size_t s = 1; s < N; s++) {
                float sum = 0; const u32 *pp = &pol[s * ND];
                for (int j = 0; j < ND; j++) if (pp[j] != NONE) sum += vf[pp[j]];
                nvf[s] = sum; mx = std::max(mx, sum);
            }
            nvf[0] = (float)(2 * D) * vf[oneStep];
            double hi = 0, lo = 1e300;
#pragma omp parallel for schedule(static) reduction(max : hi) reduction(min : lo)
            for (size_t s = 1; s < N; s++) if (vf[s] > 0) { double r = (double)nvf[s] / vf[s]; hi = std::max(hi, r); if (nvf[s] > 0) lo = std::min(lo, r); }
#pragma omp parallel for schedule(static)
            for (size_t s = 0; s < N; s++) vf[s] = nvf[s] / mx;
            lamLo = lo; lamHi = hi;
            if (hi - lo < 8e-6 || it >= 400) break;
        }
        printf("pass %d: evaluation %d its %.1fs, lambda(policy) in [%.8f, %.8f]\n", pass, it, stamp() - t2, lamLo, lamHi); fflush(stdout);
        if (pass > 0 && changed == 0) break;
    }
    // ---------- final evaluation in double ----------
    for (int th = 0; th < nth; th++) { fclose(cf[th]); remove((outp + ".cand" + std::to_string(th)).c_str()); }
    std::vector<float>().swap(nvf);
    std::vector<double> v(N), nv(N);
    for (size_t s = 0; s < N; s++) v[s] = vf[s];
    std::vector<float>().swap(vf);
    double lo = 0, hi = 0; int it;
    for (it = 1; it <= 5000; it++) {
        double mx = 0;
#pragma omp parallel for schedule(static) reduction(max : mx)
        for (size_t s = 1; s < N; s++) {
            double sum = 0; const u32 *pp = &pol[s * ND];
            for (int j = 0; j < ND; j++) if (pp[j] != NONE) sum += v[pp[j]];
            nv[s] = sum; mx = std::max(mx, sum);
        }
        nv[0] = (double)(2 * D) * v[oneStep];
        hi = 0; lo = 1e300;
#pragma omp parallel for schedule(static) reduction(max : hi) reduction(min : lo)
        for (size_t s = 1; s < N; s++) if (v[s] > 0) { double r = nv[s] / v[s]; hi = std::max(hi, r); if (nv[s] > 0) lo = std::min(lo, r); }
#pragma omp parallel for schedule(static)
        for (size_t s = 0; s < N; s++) v[s] = nv[s] / mx;
        if (it % 20 == 0) { printf("  final it %d [%.11f, %.11f]\n", it, lo, hi); fflush(stdout); }
        if (hi - lo < 2e-10) break;
    }
    printf("final evaluation: %d its, lambda in [%.11f, %.11f]  %.1fs\n", it, lo, hi, stamp()); fflush(stdout);
    // ---------- integer certificate ----------
    std::vector<double>().swap(nv);
    double vmax = 0; for (size_t s = 0; s < N; s++) vmax = std::max(vmax, v[s]);
    std::vector<u64> V(N);
    for (size_t s = 0; s < N; s++) { double x = std::ceil(v[s] / vmax * 1125899906842624.0); V[s] = (u64)std::max(1.0, x); }
    std::vector<double>().swap(v);
    // start state: V[0] must satisfy 2D*V[oneStep] <= lambda*V[0]; set V[0] = 2D*V[oneStep] (ratio 1 < lambda)
    V[0] = (u64)(2 * D) * V[oneStep];
    u64 bestNum = 0, bestDen = 1;
#pragma omp parallel
    {
        u64 lNum = 0, lDen = 1;
#pragma omp for schedule(static)
        for (size_t s = 1; s < N; s++) {
            u64 sum = 0; const u32 *pp = &pol[s * ND];
            for (int j = 0; j < ND; j++) if (pp[j] != NONE) sum += V[pp[j]];
            if ((u128)sum * lDen > (u128)lNum * V[s]) { lNum = sum; lDen = V[s]; }
        }
#pragma omp critical
        { if ((u128)lNum * bestDen > (u128)bestNum * lDen) { bestNum = lNum; bestDen = lDen; } }
    }
    u64 Q = 1000000ull;
    u64 P = (u64)(((u128)bestNum * Q + bestDen - 1) / bestDen);
    printf("integer CW: max ratio = %llu/%llu = %.12f  => certified lambda_bar = %llu/%llu\n", (unsigned long long)bestNum,
           (unsigned long long)bestDen, (double)bestNum / bestDen, (unsigned long long)P, (unsigned long long)Q);
    fflush(stdout);
    {
        std::string cp = outp + ".cert";
        FILE *f = fopen(cp.c_str(), "wb");
        u32 dd = (u32)D, kk = (u32)K; u64 NN = N;
        fwrite("SAWCERT1", 1, 8, f); fwrite(&dd, 4, 1, f); fwrite(&kk, 4, 1, f); fwrite(&NN, 8, 1, f);
        fwrite(&P, 8, 1, f); fwrite(&Q, 8, 1, f);
        fwrite(codes.data(), 8, N, f); fwrite(V.data(), 8, N, f);
        fclose(f);
        // SAWCERT2: same header and arrays, followed by the explicit automaton: for every state and every one of the
        // 2D directions (order 0..2D-1), the index of the chosen child, or 0xffffffff when no child is given.
        FILE *f2 = fopen((outp + ".cert2").c_str(), "wb");
        fwrite("SAWCERT2", 1, 8, f2); fwrite(&dd, 4, 1, f2); fwrite(&kk, 4, 1, f2); fwrite(&NN, 8, 1, f2);
        fwrite(&P, 8, 1, f2); fwrite(&Q, 8, 1, f2);
        fwrite(codes.data(), 8, N, f2); fwrite(V.data(), 8, N, f2);
        {
            std::vector<u32> row(2 * D);
            const size_t BLK = 1 << 16;
            std::vector<u32> outb; outb.reserve(BLK * 2 * D);
            for (size_t s = 0; s < N; s++) {
                for (int dir = 0; dir < 2 * D; dir++) {
                    u32 c;
                    if (s == 0) c = (u32)oneStep;
                    else if (dir == 1) c = NONE;
                    else c = pol[s * ND + (dir == 0 ? 0 : dir - 1)];
                    outb.push_back(c);
                }
                if (outb.size() >= BLK * 2 * D || s + 1 == N) { fwrite(outb.data(), 4, outb.size(), f2); outb.clear(); }
            }
        }
        fclose(f2);
        FILE *h = fopen((outp + ".hdr").c_str(), "w");
        fprintf(h, "d %d\nK %d\nstates %zu\nlambda_num %llu\nlambda_den %llu\nmax_ratio %llu/%llu\nfloat_lambda [%.11f, %.11f]\n", D, K, N,
                (unsigned long long)P, (unsigned long long)Q, (unsigned long long)bestNum, (unsigned long long)bestDen, lo, hi);
        fclose(h);
    }
    printf("certificate written %.1fs\n", stamp()); fflush(stdout);
    // ---------- path counts (u64; n <= 22) ----------
    std::vector<u64>().swap(V); std::vector<u64>().swap(codes);
    {
        std::vector<u64> cnt(N, 0), nc(N, 0);
        cnt[0] = 1;
        for (int n = 1; n <= 22; n++) {
            std::fill(nc.begin(), nc.end(), 0ull);
            nc[oneStep] += (u64)(2 * D) * cnt[0];
            for (size_t s = 1; s < N; s++) if (cnt[s]) { const u32 *pp = &pol[s * ND]; for (int j = 0; j < ND; j++) if (pp[j] != NONE) nc[pp[j]] += cnt[s]; }
            cnt.swap(nc);
            u128 tot = 0; for (size_t s = 0; s < N; s++) tot += cnt[s];
            char buf[64]; u128str(tot, buf);
            printf("  n=%2d a_n=%s\n", n, buf); fflush(stdout);
        }
    }
    printf("RESULT d=%d K=%d states=%zu lambda_float=[%.11f,%.11f] certified<=%llu/%llu total %.1fs\n", D, K, N, lo, hi,
           (unsigned long long)P, (unsigned long long)Q, stamp());
    return 0;
}
