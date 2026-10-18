// verify_cert.cpp -- independent checker for the SAW connective-constant certificates.
//
// Build:  g++ -O3 -march=native -fopenmp -o verify_cert verify_cert.cpp
// Usage:  verify_cert <file.cert> [threads] [--full]
//
// Certificate file (little endian):
//   char[8]  "SAWCERT1"
//   u32 d, u32 K, u64 N, u64 P, u64 Q          (claimed bound lambda_bar = P/Q)
//   N x u64  state codes (strictly increasing)
//   N x u64  weights V (all >= 1)
//   SAWCERT2 only: N x 2d x u32  index of the chosen child for (state, direction e = 0..2d-1); 0xffffffff = none.
//   With SAWCERT2 no search is done: (C1) is checked for the listed child.
// A state code encodes a step sequence t[0..l-1] (newest step first, BITS bits per step, BITS = 3 if 2d <= 8
// else 4; length l in bits 57..62; bit 63 set).  Direction t = 2*axis + (1 if negative).  Positions relative to
// the current vertex A = 0:  p_0 = 0, p_{i+1} = p_i - unit(t[i]).  Forbidden set F(s) = {p_1, ..., p_l}.
//
// What is checked (this is exactly the hypothesis of Lemma 1 in notes.md):
//   (C0) header sane, codes strictly increasing, the empty walk (F = empty set) is a listed state, all V >= 1;
//   (C1) for every listed state s and every unit vector e with e not in F(s), the program exhibits a listed
//        state t and a signed permutation g with g(F(t)) subset of W(s,e) := (F(s) u {0}) - e   (checked
//        explicitly, point by point), choosing among the found children one of minimal weight;
//   (C2) Q * sum_e V(t(s,e)) <= P * V(s)  in exact 128-bit integer arithmetic.
// The way children are searched for (paths through W starting at the new position, canonical relabelling, hash
// lookup) only serves to FIND candidates; soundness rests on the explicit containment test of (C1).
// Prints PASS/FAIL per claim; exit code 0 iff everything passes.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <omp.h>
typedef uint64_t u64;
typedef uint32_t u32;
typedef unsigned __int128 u128;

static int d, K, BITS;
static bool FULL = false;   // --full: look up every path prefix (slower), default: longest listed prefix per branch
static const int ML = 40;

struct V5 { int x[6]; };   // a point of Z^d, d <= 6

static int decode(u64 c, int *t) {   // returns -1 for a malformed code (too long for the 57 payload bits)
    int l = (int)((c >> 57) & 63);
    if (l * BITS > 57 || l >= ML - 2) return -1;
    for (int i = 0; i < l; i++) t[i] = (int)((c >> (BITS * i)) & ((1u << BITS) - 1));
    return l;
}
static u64 encode(int l, const int *t) {
    u64 c = (1ull << 63) | ((u64)l << 57);
    for (int i = 0; i < l; i++) c |= (u64)t[i] << (BITS * i);
    return c;
}
static void walkpos(int l, const int *t, V5 *p) {
    for (int a = 0; a < 6; a++) p[0].x[a] = 0;
    for (int i = 0; i < l; i++) {
        p[i + 1] = p[i];
        int ax = t[i] / 2, sg = (t[i] % 2) ? -1 : 1;
        p[i + 1].x[ax] -= sg;
    }
}
static bool same(const V5 &a, const V5 &b) { for (int i = 0; i < d; i++) if (a.x[i] != b.x[i]) return false; return true; }
static int dist1(const V5 &a, const V5 &b) { int s = 0; for (int i = 0; i < d; i++) s += abs(a.x[i] - b.x[i]); return s; }

// canonical relabelling used only to find candidates in the table
static void relabel(int l, int *t) {
    int newax[6] = {-1, -1, -1, -1, -1, -1}, neg[6] = {0, 0, 0, 0, 0, 0}, k = 0;
    for (int i = 0; i < l; i++) {
        int ax = t[i] / 2, s = t[i] % 2;
        if (newax[ax] == -1) { newax[ax] = k++; neg[ax] = s; }
        t[i] = 2 * newax[ax] + (s ^ neg[ax]);
    }
}

// all signed permutations of {0..d-1}
static std::vector<std::pair<std::vector<int>, std::vector<int>>> group;
static void make_group() {
    std::vector<int> perm(d); for (int i = 0; i < d; i++) perm[i] = i;
    do {
        for (int m = 0; m < (1 << d); m++) {
            std::vector<int> sg(d); for (int i = 0; i < d; i++) sg[i] = (m >> i & 1) ? -1 : 1;
            group.push_back({perm, sg});
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
}

struct Table {
    // open addressing; entry = (high bits of a bijective hash of the code) | index; the code array is read only
    // when the high bits match.  Load factor <= 0.8.
    std::vector<u64> codes; std::vector<u64> w; std::vector<u64> h; u64 mask;
    static u64 mix(u64 x) { x += 0x9e3779b97f4a7c15ull; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull; return x ^ (x >> 31); }
    void build() { u64 n = 2; while ((double)n * 0.8 < (double)codes.size() + 16) n <<= 1; h.assign(n, ~0ull); mask = n - 1;
        for (size_t i = 0; i < codes.size(); i++) { u64 x = mix(codes[i]), j = x & mask; while (h[j] != ~0ull) j = (j + 1) & mask; h[j] = (x & ~mask) | i; } }
    long get(u64 c) const { u64 x = mix(c), j = x & mask, tag = x & ~mask;
        for (;;) { u64 e = h[j]; if (e == ~0ull) return -1; if ((e & ~mask) == tag && codes[e & mask] == c) return (long)(e & mask); j = (j + 1) & mask; } }
};
static Table T;

// containment test: exists g with g(F(t)) subset of Wset (Wset given as list of points)
static bool contained(u64 childcode, const V5 *W, int nW) {
    int tt[ML]; int m = decode(childcode, tt);
    if (m < 0) return false;
    for (int i = 0; i < m; i++) if (tt[i] >= 2 * d) return false;
    V5 r[ML + 1]; walkpos(m, tt, r);
    if (d <= 3) {
        // small groups: plain enumeration of all 2^d d! signed permutations
        for (auto &g : group) {
            bool ok = true;
            for (int i = 1; i <= m && ok; i++) {
                V5 im; for (int a = 0; a < 6; a++) im.x[a] = 0;
                for (int a = 0; a < d; a++) im.x[g.first[a]] = g.second[a] * r[i].x[a];
                bool found = false;
                for (int j = 0; j < nW; j++) if (same(im, W[j])) { found = true; break; }
                ok = found;
            }
            if (ok) return true;
        }
        return false;
    }
    // larger groups: exhaustive backtracking over partial signed permutations (axis a -> sign * axis tgt[a]);
    // the image of an axis is chosen when a vertex first uses it, each vertex is tested as soon as all its axes
    // are assigned.  Explores every signed permutation compatible with the vertices tested so far, hence finds a
    // valid g if and only if one exists.
    int tgt[6], sg[6]; bool usedT[6];
    for (int a = 0; a < 6; a++) { tgt[a] = -1; sg[a] = 1; usedT[a] = false; }
    struct Rec {
        const V5 *r; const V5 *W; int nW, m; int *tgt, *sg; bool *usedT;
        bool go(int i) {
            if (i > m) return true;
            int na = -1;
            for (int a = 0; a < d; a++) if (r[i].x[a] != 0 && tgt[a] < 0) { na = a; break; }
            if (na >= 0) {
                for (int b = 0; b < d; b++) {
                    if (usedT[b]) continue;
                    for (int s = -1; s <= 1; s += 2) {
                        tgt[na] = b; sg[na] = s; usedT[b] = true;
                        if (go(i)) return true;
                        usedT[b] = false; tgt[na] = -1;
                    }
                }
                return false;
            }
            V5 im; for (int a = 0; a < 6; a++) im.x[a] = 0;
            for (int a = 0; a < d; a++) if (r[i].x[a] != 0) im.x[tgt[a]] = sg[a] * r[i].x[a];
            bool found = false;
            for (int j = 0; j < nW; j++) if (same(im, W[j])) { found = true; break; }
            return found && go(i + 1);
        }
    } R{r, W, nW, m, tgt, sg, usedT};
    return R.go(1);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s file.cert [threads]\n", argv[0]); return 2; }
    if (argc > 2) omp_set_num_threads(atoi(argv[2]));
    if (argc > 3 && strcmp(argv[3], "--full") == 0) FULL = true;
    FILE *f = fopen(argv[1], "rb");
    if (!f) { printf("FAIL cannot open certificate\n"); return 1; }
    char magic[8]; u32 dd, KK; u64 N, P, Q;
    if (fread(magic, 1, 8, f) != 8) { printf("FAIL bad magic\n"); return 1; }
    bool EXPLICIT = memcmp(magic, "SAWCERT2", 8) == 0;   // SAWCERT2: the chosen children are listed explicitly
    if (!EXPLICIT && memcmp(magic, "SAWCERT1", 8) != 0) { printf("FAIL bad magic\n"); return 1; }
    if (fread(&dd, 4, 1, f) != 1 || fread(&KK, 4, 1, f) != 1 || fread(&N, 8, 1, f) != 1 || fread(&P, 8, 1, f) != 1 || fread(&Q, 8, 1, f) != 1) { printf("FAIL header\n"); return 1; }
    d = (int)dd; K = (int)KK; BITS = (2 * d <= 8) ? 3 : 4;
    printf("certificate: d=%d K=%d states=%llu claimed lambda_bar=%llu/%llu = %.9f  (search mode: %s)\n", d, K, (unsigned long long)N,
           (unsigned long long)P, (unsigned long long)Q, (double)P / (double)Q, EXPLICIT ? "explicit children (SAWCERT2)" : (FULL ? "full" : "leaf"));
    bool allok = true;
    if (d < 2 || d > 6 || Q == 0 || N == 0 || N > 0xfffffff0ull) { printf("FAIL header values\n"); return 1; }
    T.codes.resize(N); T.w.resize(N);
    if (fread(T.codes.data(), 8, N, f) != N || fread(T.w.data(), 8, N, f) != N) { printf("FAIL truncated file\n"); return 1; }
    std::vector<u32> kids;
    if (EXPLICIT) {
        kids.resize((size_t)N * 2 * d);
        if (fread(kids.data(), 4, kids.size(), f) != kids.size()) { printf("FAIL truncated child table\n"); return 1; }
    }
    fclose(f);
    // (C0)
    bool inc = true, pos = true;
    for (u64 i = 0; i < N; i++) { if (i && T.codes[i] <= T.codes[i - 1]) inc = false; if (T.w[i] == 0 || T.w[i] > (1ull << 60)) pos = false; }
    int t0[1]; long start = -1;
    if (!EXPLICIT) { T.build(); start = T.get(encode(0, t0)); }
    else { auto it = std::lower_bound(T.codes.begin(), T.codes.end(), encode(0, t0)); if (it != T.codes.end() && *it == encode(0, t0)) start = it - T.codes.begin(); }
    printf("%s (C0a) codes strictly increasing\n", inc ? "PASS" : "FAIL");
    printf("%s (C0b) all weights in [1, 2^60]\n", pos ? "PASS" : "FAIL");
    printf("%s (C0c) the empty walk is a listed state\n", start >= 0 ? "PASS" : "FAIL");
    allok = allok && inc && pos && start >= 0;
    make_group();
    // (C1),(C2)
    long nochild = 0, notcont = 0, ineq = 0, badcode = 0;
    double worst = 0; u64 worstNum = 0, worstDen = 1;
#pragma omp parallel
    {
        double lworst = 0; u64 lNum = 0, lDen = 1;
#pragma omp for schedule(dynamic, 4096) reduction(+ : nochild, notcont, ineq, badcode)
        for (u64 s = 0; s < N; s++) {
            int t[ML]; int l = decode(T.codes[s], t);
            if (l < 0) { badcode++; continue; }
            bool okc = true; for (int i = 0; i < l; i++) if (t[i] >= 2 * d) okc = false;
            if (!okc) { badcode++; continue; }
            V5 p[ML + 1]; walkpos(l, t, p);
            u64 sum = 0;
            for (int e = 0; e < 2 * d; e++) {
                V5 ev; for (int a = 0; a < 6; a++) ev.x[a] = 0; ev.x[e / 2] = (e % 2) ? -1 : 1;
                bool forb = false;
                for (int i = 1; i <= l; i++) if (same(p[i], ev)) forb = true;
                if (forb) continue;                      // a self-avoiding walk cannot take this step
                // W = (F(s) u {0}) - e, listed with the new position (origin) first
                V5 W[ML + 2]; int nW = 0;
                for (int i = 0; i <= l; i++) { V5 w; for (int a = 0; a < 6; a++) w.x[a] = p[i].x[a] - ev.x[a]; W[nW++] = w; }
                V5 O; for (int a = 0; a < 6; a++) O.x[a] = 0;
                if (EXPLICIT) {
                    // explicit automaton: the child is given, only (C1b) has to be checked
                    u32 c = kids[(size_t)s * 2 * d + e];
                    if (c >= N) { nochild++; sum = ~0ull >> 8; break; }
                    if (!contained(T.codes[c], W, nW)) notcont++;
                    sum += T.w[c];
                    continue;
                }
                // search: self-avoiding paths O = P0, P1, ... through points of W (prefixes with m + |P_m| <= K)
                u64 best = ~0ull, bestcode = 0;
                int stackIdx[ML + 2], nextTry[ML + 2]; V5 P[ML + 2]; bool usedW[ML + 2] = {false};
                int depth = 0; P[0] = O; nextTry[0] = 0;
                {   // the empty walk itself (m = 0)
                    long id = T.get(encode(0, t0));
                    if (id >= 0 && T.w[id] < best) { best = T.w[id]; bestcode = T.codes[id]; }
                }
                long known[ML + 2];   // leaf mode: lookup result for prefix of length i (-2 = not computed yet)
                known[0] = -2;
                while (depth >= 0) {
                    if (nextTry[depth] >= nW) { if (depth > 0) usedW[stackIdx[depth]] = false; depth--; continue; }
                    int j = nextTry[depth]++;
                    if (usedW[j] || dist1(P[depth], W[j]) != 1) continue;
                    int m = depth + 1;
                    if (m + dist1(O, W[j]) > K) continue;
                    usedW[j] = true; stackIdx[m] = j; P[m] = W[j]; nextTry[m] = 0; depth = m; known[m] = -2;
                    // candidate: walk P_m (oldest) ... P_0 (newest); newest-first steps
                    int tt[ML];
                    for (int i = 0; i < m; i++) {
                        int dir = -1;
                        for (int a = 0; a < d; a++) { int df = P[i].x[a] - P[i + 1].x[a]; if (df == 1) dir = 2 * a; else if (df == -1) dir = 2 * a + 1; }
                        tt[i] = dir;
                    }
                    relabel(m, tt);   // the relabelling of a prefix is the prefix of the relabelling
                    if (FULL) {
                        long id = T.get(encode(m, tt));
                        if (id >= 0 && T.w[id] < best) { best = T.w[id]; bestcode = T.codes[id]; }
                    } else {
                        // leaf mode: only at branch ends, the longest listed prefix of the branch
                        bool ext = false;
                        for (int j2 = 0; j2 < nW && !ext; j2++)
                            if (!usedW[j2] && dist1(P[m], W[j2]) == 1 && m + 1 + dist1(O, W[j2]) <= K) ext = true;
                        if (ext) continue;
                        for (int d2 = m; d2 >= 0; d2--) {
                            if (known[d2] == -2) known[d2] = T.get(encode(d2, tt));
                            if (known[d2] >= 0) {
                                long id = known[d2];
                                if (T.w[id] < best) { best = T.w[id]; bestcode = T.codes[id]; }
                                break;
                            }
                        }
                    }
                }
                if (best == ~0ull) { nochild++; sum = ~0ull >> 8; break; }
                if (!contained(bestcode, W, nW)) { notcont++; }
                sum += best;
            }
            if ((u128)sum * Q > (u128)P * T.w[s]) ineq++;
            if ((u128)sum * lDen > (u128)lNum * T.w[s]) { lNum = sum; lDen = T.w[s]; }
        }
#pragma omp critical
        { if ((u128)lNum * worstDen > (u128)worstNum * lDen) { worstNum = lNum; worstDen = lDen; } }
        (void)lworst;
    }
    worst = (double)worstNum / (double)worstDen;
    printf("%s (C0d) every code decodes to a valid step sequence (bad: %ld)\n", badcode == 0 ? "PASS" : "FAIL", badcode);
    printf("%s (C1a) every allowed step has a listed child (missing: %ld)\n", nochild == 0 ? "PASS" : "FAIL", nochild);
    printf("%s (C1b) every chosen child satisfies g(F(t)) in (F(s) u {0}) - e for some signed permutation g (violations: %ld)\n",
           notcont == 0 ? "PASS" : "FAIL", notcont);
    printf("%s (C2) Q*sum_e V(child) <= P*V(s) for all states (violations: %ld); max ratio = %llu/%llu = %.12f\n",
           ineq == 0 ? "PASS" : "FAIL", ineq, (unsigned long long)worstNum, (unsigned long long)worstDen, worst);
    allok = allok && nochild == 0 && notcont == 0 && ineq == 0 && badcode == 0;
    if (allok) printf("PASS  => mu(Z^%d) <= lambda_bar = %llu/%llu = %.9f\n", d, (unsigned long long)P, (unsigned long long)Q, (double)P / (double)Q);
    else printf("FAIL\n");
    return allok ? 0 : 1;
}
