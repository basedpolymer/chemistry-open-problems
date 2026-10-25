// latsaw_check.cpp -- checker for LATCERT1 certificates (diamond and face-centred cubic lattices).
// Written from the mathematical specification (paper.tex, Sections 2-4), separately from the generator
// latsaw_gen.cpp: it shares no code with it (no canonical forms, no hashing of walks, no search for children).
//
// Build:  g++ -O3 -march=native -std=c++17 -pthread -o latsaw_check latsaw_check.cpp
// Usage:  latsaw_check <file.lcert> [threads=4] [counts.txt nmax]
//   counts.txt (optional): lines "n c_n" (exact numbers of n-step self-avoiding walks); the checker then also
//   computes a_n = number of n-step label sequences accepted by the certified automaton and reports a_n >= c_n
//   (a consistency test of the geometry, not part of the proof).
//
// FILE FORMAT (little endian):
//   char[8] "LATCERT1"; u32 lattice (0 = diamond, 1 = fcc); u32 K (informative only);
//   u64 N; u64 P; u64 Q      (claimed bound mu <= P/Q)
//   N x u64 codes (strictly increasing); N x u64 weights V; N x NDIR x u32 child index (0xffffffff = none)
// CODE of a state: bit 63 set; walk length l in bits 57..62; step label t_i in bits [b*i, b*i+b), b = 2 (diamond)
//   or 4 (fcc), all other bits zero.
// GEOMETRY: fcc vertices = {x in Z^3: x1+x2+x3 even}, neighbours of x: x + u_t, t = 0..11, in the order
//   u = (-1,-1,0),(-1,1,0),(1,-1,0),(1,1,0),(-1,0,-1),(-1,0,1),(1,0,-1),(1,0,1),(0,-1,-1),(0,-1,1),(0,1,-1),(0,1,1).
//   diamond vertices: type A = {x in (2Z)^3: x1+x2+x3 = 0 mod 4}, type B = A + (1,1,1); a_0..a_3 =
//   (1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1); neighbours: x + a_t (x of type A), x - a_t (x of type B).
//   The origin o is the current vertex (type A for diamond).  Positions of a state: p_0 = o,
//   p_{i+1} = p_i - s_i, s_i = u_{t_i} (fcc), s_i = (-1)^{i+1} a_{t_i} (diamond).  F(s) = {p_1, ..., p_l}.
// CLAIMS CHECKED (exactly the hypotheses of Lemmas 1 and 2 of the paper):
//   (C0) format; codes strictly increasing and decodable; 1 <= V <= 2^60; P, Q >= 1; the empty walk is listed;
//   (C1) for every listed state s and every neighbour e of o with e not in F(s): the listed child c exists, and
//        some gamma(x) = e + g x, g in the group Gamma_e of signed permutation matrices for which gamma is a
//        lattice automorphism with gamma(o) = e, satisfies gamma(F(c)) subset of F(s) u {o}  (point by point);
//   (C2) Q * sum_e V(child(s,e)) <= P * V(s) for every listed state, exact unsigned 128-bit arithmetic.
// Gamma_e is built here by brute force: all 48 signed permutation matrices g are tried and g is kept iff
// x -> e + g x maps every vertex of a box of radius 6 to a vertex and every edge inside the box to an edge
// (and is checked to contain 48 elements for fcc and 24 for diamond); that such a map is a global automorphism
// is proved in the paper (it is affine and preserves the vertex set and the neighbour vectors).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <array>
#include <thread>
#include <atomic>
#include <mutex>
#include <string>
#include <chrono>
typedef uint64_t u64;
typedef uint32_t u32;
typedef unsigned __int128 u128;

static int LATT, NSTEP, SB;   // lattice id, number of step labels, bits per label
static int STEPV[12][3];
struct Mat { int m[3][3]; };
static std::vector<Mat> GAM[12];   // Gamma_e for every neighbour e of the origin

static bool vertexQ(const int *x) {
    if (LATT == 1) return ((x[0] + x[1] + x[2]) % 2 + 2) % 2 == 0;
    int r = ((x[0] + x[1] + x[2]) % 4 + 4) % 4;
    bool allEven = x[0] % 2 == 0 && x[1] % 2 == 0 && x[2] % 2 == 0;
    bool allOdd = x[0] % 2 != 0 && x[1] % 2 != 0 && x[2] % 2 != 0;
    return (allEven && r == 0) || (allOdd && r == 3);
}
static bool typeA(const int *x) { return x[0] % 2 == 0; }   // diamond only, x a vertex
static bool edgeQ(const int *x, const int *y) {             // x, y vertices
    int d[3] = {y[0] - x[0], y[1] - x[1], y[2] - x[2]};
    for (int t = 0; t < NSTEP; t++) {
        int s = (LATT == 0 && !typeA(x)) ? -1 : 1;
        if (d[0] == s * STEPV[t][0] && d[1] == s * STEPV[t][1] && d[2] == s * STEPV[t][2]) return true;
    }
    return false;
}
static void apply(const Mat &g, const int *e, const int *x, int *y) {
    for (int r = 0; r < 3; r++) y[r] = e[r] + g.m[r][0] * x[0] + g.m[r][1] * x[1] + g.m[r][2] * x[2];
}
static void setup() {
    if (LATT == 1) {
        int U[12][3] = {{-1,-1,0},{-1,1,0},{1,-1,0},{1,1,0},{-1,0,-1},{-1,0,1},{1,0,-1},{1,0,1},{0,-1,-1},{0,-1,1},{0,1,-1},{0,1,1}};
        memcpy(STEPV, U, sizeof(U)); NSTEP = 12; SB = 4;
    } else {
        int A[4][3] = {{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
        memcpy(STEPV, A, sizeof(A)); NSTEP = 4; SB = 2;
    }
    // all 48 signed permutation matrices
    std::vector<Mat> all;
    int pr[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (int p = 0; p < 6; p++) for (int sg = 0; sg < 8; sg++) {
        Mat g; memset(g.m, 0, sizeof(g.m));
        for (int c = 0; c < 3; c++) g.m[pr[p][c]][c] = (sg >> c & 1) ? -1 : 1;
        all.push_back(g);
    }
    const int R = 6;
    for (int j = 0; j < NSTEP; j++) {
        const int *e = STEPV[j];   // neighbour of the origin (the origin is of type A for diamond)
        for (auto &g : all) {
            int o[3] = {0, 0, 0}, img[3];
            apply(g, e, o, img);
            bool ok = img[0] == e[0] && img[1] == e[1] && img[2] == e[2];
            for (int x0 = -R; x0 <= R && ok; x0++) for (int x1 = -R; x1 <= R && ok; x1++) for (int x2 = -R; x2 <= R && ok; x2++) {
                int x[3] = {x0, x1, x2};
                if (!vertexQ(x)) continue;
                int y[3]; apply(g, e, x, y);
                if (!vertexQ(y)) { ok = false; break; }
                for (int t = 0; t < NSTEP && ok; t++) {      // every edge {x, x'} with x' in the box
                    int s = (LATT == 0 && !typeA(x)) ? -1 : 1;
                    int x2v[3] = {x[0] + s * STEPV[t][0], x[1] + s * STEPV[t][1], x[2] + s * STEPV[t][2]};
                    int y2[3]; apply(g, e, x2v, y2);
                    if (!edgeQ(y, y2)) ok = false;
                }
            }
            if (ok) GAM[j].push_back(g);
        }
    }
}

// decode a code into positions; returns length or -1 if malformed
static int decodePos(u64 c, int pos[][3]) {
    if (!(c >> 63)) return -1;
    int l = (int)((c >> 57) & 63);
    if (l * SB > 57) return -1;
    u64 payload = c & ((1ULL << 57) - 1);
    if (l * SB < 57 && (payload >> (l * SB)) != 0) return -1;
    pos[0][0] = pos[0][1] = pos[0][2] = 0;
    for (int i = 0; i < l; i++) {
        int t = (int)((c >> (SB * i)) & ((1u << SB) - 1));
        if (t >= NSTEP) return -1;
        int sgn = (LATT == 1) ? 1 : ((i % 2 == 1) ? 1 : -1);   // s_i = (-1)^{i+1} a_t for diamond
        for (int a = 0; a < 3; a++) pos[i + 1][a] = pos[i][a] - sgn * STEPV[t][a];
    }
    return l;
}

static const int BR = 40, BW = 2 * BR + 1;   // membership box
static inline int bidx(const int *x) { return ((x[0] + BR) * BW + (x[1] + BR)) * BW + (x[2] + BR); }
static inline bool inBox(const int *x) { for (int a = 0; a < 3; a++) if (x[a] < -BR || x[a] > BR) return false; return true; }

struct Res { u64 bad_code = 0, missing = 0, contain_fail = 0, ineq_fail = 0, extra = 0, nonsaw = 0, pairs = 0; u64 bn = 0, bd = 1; u64 worst = 0; };

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s cert.lcert [threads] [counts.txt nmax]\n", argv[0]); return 2; }
    int threads = argc > 2 ? atoi(argv[2]) : 4; if (threads < 1) threads = 1;
    auto T0 = std::chrono::steady_clock::now();
    FILE *f = fopen(argv[1], "rb");
    if (!f) { printf("FAIL cannot open %s\n", argv[1]); return 1; }
    char mg[9] = {0}; u32 hl[2]; u64 hq[3];
    if (fread(mg, 1, 8, f) != 8 || fread(hl, 4, 2, f) != 2 || fread(hq, 8, 3, f) != 3) { printf("FAIL header\n"); return 1; }
    if (std::string(mg) != "LATCERT1") { printf("FAIL magic\n"); return 1; }
    LATT = (int)hl[0]; u64 N = hq[0], P = hq[1], Q = hq[2];
    if (LATT != 0 && LATT != 1) { printf("FAIL lattice id\n"); return 1; }
    setup();
    printf("certificate %s: lattice %s, K = %u, N = %llu, claimed mu <= %llu/%llu = %.9f\n", argv[1], LATT ? "fcc" : "diamond", hl[1],
           (unsigned long long)N, (unsigned long long)P, (unsigned long long)Q, (double)P / (double)Q);
    bool gamOK = true;
    for (int j = 0; j < NSTEP; j++) if ((int)GAM[j].size() != (LATT ? 48 : 24)) gamOK = false;
    printf("%s group Gamma_e has %zu elements for every neighbour e (expected %d)\n", gamOK ? "PASS" : "FAIL", GAM[0].size(), LATT ? 48 : 24);
    if (!gamOK) return 1;
    if (P == 0 || Q == 0 || P >= (1ULL << 60) || Q >= (1ULL << 60) || N == 0 || N >= 0xffffffffULL) { printf("FAIL header values\n"); return 1; }
    fseeko(f, 0, SEEK_END); long long fsz = ftello(f); fseeko(f, 40, SEEK_SET);
    long long expect = 40 + 16LL * (long long)N + 4LL * NSTEP * (long long)N;
    if (fsz != expect) { printf("FAIL file size %lld, expected %lld\n", fsz, expect); return 1; }
    std::vector<u64> code(N), V(N); std::vector<u32> ch((size_t)N * NSTEP);
    if (fread(code.data(), 8, N, f) != N || fread(V.data(), 8, N, f) != N || fread(ch.data(), 4, ch.size(), f) != ch.size()) { printf("FAIL read\n"); return 1; }
    fclose(f);
    // (C0)
    bool inc = true, wok = true; long long start = -1;
    for (u64 i = 0; i < N; i++) {
        if (i > 0 && code[i] <= code[i - 1]) inc = false;
        if (V[i] < 1 || V[i] > (1ULL << 60)) wok = false;
        if (code[i] == (1ULL << 63)) start = (long long)i;
    }
    printf("%s (C0a) codes strictly increasing\n", inc ? "PASS" : "FAIL");
    printf("%s (C0b) all weights in [1, 2^60]\n", wok ? "PASS" : "FAIL");
    printf("%s (C0c) the empty walk is listed (index %lld)\n", start >= 0 ? "PASS" : "FAIL", start);
    bool all = inc && wok && start >= 0;
    // (C1), (C2)
    std::vector<Res> res(threads);
    std::atomic<u64> next(0);
    auto work = [&](int w) {
        Res &r = res[w];
        std::vector<uint8_t> box((size_t)BW * BW * BW, 0);
        int ps[64][3], pc[64][3];
        for (;;) {
            u64 lo = next.fetch_add(4096); if (lo >= N) break;
            u64 hi = std::min<u64>(N, lo + 4096);
            for (u64 s = lo; s < hi; s++) {
                int l = decodePos(code[s], ps);
                if (l < 0) { r.bad_code++; continue; }
                for (int i = 0; i <= l; i++) for (int k = i + 1; k <= l; k++)
                    if (ps[i][0] == ps[k][0] && ps[i][1] == ps[k][1] && ps[i][2] == ps[k][2]) { r.nonsaw++; i = l + 1; break; }
                for (int i = 0; i <= l; i++) box[bidx(ps[i])] = 1;     // F(s) u {o}
                u128 sum = 0;
                for (int j = 0; j < NSTEP; j++) {
                    const int *e = STEPV[j];
                    bool forb = false;
                    for (int i = 1; i <= l; i++) if (ps[i][0] == e[0] && ps[i][1] == e[1] && ps[i][2] == e[2]) forb = true;
                    u32 c = ch[(size_t)s * NSTEP + j];
                    if (forb) { if (c != 0xffffffffu) r.extra++; continue; }
                    r.pairs++;
                    if (c >= N) { r.missing++; continue; }
                    int lc = decodePos(code[c], pc);
                    if (lc < 0) { r.missing++; continue; }
                    bool found = false;
                    for (const Mat &g : GAM[j]) {
                        bool ok = true;
                        for (int i = 1; i <= lc && ok; i++) {
                            int y[3]; apply(g, e, pc[i], y);
                            ok = inBox(y) && box[bidx(y)];
                        }
                        if (ok) { found = true; break; }
                    }
                    if (!found) { r.contain_fail++; continue; }
                    sum += V[c];
                }
                for (int i = 0; i <= l; i++) box[bidx(ps[i])] = 0;
                if ((u128)Q * sum > (u128)P * V[s]) r.ineq_fail++;
                u64 sm = (u64)sum;
                if ((u128)sm * r.bd > (u128)r.bn * V[s]) { r.bn = sm; r.bd = V[s]; r.worst = s; }
            }
        }
    };
    std::vector<std::thread> th;
    for (int w = 0; w < threads; w++) th.emplace_back(work, w);
    for (auto &t : th) t.join();
    Res tot;
    for (auto &r : res) {
        tot.bad_code += r.bad_code; tot.missing += r.missing; tot.contain_fail += r.contain_fail; tot.ineq_fail += r.ineq_fail;
        tot.extra += r.extra; tot.nonsaw += r.nonsaw; tot.pairs += r.pairs;
        if ((u128)r.bn * tot.bd > (u128)tot.bn * r.bd) { tot.bn = r.bn; tot.bd = r.bd; tot.worst = r.worst; }
    }
    printf("%s (C0d) every code decodes (malformed: %llu)\n", tot.bad_code == 0 ? "PASS" : "FAIL", (unsigned long long)tot.bad_code);
    printf("%s (C1a) every allowed (state, step) has a listed child (missing: %llu; pairs %llu)\n", tot.missing == 0 ? "PASS" : "FAIL",
           (unsigned long long)tot.missing, (unsigned long long)tot.pairs);
    printf("%s (C1b) gamma(F(child)) in F(s) u {o} for some automorphism gamma with gamma(o) = e (failures: %llu)\n",
           tot.contain_fail == 0 ? "PASS" : "FAIL", (unsigned long long)tot.contain_fail);
    u64 a = tot.bn, b = tot.bd; while (b) { u64 t = a % b; a = b; b = t; }
    u64 gn = tot.bn / a, gd = tot.bd / a;
    u128 ip = gn / gd, rm = gn % gd; std::string dec;
    for (int i = 0; i < 12; i++) { rm *= 10; dec += char('0' + (int)(rm / gd)); rm %= gd; }
    printf("%s (C2) Q*sum_e V(child) <= P*V(s) for all states (failures: %llu); exact max ratio = %llu/%llu = %llu.%s...\n",
           tot.ineq_fail == 0 ? "PASS" : "FAIL", (unsigned long long)tot.ineq_fail, (unsigned long long)gn, (unsigned long long)gd,
           (unsigned long long)ip, dec.c_str());
    printf("info: listed walks that are not self-avoiding: %llu; child entries on forbidden steps: %llu\n",
           (unsigned long long)tot.nonsaw, (unsigned long long)tot.extra);
    all = all && tot.bad_code == 0 && tot.missing == 0 && tot.contain_fail == 0 && tot.ineq_fail == 0;
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    printf("time %.1f s, %d threads\n", secs, threads);
    // optional: a_n >= c_n
    if (argc > 4 && all) {
        std::vector<u128> cn(64, 0); int nmax = atoi(argv[4]);
        FILE *fc = fopen(argv[3], "r"); char line[256];
        while (fc && fgets(line, sizeof line, fc)) {
            if (line[0] == '#') continue;
            int n; char num[200];
            if (sscanf(line, "%d %199s", &n, num) == 2 && n >= 0 && n < 64) { u128 v = 0; for (char *p = num; *p; p++) v = v * 10 + (*p - '0'); cn[n] = v; }
        }
        if (fc) fclose(fc);
        std::vector<u128> x(N, 1), y(N);
        bool cok = true; int neq = 0;
        for (int n = 1; n <= nmax; n++) {
            for (u64 s = 0; s < N; s++) { u128 acc = 0; for (int j = 0; j < NSTEP; j++) { u32 c = ch[(size_t)s * NSTEP + j]; if (c != 0xffffffffu) acc += x[c]; } y[s] = acc; }
            x.swap(y);
            if (cn[n] == 0) continue;
            if (x[start] < cn[n]) cok = false;
            if (x[start] == cn[n] && neq == n - 1) neq = n;
        }
        printf("%s (info) a_n >= c_n for n <= %d; a_n = c_n for n <= %d\n", cok ? "PASS" : "FAIL", nmax, neq);
        all = all && cok;
    }
    if (all) printf("PASS  => mu(%s) <= %llu/%llu = %.9f\n", LATT ? "fcc" : "diamond", (unsigned long long)P, (unsigned long long)Q, (double)P / (double)Q);
    else printf("FAIL\n");
    return all ? 0 : 1;
}
