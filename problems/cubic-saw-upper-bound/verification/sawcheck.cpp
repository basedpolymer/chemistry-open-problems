// sawcheck.cpp -- independent verifier for the SAWCERT1 / SAWCERT2 certificates of notes.md Sec. 5.
// Written by the independent verifier from the mathematical specification (Lemmas 1-2) and the format
// description only; it shares no code with src/verify_cert.cpp (which was not read before this program ran).
//
// What it checks (all integer arithmetic, no floating point in any decision):
//  (F) format: magic, header, 40-byte header + 16N (+ 8dN) bytes, codes strictly increasing, every code
//      decodable (bit 63 set, length*bits <= 57, no stray payload bits, step values < 2d), 1 <= V <= 2^60,
//      P, Q in [1, 2^60);
//  (S) the start state (empty walk, F = {}) is listed;
//  (I) for every listed state s and every direction e not in F(s): a listed child t and an explicit signed
//      permutation g with g F(t) subset of (F(s) u {0}) - e, found by exhaustive backtracking over partial
//      signed permutations (complete search over the hyperoctahedral group, restricted to the axes t uses);
//      SAWCERT2: t is the stored child; SAWCERT1: t is searched among ALL self-avoiding paths through
//      W = (F(s) u {0}) - e (optionally pruned), looked up by an independently written canonical form, and the
//      one of minimal weight is kept; the inclusion is then re-tested explicitly with the backtracking test;
//  (C) Q * sum_e V[t(s,e)] <= P * V[s] for every listed s, exact (unsigned __int128); the exact maximal ratio
//      max_s sum_e V[t]/V[s] is reported as a reduced fraction.
//
// Build: g++ -O2 -march=native -std=c++17 -pthread sawcheck.cpp -o sawcheck
// Usage: sawcheck <cert> [threads=4] [mode=full|prefix] [maxstates] [nmax for a_n counts]
//  mode "full": SAWCERT1 child search over all self-avoiding paths through W (only the safe pruning
//  k+|P_k|_1 <= max over listed states of l+|p_l|_1); "prefix": stop extending a path whose prefix is not listed.
//  Any pruning can only make the check fail, never pass wrongly. Exit code 0 iff PASS.
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64;
typedef long long i64;
typedef unsigned __int128 u128;
typedef uint32_t u32;

static int D = 0, B = 0;          // dimension, bits per step
static i64 DELTA[16];             // packed unit vector of direction j (2a -> +e_a, 2a+1 -> -e_a)
static const int SH = 7;          // bits per packed coordinate
static i64 ORIGIN = 0;            // packed origin (each coordinate offset by 64)
static const u64 LIMIT60 = 1ULL << 60;

static inline int lenOf(u64 c) { return (int)((c >> 57) & 63); }
static inline int stepOf(u64 c, int i) { return (int)((c >> (B * i)) & ((1ULL << B) - 1)); }

// decode a code into positions pos[0..L] (pos[0] = origin). Returns false if not decodable.
static bool decode(u64 c, i64* pos, int& L) {
    if (!(c >> 63)) return false;
    L = lenOf(c);
    if (L * B > 57) return false;
    u64 payload = c & ((1ULL << 57) - 1);
    if (L * B < 57 && (payload >> (L * B))) return false;
    i64 p = ORIGIN;
    pos[0] = p;
    for (int i = 0; i < L; i++) {
        int t = stepOf(c, i);
        if (t >= 2 * D) return false;
        p -= DELTA[t];
        pos[i + 1] = p;
    }
    return true;
}

static int l1norm(i64 p) {
    int s = 0;
    for (int a = 0; a < D; a++) {
        int x = (int)((p >> (SH * a)) & 127) - 64;
        s += x < 0 ? -x : x;
    }
    return s;
}

// ---------------- inclusion test: exists g in G with g F(t) subset of W ----------------
struct GSearch {
    const int* u; int m;            // child's steps
    const i64* W; int nW;           // target set
    int img[8], sg[8]; bool used[8];
    bool inW(i64 p) const { for (int i = 0; i < nW; i++) if (W[i] == p) return true; return false; }
    bool rec(int k, i64 pos) {
        if (k == m) return true;
        int a = u[k] >> 1, b = u[k] & 1;
        if (img[a] >= 0) {
            int dir = 2 * img[a] + (b ^ sg[a]);
            i64 np = pos - DELTA[dir];
            return inW(np) && rec(k + 1, np);
        }
        for (int c = 0; c < D; c++) {
            if (used[c]) continue;
            for (int s = 0; s < 2; s++) {
                int dir = 2 * c + (b ^ s);
                i64 np = pos - DELTA[dir];
                if (!inW(np)) continue;
                img[a] = c; sg[a] = s; used[c] = true;
                if (rec(k + 1, np)) return true;
                img[a] = -1; used[c] = false;
            }
        }
        return false;
    }
    bool run(const int* steps, int mm, const i64* w, int nw) {
        u = steps; m = mm; W = w; nW = nw;
        for (int a = 0; a < 8; a++) { img[a] = -1; sg[a] = 0; used[a] = false; }
        return rec(0, ORIGIN);
    }
};

// ---------------- hash table code -> index (SAWCERT1 search mode) ----------------
struct Table {
    vector<u64> key; vector<u32> val; u64 mask = 0;
    static inline u64 h(u64 x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }
    void build(const vector<u64>& codes) {
        u64 sz = 1; while (sz < 2 * codes.size() + 16) sz <<= 1;
        key.assign(sz, 0); val.assign(sz, 0); mask = sz - 1;
        for (u64 i = 0; i < codes.size(); i++) {
            u64 p = h(codes[i]) & mask;
            while (key[p]) p = (p + 1) & mask;
            key[p] = codes[i]; val[p] = (u32)i;
        }
    }
    inline long long find(u64 c) const {
        u64 p = h(c) & mask;
        while (key[p]) { if (key[p] == c) return val[p]; p = (p + 1) & mask; }
        return -1;
    }
};

struct Global {
    vector<u64> codes, V;
    u64 N = 0, P = 0, Q = 0; int K = 0;
    bool explicitCh = false;
    Table table;
    int Kmax = 0, Lmax = 0;  // max over listed states of l+|p_l|_1 and of l
    bool prefixMode = false;
};
static Global G;

struct Stats {
    u64 pairs = 0, fails_child = 0, fails_incl = 0, fails_ineq = 0, nonSAW = 0, extraPtr = 0, candidates = 0;
    u64 bestNum = 0, bestDen = 1; u64 worstState = 0;
    double searchNodes = 0;
    vector<string> msgs;
    void note(const string& s) { if (msgs.size() < 20) msgs.push_back(s); }
};

// ---------------- child search through W (SAWCERT1) ----------------
struct Searcher {
    const i64* W; int nW; int adj[24][16]; int norm[24];
    int bestIdx; u64 bestV; u64 cands; double nodes;
    int axisImg[8], axisSg[8];
    void setup(const i64* w, int nw) {
        W = w; nW = nw;
        for (int i = 0; i < nW; i++) {
            norm[i] = l1norm(W[i]);
            for (int dir = 0; dir < 2 * D; dir++) {
                i64 np = W[i] - DELTA[dir]; adj[i][dir] = -1;
                for (int j = 0; j < nW; j++) if (W[j] == np) { adj[i][dir] = j; break; }
            }
        }
    }
    void dfs(int k, i64 pos, int cur, unsigned mask, u64 code, int nextLabel) {
        nodes += 1;
        if (k >= G.Lmax) return;
        for (int dir = 0; dir < 2 * D; dir++) {
            int j;
            if (cur < 0) { i64 np = pos - DELTA[dir]; j = -1; for (int q = 0; q < nW; q++) if (W[q] == np) { j = q; break; } }
            else j = adj[cur][dir];
            if (j < 0 || (mask >> j & 1)) continue;
            if (k + 1 + norm[j] > G.Kmax) continue;  // safe: k+|p_k|_1 is non-decreasing along any walk
            int a = dir >> 1, b = dir & 1;
            bool fresh = axisImg[a] < 0; int nl = nextLabel;
            if (fresh) { axisImg[a] = nl++; axisSg[a] = b; }
            int cdir = 2 * axisImg[a] + (b ^ axisSg[a]);
            u64 nc = (code & ~(63ULL << 57)) | ((u64)(k + 1) << 57) | ((u64)cdir << (B * k));
            long long idx = G.table.find(nc);
            if (idx >= 0) {
                cands++;
                if (G.V[idx] < bestV) { bestV = G.V[idx]; bestIdx = (int)idx; }
            }
            if (idx >= 0 || !G.prefixMode) dfs(k + 1, W[j], j, mask | (1u << j), nc, nl);
            if (fresh) axisImg[a] = -1;
        }
    }
    int run(u64 startIdx) {
        bestIdx = (int)startIdx; bestV = G.V[startIdx]; cands = 1; nodes = 0;
        for (int a = 0; a < 8; a++) { axisImg[a] = -1; axisSg[a] = 0; }
        dfs(0, ORIGIN, -1, 0, 1ULL << 63, 0);
        return bestIdx;
    }
};

static u64 START = 0;  // index of start state
static int NCOUNT = 0;             // if > 0: also compute a_n = (A^n 1)_{start}, n <= NCOUNT, for the checked automaton
static vector<u32> CHOSEN;         // child actually checked for (s, dir), 0xFFFFFFFF if dir is forbidden/failed

static void processRange(u64 lo, u64 hi, const u32* ch /* may be null; ch[(s-lo)*2D+j] */, Stats& st) {
    i64 pos[64], W[64], cpos[64]; int L, cL; int csteps[64];
    GSearch gs; Searcher se;
    for (u64 s = lo; s < hi; s++) {
        u64 c = G.codes[s];
        decode(c, pos, L);  // validated beforehand
        // self-avoidance (informative only)
        { bool sa = true; for (int i = 0; i <= L && sa; i++) for (int j = i + 1; j <= L; j++) if (pos[i] == pos[j]) { sa = false; break; }
          if (!sa) st.nonSAW++; }
        u128 sum = 0;
        for (int dir = 0; dir < 2 * D; dir++) {
            i64 e = DELTA[dir];
            bool forb = false;
            for (int i = 1; i <= L; i++) if (pos[i] == ORIGIN + e) { forb = true; break; }
            if (forb) { if (ch && ch[(s - lo) * 2 * D + dir] != 0xFFFFFFFFu) st.extraPtr++; continue; }
            st.pairs++;
            // W = (F(s) u {0}) - e, deduplicated
            int nW = 0;
            for (int i = 0; i <= L; i++) { i64 w = pos[i] - e; bool dup = false; for (int q = 0; q < nW; q++) if (W[q] == w) { dup = true; break; } if (!dup) W[nW++] = w; }
            long long t;
            if (ch) {
                u32 x = ch[(s - lo) * 2 * D + dir];
                if (x == 0xFFFFFFFFu || x >= G.N) { st.fails_child++; st.note("missing/invalid child s=" + to_string(s) + " dir=" + to_string(dir)); continue; }
                t = x;
            } else {
                se.setup(W, nW);
                t = se.run(START);
                st.candidates += se.cands; st.searchNodes += se.nodes;
            }
            // explicit inclusion test
            decode(G.codes[t], cpos, cL);
            for (int i = 0; i < cL; i++) csteps[i] = stepOf(G.codes[t], i);
            if (!gs.run(csteps, cL, W, nW)) { st.fails_incl++; st.note("inclusion fails s=" + to_string(s) + " dir=" + to_string(dir) + " child=" + to_string(t)); continue; }
            sum += G.V[t];
            if (NCOUNT) CHOSEN[s * 2 * D + dir] = (u32)t;
        }
        // (C): Q*sum <= P*V[s]
        u128 lhs = (u128)G.Q * sum, rhs = (u128)G.P * G.V[s];
        if (lhs > rhs) { st.fails_ineq++; st.note("inequality fails s=" + to_string(s)); }
        u64 sm = (u64)sum;  // sum <= 12*2^60 < 2^64
        if ((u128)sm * st.bestDen > (u128)st.bestNum * G.V[s]) { st.bestNum = sm; st.bestDen = G.V[s]; st.worstState = s; }
    }
}

static u64 gcdu(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s cert [threads] [full|prefix] [maxstates]\n", argv[0]); return 2; }
    int threads = argc > 2 ? atoi(argv[2]) : 4;
    threads = max(1, min(4, threads));
    G.prefixMode = argc > 3 && string(argv[3]) == "prefix";
    auto T0 = chrono::steady_clock::now();
    FILE* f = fopen(argv[1], "rb");
    if (!f) { printf("FAIL cannot open\n"); return 1; }
    char magic[9] = {0};
    u32 dK[2]; u64 hdr[3];
    if (fread(magic, 1, 8, f) != 8 || fread(dK, 4, 2, f) != 2 || fread(hdr, 8, 3, f) != 3) { printf("FAIL header\n"); return 1; }
    string mg(magic);
    if (mg != "SAWCERT1" && mg != "SAWCERT2") { printf("FAIL magic %s\n", magic); return 1; }
    G.explicitCh = (mg == "SAWCERT2");
    D = (int)dK[0]; G.K = (int)dK[1]; G.N = hdr[0]; G.P = hdr[1]; G.Q = hdr[2];
    if (D < 2 || D > 6) { printf("FAIL d=%d\n", D); return 1; }
    B = (D <= 4) ? 3 : 4;
    if (G.P == 0 || G.Q == 0 || G.P >= LIMIT60 || G.Q >= LIMIT60) { printf("FAIL P,Q range\n"); return 1; }
    if (G.N == 0 || G.N >= 0xFFFFFFFFULL) { printf("FAIL N range\n"); return 1; }
    // file size
    fseeko(f, 0, SEEK_END); off_t fsz = ftello(f); fseeko(f, 40, SEEK_SET);
    off_t expect = 40 + (off_t)16 * G.N + (G.explicitCh ? (off_t)8 * D * G.N : 0);
    if (fsz != expect) { printf("FAIL file size %lld expected %lld\n", (long long)fsz, (long long)expect); return 1; }
    for (int a = 0; a < D; a++) { DELTA[2 * a] = (i64)1 << (SH * a); DELTA[2 * a + 1] = -((i64)1 << (SH * a)); ORIGIN += (i64)64 << (SH * a); }
    printf("file %s magic %s d=%d K=%d N=%llu P=%llu Q=%llu (P/Q=%.9f) threads=%d mode=%s\n", argv[1], magic, D, G.K,
           (unsigned long long)G.N, (unsigned long long)G.P, (unsigned long long)G.Q, (double)G.P / G.Q, threads,
           G.explicitCh ? "explicit" : (G.prefixMode ? "search-prefix" : "search-full"));
    G.codes.resize(G.N); G.V.resize(G.N);
    if (fread(G.codes.data(), 8, G.N, f) != G.N || fread(G.V.data(), 8, G.N, f) != G.N) { printf("FAIL read\n"); return 1; }
    // (F) codes, weights; (S) start state
    i64 pos[64]; int L;
    bool fmtOK = true; long long startIdx = -1;
    for (u64 i = 0; i < G.N; i++) {
        if (i && G.codes[i] <= G.codes[i - 1]) { fmtOK = false; printf("codes not strictly increasing at %llu\n", (unsigned long long)i); break; }
        if (!decode(G.codes[i], pos, L)) { fmtOK = false; printf("undecodable code at %llu\n", (unsigned long long)i); break; }
        if (G.V[i] < 1 || G.V[i] > LIMIT60) { fmtOK = false; printf("weight out of range at %llu\n", (unsigned long long)i); break; }
        if (L == 0) startIdx = (long long)i;
        G.Lmax = max(G.Lmax, L); G.Kmax = max(G.Kmax, L + l1norm(pos[L]));
    }
    if (!fmtOK) { printf("FAIL format\n"); return 1; }
    if (startIdx < 0) { printf("FAIL start state (empty walk) not listed\n"); return 1; }
    START = (u64)startIdx;
    printf("format OK; start state index %lld; Lmax=%d, max(l+|p_l|_1)=%d; minV=%llu V[start]=%llu\n", startIdx, G.Lmax, G.Kmax,
           (unsigned long long)*min_element(G.V.begin(), G.V.end()), (unsigned long long)G.V[START]);
    // canonical-form self-consistency (informative; only used by the SAWCERT1 lookup)
    if (!G.explicitCh) {
        u64 noncanon = 0;
        for (u64 i = 0; i < G.N; i++) {
            u64 c = G.codes[i]; int l = lenOf(c); int img[8], sg[8]; for (int a = 0; a < 8; a++) img[a] = -1; int nl = 0;
            u64 nc = (1ULL << 63) | ((u64)l << 57);
            for (int k = 0; k < l; k++) { int t = stepOf(c, k), a = t >> 1, b = t & 1; if (img[a] < 0) { img[a] = nl++; sg[a] = b; }
                nc |= (u64)(2 * img[a] + (b ^ sg[a])) << (B * k); }
            if (nc != c) noncanon++;
        }
        printf("codes not in my canonical form: %llu\n", (unsigned long long)noncanon);
        G.table.build(G.codes);
    }
    vector<Stats> st(threads);
    NCOUNT = argc > 5 ? atoi(argv[5]) : 0;
    if (NCOUNT) CHOSEN.assign(G.N * 2 * D, 0xFFFFFFFFu);
    if (!G.explicitCh) {
        u64 maxStates = argc > 4 ? strtoull(argv[4], 0, 10) : G.N;
        u64 n = min(maxStates, G.N);
        atomic<u64> next(0); const u64 chunk = 4096;
        vector<thread> th;
        for (int w = 0; w < threads; w++) th.emplace_back([&, w] {
            for (;;) { u64 lo = next.fetch_add(chunk); if (lo >= n) break; processRange(lo, min(n, lo + chunk), nullptr, st[w]); } });
        for (auto& t : th) t.join();
        if (n < G.N) printf("NOTE: only %llu of %llu states processed\n", (unsigned long long)n, (unsigned long long)G.N);
    } else {
        // stream the child table
        const u64 CH = 1ULL << 22;  // states per chunk
        vector<u32> buf(CH * 2 * D);
        for (u64 lo = 0; lo < G.N; lo += CH) {
            u64 hi = min(G.N, lo + CH);
            u64 cnt = (hi - lo) * 2 * D;
            if (fread(buf.data(), 4, cnt, f) != cnt) { printf("FAIL read children\n"); return 1; }
            vector<thread> th; u64 per = (hi - lo + threads - 1) / threads;
            for (int w = 0; w < threads; w++) {
                u64 a = lo + w * per, b = min(hi, a + per);
                if (a >= b) continue;
                th.emplace_back([&, a, b, w] { processRange(a, b, buf.data() + (a - lo) * 2 * D, st[w]); });
            }
            for (auto& t : th) t.join();
        }
    }
    fclose(f);
    Stats tot;
    for (auto& s : st) {
        tot.pairs += s.pairs; tot.fails_child += s.fails_child; tot.fails_incl += s.fails_incl; tot.fails_ineq += s.fails_ineq;
        tot.nonSAW += s.nonSAW; tot.extraPtr += s.extraPtr; tot.candidates += s.candidates; tot.searchNodes += s.searchNodes;
        if ((u128)s.bestNum * tot.bestDen > (u128)tot.bestNum * s.bestDen) { tot.bestNum = s.bestNum; tot.bestDen = s.bestDen; tot.worstState = s.worstState; }
        for (auto& m : s.msgs) if (tot.msgs.size() < 20) tot.msgs.push_back(m);
    }
    for (auto& m : tot.msgs) printf("  %s\n", m.c_str());
    u64 g = gcdu(tot.bestNum, tot.bestDen);
    double secs = chrono::duration<double>(chrono::steady_clock::now() - T0).count();
    printf("pairs (s,e) checked: %llu; non-SAW listed walks: %llu; pointers on forbidden directions: %llu\n",
           (unsigned long long)tot.pairs, (unsigned long long)tot.nonSAW, (unsigned long long)tot.extraPtr);
    if (!G.explicitCh) printf("candidates found: %llu, search nodes: %.0f\n", (unsigned long long)tot.candidates, tot.searchNodes);
    printf("failures: child %llu, inclusion %llu, inequality %llu\n", (unsigned long long)tot.fails_child,
           (unsigned long long)tot.fails_incl, (unsigned long long)tot.fails_ineq);
    // print the exact max ratio with 12 decimals using integer division
    u64 num = tot.bestNum / g, den = tot.bestDen / g;
    u128 ip = num / den, rem = num % den; string dec;
    for (int i = 0; i < 12; i++) { rem *= 10; dec += char('0' + (int)(rem / den)); rem %= den; }
    printf("exact max ratio = %llu/%llu = %llu.%s... (state %llu)\n", (unsigned long long)num, (unsigned long long)den,
           (unsigned long long)ip, dec.c_str(), (unsigned long long)tot.worstState);
    bool ok = tot.fails_child == 0 && tot.fails_incl == 0 && tot.fails_ineq == 0;
    printf("time %.1f s\n", secs);
    if (NCOUNT) {  // a_n = (A^n 1)_{start}; x_{n+1}(s) = sum over allowed e of x_n(child(s,e))
        vector<u128> x(G.N, 1), y(G.N);
        printf("a_n (paths of the checked automaton from the start state):\n0 1\n");
        for (int n = 1; n <= NCOUNT; n++) {
            for (u64 s = 0; s < G.N; s++) { u128 acc = 0; for (int dir = 0; dir < 2 * D; dir++) { u32 t = CHOSEN[s * 2 * D + dir]; if (t != 0xFFFFFFFFu) acc += x[t]; } y[s] = acc; }
            x.swap(y);
            u128 a = x[START]; char buf[64]; int k = 63; buf[k] = 0; if (!a) buf[--k] = '0'; while (a) { buf[--k] = char('0' + (int)(a % 10)); a /= 10; }
            printf("%d %s\n", n, buf + k);
        }
    }
    printf("%s: mu_%d <= %llu/%llu\n", ok ? "PASS" : "FAIL", D, (unsigned long long)G.P, (unsigned long long)G.Q);
    return ok ? 0 : 1;
}
