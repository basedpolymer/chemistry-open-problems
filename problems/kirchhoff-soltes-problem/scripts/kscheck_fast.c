/*
 * kscheck_fast.c -- same exact modular criterion as kscheck.c (see its header and result.md section 11),
 * re-implemented for speed: primes p < 2^28 (p1 = 268435399, p2 = 268435367), residues stored as uint32,
 * dot products accumulated in uint64 without intermediate reduction (at most 63 products < 2^56 each,
 * so the sum is < 2^62), contiguous memory access (the inverse of the unit triangular factor is stored
 * transposed).  Loops are written so that gcc -O3 -march=native vectorises them.
 *
 * Criterion: for connected H on N vertices with grounded Laplacian A (last vertex of the vertex set
 * removed), k_p(H) = N tr(A^{-1}) - 1^T A^{-1} 1 mod p = alpha(H) tau(H)^{-1} mod p.  k_p(G) != k_p(G-v)
 * proves Kf(G) != Kf(G-v).  Cut vertices are never good.  A graph passing every tested vertex for p1 is
 * re-tested on ALL vertices for p1 and p2 and printed as SURVIVOR (or INCONCLUSIVE if some p | tau for
 * both primes).  Options: -t (vertex-transitive input: test vertex 0 only before the full re-test),
 * -k (debug: print k_p1(G) and k_p1(G - v) for all v).
 * Build: gcc -O3 -march=native -o kscheck_fast kscheck_fast.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 64
typedef uint64_t u64;
typedef uint32_t u32;

static const u64 P1 = 268435399ULL;
static const u64 P2 = 268435367ULL;

static int n;
static u64 adj[MAXN];

static int parse_g6(const char *s) {
    int len = (int)strlen(s);
    while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r' || s[len-1] == ' ')) len--;
    if (len < 1) return -1;
    const unsigned char *p = (const unsigned char *)s;
    if (p[0] == '>') return -1;
    if (p[0] == '~') {
        if (len < 4 || p[1] == '~') return -1;
        n = ((p[1]-63) << 12) | ((p[2]-63) << 6) | (p[3]-63);
        p += 4; len -= 4;
    } else { n = p[0]-63; p += 1; len -= 1; }
    if (n < 1 || n > MAXN) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int need = n*(n-1)/2;
    if (len*6 < need) return -1;
    int k = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int byte = k / 6, bit = 5 - (k % 6);
            if (((p[byte]-63) >> bit) & 1) { adj[i] |= 1ULL << j; adj[j] |= 1ULL << i; }
            k++;
        }
    return n;
}

static int popc(u64 x) { return __builtin_popcountll(x); }

static int connected(u64 S) {
    if (!S) return 1;
    u64 seen = S & (~S + 1), frontier = seen;
    while (frontier) {
        u64 nf = 0, f = frontier;
        while (f) { int u = __builtin_ctzll(f); f &= f-1; nf |= adj[u] & S; }
        nf &= ~seen; seen |= nf; frontier = nf;
    }
    return seen == S;
}

static u64 inv_mod(u64 a, u64 p) {
    int64_t t = 0, nt = 1, r = (int64_t)p, nr = (int64_t)(a % p);
    while (nr) { int64_t q = r / nr, tmp;
        tmp = t - q*nt; t = nt; nt = tmp; tmp = r - q*nr; r = nr; nr = tmp; }
    if (t < 0) t += (int64_t)p;
    return (u64)t;
}

/* Fermat inverse a^(p-2) mod p (p a compile-time constant after inlining) */
static inline __attribute__((always_inline)) u64 inv_fermat(u64 a, const u64 p) {
    u64 r = 1, b = a % p, e = p - 2;
    while (e) { if (e & 1) r = (r * b) % p; b = (b * b) % p; e >>= 1; }
    return r;
}

static inline u64 dot32(const u32 *restrict x, const u32 *restrict y, int len) {
    u64 acc = 0;
    for (int k = 0; k < len; k++) acc += (u64)x[k] * (u64)y[k];
    return acc;            /* < 63 * 2^56 < 2^62 */
}

static u32 Am[MAXN][MAXN], Lm[MAXN][MAXN], W[MAXN][MAXN], YT[MAXN][MAXN];
static u64 Dinv[MAXN];

/* Kf of induced subgraph G[S] modulo p; 0 on success, -1 if p | tau.  Always inlined into the two
   wrappers below so that p is a compile-time constant (reduction by multiplication, not division). */
static inline __attribute__((always_inline)) int kf_mod_impl(u64 S, const u64 p, u64 *out) {
    int idx[MAXN], nv = 0;
    u64 s = S;
    while (s) { idx[nv++] = __builtin_ctzll(s); s &= s-1; }
    int N = nv - 1;
    if (N == 0) { *out = 0; return 0; }
    for (int i = 0; i < N; i++) {
        int u = idx[i];
        for (int j = 0; j < N; j++) {
            int w = idx[j];
            if (i == j) Am[i][j] = (u32)popc(adj[u] & S);
            else Am[i][j] = ((adj[u] >> w) & 1) ? (u32)(p - 1) : 0;
        }
    }
    /* LDL^T, left-looking: W[i][k] = L[i][k]*D[k] (k < i) */
    for (int j = 0; j < N; j++) {
        u64 dj = ((u64)Am[j][j] + p - dot32(Lm[j], W[j], j) % p) % p;
        if (dj == 0) return -1;
        u64 dji = inv_fermat(dj, p);
        Dinv[j] = dji;
        for (int i = j+1; i < N; i++) {
            u64 wij = ((u64)Am[i][j] + p - dot32(Lm[i], W[j], j) % p) % p;
            W[i][j] = (u32)wij;
            Lm[i][j] = (u32)((wij * dji) % p);
        }
    }
    /* Y = L^{-1} (unit lower); YT[k][i] = Y[i][k].  Y[i][k] = -sum_{m=k}^{i-1} L[i][m] Y[m][k]
       = -sum_{m=k}^{i-1} L[i][m] YT[k][m]  (contiguous in m) */
    for (int k = 0; k < N; k++) {
        YT[k][k] = 1;
        for (int i = k+1; i < N; i++) {
            u64 acc = dot32(&Lm[i][k], &YT[k][k], i - k) % p;
            YT[k][i] = (u32)(acc ? p - acc : 0);
        }
    }
    /* tr(A^{-1}) = sum_k Dinv[k] sum_{i<=k} Y[k][i]^2 = sum_k Dinv[k] sum_{i<=k} YT[i][k]^2 ;
       1^T A^{-1} 1 = sum_k Dinv[k] (sum_{i<=k} YT[i][k])^2 */
    static u64 sq[MAXN], rs[MAXN];
    for (int k = 0; k < N; k++) { sq[k] = 0; rs[k] = 0; }
    for (int i = 0; i < N; i++)
        for (int k = i; k < N; k++) { u64 y = YT[i][k]; sq[k] += y * y; rs[k] += y; }
    u64 T = 0, Ssum = 0;
    for (int k = 0; k < N; k++) {
        u64 a = sq[k] % p, b = rs[k] % p;
        T = (T + a * Dinv[k]) % p;
        Ssum = (Ssum + ((b * b) % p) * Dinv[k]) % p;
    }
    *out = (((u64)nv % p) * T % p + p - Ssum) % p;
    return 0;
}

static int kf_mod_p1(u64 S, u64 *out) { return kf_mod_impl(S, 268435399ULL, out); }
static int kf_mod_p2(u64 S, u64 *out) { return kf_mod_impl(S, 268435367ULL, out); }
static int kf_mod(u64 S, u64 p, u64 *out) { return p == P1 ? kf_mod_p1(S, out) : kf_mod_p2(S, out); }

static int test_vertex(int v, u64 kG, u64 p) {
    u64 all = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    u64 S = all & ~(1ULL << v);
    if (!connected(S)) return 0;
    u64 kH;
    if (kf_mod(S, p, &kH) < 0) return 2;
    return kH == kG ? 1 : 0;
}

int main(int argc, char **argv) {
    int transitive = 0, debug = 0;
    for (int a = 1; a < argc; a++) { if (!strcmp(argv[a], "-t")) transitive = 1;
                                     if (!strcmp(argv[a], "-k")) debug = 1; }
    static char line[1 << 16];
    long long cnt_total = 0, cnt_bad = 0, cnt_disc = 0, cnt_rej = 0, cnt_cut = 0,
              cnt_surv = 0, cnt_inc = 0, cnt_small = 0;
    long long byn[MAXN+1]; memset(byn, 0, sizeof byn);
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n') continue;
        if (parse_g6(line) < 0) { cnt_bad++; continue; }
        cnt_total++; byn[n]++;
        if (n < 3) { cnt_small++; continue; }
        u64 all = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
        if (!connected(all)) { cnt_disc++; continue; }
        if (debug) {
            u64 kq; int r = kf_mod(all, P1, &kq);
            printf("%d %lld", n, r < 0 ? -1LL : (long long)kq);
            for (int v = 0; v < n; v++) { u64 S = all & ~(1ULL << v), kh;
                if (!connected(S) || kf_mod(S, P1, &kh) < 0) printf(" -1"); else printf(" %llu", (unsigned long long)kh); }
            printf("\n"); continue;
        }
        u64 primes[2] = {P1, P2};
        int verdict;
        int pi = 0; u64 kG = 0;
        while (pi < 2 && kf_mod(all, primes[pi], &kG) < 0) pi++;
        if (pi == 2) verdict = 2;
        else {
            int vmax = transitive ? 1 : n;
            int rejected = 0, inconclusive = 0;
            for (int v = 0; v < vmax && !rejected; v++) {
                int r = test_vertex(v, kG, primes[pi]);
                if (r == 2) {
                    int q = 1 - pi; u64 kG2;
                    if (kf_mod(all, primes[q], &kG2) < 0) { inconclusive = 1; continue; }
                    r = test_vertex(v, kG2, primes[q]);
                    if (r == 2) { inconclusive = 1; continue; }
                }
                if (r == 0) {
                    rejected = 1;
                    if (!connected(all & ~(1ULL << v))) cnt_cut++;
                }
            }
            if (rejected) verdict = 0;
            else {
                int ok = 1, inc = inconclusive;
                for (int q = 0; q < 2 && ok; q++) {
                    u64 kq;
                    if (kf_mod(all, primes[q], &kq) < 0) { inc = 1; continue; }
                    for (int v = 0; v < n && ok; v++) {
                        int r = test_vertex(v, kq, primes[q]);
                        if (r == 0) ok = 0; else if (r == 2) inc = 1;
                    }
                }
                verdict = ok ? (inc ? 2 : 1) : 0;
            }
        }
        if (verdict == 0) cnt_rej++;
        else {
            int len = (int)strlen(line);
            while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
            if (verdict == 1) { cnt_surv++; printf("%s SURVIVOR\n", line); }
            else { cnt_inc++; printf("%s INCONCLUSIVE\n", line); }
            fflush(stdout);
        }
    }
    fprintf(stderr, "graphs_read %lld unparsable %lld n<3 %lld disconnected %lld rejected %lld "
            "(rejected_by_cut_vertex %lld) survivors %lld inconclusive %lld\n",
            cnt_total, cnt_bad, cnt_small, cnt_disc, cnt_rej, cnt_cut, cnt_surv, cnt_inc);
    for (int i = 0; i <= MAXN; i++) if (byn[i]) fprintf(stderr, "n=%d count %lld\n", i, byn[i]);
    return 0;
}
