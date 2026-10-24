/*
 * kscheck.c -- exact modular filter for the Kirchhoff Soltes problem.
 *
 * Reads graphs in graph6 format (one per line, n <= 64) from stdin and decides, for each graph G,
 * whether Kf(G - v) = Kf(G) can hold for every vertex v.
 *
 * Exact criterion used for REJECTION (rigorous):
 *   For a connected graph H on N vertices, grounded Laplacian A = L(H) with the row/column of one
 *   vertex removed, tau(H) = det A (matrix-tree theorem) and
 *        Kf(H) = N tr(A^{-1}) - 1^T A^{-1} 1        (Kf(H) = alpha(H)/tau(H), alpha(H) integer).
 *   Over F_p (p prime, p not dividing tau(H)) we compute k_p(H) = N tr(A^{-1}) - 1^T A^{-1} 1 mod p,
 *   i.e. k_p(H) = alpha(H) * tau(H)^{-1} mod p.  If k_p(G) != k_p(G - v) then the integer
 *   alpha(G) tau(G-v) - alpha(G-v) tau(G) is not divisible by p, hence nonzero, hence
 *   Kf(G) != Kf(G-v): v is not good.  A cut vertex v (G - v disconnected) is never good.
 *   If p divides tau (zero pivot mod p) the second prime is used; if both fail the graph is
 *   reported as INCONCLUSIVE (never silently rejected).
 *
 * A graph is reported as a SURVIVOR when every vertex passes the test modulo p1 AND modulo p2;
 * survivors must then be checked with exact rational arithmetic (scripts/exact_kf.py).
 *
 * Modes:  -t   ("transitive") test only vertex 0 (valid for vertex-transitive inputs; a survivor
 *              is then re-tested on all vertices anyway).
 *         default: vertices are tested in order 0,1,...,n-1 until one fails.
 * Output: survivors / inconclusive graphs on stdout (graph6 + tag); summary on stderr.
 *
 * Build:  gcc -O3 -march=native -o kscheck kscheck.c
 * Environment of the runs reported in result.md: gcc (Linux x86_64), 1 thread, nice -n 10.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 64
typedef uint64_t u64;
typedef uint32_t u32;

static const u64 P1 = 1073741789ULL; /* largest prime < 2^30 */
static const u64 P2 = 1073741783ULL; /* second largest prime < 2^30 */

static int n;
static u64 adj[MAXN];

static int parse_g6(const char *s) {
    int len = (int)strlen(s);
    while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r' || s[len-1] == ' ')) len--;
    if (len < 1) return -1;
    const unsigned char *p = (const unsigned char *)s;
    if (p[0] == '>') return -1; /* header */
    if (p[0] == '~') {
        if (len < 4 || p[1] == '~') return -1;
        n = ((p[1]-63) << 12) | ((p[2]-63) << 6) | (p[3]-63);
        p += 4; len -= 4;
    } else { n = p[0]-63; p += 1; len -= 1; }
    if (n < 1 || n > MAXN) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int k = 0; /* bit index */
    int need = n*(n-1)/2;
    if (len*6 < need) return -1;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int byte = k / 6, bit = 5 - (k % 6);
            if (((p[byte]-63) >> bit) & 1) { adj[i] |= 1ULL << j; adj[j] |= 1ULL << i; }
            k++;
        }
    return n;
}

static int popc(u64 x) { return __builtin_popcountll(x); }

/* is induced subgraph on vertex set S connected? */
static int connected(u64 S) {
    if (!S) return 1;
    u64 seen = S & (~S + 1), frontier = seen;
    while (frontier) {
        u64 nf = 0;
        u64 f = frontier;
        while (f) { int u = __builtin_ctzll(f); f &= f-1; nf |= adj[u] & S; }
        nf &= ~seen; seen |= nf; frontier = nf;
    }
    return seen == S;
}

static u64 inv_mod(u64 a, u64 p) { /* a != 0 mod p */
    int64_t t = 0, nt = 1, r = (int64_t)p, nr = (int64_t)(a % p);
    while (nr) { int64_t q = r / nr, tmp;
        tmp = t - q*nt; t = nt; nt = tmp; tmp = r - q*nr; r = nr; nr = tmp; }
    if (t < 0) t += (int64_t)p;
    return (u64)t;
}

/* dot product sum_{k<len} x[k]*y[k] mod p with lazy reduction (p < 2^30: 15 products < 2^64) */
static inline u64 dotmod(const u64 *x, const u64 *y, int len, u64 p) {
    u64 acc = 0; int k = 0;
    while (k + 15 <= len) {
        for (int t = 0; t < 15; t++) acc += x[k+t]*y[k+t];
        acc %= p; k += 15;
    }
    for (; k < len; k++) acc += x[k]*y[k];
    return acc % p;
}

static u64 A[MAXN][MAXN], Lm[MAXN][MAXN], W[MAXN][MAXN], Y[MAXN][MAXN], D[MAXN], Dinv[MAXN];

/* Kf of induced subgraph G[S] modulo p; returns 0 and sets *out on success, -1 if p | tau. */
static int kf_mod(u64 S, u64 p, u64 *out) {
    int idx[MAXN], N = 0;
    u64 s = S;
    while (s) { idx[N++] = __builtin_ctzll(s); s &= s-1; }
    int nv = N; N = nv - 1;              /* ground at the last vertex of S */
    if (N == 0) { *out = 0; return 0; }
    for (int i = 0; i < N; i++) {
        int u = idx[i];
        for (int j = 0; j < N; j++) {
            int w = idx[j];
            if (i == j) A[i][j] = (u64)popc(adj[u] & S);
            else A[i][j] = ((adj[u] >> w) & 1) ? p - 1 : 0;
        }
    }
    /* LDL^T: W[i][j] = L[i][j]*D[j] (i>j), stored row-wise for contiguous dot products */
    for (int j = 0; j < N; j++) {
        u64 dj = (A[j][j] + p - dotmod(Lm[j], W[j], j, p)) % p;
        if (dj == 0) return -1;
        D[j] = dj; Dinv[j] = inv_mod(dj, p);
        for (int i = j+1; i < N; i++) {
            u64 wij = (A[i][j] + p - dotmod(Lm[i], W[j], j, p)) % p;
            W[i][j] = wij;              /* = L[i][j] * D[j] */
            Lm[i][j] = (wij * Dinv[j]) % p;
        }
        /* W[j][k] for k<j must hold L[j][k]*D[k]: set when row j was processed (as W[j][k]) */
    }
    /* Y = L^{-1}, unit lower triangular. Y[i][k] = - sum_{m=k}^{i-1} L[i][m] Y[m][k] */
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < i; k++) {
            u64 acc = 0; int cnt = 0;
            for (int m = k; m < i; m++) {
                acc += Lm[i][m] * Y[m][k];
                if (++cnt == 15) { acc %= p; cnt = 0; }
            }
            acc %= p;
            Y[i][k] = acc ? p - acc : 0;
        }
        Y[i][i] = 1;
    }
    u64 T = 0, Ssum = 0;
    for (int k = 0; k < N; k++) {
        u64 sq = 0, rs = 0;
        for (int i = 0; i <= k; i++) { sq = (sq + Y[k][i]*Y[k][i]) % p; rs += Y[k][i]; }
        rs %= p;
        T = (T + sq * Dinv[k]) % p;
        Ssum = (Ssum + ((rs*rs) % p) * Dinv[k]) % p;
    }
    *out = (((u64)nv % p) * T % p + p - Ssum) % p;
    return 0;
}

/* returns 1 if vertex v passes mod p (or 2 if undecidable for this p), 0 if certainly not good */
static int test_vertex(int v, u64 kG, u64 p) {
    u64 all = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    u64 S = all & ~(1ULL << v);
    if (!connected(S)) return 0;          /* cut vertex: never good */
    u64 kH;
    if (kf_mod(S, p, &kH) < 0) return 2;
    return kH == kG ? 1 : 0;
}

int main(int argc, char **argv) {
    int transitive = 0;
    for (int a = 1; a < argc; a++) if (!strcmp(argv[a], "-t")) transitive = 1;
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
        u64 primes[2] = {P1, P2};
        int verdict = -1; /* 0 rejected, 1 survivor, 2 inconclusive */
        int pi = 0; u64 kG = 0;
        /* choose a prime not dividing tau(G) */
        while (pi < 2 && kf_mod(all, primes[pi], &kG) < 0) pi++;
        if (pi == 2) { verdict = 2; }
        else {
            int vmax = transitive ? 1 : n;
            int rejected = 0, inconclusive = 0;
            for (int v = 0; v < vmax && !rejected; v++) {
                int r = test_vertex(v, kG, primes[pi]);
                if (r == 2) { /* try the other prime for this vertex */
                    int q = 1 - pi; u64 kG2;
                    if (kf_mod(all, primes[q], &kG2) < 0) { inconclusive = 1; continue; }
                    r = test_vertex(v, kG2, primes[q]);
                    if (r == 2) { inconclusive = 1; continue; }
                }
                if (r == 0) {
                    rejected = 1;
                    u64 S = all & ~(1ULL << v);
                    if (!connected(S)) cnt_cut++;
                }
            }
            if (rejected) verdict = 0;
            else {
                /* full re-test on all vertices with both primes */
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
