/*
 * goodcount.c -- rigorous upper bound on the number of good vertices of every graph read (graph6, n <= 64).
 *
 * A vertex v of a connected graph G is good if G - v is connected and Kf(G - v) = Kf(G).
 * For a connected graph H on N vertices with grounded Laplacian A (one row/column deleted),
 * tau(H) = det A and Kf(H) = N tr(A^{-1}) - 1^T A^{-1} 1 = alpha(H)/tau(H) with alpha(H) an integer.
 * Modulo a prime p not dividing tau(H) we compute k_p(H) = alpha(H) tau(H)^{-1} mod p.
 * If k_p(G) != k_p(G - v) then alpha(G) tau(G-v) - alpha(G-v) tau(G) is not divisible by p, hence
 * nonzero, hence v is NOT good.  A cut vertex is never good.  A vertex is counted as "possibly good"
 * (pg) unless it is rejected for at least one of the two primes p1, p2 (if p | tau the test with that
 * prime is skipped, which can only increase the count).  Therefore
 *        number of good vertices of G  <=  pg(G).
 * The modular arithmetic (LDL^T factorisation and inverse over F_p) is taken from
 * ../kirchhoff-soltes-problem/scripts/kscheck.c.
 *
 * Output (stdout): one line per graph with pg(G) >= threshold:  graph6  n  pg
 * Summary (stderr): for each n, histogram of pg over connected graphs.
 * Usage: geng -c -q N | ./goodcount THRESHOLD_NUM THRESHOLD_DEN
 *        reports graphs with pg * DEN >= NUM * n   (e.g. 2 3 : at least two thirds possibly good).
 * Build: gcc -O3 -march=native -o goodcount goodcount.c
 * Environment: Fedora 44 (WSL), gcc; runtime about 1 min for all connected graphs with n <= 10.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 64
typedef uint64_t u64;
static const u64 P1 = 1073741789ULL; /* largest prime < 2^30 */
static const u64 P2 = 1073741783ULL; /* second largest prime < 2^30 */
static int n;
static u64 adj[MAXN];
static int parse_g6(const char *s) {
    int len = (int)strlen(s);
    while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r' || s[len-1] == ' ')) len--;
    if (len < 1) return -1;
    const unsigned char *p = (const unsigned char *)s;
    if (p[0] == '>') return -1;
    if (p[0] == '~') { if (len < 4 || p[1] == '~') return -1;
        n = ((p[1]-63) << 12) | ((p[2]-63) << 6) | (p[3]-63); p += 4; len -= 4; }
    else { n = p[0]-63; p += 1; len -= 1; }
    if (n < 1 || n > MAXN) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int k = 0, need = n*(n-1)/2;
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
static int connected(u64 S) {
    if (!S) return 1;
    u64 seen = S & (~S + 1), frontier = seen;
    while (frontier) { u64 nf = 0, f = frontier;
        while (f) { int u = __builtin_ctzll(f); f &= f-1; nf |= adj[u] & S; }
        nf &= ~seen; seen |= nf; frontier = nf; }
    return seen == S;
}
static u64 inv_mod(u64 a, u64 p) {
    int64_t t = 0, nt = 1, r = (int64_t)p, nr = (int64_t)(a % p);
    while (nr) { int64_t q = r / nr, tmp; tmp = t - q*nt; t = nt; nt = tmp; tmp = r - q*nr; r = nr; nr = tmp; }
    if (t < 0) t += (int64_t)p;
    return (u64)t;
}
static inline u64 dotmod(const u64 *x, const u64 *y, int len, u64 p) {
    u64 acc = 0; int k = 0;
    while (k + 15 <= len) { for (int t = 0; t < 15; t++) acc += x[k+t]*y[k+t]; acc %= p; k += 15; }
    for (; k < len; k++) acc += x[k]*y[k];
    return acc % p;
}
static u64 A[MAXN][MAXN], Lm[MAXN][MAXN], W[MAXN][MAXN], Y[MAXN][MAXN], D[MAXN], Dinv[MAXN];
static int kf_mod(u64 S, u64 p, u64 *out) {
    int idx[MAXN], N = 0; u64 s = S;
    while (s) { idx[N++] = __builtin_ctzll(s); s &= s-1; }
    int nv = N; N = nv - 1;
    if (N == 0) { *out = 0; return 0; }
    for (int i = 0; i < N; i++) { int u = idx[i];
        for (int j = 0; j < N; j++) { int w = idx[j];
            if (i == j) A[i][j] = (u64)popc(adj[u] & S);
            else A[i][j] = ((adj[u] >> w) & 1) ? p - 1 : 0; } }
    for (int j = 0; j < N; j++) {
        u64 dj = (A[j][j] + p - dotmod(Lm[j], W[j], j, p)) % p;
        if (dj == 0) return -1;
        D[j] = dj; Dinv[j] = inv_mod(dj, p);
        for (int i = j+1; i < N; i++) {
            u64 wij = (A[i][j] + p - dotmod(Lm[i], W[j], j, p)) % p;
            W[i][j] = wij; Lm[i][j] = (wij * Dinv[j]) % p; }
    }
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < i; k++) { u64 acc = 0; int cnt = 0;
            for (int m = k; m < i; m++) { acc += Lm[i][m] * Y[m][k]; if (++cnt == 15) { acc %= p; cnt = 0; } }
            acc %= p; Y[i][k] = acc ? p - acc : 0; }
        Y[i][i] = 1;
    }
    u64 T = 0, Ssum = 0;
    for (int k = 0; k < N; k++) { u64 sq = 0, rs = 0;
        for (int i = 0; i <= k; i++) { sq = (sq + Y[k][i]*Y[k][i]) % p; rs += Y[k][i]; }
        rs %= p; T = (T + sq * Dinv[k]) % p; Ssum = (Ssum + ((rs*rs) % p) * Dinv[k]) % p; }
    *out = (((u64)nv % p) * T % p + p - Ssum) % p;
    return 0;
}
int main(int argc, char **argv) {
    long num = argc > 2 ? atol(argv[1]) : 2, den = argc > 2 ? atol(argv[2]) : 3;
    static char line[1 << 16];
    static long long hist[MAXN+1][MAXN+1]; long long ntot = 0, ndisc = 0;
    int nmax = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>' || line[0] == '\n') continue;
        if (parse_g6(line) < 0) continue;
        u64 all = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
        if (!connected(all)) { ndisc++; continue; }
        ntot++; if (n > nmax) nmax = n;
        u64 kG[2]; int okG[2];
        u64 pr[2] = {P1, P2};
        for (int t = 0; t < 2; t++) okG[t] = (kf_mod(all, pr[t], &kG[t]) == 0);
        int pg = 0;
        for (int v = 0; v < n; v++) {
            u64 S = all & ~(1ULL << v);
            if (n >= 2 && !connected(S)) continue;      /* cut vertex: not good */
            if (n == 1) continue;
            int rejected = 0;
            for (int t = 0; t < 2 && !rejected; t++) {
                u64 kH;
                if (!okG[t] || kf_mod(S, pr[t], &kH) < 0) continue;   /* test skipped: conservative */
                if (kH != kG[t]) rejected = 1;
            }
            if (!rejected) pg++;
        }
        hist[n][pg]++;
        if ((long)pg * den >= num * (long)n && pg > 0) {
            int len = (int)strlen(line); while (len && (line[len-1]=='\n'||line[len-1]=='\r')) line[--len] = 0;
            printf("%s %d %d\n", line, n, pg);
        }
    }
    fprintf(stderr, "connected graphs: %lld (skipped disconnected: %lld)\n", ntot, ndisc);
    for (int m = 1; m <= nmax; m++) {
        long long tot = 0; for (int g = 0; g <= m; g++) tot += hist[m][g];
        if (!tot) continue;
        fprintf(stderr, "n=%d graphs=%lld pg-histogram:", m, tot);
        for (int g = 0; g <= m; g++) if (hist[m][g]) fprintf(stderr, " %d:%lld", g, hist[m][g]);
        fprintf(stderr, "\n");
    }
    return 0;
}
