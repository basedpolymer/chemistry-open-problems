/*
 * gadgetbound.c -- exact computation of the splice efficiency of every rooted graph read (graph6, n <= 12).
 *
 * For a connected graph B on n vertices and a root x, let m = n - 1, R = Rt_B(x) and, for every vertex v != x with
 * B - v connected, rho(v) = Rt_{B-v}(x) - Rt_B(x) and delta(v) = Kf(B) - Kf(B - v).  A class is a maximal set of
 * such vertices with equal (rho, delta).  The efficiency of a class O with rho_O > 0 is
 *      eff(O) = (|O|/m) * min(1, rho_O * m / R),         eff(O) = 0 if rho_O <= 0,
 * and eff(B, x) is the maximum over the classes (Proposition 2.4 of the manuscript).
 *
 * Arithmetic is exact: for a connected graph H with grounded Laplacian A (last vertex deleted) we compute
 * tau = det A and adj(A) by fraction-free Gauss-Jordan elimination (Bareiss) in __int128; the pivots are leading
 * principal minors of a positive definite matrix, hence positive, and every intermediate value is a minor of
 * [A | I], bounded by Hadamard's inequality (< 10^16 for n <= 12).  Then tau*r(u,w) = adj_uu + adj_ww - 2 adj_uw,
 * tau*r(u,g) = adj_uu, so tau*Rt(x) and tau*Kf(H) are integers; all comparisons are done by cross-multiplication.
 *
 * Output (stdout): every (graph6, root, class) with eff > NUM/DEN (default 2/3), and the summary (stderr):
 * number of graphs and roots, the largest eff found among classes with 3|O| > 2m (exact fraction printed as
 * |O|/m * rho m / R), and the number of (root, class) pairs with eff exactly equal to NUM/DEN.
 * Usage: geng -c -q N | ./gadgetbound [NUM DEN]     (./gadgetbound -d prints every class, used by the tests)
 * Build: gcc -O3 -march=native -o gadgetbound gadgetbound.c
 * Environment: Fedora 44 (WSL), gcc.  Runtime: n = 9 about 1 min, n = 10 about 30 min (one core).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef __int128 i128;
typedef int64_t i64;
#define MAXN 12
static int n;
static uint32_t adj[MAXN];

static int parse_g6(const char *s) {
    int len = (int)strlen(s);
    while (len > 0 && (s[len-1] == '\n' || s[len-1] == '\r' || s[len-1] == ' ')) len--;
    if (len < 1 || s[0] == '>') return -1;
    const unsigned char *p = (const unsigned char *)s;
    n = p[0] - 63; p++; len--;
    if (n < 1 || n > MAXN) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int k = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int byte = k / 6, bit = 5 - (k % 6);
            if (((p[byte]-63) >> bit) & 1) { adj[i] |= 1u << j; adj[j] |= 1u << i; }
            k++;
        }
    return n;
}
static int connected(uint32_t S) {
    if (!S) return 1;
    uint32_t seen = S & (~S + 1), fr = seen;
    while (fr) { uint32_t nf = 0, f = fr;
        while (f) { int u = __builtin_ctz(f); f &= f-1; nf |= adj[u] & S; }
        nf &= ~seen; seen |= nf; fr = nf; }
    return seen == S;
}
/* For H = B[S] connected: tau, and tau*Rt(u) for every u in S (indexed by vertex), tau*Kf. */
static void kf_exact(uint32_t S, i128 *tau, i128 *tRt, i128 *tKf) {
    int idx[MAXN], N = 0; uint32_t s = S;
    while (s) { idx[N++] = __builtin_ctz(s); s &= s-1; }
    int nv = N; N = nv - 1;
    for (int u = 0; u < MAXN; u++) tRt[u] = 0;
    if (N == 0) { *tau = 1; *tKf = 0; return; }
    i128 M[MAXN][2*MAXN];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 2*N; j++) M[i][j] = 0;
        M[i][i] = __builtin_popcount(adj[idx[i]] & S);
        for (int j = 0; j < N; j++) if (j != i && ((adj[idx[i]] >> idx[j]) & 1)) M[i][j] = -1;
        M[i][N+i] = 1;
    }
    i128 prev = 1;
    for (int k = 0; k < N; k++) {
        i128 piv = M[k][k];
        if (piv <= 0) { fprintf(stderr, "nonpositive pivot\n"); exit(2); }
        for (int i = 0; i < N; i++) {
            if (i == k) continue;
            i128 f = M[i][k];
            for (int j = 0; j < 2*N; j++) {
                i128 num = piv * M[i][j] - f * M[k][j];
                if (num % prev != 0) { fprintf(stderr, "inexact division\n"); exit(2); }
                M[i][j] = num / prev;
            }
        }
        prev = piv;
    }
    /* now left block = det * I (rows other than the last carry det after full elimination) */
    i128 det = M[N-1][N-1];
    /* every row i: M[i][i] = det, right block = adj */
    for (int i = 0; i < N; i++) if (M[i][i] != det) { fprintf(stderr, "Gauss-Jordan check failed\n"); exit(2); }
    *tau = det;
    /* tau * r(u,w) */
    i128 kfs = 0;
    for (int i = 0; i < N; i++) {
        i128 rig = M[i][N+i];                       /* tau r(idx[i], g) */
        tRt[idx[i]] += rig; tRt[idx[nv-1]] += rig; kfs += rig;
        for (int j = i+1; j < N; j++) {
            i128 r = M[i][N+i] + M[j][N+j] - 2*M[i][N+j];
            tRt[idx[i]] += r; tRt[idx[j]] += r; kfs += r;
        }
    }
    *tKf = kfs;
}
static void print128(i128 x) {
    char buf[64]; int k = 63; buf[k] = 0; int neg = x < 0; if (neg) x = -x;
    do { buf[--k] = '0' + (int)(x % 10); x /= 10; } while (x);
    if (neg) buf[--k] = '-';
    fputs(buf + k, stdout);
}
int main(int argc, char **argv) {
    int debug = (argc > 1 && !strcmp(argv[1], "-d"));      /* -d: print every class (graph6 x g eff_num eff_den) */
    i64 NUM = (!debug && argc > 2) ? atoll(argv[1]) : 2, DEN = (!debug && argc > 2) ? atoll(argv[2]) : 3;
    static char line[4096];
    long long ngraphs = 0, nroots = 0, nequal = 0, nabove = 0;
    /* best eff among classes with 3g > 2m, stored as fraction bn/bd */
    long double best = -1; char bestg6[64] = ""; int bestx = -1, bestg = 0, bestm = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (parse_g6(line) < 0) continue;
        uint32_t all = (1u << n) - 1;
        if (!connected(all)) continue;
        ngraphs++;
        i128 tau0, tRt0[MAXN], tKf0;
        kf_exact(all, &tau0, tRt0, &tKf0);
        i128 tauv[MAXN], tRtv[MAXN][MAXN], tKfv[MAXN]; int ok[MAXN];
        for (int v = 0; v < n; v++) {
            uint32_t S = all & ~(1u << v);
            ok[v] = (n >= 2) && connected(S);
            if (ok[v]) kf_exact(S, &tauv[v], tRtv[v], &tKfv[v]);
        }
        int m = n - 1;
        for (int x = 0; x < n; x++) {
            nroots++;
            /* rho(v) = (tRtv[v][x]*tau0 - tRt0[x]*tauv[v]) / (tau0*tauv[v]);  delta(v) = (tKf0*tauv - tKfv*tau0)/(tau0 tauv) */
            int cls[MAXN]; for (int v = 0; v < n; v++) cls[v] = -1;
            for (int v = 0; v < n; v++) {
                if (v == x || !ok[v] || cls[v] >= 0) continue;
                cls[v] = v;
                /* rho(v) = rnv/(tau0 tauv), delta(v) = dnv/(tau0 tauv): equal iff cross products agree (tau0 cancels) */
                i128 rnv = tRtv[v][x]*tau0 - tRt0[x]*tauv[v], dnv = tKf0*tauv[v] - tKfv[v]*tau0;
                for (int w = v+1; w < n; w++) {
                    if (w == x || !ok[w] || cls[w] >= 0) continue;
                    i128 rnw = tRtv[w][x]*tau0 - tRt0[x]*tauv[w], dnw = tKf0*tauv[w] - tKfv[w]*tau0;
                    if (rnv * tauv[w] == rnw * tauv[v] && dnv * tauv[w] == dnw * tauv[v]) cls[w] = v;
                }
            }
            for (int v = 0; v < n; v++) {
                if (cls[v] != v) continue;
                int g = 0; for (int w = 0; w < n; w++) if (cls[w] == v) g++;
                i128 rn = tRtv[v][x]*tau0 - tRt0[x]*tauv[v];      /* rho = rn / (tau0 tauv) */
                /* rho*m/R = rn*m / (tauv * tRt0[x]);  eff = (g/m) * min(1, that);  eff = 0 if rho <= 0 */
                i128 qn = rn * m, qd = tauv[v] * tRt0[x];
                i128 en, ed;                                          /* eff = en/ed */
                if (rn <= 0) { en = 0; ed = 1; }
                else if (qn >= qd) { en = g; ed = m; } else { en = (i128)g * qn; ed = (i128)m * qd; }
                if (debug) {   /* D graph6 x g eff_num eff_den rho_num a_num common_den  (rho = rho_num/den, a = a_num/den) */
                    line[strcspn(line, "\n")] = 0;
                    i128 den = tau0 * tauv[v], an = tRt0[x] * tauv[v] - (i128)m * rn;
                    printf("D %s %d %d ", line, x, g); print128(en); printf(" "); print128(ed);
                    printf(" "); print128(rn); printf(" "); print128(an); printf(" "); print128(den); printf("\n");
                }
                if (rn <= 0) continue;
                /* compare with NUM/DEN */
                i128 lhs = en * DEN, rhs = ed * NUM;
                if (lhs > rhs) {
                    nabove++;
                    line[strcspn(line, "\n")] = 0;
                    printf("ABOVE %s x=%d g=%d m=%d eff=", line, x, g, m);
                    print128(en); printf("/"); print128(ed); printf("\n");
                } else if (lhs == rhs) nequal++;
                /* informational maximum (long double); the rigorous test is the exact comparison above */
                if (3*g > 2*m && (long double)en / (long double)ed > best) {
                    best = (long double)en / (long double)ed; bestx = x; bestg = g; bestm = m;
                    snprintf(bestg6, sizeof bestg6, "%.60s", line); bestg6[strcspn(bestg6, "\n")] = 0;
                }
            }
        }
    }
    fprintf(stderr, "graphs=%lld roots=%lld classes_with_eff_above_%lld/%lld=%lld classes_with_eff_equal=%lld\n",
            ngraphs, nroots, (long long)NUM, (long long)DEN, nabove, nequal);
    fflush(stdout);
    if (bestx >= 0) {
        double e = (double)best;
        fprintf(stderr, "max eff among classes with 3g > 2m: %.6f at %s x=%d g=%d m=%d\n", e, bestg6, bestx, bestg, bestm);
    } else fprintf(stderr, "no class with 3g > 2m and rho > 0\n");
    return 0;
}
