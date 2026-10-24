/*
 * dualspec.c -- exhaustive spectral-collision scan of fullerene duals (pipeline: buckygen -q Nd | dualspec ...).
 *
 * Reads planar_code (triangulations = duals T_n of the fullerenes C_n, as written by buckygen WITHOUT -d)
 * from stdin.  For every graph (index i = 0,1,2,... in buckygen's deterministic output order) and for every
 * requested matrix type it computes the closed-walk numbers  t_k = tr(M^k) mod 2^64, k = 1..K, of an integer
 * symmetric matrix M built from the graph, hashes the vector (t_1..t_K) to 128 bits, and inserts it into a
 * hash table.  Two graphs with the same spectrum of M have the same integer traces for every k, hence the
 * same residues mod 2^64 and the same key: the scan has NO false negatives.  Every key collision is printed
 * (type, first index, new index); the collisions are then settled EXACTLY (integer characteristic polynomial,
 * python-flint) by verify/settle scripts.  Hash-equal but non-cospectral pairs (possible because K may be < the
 * number of vertices, or by a 128-bit hash accident) are eliminated there.
 *
 * Matrix types (bit mask argument):
 *   1  : A(T)            adjacency matrix of the whole dual T_n           (m = n/2+2 vertices)
 *   2  : A(T6)           adjacency matrix of the subgraph induced by degree-6 vertices (hexagons)
 *   4  : A(T5)           adjacency matrix of the subgraph induced by degree-5 vertices (pentagons)
 *   8  : 2A + D  on T    (spectrum determines the (1/2,1/4)-character tr exp(A/2 + D/4))
 *  16  : A + D   on T    (signless Laplacian; (1/2,1/2)-character)
 *  32  : A + tP  on T    P = diag(1 on pentagons), t = fixed odd 64-bit constant: detects pairs with equal
 *                        bivariate polynomial det(xI - A - sP) (=> A + cD cospectral for every real c)
 *
 * Usage: buckygen -q 60d | dualspec MASK [K [print5]]
 *   K = number of traces (default: m, the number of vertices, i.e. the full spectrum information).
 * Output (stdout): lines "C type i j" for collisions, and a final "SUMMARY ..." line per type:
 *   count of graphs, number of distinct keys, number of keys hit by >= 2 graphs.
 * Environment: gcc -O3 -march=native; single thread; memory ~ 24 bytes per graph per type.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXM 128
typedef uint64_t u64;

static int nv;
static int adj[MAXM][8];
static int deg[MAXM];

/* ---------- planar code reading ---------- */
static int read_graph(FILE *f) {
    int c = getc(f);
    if (c == EOF) return 0;
    if (c == 0) { fprintf(stderr, "nv>255 unsupported\n"); exit(2); }
    nv = c;
    if (nv > MAXM) { fprintf(stderr, "nv > MAXM\n"); exit(2); }
    for (int v = 0; v < nv; v++) {
        deg[v] = 0;
        while (1) {
            int x = getc(f);
            if (x == EOF) { fprintf(stderr, "truncated input\n"); exit(2); }
            if (x == 0) break;
            if (deg[v] >= 8) { fprintf(stderr, "degree > 8\n"); exit(2); }
            adj[v][deg[v]++] = x - 1;
        }
    }
    return 1;
}

/* ---------- traces of M^k mod 2^64 for M = a*A_sub + diag(dd) on vertex subset ---------- */
#ifdef WORD32
typedef uint32_t wt;   /* traces mod 2^32 (still a ring homomorphism Z -> Z/2^32: no false negatives) */
#else
typedef uint64_t wt;   /* traces mod 2^64 */
#endif
static wt X[MAXM][MAXM] __attribute__((aligned(64))), Y[MAXM][MAXM] __attribute__((aligned(64)));

static void traces(int k, const int *verts, const u64 *dd, u64 a, int K, u64 *t) {
    /* local relabelling of the induced subgraph */
    static int loc[MAXM];
    static int ladj[MAXM][8], ldeg[MAXM];
    for (int v = 0; v < nv; v++) loc[v] = -1;
    for (int i = 0; i < k; i++) loc[verts[i]] = i;
    for (int i = 0; i < k; i++) {
        int v = verts[i];
        ldeg[i] = 0;
        for (int j = 0; j < deg[v]; j++) {
            int w = loc[adj[v][j]];
            if (w >= 0) ladj[i][ldeg[i]++] = w;
        }
    }
    for (int i = 0; i < k; i++) {
        memset(X[i], 0, sizeof(wt) * k);
        X[i][i] = 1;
    }
    /* t[1..K]; t[0] = k */
    t[0] = (u64)k;
    int steps = (K + 1) / 2;
    for (int s = 0; s < steps; s++) {
        /* Y = M X */
        for (int u = 0; u < k; u++) {
            wt *yu = Y[u];
            const wt *xu = X[u];
            wt d = (wt)dd[u];
            wt aa = (wt)a;
            for (int c = 0; c < k; c++) yu[c] = d * xu[c];
            for (int j = 0; j < ldeg[u]; j++) {
                const wt *xw = X[ladj[u][j]];
                if (a == 1) { for (int c = 0; c < k; c++) yu[c] += xw[c]; }
                else { for (int c = 0; c < k; c++) yu[c] += aa * xw[c]; }
            }
        }
        /* X = M^s, Y = M^{s+1}: tr(M^{2s}) = sum X^2 (already have if s>0), tr(M^{2s+1}) = sum X.*Y,
           tr(M^{2s+2}) = sum Y.*Y */
        wt s1 = 0, s2 = 0;
        for (int u = 0; u < k; u++) {
            const wt *xu = X[u], *yu = Y[u];
            for (int c = 0; c < k; c++) { s1 += xu[c] * yu[c]; s2 += yu[c] * yu[c]; }
        }
        if (2 * s + 1 <= K) t[2 * s + 1] = (u64)s1;
        if (2 * s + 2 <= K) t[2 * s + 2] = (u64)s2;
        for (int u = 0; u < k; u++) memcpy(X[u], Y[u], sizeof(wt) * k);
    }
}

/* ---------- hashing ---------- */
static inline u64 mix64(u64 z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void hash128(const u64 *t, int K, u64 *h1, u64 *h2) {
    u64 a = 0x243f6a8885a308d3ULL, b = 0x13198a2e03707344ULL;
    for (int i = 0; i <= K; i++) {
        a = mix64(a ^ (t[i] + 0x9e3779b97f4a7c15ULL * (u64)(i + 1)));
        b = mix64(b + mix64(t[i] ^ (0xa4093822299f31d0ULL + (u64)i)));
    }
    *h1 = a; *h2 = b;
}

typedef struct { u64 h1, h2; uint32_t idx; uint32_t cnt; } entry;
typedef struct { entry *tab; u64 cap, used, multi; } table;

static void tinit(table *T, u64 cap) {
    T->cap = cap; T->used = 0; T->multi = 0;
    T->tab = calloc(cap, sizeof(entry));
    if (!T->tab) { fprintf(stderr, "out of memory\n"); exit(3); }
}
static entry *tfind(table *T, u64 h1, u64 h2) {
    u64 p = h1 & (T->cap - 1);
    while (1) {
        entry *e = &T->tab[p];
        if (e->cnt == 0) return e;
        if (e->h1 == h1 && e->h2 == h2) return e;
        p = (p + 1) & (T->cap - 1);
    }
}
static void tgrow(table *T) {
    table N; tinit(&N, T->cap * 2);
    for (u64 i = 0; i < T->cap; i++) if (T->tab[i].cnt) {
        entry *e = tfind(&N, T->tab[i].h1, T->tab[i].h2);
        *e = T->tab[i];
    }
    N.used = T->used; N.multi = T->multi;
    free(T->tab); *T = N;
}
/* returns index of previous graph with the same key, or -1 */
static long tinsert(table *T, u64 h1, u64 h2, uint32_t idx) {
    if (2 * (T->used + 1) > T->cap) tgrow(T);
    entry *e = tfind(T, h1, h2);
    if (e->cnt == 0) { e->h1 = h1; e->h2 = h2; e->idx = idx; e->cnt = 1; T->used++; return -1; }
    if (e->cnt == 1) T->multi++;
    e->cnt++;
    return (long)e->idx;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: dualspec MASK [K]\n"); return 1; }
    int mask = atoi(argv[1]);
    int Karg = argc > 2 ? atoi(argv[2]) : 0;
    int print5 = argc > 3 ? atoi(argv[3]) : 0; /* print T5 collisions too */
    char hdr[16];
    if (fread(hdr, 1, 15, stdin) != 15 || memcmp(hdr, ">>planar_code<<", 15) != 0) {
        fprintf(stderr, "bad header\n"); return 2;
    }
    table T[6];
    for (int ty = 0; ty < 6; ty++) tinit(&T[ty], 1 << 16);
    u64 count = 0;
    u64 t[2 * MAXM + 2];
    int all[MAXM], v6[MAXM], v5[MAXM];
    u64 dzero[MAXM];
    const u64 TT = 0x9e3779b97f4a7c15ULL; /* odd constant for A + tP */
    long long nsum5 = 0, nsum6 = 0;
    int m_seen = -1;
    while (read_graph(stdin)) {
        if (m_seen < 0) m_seen = nv;
        if (nv != m_seen) { fprintf(stderr, "mixed sizes\n"); return 2; }
        int n5 = 0, n6 = 0;
        for (int v = 0; v < nv; v++) {
            all[v] = v; dzero[v] = 0;
            if (deg[v] == 5) { v5[n5++] = v; }
            else if (deg[v] == 6) { v6[n6++] = v; }
            else { fprintf(stderr, "degree not 5/6 in graph %llu\n", (unsigned long long)count); return 2; }
        }
        if (n5 != 12) { fprintf(stderr, "n5 != 12\n"); return 2; }
        nsum5 += n5; nsum6 += n6;
        uint32_t idx = (uint32_t)count;
        for (int ty = 0; ty < 6; ty++) {
            if (!(mask & (1 << ty))) continue;
            int k; const int *verts; u64 a = 1; u64 dd[MAXM];
            if (ty == 1) { k = n6; verts = v6; for (int i = 0; i < k; i++) dd[i] = 0; }
            else if (ty == 2) { k = n5; verts = v5; for (int i = 0; i < k; i++) dd[i] = 0; }
            else {
                k = nv; verts = all;
                for (int i = 0; i < k; i++) {
                    if (ty == 0) dd[i] = 0;
                    else if (ty == 3) { dd[i] = (u64)deg[i]; }
                    else if (ty == 4) { dd[i] = (u64)deg[i]; }
                    else dd[i] = (deg[i] == 5) ? TT : 0;
                }
                if (ty == 3) a = 2;
            }
            int K = Karg > 0 ? (Karg < k ? Karg : k) : k;
            if (k == 0) { t[0] = 0; K = 0; }
            else traces(k, verts, dd, a, K, t);
            u64 h1, h2;
            hash128(t, K, &h1, &h2);
            long prev = tinsert(&T[ty], h1, h2, idx);
            if (prev >= 0 && (ty != 2 || print5)) printf("C %d %ld %u\n", ty, prev, idx);
        }
        count++;
    }
    (void)dzero;
    const char *names[6] = {"A(T)", "A(T6)", "A(T5)", "2A+D(T)", "A+D(T)", "A+tP(T)"};
    for (int ty = 0; ty < 6; ty++) {
        if (!(mask & (1 << ty))) continue;
        u64 ingraphs = 0;
        for (u64 q = 0; q < T[ty].cap; q++) if (T[ty].tab[q].cnt > 1) ingraphs += T[ty].tab[q].cnt;
        printf("SUMMARY graphs_in_multiple_keys=%llu ", (unsigned long long)ingraphs);
        printf("SUMMARY type=%d name=%s m=%d graphs=%llu distinct_keys=%llu keys_with_multiple_graphs=%llu K=%d\n",
               ty, names[ty], m_seen, (unsigned long long)count, (unsigned long long)T[ty].used,
               (unsigned long long)T[ty].multi, Karg);
    }
    return 0;
}
