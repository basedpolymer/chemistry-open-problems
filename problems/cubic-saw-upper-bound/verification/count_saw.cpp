// count_saw.cpp -- brute-force count of n-step self-avoiding walks on Z^d (independent of every certificate).
// c_n = 2d * (number of SAWs whose first step is +e_1). Single thread, plain depth-first search on a hash-free
// occupancy grid. Build: g++ -O2 -march=native count_saw.cpp -o count_saw ; usage: count_saw d nmax
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long u64;
static int D, NMAX, W;
static vector<char> occ;
static long long stride[8];
static vector<u64> cnt;
static void dfs(long long p, int n) {
    cnt[n]++;
    if (n == NMAX) return;
    for (int a = 0; a < D; a++)
        for (int s = -1; s <= 1; s += 2) {
            long long q = p + s * stride[a];
            if (occ[q]) continue;
            occ[q] = 1; dfs(q, n + 1); occ[q] = 0;
        }
}
int main(int argc, char** argv) {
    D = atoi(argv[1]); NMAX = atoi(argv[2]);
    W = 2 * NMAX + 3;
    long long tot = 1;
    for (int a = 0; a < D; a++) { stride[a] = tot; tot *= W; }
    occ.assign(tot, 0); cnt.assign(NMAX + 1, 0);
    long long o = 0; for (int a = 0; a < D; a++) o += (NMAX + 1) * stride[a];
    occ[o] = 1; occ[o + stride[0]] = 1;
    dfs(o + stride[0], 1);  // walks starting with +e_1
    printf("0 1\n");
    for (int n = 1; n <= NMAX; n++) printf("%d %llu\n", n, (unsigned long long)(cnt[n] * 2 * D));
}
