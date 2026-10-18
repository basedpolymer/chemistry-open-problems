// Exploration tool: count Ponitz-Tittmann-type states for Z^d (d<=3 here).
// A state candidate is a SAW w = (B=v_0, ..., v_l=A) such that every prefix satisfies j + |v_j - B|_1 <= k
// and l + |A-B|_1 <= k ("loop size" with Manhattan closing distance).
// Usage: count_states d k
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
static int D, K;
static const int R = 40;
static unsigned char occ[2 * R + 1][2 * R + 1][2 * R + 1];
static long long cnt[64];
static long long cntd[64][64];
static int pos[3];
static void dfs(int l) {
    int d = abs(pos[0]) + abs(pos[1]) + abs(pos[2]);
    if (l + d > K) return;
    cnt[l]++;
    cntd[l][d]++;
    if (l == K) return;
    for (int a = 0; a < D; a++)
        for (int s = -1; s <= 1; s += 2) {
            pos[a] += s;
            unsigned char &o = occ[pos[0] + R][pos[1] + R][pos[2] + R];
            if (!o) { o = 1; dfs(l + 1); o = 0; }
            pos[a] -= s;
        }
}
int main(int argc, char **argv) {
    D = atoi(argv[1]);
    K = atoi(argv[2]);
    occ[R][R][R] = 1;
    dfs(0);
    long long tot = 0;
    for (int l = 0; l <= K; l++) {
        printf("l=%2d count=%lld  dist:", l, cnt[l]);
        for (int d = 0; d <= K; d++)
            if (cntd[l][d]) printf(" %d:%lld", d, cntd[l][d]);
        printf("\n");
        tot += cnt[l];
    }
    printf("total=%lld  /(2^d d!) approx %.4g\n", tot, tot / (D == 3 ? 48.0 : D == 2 ? 8.0 : 2.0));
}
