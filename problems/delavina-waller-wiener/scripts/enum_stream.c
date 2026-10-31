/* enum_stream.c  --  read connected graphs in graph6 (n <= 16) from stdin and, for each
 * diameter, keep the K largest Wiener indices with a histogram of the top values.
 *
 * Usage: geng -c N [mine:maxe] [res/mod] | ./enum_stream N OUTFILE
 * Output: for each diameter, "diam max W, #graphs, count of top values, and one graph6 per top value".
 * Deterministic. Compile: gcc -O3 -march=native -o enum_stream enum_stream.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 16
#define MAXW 2000
static long long hist[MAXN][MAXW];
static char sample[MAXN][MAXW][32];
static char sample2[MAXN][MAXW][8][32];
static long long total;

int main(int argc, char **argv)
{
    int N = atoi(argv[1]);
    FILE *out = fopen(argv[2], "w");
    char line[128];
    uint32_t adj[MAXN];
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line);
        if (len && line[len - 1] == '\n') line[--len] = 0;
        int n = line[0] - 63;
        memset(adj, 0, sizeof adj);
        int pos = 1, bit = 0, x = 0, cnt = 0;
        for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
            if (bit == 0) { x = line[pos++] - 63; bit = 6; }
            bit--;
            if ((x >> bit) & 1) { adj[i] |= 1u << j; adj[j] |= 1u << i; }
        }
        uint32_t all = (1u << n) - 1;
        int W2 = 0, diam = 0, conn = 1;
        for (int s = 0; s < n; s++) {
            uint32_t seen = 1u << s, fr = seen;
            int k = 0;
            while (seen != all) {
                uint32_t nx = 0, f = fr;
                while (f) { int v = __builtin_ctz(f); f &= f - 1; nx |= adj[v]; }
                nx &= ~seen;
                if (!nx) { conn = 0; break; }
                k++;
                W2 += k * __builtin_popcount(nx);
                seen |= nx; fr = nx;
            }
            if (!conn) break;
            if (k > diam) diam = k;
        }
        if (!conn) continue;
        total++;
        int W = W2 / 2;
        if (hist[diam][W] < 8) strncpy(sample2[diam][W][hist[diam][W]], line, 31);
        hist[diam][W]++;
    }
    fprintf(out, "n=%d total_connected_read=%lld\n", N, total);
    for (int d = 1; d < MAXN; d++) {
        int top = -1;
        long long c = 0;
        for (int w = 0; w < MAXW; w++) if (hist[d][w]) { top = w; c += hist[d][w]; }
        if (top < 0) continue;
        fprintf(out, "diam=%d count=%lld maxW=%d\n", d, c, top);
        int shown = 0;
        for (int w = top; w >= 0 && shown < 8; w--) if (hist[d][w]) {
            fprintf(out, "  W=%d count=%lld :", w, hist[d][w]);
            for (int i = 0; i < 8 && i < hist[d][w]; i++) fprintf(out, " %s", sample2[d][w][i]);
            fprintf(out, "\n");
            shown++;
        }
    }
    fclose(out);
    return 0;
}
