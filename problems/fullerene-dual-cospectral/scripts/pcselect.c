/*
 * pcselect.c -- copy selected graphs (0-based indices, in input order) from a planar_code stream.
 * Usage: buckygen -q 100d | pcselect indexfile > selected.pc
 *   indexfile: whitespace-separated non-negative integers (any order, duplicates allowed).
 * Output: planar_code with header, the selected graphs in increasing index order; a text line per graph
 * "index position" is written to stderr so that the caller can map positions back to indices.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmp(const void *a, const void *b) {
    long x = *(const long *)a, y = *(const long *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: pcselect indexfile\n"); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("indexfile"); return 1; }
    long cap = 1024, n = 0, *idx = malloc(cap * sizeof(long)), x;
    while (fscanf(f, "%ld", &x) == 1) {
        if (n == cap) { cap *= 2; idx = realloc(idx, cap * sizeof(long)); }
        idx[n++] = x;
    }
    fclose(f);
    qsort(idx, n, sizeof(long), cmp);
    long u = 0;
    for (long i = 0; i < n; i++) if (u == 0 || idx[i] != idx[u - 1]) idx[u++] = idx[i];
    n = u;
    char hdr[16];
    if (fread(hdr, 1, 15, stdin) != 15 || memcmp(hdr, ">>planar_code<<", 15) != 0) {
        fprintf(stderr, "bad header\n"); return 2;
    }
    fwrite(hdr, 1, 15, stdout);
    static unsigned char buf[1 << 16];
    long gi = 0, p = 0, pos = 0;
    int c;
    while (p < n && (c = getc(stdin)) != EOF) {
        int nv = c, len = 0;
        buf[len++] = (unsigned char)c;
        for (int v = 0; v < nv; v++) {
            while (1) {
                int y = getc(stdin);
                if (y == EOF) { fprintf(stderr, "truncated\n"); return 2; }
                buf[len++] = (unsigned char)y;
                if (y == 0) break;
            }
        }
        if (gi == idx[p]) {
            fwrite(buf, 1, len, stdout);
            fprintf(stderr, "%ld %ld\n", gi, pos++);
            p++;
        }
        gi++;
    }
    if (p < n) { fprintf(stderr, "index beyond end of input\n"); return 3; }
    return 0;
}
