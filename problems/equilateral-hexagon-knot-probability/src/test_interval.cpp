// Reads lines "op alo ahi blo bhi" (hex-float doubles) and prints the interval result (hex floats).
// Used by scripts/test_interval.py, which checks enclosure against exact rational arithmetic.
#include <cstdio>
#include <cstring>
#include "interval.hpp"
int main() {
    hex_init_rounding();
    char op[8]; double a, b, c, d;
    while (scanf("%7s %la %la %la %la", op, &a, &b, &c, &d) == 5) {
        I x(a, b), y(c, d), r;
        if (!strcmp(op, "add")) r = x + y;
        else if (!strcmp(op, "sub")) r = x - y;
        else if (!strcmp(op, "mul")) r = x * y;
        else if (!strcmp(op, "div")) r = x / y;
        else if (!strcmp(op, "sqrt")) r = isqrt(x);
        printf("%a %a\n", r.lo, r.hi);
    }
    return 0;
}
