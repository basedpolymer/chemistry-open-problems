// Minimal rigorous interval arithmetic (outward rounding by one ulp after each IEEE operation)
// and first-order forward automatic differentiation over a box (value and gradient enclosures).
//
// Rigour: IEEE-754 binary64 +, -, *, /, sqrt are correctly rounded (round to nearest), so the exact
// result lies within one ulp of the computed one; pushing the lower end one ulp down and the upper end
// one ulp up yields an enclosure.  Must be compiled WITHOUT -ffast-math and with -ffp-contract=off
// (no fused multiply-add contraction) on an SSE2 target (x86-64 default).  No libm transcendental
// function is used anywhere: trigonometric values come from a table certified with mpmath.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>

static inline double nup(double x) {
    if (std::isnan(x)) return x;
    if (x == std::numeric_limits<double>::infinity()) return x;
    if (x == 0.0) return 4.9406564584124654e-324;
    uint64_t u; std::memcpy(&u, &x, 8);
    if (x > 0) u++; else u--;
    std::memcpy(&x, &u, 8);
    return x;
}
static inline double ndn(double x) { return -nup(-x); }

struct I {
    double lo, hi;
    I() : lo(0), hi(0) {}
    I(double a) : lo(a), hi(a) {}
    I(double a, double b) : lo(a), hi(b) {}
    bool pos() const { return lo > 0; }
    bool neg() const { return hi < 0; }
    double mid() const { return 0.5 * (lo + hi); }
};
static const double INF = std::numeric_limits<double>::infinity();
static inline I entire() { return I(-INF, INF); }
static inline I operator+(I a, I b) { return I(ndn(a.lo + b.lo), nup(a.hi + b.hi)); }
static inline I operator-(I a, I b) { return I(ndn(a.lo - b.hi), nup(a.hi - b.lo)); }
static inline I operator-(I a) { return I(-a.hi, -a.lo); }
static inline I operator*(I a, I b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    double lo = std::min(std::min(p1, p2), std::min(p3, p4));
    double hi = std::max(std::max(p1, p2), std::max(p3, p4));
    if (std::isnan(lo) || std::isnan(hi)) return entire();   // 0*inf
    return I(ndn(lo), nup(hi));
}
static inline I operator/(I a, I b) {
    if (b.lo <= 0 && b.hi >= 0) return entire();
    double p1 = a.lo / b.lo, p2 = a.lo / b.hi, p3 = a.hi / b.lo, p4 = a.hi / b.hi;
    double lo = std::min(std::min(p1, p2), std::min(p3, p4));
    double hi = std::max(std::max(p1, p2), std::max(p3, p4));
    if (std::isnan(lo) || std::isnan(hi)) return entire();
    return I(ndn(lo), nup(hi));
}
// sqrt of an interval whose exact (mathematical) argument is known to be >= 0
static inline I isqrt(I a) {
    double lo = a.lo <= 0 ? 0.0 : std::max(0.0, ndn(std::sqrt(a.lo)));
    double hi = a.hi <= 0 ? 0.0 : nup(std::sqrt(a.hi));
    if (a.hi == INF) hi = INF;
    return I(lo, hi);
}
static inline I hull(I a, I b) { return I(std::min(a.lo, b.lo), std::max(a.hi, b.hi)); }
static inline I meet(I a, I b) { return I(std::max(a.lo, b.lo), std::min(a.hi, b.hi)); }

// first-order AD with interval value and interval gradient (3 variables)
struct A {
    I v; I g[3];
    A() {}
    A(double c) : v(c) { g[0] = g[1] = g[2] = I(0.0); }
    A(I c) : v(c) { g[0] = g[1] = g[2] = I(0.0); }
};
static inline A operator+(const A& a, const A& b) { A r; r.v = a.v + b.v; for (int k = 0; k < 3; k++) r.g[k] = a.g[k] + b.g[k]; return r; }
static inline A operator-(const A& a, const A& b) { A r; r.v = a.v - b.v; for (int k = 0; k < 3; k++) r.g[k] = a.g[k] - b.g[k]; return r; }
static inline A operator-(const A& a) { A r; r.v = -a.v; for (int k = 0; k < 3; k++) r.g[k] = -a.g[k]; return r; }
static inline A operator*(const A& a, const A& b) { A r; r.v = a.v * b.v; for (int k = 0; k < 3; k++) r.g[k] = a.g[k] * b.v + a.v * b.g[k]; return r; }
static inline A operator/(const A& a, const A& b) { A r; r.v = a.v / b.v; for (int k = 0; k < 3; k++) r.g[k] = (a.g[k] - r.v * b.g[k]) / b.v; return r; }
static inline A asqrt(const A& a) { A r; r.v = isqrt(a.v); I two_r = I(2.0) * r.v; for (int k = 0; k < 3; k++) r.g[k] = a.g[k] / two_r; return r; }
