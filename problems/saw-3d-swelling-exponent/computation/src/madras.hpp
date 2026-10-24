// madras.hpp -- Madras's reflection map T (Madras, Canad. Math. Bull. 57 (2014), Section 2, as restated in
// lower-bound/result.md Section 1) and the exact multiplicity K(T omega) = #{theta : ||theta(N)||_inf <= m,
// T theta = T omega}, computed from the definitions (written independently of lower-bound/scripts/).
//
// Definitions (omega(0) = 0): span = max_k ||omega(k)||_inf;  L = min{k > 0 : ||omega(k)||_inf = span};
// I = min{i : |omega_i(L)| = span};  H = omega_I(L);  (T omega)(k) = omega(k) for k <= L and the reflection of
// omega(k) through the hyperplane {x_I = H} for k > L.
//
// Candidates for theta with T theta = psi: theta is determined by its Madras data (L', I', H'):
// theta = psi[0, L'] followed by the reflection of psi[L', N] through {x_I' = H'}, where L' is the first k with
// psi_I'(k) = H' (theta and psi agree up to L', and ||theta(k)||_inf < |H'| for k < L').  The endpoint condition
// ||theta(N)||_inf <= m gives |psi_I'(N) - 2H'| <= m and |psi_i(N)| <= m for i != I'.  For each such (I', H') we build
// theta and accept it iff theta is self-avoiding, ||theta(N)||_inf <= m and Madras(theta) = (L', I', H')
// (then T theta = psi automatically, the reflection being an involution).  Partial counts:
//   cand  = number of (I', H') satisfying the endpoint condition;
//   ren   = those for which L' exists and psi[L', N] lies on the far side (sign(H') psi_I'(k) >= |H'| for k >= L'),
//           equivalently theta[L', N] does not go beyond the level H';
//   saw   = those for which, in addition, theta is self-avoiding.
#pragma once
#include <cstdint>
#include <cstdlib>
#include <vector>

struct P3 {
  int32_t c[3];
};
inline bool eqP(const P3& a, const P3& b) { return a.c[0] == b.c[0] && a.c[1] == b.c[1] && a.c[2] == b.c[2]; }
inline int ninf(const P3& p) {
  int a = std::abs(p.c[0]), b = std::abs(p.c[1]), c = std::abs(p.c[2]);
  return a > b ? (a > c ? a : c) : (b > c ? b : c);
}
inline uint64_t keyP(const P3& p) {
  const int64_t O = 1 << 20;
  return ((uint64_t)(p.c[0] + O) << 42) | ((uint64_t)(p.c[1] + O) << 21) | (uint64_t)(p.c[2] + O);
}
inline uint64_t hmix(uint64_t k) {
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  return k;
}

// hash set with generation stamps (clear in O(1))
struct GenSet {
  std::vector<uint64_t> K;
  std::vector<uint32_t> G;
  uint32_t g = 1;
  uint64_t mask = 0;
  void init(size_t n) {
    size_t cap = 1;
    while (cap < 4 * n + 16) cap <<= 1;
    K.assign(cap, 0);
    G.assign(cap, 0);
    mask = cap - 1;
    g = 1;
  }
  void clear() {
    if (++g == 0) {
      std::fill(G.begin(), G.end(), 0);
      g = 1;
    }
  }
  bool insert(uint64_t k) {  // false if already present
    uint64_t i = hmix(k) & mask;
    while (G[i] == g) {
      if (K[i] == k) return false;
      i = (i + 1) & mask;
    }
    G[i] = g;
    K[i] = k;
    return true;
  }
  bool has(uint64_t k) const {
    uint64_t i = hmix(k) & mask;
    while (G[i] == g) {
      if (K[i] == k) return true;
      i = (i + 1) & mask;
    }
    return false;
  }
};

struct MData {
  int span, L, I, H;
};

inline MData madrasData(const std::vector<P3>& w) {
  int n = (int)w.size() - 1;
  MData d{0, -1, -1, 0};
  for (int k = 0; k <= n; k++) {
    int v = ninf(w[k]);
    if (v > d.span) d.span = v;
  }
  for (int k = 1; k <= n; k++)
    if (ninf(w[k]) == d.span) {
      d.L = k;
      break;
    }
  for (int i = 0; i < 3; i++)
    if (std::abs(w[d.L].c[i]) == d.span) {
      d.I = i;
      break;
    }
  d.H = w[d.L].c[d.I];
  return d;
}

inline void reflectAfter(const std::vector<P3>& w, int L, int I, int H, std::vector<P3>& out) {
  out = w;
  for (size_t k = L + 1; k < w.size(); k++) out[k].c[I] = 2 * H - w[k].c[I];
}

struct KResult {
  int K, cand, ren, saw;
};

// K(T omega) for the near set {||theta(N)||_inf <= m}; omega must satisfy ||omega(N)||_inf <= m.
inline KResult multiplicity(const std::vector<P3>& omega, int m, GenSet& hs, std::vector<P3>& psi,
                            std::vector<P3>& th) {
  int n = (int)omega.size() - 1;
  MData d = madrasData(omega);
  reflectAfter(omega, d.L, d.I, d.H, psi);
  KResult r{0, 0, 0, 0};
  static thread_local std::vector<int> first;
  for (int Ip = 0; Ip < 3; Ip++) {
    bool okOther = true;
    for (int i = 0; i < 3; i++)
      if (i != Ip && std::abs(psi[n].c[i]) > m) okOther = false;
    if (!okOther) continue;
    // first hitting time of every level of coordinate Ip; |psi_i(k)| <= 3 span <= 3n (psi_I = 2H - omega_I)
    first.assign(6 * n + 3, -1);
    for (int k = n; k >= 0; k--) first[psi[k].c[Ip] + 3 * n + 1] = k;
    int b = psi[n].c[Ip];
    // integers H' != 0 with |b - 2H'| <= m
    int lo = (b - m) >= 0 ? (b - m + 1) / 2 : -((m - b) / 2);  // ceil((b-m)/2)
    int hi = (b + m) >= 0 ? (b + m) / 2 : -((-(b + m) + 1) / 2);  // floor((b+m)/2)
    for (int Hp = lo; Hp <= hi; Hp++) {
      if (Hp == 0) continue;
      r.cand++;
      int Lp = (Hp >= -3 * n - 1 && Hp <= 3 * n + 1) ? first[Hp + 3 * n + 1] : -1;
      if (Lp < 1) continue;
      int sg = Hp > 0 ? 1 : -1, Ha = Hp > 0 ? Hp : -Hp;
      bool far = true;
      for (int k = Lp; k <= n; k++)
        if (sg * psi[k].c[Ip] < Ha) {
          far = false;
          break;
        }
      if (!far) continue;
      r.ren++;
      reflectAfter(psi, Lp, Ip, Hp, th);
      hs.clear();
      bool saw = true;
      for (int k = 0; k <= n; k++)
        if (!hs.insert(keyP(th[k]))) {
          saw = false;
          break;
        }
      if (!saw) continue;
      r.saw++;
      if (ninf(th[n]) > m) continue;
      MData d2 = madrasData(th);
      if (d2.L == Lp && d2.I == Ip && d2.H == Hp) r.K++;
    }
  }
  return r;
}
