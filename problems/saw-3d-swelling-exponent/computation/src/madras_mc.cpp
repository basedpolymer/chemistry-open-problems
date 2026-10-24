// madras_mc.cpp -- HEURISTIC Monte Carlo (pivot algorithm) for the lower-bound solver's candidate statements
// CL1 (average multiplicity of Madras's map on the near set), CL2 (renewal depths at the top) and CL5 (span).
// NUMERICS ONLY: nothing here is a proof.  Written independently of lower-bound/scripts/ from the definitions in
// lower-bound/result.md (Madras data, renewal depths) and lower-bound/candidate-lemmas.md; K(T omega) is computed by
// src/madras.hpp, validated exhaustively against direct image counts (src/madras_exact.cpp, n <= 9).
//
// Chain: pivot algorithm as in src/pivot.cpp (shorter side pivoted; the walk is re-anchored at omega(0) = 0 before
// each measurement; the anchored walk is uniformly distributed on S_N in equilibrium, orientation included).
// Measurements:
//   every SN attempts: if ||omega(N)||_inf <= m for some m in MS, compute K, cand, ren, saw (madras.hpp) for every
//     such m; ratio estimators kappa_N(m) = sum K / #near etc., errors by jackknife over 64 consecutive blocks;
//   every SR attempts: Madras data, span, tail/head renewal indicators at depths j = 1..JMAX (unconditional).
// Usage: madras_mc N attempts seed [SN=10] [SR=100]
// Build: g++ -O3 -march=native -std=c++17 -Isrc madras_mc.cpp -o madras_mc
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "madras.hpp"

struct Rng {  // xoshiro256**
  uint64_t s[4];
  static uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  explicit Rng(uint64_t seed) {
    for (int i = 0; i < 4; i++) s[i] = splitmix(seed);
  }
  static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  inline uint64_t next() {
    uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return r;
  }
  inline uint64_t below(uint64_t m) { return (uint64_t)(((unsigned __int128)next() * m) >> 64); }
};

struct HashMap {  // site -> index; open addressing, linear probing, backward-shift deletion
  std::vector<uint64_t> K;
  std::vector<int32_t> V;
  uint64_t mask;
  void init(size_t n) {
    size_t cap = 1;
    while (cap < 4 * n) cap <<= 1;
    K.assign(cap, 0);
    V.assign(cap, -1);
    mask = cap - 1;
  }
  int32_t find(uint64_t k) const {
    uint64_t i = hmix(k) & mask;
    while (K[i]) {
      if (K[i] == k) return V[i];
      i = (i + 1) & mask;
    }
    return -1;
  }
  void insert(uint64_t k, int32_t v) {
    uint64_t i = hmix(k) & mask;
    while (K[i] && K[i] != k) i = (i + 1) & mask;
    K[i] = k;
    V[i] = v;
  }
  void erase(uint64_t k) {
    uint64_t i = hmix(k) & mask;
    while (K[i] != k) {
      if (!K[i]) return;
      i = (i + 1) & mask;
    }
    uint64_t j = i;
    for (;;) {
      j = (j + 1) & mask;
      if (!K[j]) break;
      uint64_t h = hmix(K[j]) & mask;
      bool inRange = (i <= j) ? (h > i && h <= j) : (h > i || h <= j);
      if (!inRange) {
        K[i] = K[j];
        V[i] = V[j];
        i = j;
      }
    }
    K[i] = 0;
    V[i] = -1;
  }
};

struct Sym {
  int perm[3], sign[3];
};
static Sym SYM[48];
static void initSym() {
  int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  int k = 0;
  for (int p = 0; p < 6; p++)
    for (int s = 0; s < 8; s++) {
      for (int i = 0; i < 3; i++) {
        SYM[k].perm[i] = perms[p][i];
        SYM[k].sign[i] = (s >> i) & 1 ? -1 : 1;
      }
      k++;
    }
}

int main(int argc, char** argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s N attempts seed [SN] [SR]\n", argv[0]);
    return 2;
  }
  int n = atoi(argv[1]);
  long long T = atoll(argv[2]);
  uint64_t seed = strtoull(argv[3], 0, 10);
  int SN = argc > 4 ? atoi(argv[4]) : 10, SR = argc > 5 ? atoi(argv[5]) : 100;
  const int NM = 8;
  const int MS[NM] = {2, 4, 8, 16, 32, 64, 128, 256};
  const int JMAX = 512;
  const int NBLK = 64;
  initSym();
  Rng rng(seed * 7919ULL + (uint64_t)n);
  std::vector<P3> w(n + 1), nw(n + 1), anc(n + 1), psi, th;
  for (int i = 0; i <= n; i++) w[i] = P3{{i, 0, 0}};
  HashMap H;
  H.init(n + 1);
  for (int i = 0; i <= n; i++) H.insert(keyP(w[i]), i);
  GenSet hs;
  hs.init(n + 1);
  long long B = 20LL * n;  // burn-in attempts
  long long acc = 0;
  // accumulators per block
  std::vector<double> bNear(NBLK * NM, 0), bK(NBLK * NM, 0), bCand(NBLK * NM, 0), bRen(NBLK * NM, 0),
      bSaw(NBLK * NM, 0), bSpanLe(NBLK * NM, 0);
  std::vector<long long> histK(NM * 10, 0);
  long long nNearSamples = 0, nSamplesNear = 0;
  std::vector<double> rt(JMAX + 1, 0), rh(JMAX + 1, 0);
  std::vector<double> bRt(NBLK * 8, 0), bRh(NBLK * 8, 0);  // block sums at j = 1,2,4,...,128 for errors
  long long nR = 0;
  std::vector<double> spans;
  std::vector<int> ys(n + 1);
  auto t0 = std::chrono::steady_clock::now();
  for (long long it = 0; it < B + T; it++) {
    int k = 1 + (int)rng.below(n - 1);
    const Sym& g = SYM[1 + rng.below(47)];
    bool after = (k >= n - k);
    int lo = after ? k + 1 : 0, hi = after ? n : k - 1;
    P3 c = w[k];
    bool ok = true;
    int len = hi - lo + 1;
    for (int t = 0; t < len; t++) {
      int j = after ? (k + 1 + t) : (k - 1 - t);
      int d[3] = {w[j].c[0] - c.c[0], w[j].c[1] - c.c[1], w[j].c[2] - c.c[2]};
      P3 q;
      for (int i = 0; i < 3; i++) q.c[i] = c.c[i] + g.sign[i] * d[g.perm[i]];
      int32_t o = H.find(keyP(q));
      if (o >= 0 && (after ? (o <= k) : (o >= k))) {
        ok = false;
        break;
      }
      nw[j] = q;
    }
    if (ok) {
      acc++;
      for (int j = lo; j <= hi; j++) H.erase(keyP(w[j]));
      for (int j = lo; j <= hi; j++) {
        w[j] = nw[j];
        H.insert(keyP(w[j]), j);
      }
      if ((acc & 4095) == 0) {  // re-centre to keep coordinates bounded
        P3 s = w[0];
        H.init(n + 1);
        for (int i = 0; i <= n; i++) {
          for (int a = 0; a < 3; a++) w[i].c[a] -= s.c[a];
          H.insert(keyP(w[i]), i);
        }
      }
    }
    if (it < B) continue;
    long long s = it - B;
    int blk = (int)(s * NBLK / T);
    if (s % SN == 0) {
      int en = 0;
      for (int a = 0; a < 3; a++) {
        int v = std::abs(w[n].c[a] - w[0].c[a]);
        if (v > en) en = v;
      }
      if (en <= MS[NM - 1]) {
        nSamplesNear++;
        for (int i = 0; i <= n; i++)
          for (int a = 0; a < 3; a++) anc[i].c[a] = w[i].c[a] - w[0].c[a];
        MData d = madrasData(anc);
        for (int mi = 0; mi < NM; mi++) {
          int m = MS[mi];
          if (en > m) continue;
          KResult r = multiplicity(anc, m, hs, psi, th);
          bNear[blk * NM + mi] += 1;
          bK[blk * NM + mi] += r.K;
          bCand[blk * NM + mi] += r.cand;
          bRen[blk * NM + mi] += r.ren;
          bSaw[blk * NM + mi] += r.saw;
          bSpanLe[blk * NM + mi] += (d.span <= m);
          histK[mi * 10 + (r.K < 9 ? r.K : 9)]++;
        }
      }
      nNearSamples++;
    }
    if (s % SR == 0) {
      for (int i = 0; i <= n; i++)
        for (int a = 0; a < 3; a++) anc[i].c[a] = w[i].c[a] - w[0].c[a];
      MData d = madrasData(anc);
      spans.push_back(d.span);
      int sg = d.H > 0 ? 1 : -1, Ha = std::abs(d.H);
      for (int i = 0; i <= n; i++) ys[i] = sg * anc[i].c[d.I];
      // tail: tau_j = first k >= L with y = Ha - j; renewal iff max_{k >= tau_j} y <= Ha - j
      // suffix max
      std::vector<int> smax(n + 2, -1000000000);
      for (int i = n; i >= 0; i--) smax[i] = std::max(ys[i], smax[i + 1]);
      int cur = Ha;
      for (int i = d.L; i <= n; i++) {
        while (ys[i] < cur) {  // levels cur-1 .. ys[i] reached for the first time (steps are 0 or +-1)
          cur--;
          int j = Ha - cur;
          if (j <= JMAX && smax[i] <= cur) {
            rt[j] += 1;
            for (int b = 0; b < 8; b++)
              if (j == (1 << b)) bRt[blk * 8 + b] += 1;
          }
        }
      }
      // head: sigma_j = first k with y = Ha - j (Ha - j >= 0); renewal iff min_{sigma_j <= k <= L} y >= Ha - j
      std::vector<int> smin(d.L + 2, 1000000000);
      for (int i = d.L; i >= 0; i--) smin[i] = std::min(ys[i], smin[i + 1]);
      int top = 0;  // y(0) = 0: level 0 first reached at k = 0
      {
        int lev = 0;
        int j = Ha - lev;
        if (j >= 1 && j <= JMAX && smin[0] >= lev) {
          rh[j] += 1;
          for (int b = 0; b < 8; b++)
            if (j == (1 << b)) bRh[blk * 8 + b] += 1;
        }
      }
      for (int i = 1; i <= d.L; i++) {
        while (ys[i] > top) {
          top++;
          int j = Ha - top;
          if (j >= 1 && j <= JMAX && smin[i] >= top) {
            rh[j] += 1;
            for (int b = 0; b < 8; b++)
              if (j == (1 << b)) bRh[blk * 8 + b] += 1;
          }
        }
      }
      nR++;
    }
  }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  printf("# madras_mc N=%d attempts=%lld burnin=%lld seed=%llu SN=%d SR=%d acc=%.4f time_s=%.0f\n", n, T, B,
         (unsigned long long)seed, SN, SR, (double)acc / (B + T), secs);
  // CL1: jackknife over blocks for ratios
  printf("# CL1 columns: m  n_near  P(near)  kappa err  cand err  ren err  saw err  frac(span<=m)  hist K=1..8,9+\n");
  for (int mi = 0; mi < NM; mi++) {
    double totN = 0, tot[5] = {0, 0, 0, 0, 0};
    for (int b = 0; b < NBLK; b++) {
      totN += bNear[b * NM + mi];
      tot[0] += bK[b * NM + mi];
      tot[1] += bCand[b * NM + mi];
      tot[2] += bRen[b * NM + mi];
      tot[3] += bSaw[b * NM + mi];
      tot[4] += bSpanLe[b * NM + mi];
    }
    if (totN == 0) {
      printf("M %d 0\n", MS[mi]);
      continue;
    }
    double est[5], err[5];
    for (int q = 0; q < 5; q++) {
      est[q] = tot[q] / totN;
      double sj = 0, sj2 = 0;
      int nb = 0;
      for (int b = 0; b < NBLK; b++) {
        double nN = totN - bNear[b * NM + mi];
        if (nN <= 0) continue;
        double v[5] = {bK[b * NM + mi], bCand[b * NM + mi], bRen[b * NM + mi], bSaw[b * NM + mi], bSpanLe[b * NM + mi]};
        double e = (tot[q] - v[q]) / nN;
        sj += e;
        sj2 += e * e;
        nb++;
      }
      double mean = sj / nb;
      err[q] = std::sqrt((double)(nb - 1) / nb * (sj2 - nb * mean * mean));
    }
    printf("M %d %.0f %.5f  %.4f %.4f  %.3f %.3f  %.4f %.4f  %.4f %.4f  %.4f ", MS[mi], totN, totN / nNearSamples,
           est[0], err[0], est[1], err[1], est[2], err[2], est[3], err[3], est[4]);
    for (int k = 1; k <= 9; k++) printf(" %lld", histK[mi * 10 + k]);
    printf("\n");
  }
  // CL2: renewal densities
  printf("# CL2 columns: j  r_N(j) (tail)  r'_N(j) (head)   [samples %lld]\n", nR);
  for (int j = 1; j <= JMAX; j++)
    if (j <= 16 || (j & (j - 1)) == 0 || j % 32 == 0) printf("J %d %.5f %.5f\n", j, rt[j] / nR, rh[j] / nR);
  printf("# CL2 jackknife errors at j = 1,2,4,...,128: j  r err  r' err\n");
  for (int b2 = 0; b2 < 8; b2++) {
    int j = 1 << b2;
    double tr = rt[j], th2 = rh[j];
    double s1 = 0, s12 = 0, s2 = 0, s22 = 0;
    double nPer = (double)nR / NBLK;
    for (int b = 0; b < NBLK; b++) {
      double e1 = (tr - bRt[b * 8 + b2]) / (nR - nPer), e2 = (th2 - bRh[b * 8 + b2]) / (nR - nPer);
      s1 += e1;
      s12 += e1 * e1;
      s2 += e2;
      s22 += e2 * e2;
    }
    double m1 = s1 / NBLK, m2 = s2 / NBLK;
    printf("JE %d %.5f %.5f %.5f %.5f\n", j, tr / nR, std::sqrt((NBLK - 1.0) / NBLK * (s12 - NBLK * m1 * m1)), th2 / nR,
           std::sqrt((NBLK - 1.0) / NBLK * (s22 - NBLK * m2 * m2)));
  }
  // CL5: span distribution
  std::vector<double> sp = spans;
  std::sort(sp.begin(), sp.end());
  double c3 = std::cbrt((double)n);
  printf("# CL5: span/N^{1/3} quantiles (N^{1/3} = %.3f), samples %zu\n", c3, sp.size());
  printf("Q");
  for (double q : {0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95})
    printf(" q%.3f=%.3f", q, sp[(size_t)(q * (sp.size() - 1))] / c3);
  printf("\n");
  for (double eta : {0.0, 0.05, 0.1, 0.15, 0.2, 0.25}) {
    double thr = std::pow((double)n, 1.0 / 3 + eta);
    size_t cnt = 0;
    for (double v : sp)
      if (v < thr) cnt++;
    printf("P eta=%.2f P(span < N^{1/3+eta}) = %.5f\n", eta, (double)cnt / sp.size());
  }
  return 0;
}
