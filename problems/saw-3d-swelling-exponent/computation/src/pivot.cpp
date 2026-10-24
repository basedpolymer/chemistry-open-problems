// pivot.cpp -- pivot algorithm (Lal 1969, Madras-Sokal 1988) for n-step SAWs on Z^3, hash-map implementation.
// NUMERICS ONLY: Monte Carlo estimate of E_n |omega_n|^2; it proves nothing.
//
// Move: choose a pivot site k uniformly in {1, ..., n-1} and a lattice symmetry g uniformly among the 47
// non-identity elements of B_3; apply x -> omega_k + g(x - omega_k) to the shorter of the two sub-walks (this is
// equivalent, up to a global rigid motion which leaves |omega_n - omega_0| unchanged, to pivoting the end part with
// g^{-1}).  Self-avoidance of the proposal is checked with a hash map site -> index, testing the moved sites in order
// of increasing distance (along the walk) from the pivot; accepted moves update the map.  The chain is reversible
// w.r.t. the uniform measure on SAWs modulo rigid motions and ergodic (Madras-Sokal 1988).
// Start: straight rod, burn-in of B attempts; then T attempts, recording R^2 = |omega_n - omega_0|^2 after every
// attempt.  Error bars: batch means over 2^m batches; the integrated autocorrelation time is estimated from the
// growth of the batch-mean variance.
//
// Also records the end-step correlations s_1 . s_n and s_1 . (omega_n - omega_0) (invariant under rigid motions).
// Usage: pivot n attempts burnin seed      Output: summary lines "R ..." and "S ..." (see printf at the end).
// Build: g++ -O3 -march=native -std=c++17 pivot.cpp -o pivot
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

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

struct P3 {
  int32_t x, y, z;
};
static const int64_t OFFS = 1 << 20;
static inline uint64_t key(const P3& p) {
  return ((uint64_t)(p.x + OFFS) << 42) | ((uint64_t)(p.y + OFFS) << 21) | (uint64_t)(p.z + OFFS);
}
static inline uint64_t hsh(uint64_t k) {
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  return k;
}

struct HashMap {  // open addressing, linear probing, backward-shift deletion; key 0 = empty (never a real key)
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
  inline int32_t find(uint64_t k) const {
    uint64_t i = hsh(k) & mask;
    while (K[i]) {
      if (K[i] == k) return V[i];
      i = (i + 1) & mask;
    }
    return -1;
  }
  inline void insert(uint64_t k, int32_t v) {
    uint64_t i = hsh(k) & mask;
    while (K[i] && K[i] != k) i = (i + 1) & mask;
    K[i] = k;
    V[i] = v;
  }
  inline void erase(uint64_t k) {
    uint64_t i = hsh(k) & mask;
    while (K[i] != k) {
      if (!K[i]) return;
      i = (i + 1) & mask;
    }
    uint64_t j = i;
    for (;;) {
      j = (j + 1) & mask;
      if (!K[j]) break;
      uint64_t h = hsh(K[j]) & mask;
      // move K[j] into hole i if its home h is not cyclically in (i, j]
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

// the 48 signed permutation matrices; index 0 = identity
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
static inline P3 apply(const Sym& g, const P3& v) {
  int c[3] = {v.x, v.y, v.z};
  return P3{g.sign[0] * c[g.perm[0]], g.sign[1] * c[g.perm[1]], g.sign[2] * c[g.perm[2]]};
}

int main(int argc, char** argv) {
  if (argc < 5) {
    fprintf(stderr, "usage: %s n attempts burnin seed\n", argv[0]);
    return 2;
  }
  int n = atoi(argv[1]);
  long long T = atoll(argv[2]), B = atoll(argv[3]);
  uint64_t seed = strtoull(argv[4], 0, 10);
  initSym();
  Rng rng(seed * 1000003ULL + (uint64_t)n);
  std::vector<P3> w(n + 1), nw(n + 1);
  for (int i = 0; i <= n; i++) w[i] = P3{i, 0, 0};
  HashMap H;
  H.init(n + 1);
  for (int i = 0; i <= n; i++) H.insert(key(w[i]), i);
  auto t0 = std::chrono::steady_clock::now();
  long long acc = 0, accSinceRecenter = 0;
  // batch means
  const int NB = 1024;  // number of finest batches
  long long perBatch = T / NB;
  if (perBatch < 1) perBatch = 1;
  std::vector<double> bm;
  bm.reserve(NB + 1);
  double bsum = 0;
  long long bcnt = 0;
  long double tot = 0, tot2 = 0;
  long long rec = 0;
  // end-step correlations: s_1 . s_n and s_1 . (omega_n - omega_0), with batch means for errors
  long double ts = 0, tsr = 0;
  std::vector<double> bms, bmr;
  double bss = 0, bsr = 0;
  for (long long it = 0; it < B + T; it++) {
    int k = 1 + (int)rng.below(n - 1);
    const Sym& g = SYM[1 + rng.below(47)];
    bool after = (k >= n - k);  // move the shorter side
    int lo = after ? k + 1 : 0, hi = after ? n : k - 1;  // moved indices [lo, hi]
    P3 c = w[k];
    bool ok = true;
    int len = hi - lo + 1;
    for (int t = 0; t < len; t++) {
      int j = after ? (k + 1 + t) : (k - 1 - t);
      P3 d{w[j].x - c.x, w[j].y - c.y, w[j].z - c.z};
      P3 r = apply(g, d);
      P3 q{c.x + r.x, c.y + r.y, c.z + r.z};
      int32_t o = H.find(key(q));
      if (o >= 0 && (after ? (o <= k) : (o >= k))) {
        ok = false;
        break;
      }
      nw[j] = q;
    }
    if (ok) {
      acc++;
      for (int j = lo; j <= hi; j++) H.erase(key(w[j]));
      for (int j = lo; j <= hi; j++) {
        w[j] = nw[j];
        H.insert(key(w[j]), j);
      }
      if (++accSinceRecenter >= 4096) {  // keep coordinates bounded: translate so that w[0] = 0
        accSinceRecenter = 0;
        P3 s = w[0];
        bool need = false;
        for (int i = 0; i <= n; i += (n / 64 + 1))
          if (std::abs(w[i].x) > (1 << 18) || std::abs(w[i].y) > (1 << 18) || std::abs(w[i].z) > (1 << 18)) need = true;
        if (need) {
          H.init(n + 1);
          for (int i = 0; i <= n; i++) {
            w[i] = P3{w[i].x - s.x, w[i].y - s.y, w[i].z - s.z};
            H.insert(key(w[i]), i);
          }
        }
      }
    }
    if (it >= B) {
      double dx = w[n].x - w[0].x, dy = w[n].y - w[0].y, dz = w[n].z - w[0].z;
      double r2 = dx * dx + dy * dy + dz * dz;
      tot += r2;
      tot2 += (long double)r2 * r2;
      rec++;
      bsum += r2;
      int s1x = w[1].x - w[0].x, s1y = w[1].y - w[0].y, s1z = w[1].z - w[0].z;
      int snx = w[n].x - w[n - 1].x, sny = w[n].y - w[n - 1].y, snz = w[n].z - w[n - 1].z;
      double d11 = s1x * snx + s1y * sny + s1z * snz;
      double d1r = s1x * dx + s1y * dy + s1z * dz;
      ts += d11;
      tsr += d1r;
      bss += d11;
      bsr += d1r;
      if (++bcnt == perBatch) {
        bm.push_back(bsum / perBatch);
        bms.push_back(bss / perBatch);
        bmr.push_back(bsr / perBatch);
        bsum = 0;
        bss = 0;
        bsr = 0;
        bcnt = 0;
      }
    }
  }
  // final self-avoidance check (sanity)
  {
    HashMap C;
    C.init(n + 1);
    for (int i = 0; i <= n; i++) {
      if (C.find(key(w[i])) >= 0) {
        fprintf(stderr, "ERROR: final walk not self-avoiding\n");
        return 1;
      }
      C.insert(key(w[i]), i);
      if (i > 0) {
        int dd = std::abs(w[i].x - w[i - 1].x) + std::abs(w[i].y - w[i - 1].y) + std::abs(w[i].z - w[i - 1].z);
        if (dd != 1) {
          fprintf(stderr, "ERROR: final walk not nearest-neighbour\n");
          return 1;
        }
      }
    }
  }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  double mean = (double)(tot / rec);
  double var = (double)(tot2 / rec) - mean * mean;
  // batch-mean errors at several batch sizes (blocking)
  std::vector<double> cur = bm;
  double bestErr = 0;
  printf("# n=%d attempts=%lld burnin=%lld seed=%llu blocking:", n, T, B, (unsigned long long)seed);
  long long bs = perBatch;
  while (cur.size() >= 16) {
    double m = 0, v = 0;
    for (double x : cur) m += x;
    m /= cur.size();
    for (double x : cur) v += (x - m) * (x - m);
    v /= (cur.size() - 1);
    double err = std::sqrt(v / cur.size());
    printf(" [%lld:%.4g]", bs, err);
    if (err > bestErr) bestErr = err;
    std::vector<double> nx;
    for (size_t i = 0; i + 1 < cur.size(); i += 2) nx.push_back(0.5 * (cur[i] + cur[i + 1]));
    cur.swap(nx);
    bs *= 2;
  }
  printf("\n");
  double naiveErr = std::sqrt(var / rec);
  double tauint = 0.5 * (bestErr / naiveErr) * (bestErr / naiveErr);
  // blocking error for the end-step correlations (largest over levels)
  auto blockErr = [](std::vector<double> v) {
    double best = 0;
    while (v.size() >= 16) {
      double m = 0, q = 0;
      for (double x : v) m += x;
      m /= v.size();
      for (double x : v) q += (x - m) * (x - m);
      q /= (v.size() - 1);
      double e = std::sqrt(q / v.size());
      if (e > best) best = e;
      std::vector<double> nx;
      for (size_t i = 0; i + 1 < v.size(); i += 2) nx.push_back(0.5 * (v[i] + v[i + 1]));
      v.swap(nx);
    }
    return best;
  };
  // n attempts acc_frac mean_R2 err(blocking max) tau_int(attempts) var/mean^2 seconds
  printf("R %d %lld %.6f %.8g %.4g %.3g %.5f %.1f\n", n, T, (double)acc / (B + T), mean, bestErr, tauint,
         var / (mean * mean), secs);
  // n  E[s_1 . s_n] err  E[s_1 . (omega_n - omega_0)] err
  printf("S %d %.6e %.2e %.6e %.2e\n", n, (double)(ts / rec), blockErr(bms), (double)(tsr / rec), blockErr(bmr));
  return 0;
}
