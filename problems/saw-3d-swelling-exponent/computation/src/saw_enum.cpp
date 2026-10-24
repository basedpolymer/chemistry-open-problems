// saw_enum.cpp -- exact enumeration of self-avoiding walks (SAWs) on Z^3, reduced by the 48 lattice symmetries.
//
// For every n <= N it computes, exactly (unsigned 64-bit counters per lattice point, 128-bit totals):
//   * the endpoint distribution c_n(x), stored per orbit of the hyperoctahedral group B_3 (48 elements);
//   * h_n, the number of half-space walks: x_1(omega_i) > 0 for 1 <= i <= n (h_0 = 1);
//   * b_n, the number of bridges, Hammersley-Welsh / Madras-Slade convention
//         0 = x_1(omega_0) < x_1(omega_i) <= x_1(omega_n)   for 1 <= i <= n     (b_0 = 1);
//     and the endpoint distributions of half-space walks and bridges, stored per orbit of the stabiliser of e_1
//     (8 elements: permutations and sign changes of x_2, x_3);
//   * for 1 <= k < n <= N, the distribution of the intermediate point omega_k under the n-step walks (per B_3-orbit),
//     computed from weighted descendant counts (|omega_k| is B_3-invariant, so canonical prefixes suffice).
//
// Symmetry reduction.  Every walk is mapped to a unique canonical representative: first step +e1; the first step
// not parallel to e1 is +e2; the first step not in the (e1,e2)-plane is +e3.  Canonical walks are of three types:
// straight (orbit size 6), planar (orbit size 24), genuinely 3D (orbit size 48).  A canonical walk w contributes its
// orbit size |G w| to the orbit class of its endpoint; c_n(x) = N_n(O(x)) / |O(x)| (the division is checked to be
// exact).  For half-space walks and bridges (first step necessarily +e1), the walks of G w with first step +e1 are
// the |G w|/6 images under Stab(e1); w contributes |G w|/6 in {1,4,8}.
//
// Usage:  saw_enum N [threads=1] [prefix_depth=auto]  > out.txt
// Output lines (all counts exact integers):
//   T n c_n h_n b_n                total counts
//   A n a b c Nsum  orbit          all walks, endpoint orbit (a>=b>=c>=0), Nsum = #walks ending in the orbit
//   H n x1 b c Nsum orbit          half-space walks, endpoint class (x1, b>=c>=0) under Stab(e1)
//   B n x1 b c Nsum orbit          bridges, same classes
//   M k n a b c Nsum orbit         n-step walks whose k-th point lies in the orbit (a,b,c), 1 <= k < n
// Build: g++ -O3 -march=native -std=c++17 -pthread saw_enum.cpp -o saw_enum
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <thread>
#include <tuple>
#include <vector>

typedef unsigned __int128 u128;
static int N, L, W, W3;
static int OFF[6];
static const int DX[6] = {1, -1, 0, 0, 0, 0};
static const int DIRY[6] = {0, 0, 1, -1, 0, 0};
static const int DIRZ[6] = {0, 0, 0, 0, 1, -1};
static const uint64_t WT[3] = {6, 24, 48};
static const uint64_t WB[3] = {1, 4, 8};
static int NC;                       // number of B_3-orbit classes with |x|_1 <= N
static std::vector<int> CLS;         // box index -> class id
static std::vector<std::tuple<int, int, int>> CLASSREP;
static int PK;                       // prefix depth K (global, used by workers)

static inline int idx(int x, int y, int z) { return (x + L) + W * ((y + L) + W * (z + L)); }

// allowed directions and resulting types
// type 0 (straight along +e1 so far): +x -> 0, +y -> 1
// type 1 (planar in e1,e2): +-x, +-y -> 1, +z -> 2
// type 2: all six -> 2
static inline int newtype(int t, int d) {
  if (t == 0) return d == 0 ? 0 : 1;
  if (t == 1) return d < 4 ? 1 : 2;
  return 2;
}
static const int NA[3] = {2, 5, 6};
static const int AD[3][6] = {{0, 2, -1, -1, -1, -1}, {0, 1, 2, 3, 4, -1}, {0, 1, 2, 3, 4, 5}};

struct Counts {
  std::vector<std::vector<uint64_t>> all, half, bridge;  // [n][box index]
  std::vector<uint64_t> tab;                              // [(k*(N+1)+n)*NC + class]
  void init() {
    tab.assign((size_t)(N + 1) * (N + 1) * NC, 0);
    all.assign(N + 1, std::vector<uint64_t>());
    half.assign(N + 1, std::vector<uint64_t>());
    bridge.assign(N + 1, std::vector<uint64_t>());
    for (int m = 0; m <= N; m++) {
      all[m].assign(W3, 0);
      half[m].assign(W3, 0);
      bridge[m].assign(W3, 0);
    }
  }
  void add(const Counts& o) {
    for (size_t i = 0; i < tab.size(); i++) tab[i] += o.tab[i];
    for (int m = 0; m <= N; m++)
      for (int i = 0; i < W3; i++) {
        all[m][i] += o.all[m][i];
        half[m][i] += o.half[m][i];
        bridge[m][i] += o.bridge[m][i];
      }
  }
};

// Descendant counts.  For a canonical node u at depth m, D[j] (0 <= j <= N-m) is the sum of orbit weights WT[type]
// of the canonical walks of length m+j extending u (D[0] = weight of u itself).  The number of (m+j)-step SAWs
// whose m-th point lies in the B_3-orbit O equals the sum of D[j] over canonical m-step u with u_m in O.
static const int MAXD = 40;

struct Worker {
  std::vector<uint8_t> occ;
  Counts C;
  uint64_t** hA;
  uint64_t** hH;
  uint64_t** hB;
  uint64_t* T;         // tab
  int minDepthRecord;  // record nodes (A/H/B) only at depth >= this
  Worker() {
    occ.assign(W3, 0);
    C.init();
    hA = new uint64_t*[N + 1];
    hH = new uint64_t*[N + 1];
    hB = new uint64_t*[N + 1];
    for (int m = 0; m <= N; m++) {
      hA[m] = C.all[m].data();
      hH[m] = C.half[m].data();
      hB[m] = C.bridge[m].data();
    }
    T = C.tab.data();
    minDepthRecord = 0;
  }
  inline uint64_t& tab(int k, int n, int cls) { return T[((size_t)k * (N + 1) + n) * NC + cls]; }
  inline void record(int m, int pos, int t, int x1, int maxx_before, bool hs) {
    if (m < minDepthRecord) return;
    hA[m][pos] += WT[t];
    if (hs) {
      hH[m][pos] += WB[t];
      if (x1 >= maxx_before) hB[m][pos] += WB[t];
    }
  }
  // node at depth m (already recorded), PK <= m < N.  Adds its descendant counts, shifted by one, into Dp.
  // leaf parent: node at depth N-1
  template <int TY>
  __attribute__((always_inline)) inline void leafParent(int m, int pos, int x1, int maxx, bool hs, uint64_t* Dp) {
    uint64_t* h = hA[N];
    uint64_t d1 = 0;
    if (!hs) {
      for (int k = 0; k < NA[TY]; k++) {
        int d = AD[TY][k];
        int q = pos + OFF[d];
        if (!occ[q]) {
          uint64_t w = WT[newtype(TY, d)];
          h[q] += w;
          d1 += w;
        }
      }
    } else {
      for (int k = 0; k < NA[TY]; k++) {
        int d = AD[TY][k];
        int q = pos + OFF[d];
        if (!occ[q]) {
          int t2 = newtype(TY, d);
          h[q] += WT[t2];
          d1 += WT[t2];
          int x2 = x1 + DX[d];
          if (x2 > 0) {
            hH[N][q] += WB[t2];
            if (x2 >= maxx) hB[N][q] += WB[t2];
          }
        }
      }
    }
    if (m > PK) tab(m, N, CLS[pos]) += d1;
    Dp[1] += WT[TY];
    Dp[2] += d1;
  }
  // node at depth m (already recorded), PK <= m < N.  Adds its descendant counts, shifted by one, into Dp.
  template <int TY>
  void dfs(int m, int pos, int x1, int maxx, bool hs, uint64_t* Dp) {
    if (m == N - 1) {
      leafParent<TY>(m, pos, x1, maxx, hs, Dp);
      return;
    }
    if (m == N - 2) {
      uint64_t D[3] = {WT[TY], 0, 0};
      for (int k = 0; k < NA[TY]; k++) {
        int d = AD[TY][k];
        int q = pos + OFF[d];
        if (occ[q]) continue;
        int x2 = x1 + DX[d];
        bool hs2 = hs && (x2 > 0);
        int t2 = newtype(TY, d);
        record(m + 1, q, t2, x2, maxx, hs2);
        int mx2 = x2 > maxx ? x2 : maxx;
        occ[q] = 1;
        if (t2 == 0)
          leafParent<0>(m + 1, q, x2, mx2, hs2, D);
        else if (t2 == 1)
          leafParent<1>(m + 1, q, x2, mx2, hs2, D);
        else
          leafParent<2>(m + 1, q, x2, mx2, hs2, D);
        occ[q] = 0;
      }
      if (m > PK) {
        int c = CLS[pos];
        tab(m, m + 1, c) += D[1];
        tab(m, m + 2, c) += D[2];
      }
      Dp[1] += D[0];
      Dp[2] += D[1];
      Dp[3] += D[2];
      return;
    }
    uint64_t D[MAXD];
    int len = N - m;
    for (int j = 0; j <= len; j++) D[j] = 0;
    D[0] = WT[TY];
    for (int k = 0; k < NA[TY]; k++) {
      int d = AD[TY][k];
      int q = pos + OFF[d];
      if (occ[q]) continue;
      int x2 = x1 + DX[d];
      bool hs2 = hs && (x2 > 0);
      int t2 = newtype(TY, d);
      record(m + 1, q, t2, x2, maxx, hs2);
      int mx2 = x2 > maxx ? x2 : maxx;
      occ[q] = 1;
      if (t2 == 0)
        dfs<0>(m + 1, q, x2, mx2, hs2, D);
      else if (t2 == 1)
        dfs<1>(m + 1, q, x2, mx2, hs2, D);
      else
        dfs<2>(m + 1, q, x2, mx2, hs2, D);
      occ[q] = 0;
    }
    if (m > PK) {
      int c = CLS[pos];
      for (int j = 1; j <= len; j++) tab(m, m + j, c) += D[j];
    }
    for (int j = 0; j <= len; j++) Dp[j + 1] += D[j];
  }
  void start(int m, int pos, int t, int x1, int maxx, bool hs, uint64_t* Dp) {
    if (t == 0)
      dfs<0>(m, pos, x1, maxx, hs, Dp);
    else if (t == 1)
      dfs<1>(m, pos, x1, maxx, hs, Dp);
    else
      dfs<2>(m, pos, x1, maxx, hs, Dp);
  }
};

// prefix generation: canonical walks of length exactly K; the master records all nodes of depth <= K
static std::vector<std::vector<uint8_t>> prefixes;
static void genPrefixes(Worker& w, int K, int m, int pos, int t, int x1, int maxx, bool hs, std::vector<uint8_t>& st,
                        uint64_t* Dp) {
  uint64_t D[MAXD];
  int len = K - m;
  for (int j = 0; j <= len; j++) D[j] = 0;
  D[0] = WT[t];
  if (m == K) {
    prefixes.push_back(st);
  } else {
    for (int k = 0; k < NA[t]; k++) {
      int d = AD[t][k];
      int q = pos + OFF[d];
      if (w.occ[q]) continue;
      int x2 = x1 + DX[d];
      bool hs2 = hs && (x2 > 0);
      int t2 = newtype(t, d);
      w.record(m + 1, q, t2, x2, maxx, hs2);
      int mx2 = x2 > maxx ? x2 : maxx;
      w.occ[q] = 1;
      st.push_back((uint8_t)d);
      genPrefixes(w, K, m + 1, q, t2, x2, mx2, hs2, st, D);
      st.pop_back();
      w.occ[q] = 0;
    }
  }
  int c = CLS[pos];
  for (int j = 1; j <= len; j++) w.tab(m, m + j, c) += D[j];
  for (int j = 0; j <= len; j++) Dp[j + 1] += D[j];
}

static std::atomic<size_t> nextTask(0);
static void runWorker(Worker* w, int K) {
  w->minDepthRecord = K + 1;
  std::vector<int> visited;
  visited.reserve(K + 1);
  for (;;) {
    size_t i = nextTask.fetch_add(1);
    if (i >= prefixes.size()) break;
    const std::vector<uint8_t>& st = prefixes[i];
    // replay
    int pos = idx(0, 0, 0), t = 0, x1 = 0, maxx = 0;
    bool hs = true;
    visited.clear();
    w->occ[pos] = 1;
    visited.push_back(pos);
    for (uint8_t d : st) {
      pos += OFF[d];
      x1 += DX[d];
      hs = hs && (x1 > 0);
      if (x1 > maxx) maxx = x1;
      t = newtype(t, d);
      w->occ[pos] = 1;
      visited.push_back(pos);
    }
    uint64_t Dr[MAXD + 1];
    for (int j = 0; j <= MAXD; j++) Dr[j] = 0;
    w->start(K, pos, t, x1, maxx, hs, Dr);  // Dr[1+j] = weighted descendants of the root at depth K+j
    for (int k = 1; k <= K; k++) {
      int c = CLS[visited[k]];
      for (int n = K + 1; n <= N; n++) w->tab(k, n, c) += Dr[1 + n - K];
    }
    for (int p : visited) w->occ[p] = 0;
  }
}

static void print128(u128 v) {
  char buf[64];
  int i = 63;
  buf[i] = 0;
  if (v == 0) buf[--i] = '0';
  while (v > 0) {
    buf[--i] = char('0' + (int)(v % 10));
    v /= 10;
  }
  fputs(buf + i, stdout);
}
static int orbitB3(int a, int b, int c) {
  int perms = (a == b && b == c) ? 1 : ((a == b || b == c) ? 3 : 6);
  int nz = (a > 0) + (b > 0) + (c > 0);
  return perms << nz;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s N [threads] [prefix_depth]\n", argv[0]);
    return 2;
  }
  N = atoi(argv[1]);
  if (N < 2 || N > MAXD - 4) {
    fprintf(stderr, "need 2 <= N <= %d\n", MAXD - 4);
    return 2;
  }
  int threads = argc > 2 ? atoi(argv[2]) : 1;
  int K = argc > 3 ? atoi(argv[3]) : (N >= 14 ? 10 : N - 2);
  if (K > N) K = N;
  if (K < 1) K = 1;
  PK = K;
  if (threads < 1) threads = 1;
  if (threads > 3) {
    fprintf(stderr, "resource cap: at most 3 threads\n");
    threads = 3;
  }
  L = N + 2;
  W = 2 * L + 1;
  W3 = W * W * W;
  OFF[0] = 1;
  OFF[1] = -1;
  OFF[2] = W;
  OFF[3] = -W;
  OFF[4] = W * W;
  OFF[5] = -W * W;
  // classes
  {
    std::map<std::tuple<int, int, int>, int> id;
    for (int a = 0; a <= N; a++)
      for (int b = 0; b <= a; b++)
        for (int c = 0; c <= b; c++)
          if (a + b + c <= N) {
            int k = (int)CLASSREP.size();
            id[std::make_tuple(a, b, c)] = k;
            CLASSREP.push_back(std::make_tuple(a, b, c));
          }
    NC = (int)CLASSREP.size();
    CLS.assign(W3, 0);
    for (int z = -L; z <= L; z++)
      for (int y = -L; y <= L; y++)
        for (int x = -L; x <= L; x++) {
          int a = abs(x), b = abs(y), c = abs(z);
          if (a < b) std::swap(a, b);
          if (b < c) std::swap(b, c);
          if (a < b) std::swap(a, b);
          auto it = id.find(std::make_tuple(a, b, c));
          CLS[idx(x, y, z)] = it == id.end() ? 0 : it->second;
        }
  }
  auto t0 = std::chrono::steady_clock::now();

  Worker master;
  int o = idx(0, 0, 0);
  master.occ[o] = 1;
  master.C.all[0][o] += 1;  // the 0-step walk (orbit size 1)
  master.C.half[0][o] += 1;
  master.C.bridge[0][o] += 1;
  // depth-1 node: step +e1, type 0, weight 6
  int p1 = o + OFF[0];
  master.record(1, p1, 0, 1, 0, true);
  master.occ[p1] = 1;
  std::vector<uint8_t> st;
  st.push_back(0);
  uint64_t Dm[MAXD + 1] = {0};
  genPrefixes(master, K, 1, p1, 0, 1, 1, true, st, Dm);
  master.occ[p1] = 0;
  master.occ[o] = 0;
  fprintf(stderr, "prefix depth K=%d, %zu prefixes, threads=%d\n", K, prefixes.size(), threads);

  std::vector<Worker*> ws;
  std::vector<std::thread> th;
  if (K < N) {
    for (int i = 0; i < threads; i++) ws.push_back(new Worker());
    for (int i = 0; i < threads; i++) th.emplace_back(runWorker, ws[i], K);
    for (auto& x : th) x.join();
    for (auto* w : ws) master.C.add(w->C);
  }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  printf("# saw_enum N=%d threads=%d prefix_depth=%d prefixes=%zu time_s=%.2f\n", N, threads, K, prefixes.size(),
         secs);

  // aggregate into classes and print
  bool ok = true;
  for (int m = 0; m <= N; m++) {
    std::map<std::tuple<int, int, int>, u128> A, H, B;
    u128 cn = 0, hn = 0, bn = 0;
    for (int z = -m; z <= m; z++)
      for (int y = -m; y <= m; y++)
        for (int x = -m; x <= m; x++) {
          int i = idx(x, y, z);
          uint64_t va = master.C.all[m][i], vh = master.C.half[m][i], vb = master.C.bridge[m][i];
          if (va) {
            int a = abs(x), b = abs(y), c = abs(z);
            if (a < b) std::swap(a, b);
            if (b < c) std::swap(b, c);
            if (a < b) std::swap(a, b);
            A[std::make_tuple(a, b, c)] += va;
            cn += va;
          }
          if (vh || vb) {
            int b = abs(y), c = abs(z);
            if (b < c) std::swap(b, c);
            if (vh) H[std::make_tuple(x, b, c)] += vh, hn += vh;
            if (vb) B[std::make_tuple(x, b, c)] += vb, bn += vb;
          }
        }
    printf("T %d ", m);
    print128(cn);
    printf(" ");
    print128(hn);
    printf(" ");
    print128(bn);
    printf("\n");
    for (auto& kv : A) {
      int a, b, c;
      std::tie(a, b, c) = kv.first;
      int orb = orbitB3(a, b, c);
      if (kv.second % orb != 0) ok = false;
      printf("A %d %d %d %d ", m, a, b, c);
      print128(kv.second);
      printf(" %d\n", orb);
    }
    for (int pass = 0; pass < 2; pass++) {
      auto& M = pass == 0 ? H : B;
      for (auto& kv : M) {
        int x, b, c;
        std::tie(x, b, c) = kv.first;
        int orb = (b == c ? 1 : 2) * (b > 0 ? 2 : 1) * (c > 0 ? 2 : 1);
        if (kv.second % orb != 0) ok = false;
        printf("%c %d %d %d %d ", pass == 0 ? 'H' : 'B', m, x, b, c);
        print128(kv.second);
        printf(" %d\n", orb);
      }
    }
  }
  // intermediate-point distributions; consistency: for each (k,n) the total must equal c_n
  bool okM = true;
  {
    std::vector<u128> cn(N + 1, 0);
    for (int m = 0; m <= N; m++)
      for (int i = 0; i < W3; i++) cn[m] += master.C.all[m][i];
    for (int k = 1; k <= N; k++)
      for (int n = k + 1; n <= N; n++) {
        u128 tot = 0;
        for (int c = 0; c < NC; c++) {
          uint64_t v = master.tab(k, n, c);
          if (!v) continue;
          tot += v;
          int a, b, cc;
          std::tie(a, b, cc) = CLASSREP[c];
          printf("M %d %d %d %d %d %llu %d\n", k, n, a, b, cc, (unsigned long long)v, orbitB3(a, b, cc));
        }
        if (tot != cn[n]) okM = false;
      }
  }
  printf("# orbit divisibility check: %s\n", ok ? "PASS" : "FAIL");
  printf("# intermediate-point totals equal c_n: %s\n", okM ? "PASS" : "FAIL");
  fprintf(stderr, "done in %.2f s, divisibility %s, M-totals %s\n", secs, ok ? "PASS" : "FAIL", okM ? "PASS" : "FAIL");
  return (ok && okM) ? 0 : 1;
}
