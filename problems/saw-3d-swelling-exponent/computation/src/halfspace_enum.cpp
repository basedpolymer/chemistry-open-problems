// halfspace_enum.cpp -- exact enumeration of half-space walks on Z^3 (x1(w_i) > 0 for 1 <= i <= n), reduced by the
// 8 symmetries fixing e1 (sign changes and exchange of x2, x3).  For n <= N it counts
//   B n T c   : n-step bridges (0 = x1(w_0) < x1(w_i) <= x1(w_n), i >= 1) of height x1(w_n) = T;
//   R n T c   : n-step bridges of height 2T having a renewal point at height T, i.e. an index i with
//               x1(w_j) <= x1(w_i) = T < x1(w_k) for all j <= i < k (Madras-Slade (4.2.1));
//   A n M c   : n-step arches: x1(w_n) = 0, x1(w_i) > 0 for 0 < i < n, max_i x1(w_i) = M;
//   H n c     : n-step half-space walks (total, for cross-checking).
// An n-step arch is an (n-1)-step half-space walk ending at height 1, not at e1, followed by the step -e1 (the only
// visited site at height 0 is the origin).  Renewal test: let i* be the last index with x1 <= T; a renewal at height
// T exists iff max_{j <= i*} x1(w_j) <= T (then x1(w_{i*}) = T automatically).
// Canonical representatives: first step +e1; first step with a nonzero x2/x3 component is +e2; first step with a
// nonzero x3 component, if any, is +e3.  Weights 1 / 4 / 8 (orbit sizes under the stabiliser of e1).
// Usage: halfspace_enum N [threads<=3]     Build: g++ -O3 -march=native -std=c++17 -pthread halfspace_enum.cpp
// Runtime (3 threads, Xeon 2.8 GHz, shared machine): N = 17 56 s; N = 19 17 min (measured).
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

static int N, L, W, W3;
static int OFF[6];
static const int DX[6] = {1, -1, 0, 0, 0, 0};
static const uint64_t WB[3] = {1, 4, 8};
static const int NA[3] = {3, 5, 6};
static const int AD[3][6] = {{0, 1, 2, -1, -1, -1}, {0, 1, 2, 3, 4, -1}, {0, 1, 2, 3, 4, 5}};
static inline int newtype(int t, int d) {
  if (t == 0) return d < 2 ? 0 : 1;
  if (t == 1) return d < 4 ? 1 : 2;
  return 2;
}

struct Acc {
  std::vector<uint64_t> B, R, A, H;  // [n*(N+2) + level]
  Acc() : B((N + 1) * (N + 2), 0), R((N + 1) * (N + 2), 0), A((N + 2) * (N + 2), 0), H(N + 1, 0) {}
  void add(const Acc& o) {
    for (size_t i = 0; i < B.size(); i++) B[i] += o.B[i], R[i] += o.R[i];
    for (size_t i = 0; i < A.size(); i++) A[i] += o.A[i];
    for (size_t i = 0; i < H.size(); i++) H[i] += o.H[i];
  }
};

struct Worker {
  std::vector<uint8_t> occ;
  Acc acc;
  int hgt[64];
  int e1;
  int minRecord;
  Worker() : occ(W3, 0), minRecord(0) {}
  // record node at depth m: position pos, height y = hgt[m], prefix max mx (including m), type t
  inline void record(int m, int pos, int t, int mx) {
    if (m < minRecord) return;
    uint64_t w = WB[t];
    int y = hgt[m];
    acc.H[m] += w;
    if (y == mx) {  // bridge of height y
      acc.B[m * (N + 2) + y] += w;
      if ((y & 1) == 0) {
        int T = y / 2;
        int ist = -1;
        for (int i = m; i >= 0; i--)
          if (hgt[i] <= T) {
            ist = i;
            break;
          }
        bool ok = true;
        for (int j = 0; j <= ist; j++)
          if (hgt[j] > T) {
            ok = false;
            break;
          }
        if (ok) acc.R[m * (N + 2) + T] += w;
      }
    }
    if (y == 1 && pos != e1 && m + 1 <= N) acc.A[(m + 1) * (N + 2) + mx] += w;
  }
  void dfs(int m, int pos, int t, int mx) {
    if (m == N) return;
    for (int k = 0; k < NA[t]; k++) {
      int d = AD[t][k];
      int q = pos + OFF[d];
      int y = hgt[m] + DX[d];
      if (y < 1 || occ[q]) continue;
      int t2 = newtype(t, d);
      hgt[m + 1] = y;
      int mx2 = y > mx ? y : mx;
      record(m + 1, q, t2, mx2);
      if (m + 1 < N) {
        occ[q] = 1;
        dfs(m + 1, q, t2, mx2);
        occ[q] = 0;
      }
    }
  }
};

static std::vector<std::vector<uint8_t>> prefixes;
static void gen(Worker& w, int K, int m, int pos, int t, int mx, std::vector<uint8_t>& st) {
  if (m == K) {
    prefixes.push_back(st);
    return;
  }
  for (int k = 0; k < NA[t]; k++) {
    int d = AD[t][k];
    int q = pos + OFF[d];
    int y = w.hgt[m] + DX[d];
    if (y < 1 || w.occ[q]) continue;
    int t2 = newtype(t, d);
    w.hgt[m + 1] = y;
    int mx2 = y > mx ? y : mx;
    w.record(m + 1, q, t2, mx2);
    w.occ[q] = 1;
    st.push_back((uint8_t)d);
    gen(w, K, m + 1, q, t2, mx2, st);
    st.pop_back();
    w.occ[q] = 0;
  }
}

int main(int argc, char** argv) {
  N = argc > 1 ? atoi(argv[1]) : 10;
  int threads = argc > 2 ? atoi(argv[2]) : 1;
  if (threads > 3) threads = 3;
  if (threads < 1) threads = 1;
  if (N < 3 || N > 40) return 2;
  int K = N >= 14 ? 9 : N - 2;
  L = N + 2;
  W = 2 * L + 1;
  W3 = W * W * W;
  OFF[0] = 1;
  OFF[1] = -1;
  OFF[2] = W;
  OFF[3] = -W;
  OFF[4] = W * W;
  OFF[5] = -W * W;
  int o = L + W * L + W * W * L;
  Worker master;
  master.e1 = o + 1;
  master.hgt[0] = 0;
  master.hgt[1] = 1;
  master.occ[o] = master.occ[o + 1] = 1;
  master.acc.H[0] = 1;
  master.acc.B[0] = 1;  // the 0-step bridge, height 0
  master.acc.R[0] = 1;  // height 0 = 2*0: renewal at 0 (index 0), consistent with b_0 b_0
  master.record(1, o + 1, 0, 1);
  std::vector<uint8_t> st;
  gen(master, K, 1, o + 1, 0, 1, st);
  std::atomic<size_t> next(0);
  std::vector<Worker*> ws;
  for (int i = 0; i < threads; i++) {
    ws.push_back(new Worker());
    ws.back()->e1 = o + 1;
    ws.back()->minRecord = K + 1;
  }
  std::vector<std::thread> th;
  for (int i = 0; i < threads; i++)
    th.emplace_back([&, i]() {
      Worker* w = ws[i];
      std::vector<int> vis;
      for (;;) {
        size_t j = next.fetch_add(1);
        if (j >= prefixes.size()) break;
        int pos = o + 1, t = 0, mx = 1;
        vis.clear();
        vis.push_back(o);
        vis.push_back(o + 1);
        w->occ[o] = w->occ[o + 1] = 1;
        w->hgt[0] = 0;
        w->hgt[1] = 1;
        int m = 1;
        for (uint8_t d : prefixes[j]) {
          pos += OFF[d];
          w->hgt[m + 1] = w->hgt[m] + DX[d];
          if (w->hgt[m + 1] > mx) mx = w->hgt[m + 1];
          t = newtype(t, d);
          m++;
          w->occ[pos] = 1;
          vis.push_back(pos);
        }
        w->dfs(m, pos, t, mx);
        for (int p : vis) w->occ[p] = 0;
      }
    });
  for (auto& x : th) x.join();
  for (auto* w : ws) master.acc.add(w->acc);
  printf("# halfspace_enum N=%d threads=%d prefix_depth=%d\n", N, threads, K);
  for (int n = 0; n <= N; n++) printf("H %d %llu\n", n, (unsigned long long)master.acc.H[n]);
  for (int n = 0; n <= N; n++)
    for (int T = 0; T <= N + 1; T++) {
      uint64_t v = master.acc.B[n * (N + 2) + T];
      if (v) printf("B %d %d %llu\n", n, T, (unsigned long long)v);
    }
  for (int n = 0; n <= N; n++)
    for (int T = 0; T <= N + 1; T++) {
      uint64_t v = master.acc.R[n * (N + 2) + T];
      if (v) printf("R %d %d %llu\n", n, T, (unsigned long long)v);
    }
  for (int n = 0; n <= N; n++)
    for (int M = 0; M <= N + 1; M++) {
      uint64_t v = master.acc.A[n * (N + 2) + M];
      if (v) printf("A %d %d %llu\n", n, M, (unsigned long long)v);
    }
  return 0;
}
