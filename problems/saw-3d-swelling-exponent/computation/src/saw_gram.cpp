// saw_gram.cpp -- exact Gram sums G_n(i,j) = sum over n-step SAWs on Z^3 of  omega_i . omega_j  (1 <= i <= j <= n),
// for every n <= N.  Step correlations follow:  sum_omega s_i . s_j = G(i,j) - G(i-1,j) - G(i,j-1) + G(i-1,j-1),
// s_i = omega_i - omega_{i-1}, G(0,.) = 0.
//
// Symmetry reduction as in saw_enum.cpp: canonical walks (first step +e1, first non-e1 step +e2, first step out of
// the (e1,e2)-plane +e3) with orbit weights 6 / 24 / 48.  omega_i . omega_j is invariant under the 48 lattice
// symmetries, so  G_n(i,j) = sum over canonical j-step walks u of (u_i . u_j) * D_n(u),  where D_n(u) is the sum of
// orbit weights of the canonical n-step extensions of u.  Leaves are handled at their parent through the weighted
// vector sum of the children's endpoints.
//
// Usage: saw_gram N [threads<=3] > out.txt      Output:  "T n c_n"  and  "G n i j value"  (value: signed integer)
// Build: g++ -O3 -march=native -std=c++17 -pthread saw_gram.cpp -o saw_gram
// Runtime (3 threads, Xeon 2.8 GHz): N = 16 ~ 1 min, N = 17 ~ 4 min, N = 18 ~ 20 min (estimates).
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

static int N, L, W, W3, PK;
static int OFF[6];
static const int DXX[6] = {1, -1, 0, 0, 0, 0};
static const int DYY[6] = {0, 0, 1, -1, 0, 0};
static const int DZZ[6] = {0, 0, 0, 0, 1, -1};
static const int64_t WT[3] = {6, 24, 48};
static const int NA[3] = {2, 5, 6};
static const int AD[3][6] = {{0, 2, -1, -1, -1, -1}, {0, 1, 2, 3, 4, -1}, {0, 1, 2, 3, 4, 5}};
static inline int newtype(int t, int d) {
  if (t == 0) return d == 0 ? 0 : 1;
  if (t == 1) return d < 4 ? 1 : 2;
  return 2;
}
static const int MAXD = 40;

struct Worker {
  std::vector<uint8_t> occ;
  std::vector<int64_t> G;  // index (n*(N+1)+i)*(N+1)+j
  int X[MAXD + 1], Y[MAXD + 1], Z[MAXD + 1];
  Worker() {
    occ.assign(W3, 0);
    G.assign((size_t)(N + 1) * (N + 1) * (N + 1), 0);
  }
  inline int64_t& g(int n, int i, int j) { return G[((size_t)n * (N + 1) + i) * (N + 1) + j]; }
  // contributions of node at depth m with descendant counts D[0..N-m]
  inline void addNode(int m, const int64_t* D, int len) {
    for (int i = 1; i <= m; i++) {
      int64_t d = (int64_t)X[i] * X[m] + (int64_t)Y[i] * Y[m] + (int64_t)Z[i] * Z[m];
      for (int j = 0; j <= len; j++) g(m + j, i, m) += d * D[j];
    }
  }
  template <int TY>
  void dfs(int m, int pos, int64_t* Dp) {
    if (m == N - 1) {
      int64_t d1 = 0, vx = 0, vy = 0, vz = 0, s2 = 0;
      for (int k = 0; k < NA[TY]; k++) {
        int d = AD[TY][k];
        int q = pos + OFF[d];
        if (!occ[q]) {
          int64_t w = WT[newtype(TY, d)];
          int x = X[m] + DXX[d], y = Y[m] + DYY[d], z = Z[m] + DZZ[d];
          d1 += w;
          vx += w * x;
          vy += w * y;
          vz += w * z;
          s2 += w * (int64_t)(x * x + y * y + z * z);
        }
      }
      if (m > PK) {
        int64_t D[2] = {WT[TY], d1};
        addNode(m, D, 1);
      }
      for (int i = 1; i <= m; i++) g(N, i, N) += X[i] * vx + Y[i] * vy + Z[i] * vz;
      g(N, N, N) += s2;
      Dp[1] += WT[TY];
      Dp[2] += d1;
      return;
    }
    int64_t D[MAXD];
    int len = N - m;
    for (int j = 0; j <= len; j++) D[j] = 0;
    D[0] = WT[TY];
    for (int k = 0; k < NA[TY]; k++) {
      int d = AD[TY][k];
      int q = pos + OFF[d];
      if (occ[q]) continue;
      int t2 = newtype(TY, d);
      X[m + 1] = X[m] + DXX[d];
      Y[m + 1] = Y[m] + DYY[d];
      Z[m + 1] = Z[m] + DZZ[d];
      occ[q] = 1;
      if (t2 == 0)
        dfs<0>(m + 1, q, D);
      else if (t2 == 1)
        dfs<1>(m + 1, q, D);
      else
        dfs<2>(m + 1, q, D);
      occ[q] = 0;
    }
    if (m > PK) addNode(m, D, len);
    for (int j = 0; j <= len; j++) Dp[j + 1] += D[j];
  }
  void start(int m, int pos, int t, int64_t* Dp) {
    if (t == 0)
      dfs<0>(m, pos, Dp);
    else if (t == 1)
      dfs<1>(m, pos, Dp);
    else
      dfs<2>(m, pos, Dp);
  }
};

static inline int idx(int x, int y, int z) { return (x + L) + W * ((y + L) + W * (z + L)); }
static std::vector<std::vector<uint8_t>> prefixes;

// master: canonical walks of length <= K, all their G contributions for n <= K
static void genPrefixes(Worker& w, int K, int m, int pos, int t, std::vector<uint8_t>& st, int64_t* Dp) {
  int64_t D[MAXD];
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
      int t2 = newtype(t, d);
      w.X[m + 1] = w.X[m] + DXX[d];
      w.Y[m + 1] = w.Y[m] + DYY[d];
      w.Z[m + 1] = w.Z[m] + DZZ[d];
      w.occ[q] = 1;
      st.push_back((uint8_t)d);
      genPrefixes(w, K, m + 1, q, t2, st, D);
      st.pop_back();
      w.occ[q] = 0;
    }
  }
  w.addNode(m, D, len);
  for (int j = 0; j <= len; j++) Dp[j + 1] += D[j];
}

static std::atomic<size_t> nextTask(0);
static void runWorker(Worker* w, int K) {
  std::vector<int> visited;
  for (;;) {
    size_t ti = nextTask.fetch_add(1);
    if (ti >= prefixes.size()) break;
    const std::vector<uint8_t>& st = prefixes[ti];
    int pos = idx(0, 0, 0), t = 0;
    visited.clear();
    w->occ[pos] = 1;
    visited.push_back(pos);
    w->X[0] = w->Y[0] = w->Z[0] = 0;
    int m = 0;
    for (uint8_t d : st) {
      pos += OFF[d];
      t = newtype(t, d);
      w->X[m + 1] = w->X[m] + DXX[d];
      w->Y[m + 1] = w->Y[m] + DYY[d];
      w->Z[m + 1] = w->Z[m] + DZZ[d];
      m++;
      w->occ[pos] = 1;
      visited.push_back(pos);
    }
    int64_t Dr[MAXD + 1] = {0};
    w->start(K, pos, t, Dr);  // Dr[1+j]: weighted descendants of the root at depth K+j
    // contributions of the prefix nodes (depth k <= K) for n > K
    for (int k = 1; k <= K; k++)
      for (int i = 1; i <= k; i++) {
        int64_t d = (int64_t)w->X[i] * w->X[k] + (int64_t)w->Y[i] * w->Y[k] + (int64_t)w->Z[i] * w->Z[k];
        for (int n = K + 1; n <= N; n++) w->g(n, i, k) += d * Dr[1 + n - K];
      }
    for (int p : visited) w->occ[p] = 0;
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s N [threads]\n", argv[0]);
    return 2;
  }
  N = atoi(argv[1]);
  if (N < 3 || N > 30) {
    fprintf(stderr, "need 3 <= N <= 30\n");
    return 2;
  }
  int threads = argc > 2 ? atoi(argv[2]) : 1;
  if (threads < 1) threads = 1;
  if (threads > 3) threads = 3;
  int K = argc > 3 ? atoi(argv[3]) : (N >= 14 ? 10 : N - 2);
  if (K < 1) K = 1;
  if (K > N - 1) K = N - 1;
  PK = K;
  L = N + 2;
  W = 2 * L + 1;
  W3 = W * W * W;
  OFF[0] = 1;
  OFF[1] = -1;
  OFF[2] = W;
  OFF[3] = -W;
  OFF[4] = W * W;
  OFF[5] = -W * W;
  auto t0 = std::chrono::steady_clock::now();
  Worker master;
  int o = idx(0, 0, 0), p1 = o + OFF[0];
  master.occ[o] = master.occ[p1] = 1;
  master.X[0] = master.Y[0] = master.Z[0] = 0;
  master.X[1] = 1;
  master.Y[1] = master.Z[1] = 0;
  std::vector<uint8_t> st(1, 0);
  int64_t Dm[MAXD + 1] = {0};
  genPrefixes(master, K, 1, p1, 0, st, Dm);
  master.occ[o] = master.occ[p1] = 0;
  std::vector<Worker*> ws;
  std::vector<std::thread> th;
  for (int i = 0; i < threads; i++) ws.push_back(new Worker());
  for (int i = 0; i < threads; i++) th.emplace_back(runWorker, ws[i], K);
  for (auto& x : th) x.join();
  for (auto* w : ws)
    for (size_t i = 0; i < master.G.size(); i++) master.G[i] += w->G[i];
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  printf("# saw_gram N=%d threads=%d prefix_depth=%d time_s=%.2f\n", N, threads, K, secs);
  // c_n from the diagonal identity sum s_i.s_i = c_n: c_n = G(1,1) since |omega_1| = 1
  for (int n = 1; n <= N; n++) printf("T %d %lld\n", n, (long long)master.g(n, 1, 1));
  for (int n = 1; n <= N; n++)
    for (int i = 1; i <= n; i++)
      for (int j = i; j <= n; j++) printf("G %d %d %d %lld\n", n, i, j, (long long)master.g(n, i, j));
  fprintf(stderr, "done in %.2f s\n", secs);
  return 0;
}
