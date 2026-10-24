// end_steps.cpp -- independent check of the first-step / last-step correlation of SAWs on Z^3.
// Enumerates every n-step SAW with first step +e1 (no other symmetry reduction; plain DFS with an occupancy array),
// and counts them by the direction of the last step.  For every n <= N prints
//   E n  N(+x) N(-x) N(+y) N(-y) N(+z) N(-z)
// Since the uniform measure is invariant under the 48 lattice symmetries,
//   sum over all n-step SAWs of s_1 . s_n = 6 (N(+x) - N(-x)),  and  c_n = 6 sum_d N(d).
// Parallelised over the 25 two-step continuations of the first step (threads <= 3).
// Usage: end_steps N [threads]     Build: g++ -O3 -march=native -std=c++17 -pthread end_steps.cpp -o end_steps
// Runtime (3 threads, Xeon 2.8 GHz): N = 16 about 1-2 min; N = 18 about 25-40 min.
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

static int N, W, W3;
static int OFF[6];

struct Job {
  int d2, d3;
};

struct Counter {
  std::vector<uint8_t> occ;
  std::vector<std::vector<uint64_t>> cnt;  // [n][last direction]
  Counter() : occ(W3, 0), cnt(N + 1, std::vector<uint64_t>(6, 0)) {}
  void dfs(int m, int pos) {  // walk of length m ends at pos
    for (int d = 0; d < 6; d++) {
      int q = pos + OFF[d];
      if (occ[q]) continue;
      cnt[m + 1][d]++;
      if (m + 1 < N) {
        occ[q] = 1;
        dfs(m + 1, q);
        occ[q] = 0;
      }
    }
  }
};

int main(int argc, char** argv) {
  N = argc > 1 ? atoi(argv[1]) : 10;
  int threads = argc > 2 ? atoi(argv[2]) : 1;
  if (threads > 3) threads = 3;
  if (N < 3) return 2;
  int L = N + 2;
  W = 2 * L + 1;
  W3 = W * W * W;
  OFF[0] = 1;
  OFF[1] = -1;
  OFF[2] = W;
  OFF[3] = -W;
  OFF[4] = W * W;
  OFF[5] = -W * W;
  int o = L + W * L + W * W * L;
  // jobs: steps 2 and 3 (step 1 is +x = direction 0)
  std::vector<Job> jobs;
  std::vector<std::vector<uint64_t>> head(N + 1, std::vector<uint64_t>(6, 0));
  head[1][0] = 1;
  for (int d2 = 0; d2 < 6; d2++) {
    if (d2 == 1) continue;  // immediate reversal
    head[2][d2]++;
    for (int d3 = 0; d3 < 6; d3++) {
      int p1 = o + OFF[0], p2 = p1 + OFF[d2], p3 = p2 + OFF[d3];
      if (p3 == o || p3 == p1) continue;
      head[3][d3]++;
      jobs.push_back(Job{d2, d3});
    }
  }
  std::atomic<size_t> next(0);
  std::vector<Counter*> cs;
  for (int t = 0; t < threads; t++) cs.push_back(new Counter());
  std::vector<std::thread> th;
  for (int t = 0; t < threads; t++)
    th.emplace_back([&, t]() {
      Counter* c = cs[t];
      for (;;) {
        size_t i = next.fetch_add(1);
        if (i >= jobs.size()) break;
        int p1 = o + OFF[0], p2 = p1 + OFF[jobs[i].d2], p3 = p2 + OFF[jobs[i].d3];
        c->occ[o] = c->occ[p1] = c->occ[p2] = c->occ[p3] = 1;
        if (N > 3) c->dfs(3, p3);
        c->occ[o] = c->occ[p1] = c->occ[p2] = c->occ[p3] = 0;
      }
    });
  for (auto& x : th) x.join();
  printf("# end_steps N=%d: walks with first step +e1, counted by last-step direction (+x -x +y -y +z -z)\n", N);
  for (int n = 1; n <= N; n++) {
    uint64_t v[6];
    for (int d = 0; d < 6; d++) {
      v[d] = (n <= 3 ? head[n][d] : 0);
      for (auto* c : cs) v[d] += c->cnt[n][d];
    }
    printf("E %d %llu %llu %llu %llu %llu %llu\n", n, (unsigned long long)v[0], (unsigned long long)v[1],
           (unsigned long long)v[2], (unsigned long long)v[3], (unsigned long long)v[4], (unsigned long long)v[5]);
  }
  return 0;
}
