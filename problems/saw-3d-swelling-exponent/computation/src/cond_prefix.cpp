// cond_prefix.cpp -- per-prefix version of the inequality (I2) E_n[omega_n . omega_1] >= 1 (lower-bound CL3).
// Enumerates every SAW on Z^3 with first step +e1 (so omega_n . omega_1 = x1(omega_n)), and for every k <= KMAX and
// every k-step prefix u (all of them, no symmetry reduction) computes, for each n with k <= n <= N,
//   count(u, n) = #{n-step SAWs extending u},  sumx(u, n) = sum of x1(omega_n) over them.
// Prints, for each (k, n): the number of prefixes u with count > 0, the number with E[x1(omega_n) | u] < 1, the
// fraction of walks lying in such prefixes, and the minimal conditional mean (exact fraction) with its prefix.
// Output lines:  C k n nprefixes nfail walks_in_fail total_walks min_num min_den prefix(steps as digits 0..5 =
// +x -x +y -y +z -z)
// Usage: cond_prefix N KMAX        Build: g++ -O3 -march=native -std=c++17 cond_prefix.cpp -o cond_prefix
// Runtime (single thread): N = 15 about 1-2 min.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static int N, KMAX, W, W3;
static int OFF[6];
static const int DX[6] = {1, -1, 0, 0, 0, 0};
static std::vector<uint8_t> occ;
static int64_t cnt[16][48], sx[16][48];
struct Stat {
  int64_t npre = 0, nfail = 0, wfail = 0, wtot = 0;
  int64_t mnum = 1, mden = 0;  // min conditional mean mnum/mden (mden = 0: unset)
  std::string arg;
};
static Stat st[16][48];
static int path[48];

static void finalize(int k) {
  for (int n = k; n <= N; n++) {
    int64_t c = cnt[k][n], s = sx[k][n];
    if (c == 0) continue;
    Stat& S = st[k][n];
    S.npre++;
    S.wtot += c;
    if (s < c) {
      S.nfail++;
      S.wfail += c;
    }
    // compare s/c < mnum/mden
    if (S.mden == 0 || (__int128)s * S.mden < (__int128)S.mnum * c) {
      S.mnum = s;
      S.mden = c;
      S.arg.clear();
      for (int i = 1; i <= k; i++) S.arg.push_back(char('0' + path[i]));
    }
    cnt[k][n] = sx[k][n] = 0;
  }
}

static void dfs(int m, int pos, int x) {
  int kk = m < KMAX ? m : KMAX;
  for (int k = 1; k <= kk; k++) {
    cnt[k][m]++;
    sx[k][m] += x;
  }
  if (m < N) {
    for (int d = 0; d < 6; d++) {
      int q = pos + OFF[d];
      if (occ[q]) continue;
      occ[q] = 1;
      path[m + 1] = d;
      dfs(m + 1, q, x + DX[d]);
      occ[q] = 0;
    }
  }
  if (m >= 1 && m <= KMAX) finalize(m);
}

int main(int argc, char** argv) {
  N = argc > 1 ? atoi(argv[1]) : 10;
  KMAX = argc > 2 ? atoi(argv[2]) : 4;
  if (N > 40 || KMAX > 15 || KMAX > N) return 2;
  int L = N + 2;
  W = 2 * L + 1;
  W3 = W * W * W;
  OFF[0] = 1;
  OFF[1] = -1;
  OFF[2] = W;
  OFF[3] = -W;
  OFF[4] = W * W;
  OFF[5] = -W * W;
  occ.assign(W3, 0);
  int o = L + W * L + W * W * L;
  occ[o] = occ[o + 1] = 1;
  path[1] = 0;
  dfs(1, o + 1, 1);
  printf("# cond_prefix N=%d KMAX=%d: walks with first step +e1; E[x1(w_n) | k-step prefix]\n", N, KMAX);
  for (int k = 1; k <= KMAX; k++)
    for (int n = k; n <= N; n++) {
      Stat& S = st[k][n];
      printf("C %d %d %lld %lld %lld %lld %lld %lld %s\n", k, n, (long long)S.npre, (long long)S.nfail,
             (long long)S.wfail, (long long)S.wtot, (long long)S.mnum, (long long)S.mden, S.arg.c_str());
    }
  return 0;
}
