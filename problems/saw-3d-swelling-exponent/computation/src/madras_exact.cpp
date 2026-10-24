// madras_exact.cpp -- exhaustive exact statistics of the multiplicity K(T omega) of Madras's map on the near sets
// {||omega(n)||_inf <= m}, m = 1, 2, 3, 4, for every n <= N (all n-step SAWs on Z^3, no symmetry reduction).
// Self-test (n <= NV): the routine multiplicity() of madras.hpp is compared, for every near walk, with the direct
// count #{near theta : T theta = T omega} obtained by tabulating all images T theta.
// Output lines:
//   X n m near sumK hist[K=1..8+] near_span>=2 hist_span>=2[K=1..8+] nimages
// (nimages = number of distinct images T omega of near walks; only for n <= NV, else 0).
// Usage: madras_exact N [NV=9]     Build: g++ -O3 -march=native -std=c++17 madras_exact.cpp -o madras_exact
// Runtime: N = 13 about 1 min, N = 14 about 5 min (single thread).
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "madras.hpp"

static int N, NV;
static std::vector<P3> path;
static GenSet occ;
static GenSet hs;
static std::vector<P3> psi, th;
static const int MM = 4;
struct Acc {
  long long near = 0, sumK = 0, hist[9] = {0}, near2 = 0, hist2[9] = {0};
};
static Acc acc[40][MM + 1];
static bool selfTestOK = true;
// for the self-test: images of near walks, per (n, m)
static std::map<std::string, int> imgCount[40][MM + 1];
static std::vector<std::pair<std::string, int>> walkImgK[40][MM + 1];  // (image, computed K)

static std::string enc(const std::vector<P3>& w) {
  std::string s;
  for (auto& p : w)
    for (int i = 0; i < 3; i++) s += std::to_string(p.c[i]) + ",";
  return s;
}

static void visit() {
  int n = (int)path.size() - 1;
  if (n >= 1) {
    int en = ninf(path[n]);
    for (int m = 1; m <= MM; m++) {
      if (en > m) continue;
      KResult r = multiplicity(path, m, hs, psi, th);
      if (r.K < 1) selfTestOK = false;  // omega itself must be counted
      Acc& a = acc[n][m];
      a.near++;
      a.sumK += r.K;
      a.hist[r.K < 8 ? r.K : 8]++;
      MData d = madrasData(path);
      if (d.span >= 2) {
        a.near2++;
        a.hist2[r.K < 8 ? r.K : 8]++;
      }
      if (n <= NV) {
        std::string e = enc(psi);
        imgCount[n][m][e]++;
        walkImgK[n][m].push_back({e, r.K});
      }
    }
  }
  if (n == N) return;
  static const int D[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  P3 e = path[n];
  for (int d = 0; d < 6; d++) {
    P3 q{{e.c[0] + D[d][0], e.c[1] + D[d][1], e.c[2] + D[d][2]}};
    bool used = false;
    for (auto& p : path)
      if (eqP(p, q)) {
        used = true;
        break;
      }
    if (used) continue;
    path.push_back(q);
    visit();
    path.pop_back();
  }
}

int main(int argc, char** argv) {
  N = argc > 1 ? atoi(argv[1]) : 10;
  NV = argc > 2 ? atoi(argv[2]) : 9;
  hs.init(N + 1);
  path.push_back(P3{{0, 0, 0}});
  visit();
  // self-test
  for (int n = 1; n <= NV && n <= N; n++)
    for (int m = 1; m <= MM; m++)
      for (auto& pr : walkImgK[n][m])
        if (imgCount[n][m][pr.first] != pr.second) selfTestOK = false;
  printf("# madras_exact N=%d: near sets ||w(n)||_inf <= m; self-test (direct image counts, n <= %d): %s\n", N, NV,
         selfTestOK ? "PASS" : "FAIL");
  printf("# X n m near sumK hist[K=1..8+] near_span>=2 hist_span>=2[K=1..8+] nimages\n");
  for (int n = 1; n <= N; n++)
    for (int m = 1; m <= MM; m++) {
      Acc& a = acc[n][m];
      printf("X %d %d %lld %lld", n, m, a.near, a.sumK);
      for (int k = 1; k <= 8; k++) printf(" %lld", a.hist[k]);
      printf(" %lld", a.near2);
      for (int k = 1; k <= 8; k++) printf(" %lld", a.hist2[k]);
      printf(" %zu\n", n <= NV ? imgCount[n][m].size() : (size_t)0);
    }
  return selfTestOK ? 0 : 1;
}
