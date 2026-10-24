// saw_naive.cpp -- deliberately simple, independent cross-check of saw_enum.cpp.
// No symmetry reduction: enumerates every n-step SAW on Z^3 from the origin (all 6 first steps); self-avoidance is
// checked by a linear scan over the walk so far.  Records the full endpoint distribution c_n(x) at every point x,
// checks that it is invariant under all 48 signed permutations, and prints the same class format as saw_enum:
//   T n c_n h_n b_n ;  A n a b c Nsum orbit ;  H n x1 b c Nsum orbit ;  B n x1 b c Nsum orbit ;
//   M k n a b c Nsum orbit  (number of n-walks whose k-th point lies in the orbit of (a,b,c), 1 <= k < n)
// Half-space walk: x1(w_i) > 0 for 1<=i<=n.  Bridge: 0 < x1(w_i) <= x1(w_n) for 1<=i<=n.
// Usage: saw_naive N > out.txt     (N <= 12 recommended: c_12 ~ 2e8 walks, ~1 min)
// Build: g++ -O2 -std=c++17 saw_naive.cpp -o saw_naive
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

typedef std::array<int, 3> P;
static int N;
static std::vector<P> walk;
static std::vector<std::map<P, unsigned long long>> cnt, half, bridge;
static std::map<std::array<int, 5>, unsigned long long> Mtab;  // (k, n, a, b, c) -> #n-walks with w_k in orbit (a,b,c)
static P sorted_abs(const P& x) {
  P a = {abs(x[0]), abs(x[1]), abs(x[2])};
  std::sort(a.begin(), a.end());
  std::reverse(a.begin(), a.end());
  return a;
}

static void rec(int m) {
  // record walk of length m
  const P e = walk[m];  // copy: walk may reallocate below
  cnt[m][e]++;
  for (int k = 1; k < m; k++) {
    P a = sorted_abs(walk[k]);
    Mtab[{k, m, a[0], a[1], a[2]}]++;
  }
  bool hs = true;
  int mx = 0;
  for (int i = 1; i <= m; i++) {
    if (walk[i][0] <= 0) hs = false;
    mx = std::max(mx, walk[i][0]);
  }
  if (hs) {
    half[m][e]++;
    if (e[0] >= mx) bridge[m][e]++;
  }
  if (m == N) return;
  static const int D[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  for (int d = 0; d < 6; d++) {
    P q = {e[0] + D[d][0], e[1] + D[d][1], e[2] + D[d][2]};
    bool free = true;
    for (int i = 0; i <= m; i++)
      if (walk[i] == q) {
        free = false;
        break;
      }
    if (!free) continue;
    walk.push_back(q);
    rec(m + 1);
    walk.pop_back();
  }
}

int main(int argc, char** argv) {
  N = argc > 1 ? atoi(argv[1]) : 8;
  cnt.resize(N + 1);
  half.resize(N + 1);
  bridge.resize(N + 1);
  walk.push_back(P{0, 0, 0});
  rec(0);
  bool inv = true;
  printf("# saw_naive N=%d\n", N);
  for (int m = 0; m <= N; m++) {
    // invariance check under the 48 signed permutations
    int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (auto& kv : cnt[m])
      for (int p = 0; p < 6; p++)
        for (int s = 0; s < 8; s++) {
          P q;
          for (int k = 0; k < 3; k++) q[k] = ((s >> k) & 1 ? -1 : 1) * kv.first[perm[p][k]];
          auto it = cnt[m].find(q);
          if (it == cnt[m].end() || it->second != kv.second) inv = false;
        }
    // stabiliser of e1 for half-space / bridges
    for (auto* M : {&half[m], &bridge[m]})
      for (auto& kv : *M)
        for (int sw = 0; sw < 2; sw++)
          for (int s = 0; s < 4; s++) {
            P q = kv.first;
            if (sw) std::swap(q[1], q[2]);
            if (s & 1) q[1] = -q[1];
            if (s & 2) q[2] = -q[2];
            auto it = M->find(q);
            if (it == M->end() || it->second != kv.second) inv = false;
          }
    unsigned long long cn = 0, hn = 0, bn = 0;
    std::map<P, unsigned long long> A, H, B;
    for (auto& kv : cnt[m]) {
      P a = {abs(kv.first[0]), abs(kv.first[1]), abs(kv.first[2])};
      std::sort(a.begin(), a.end());
      std::reverse(a.begin(), a.end());
      A[a] += kv.second;
      cn += kv.second;
    }
    for (int pass = 0; pass < 2; pass++)
      for (auto& kv : (pass ? bridge[m] : half[m])) {
        P a = {kv.first[0], std::max(abs(kv.first[1]), abs(kv.first[2])), std::min(abs(kv.first[1]), abs(kv.first[2]))};
        (pass ? B : H)[a] += kv.second;
        (pass ? bn : hn) += kv.second;
      }
    printf("T %d %llu %llu %llu\n", m, cn, hn, bn);
    for (auto& kv : A) {
      int a = kv.first[0], b = kv.first[1], c = kv.first[2];
      int perms = (a == b && b == c) ? 1 : ((a == b || b == c) ? 3 : 6);
      int orb = perms << ((a > 0) + (b > 0) + (c > 0));
      printf("A %d %d %d %d %llu %d\n", m, a, b, c, kv.second, orb);
    }
    for (int pass = 0; pass < 2; pass++)
      for (auto& kv : (pass ? B : H)) {
        int x = kv.first[0], b = kv.first[1], c = kv.first[2];
        int orb = (b == c ? 1 : 2) * (b > 0 ? 2 : 1) * (c > 0 ? 2 : 1);
        printf("%c %d %d %d %d %llu %d\n", pass ? 'B' : 'H', m, x, b, c, kv.second, orb);
      }
  }
  for (auto& kv : Mtab) {
    int a = kv.first[2], b = kv.first[3], c = kv.first[4];
    int perms = (a == b && b == c) ? 1 : ((a == b || b == c) ? 3 : 6);
    int orb = perms << ((a > 0) + (b > 0) + (c > 0));
    printf("M %d %d %d %d %d %llu %d\n", kv.first[0], kv.first[1], a, b, c, kv.second, orb);
  }
  printf("# symmetry invariance check: %s\n", inv ? "PASS" : "FAIL");
  return inv ? 0 : 1;
}
