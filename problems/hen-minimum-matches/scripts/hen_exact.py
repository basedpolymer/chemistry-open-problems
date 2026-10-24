"""Exact tools: (1) decide feasibility of a match set by an integer max-flow, (2) check a solution exactly.

A match set M (pairs (i,j)) is feasible iff the transshipment network
    source -> (i,t)  [cap QH[i,t]],  (i,t) -> (i,t+1) [cap inf],  (i,t) -> (j,t) for (i,j) in M [cap inf],
    (j,t) -> sink [cap QC[j,t]]
carries a flow of value sum QH = sum QC.  All data are scaled by the lcm L of the decimal denominators, so the
capacities are integers and the max-flow is an exact integer computation (Dinic's algorithm).

check_solution() is independent of the flow computation: it takes rational q_{ijt} and verifies, with
fractions.Fraction, that
  * q_{ijt} >= 0, and q_{ijt} > 0 only if (i,j) is a declared match;
  * every cold balance  sum_i q_{ijt} = QC[j,t]  holds exactly;
  * the hot residuals r_{it} = sum_{s<=t} (QH[i,s] - sum_j q_{ijs}) are >= 0 for all t and r_{i,k-1} = 0
    (heat of hot i in interval s is only used in intervals t >= s: temperature feasibility);
  * every declared match carries positive heat (so the match count is exactly the number of used pairs).
Environment: Python >= 3.10, standard library only.
"""
from fractions import Fraction
from collections import deque
import json


class Dinic:
    def __init__(self, N):
        self.N = N
        self.g = [[] for _ in range(N)]

    def add(self, u, v, c):
        self.g[u].append([v, c, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])
        return (u, len(self.g[u]) - 1)

    def bfs(self, s, t):
        self.level = [-1] * self.N
        self.level[s] = 0
        dq = deque([s])
        while dq:
            u = dq.popleft()
            for v, c, _ in self.g[u]:
                if c > 0 and self.level[v] < 0:
                    self.level[v] = self.level[u] + 1
                    dq.append(v)
        return self.level[t] >= 0

    def dfs(self, u, t, f):
        if u == t:
            return f
        while self.it[u] < len(self.g[u]):
            e = self.g[u][self.it[u]]
            v, c, rev = e
            if c > 0 and self.level[v] == self.level[u] + 1:
                d = self.dfs(v, t, min(f, c))
                if d > 0:
                    e[1] -= d
                    self.g[v][rev][1] += d
                    return d
            self.it[u] += 1
        return 0

    def maxflow(self, s, t):
        flow = 0
        INF = float("inf")
        while self.bfs(s, t):
            self.it = [0] * self.N
            while True:
                f = self.dfs(s, t, INF)
                if f == 0:
                    break
                flow += f
        return flow


def flow_for_matches(inst, matches):
    """Return (feasible, q) with q a dict (i,j,t) -> Fraction (exact), or (False, deficit)."""
    import sys
    sys.setrecursionlimit(100000)
    L = inst.denominators_lcm()
    n, m, k = inst.n, inst.m, inst.k
    S, T = 0, 1
    hid = lambda i, t: 2 + i * k + t
    cid = lambda j, t: 2 + n * k + j * k + t
    G = Dinic(2 + (n + m) * k)
    BIG = sum(int(v * L) for row in inst.QH for v in row) + 1
    for i in range(n):
        for t in range(k):
            if inst.QH[i][t] > 0:
                G.add(S, hid(i, t), int(inst.QH[i][t] * L))
            if t + 1 < k:
                G.add(hid(i, t), hid(i, t + 1), BIG)
    arcs = {}
    for (i, j) in matches:
        for t in range(k):
            if inst.QC[j][t] > 0:
                arcs[i, j, t] = G.add(hid(i, t), cid(j, t), BIG)
    for j in range(m):
        for t in range(k):
            if inst.QC[j][t] > 0:
                G.add(cid(j, t), T, int(inst.QC[j][t] * L))
    total = sum(int(v * L) for row in inst.QC for v in row)
    f = G.maxflow(S, T)
    if f != total:
        return False, Fraction(total - f, L)
    q = {}
    for key, (u, idx) in arcs.items():
        v, c, rev = G.g[u][idx]
        sent = G.g[v][rev][1]  # reverse residual = flow
        if sent > 0:
            q[key] = Fraction(sent, L)
    return True, q


def check_solution(inst, matches, q):
    """Exact verification.  matches: iterable of (i,j); q: dict (i,j,t)->Fraction.  Returns list of errors."""
    errs = []
    n, m, k = inst.n, inst.m, inst.k
    M = set(tuple(p) for p in matches)
    if len(M) != len(list(matches)):
        errs.append("duplicate matches")
    for (i, j, t), v in q.items():
        if not (0 <= i < n and 0 <= j < m and 0 <= t < k):
            errs.append(f"index out of range {(i, j, t)}")
            continue
        if v < 0:
            errs.append(f"negative q{(i, j, t)}")
        if v > 0 and (i, j) not in M:
            errs.append(f"heat on undeclared match {(i, j)}")
    for j in range(m):
        for t in range(k):
            s = sum((v for (a, b, c), v in q.items() if b == j and c == t), Fraction(0))
            if s != inst.QC[j][t]:
                errs.append(f"cold balance j={j} t={t}: {s} != {inst.QC[j][t]}")
    for i in range(n):
        r = Fraction(0)
        for t in range(k):
            used = sum((v for (a, b, c), v in q.items() if a == i and c == t), Fraction(0))
            r = r + inst.QH[i][t] - used
            if r < 0:
                errs.append(f"hot residual negative i={i} t={t}: {r} (temperature infeasible)")
        if r != 0:
            errs.append(f"hot stream {i} not fully used: residual {r}")
    for (i, j) in M:
        if sum((v for (a, b, c), v in q.items() if a == i and b == j), Fraction(0)) <= 0:
            errs.append(f"declared match {(i, j)} carries no heat")
    return errs


def save_solution(path, inst, matches, q, meta=None):
    d = {"instance": inst.name, "n": inst.n, "m": inst.m, "k": inst.k,
         "number_of_matches": len(matches),
         "matches": sorted([list(p) for p in matches]),
         "q": [[i, j, t, str(v)] for (i, j, t), v in sorted(q.items())],
         "meta": meta or {}}
    with open(path, "w") as f:
        json.dump(d, f, indent=1)


def load_solution(path):
    with open(path) as f:
        d = json.load(f)
    matches = [tuple(p) for p in d["matches"]]
    q = {(a, b, c): Fraction(v) for a, b, c, v in d["q"]}
    return d, matches, q


def prune_matches(inst, matches):
    """Greedily drop matches that are not needed (exact max-flow test).  Returns a feasible subset."""
    cur = list(matches)
    ok, _ = flow_for_matches(inst, cur)
    assert ok, "initial match set infeasible"
    changed = True
    while changed:
        changed = False
        for e in list(cur):
            trial = [f for f in cur if f != e]
            ok, _ = flow_for_matches(inst, trial)
            if ok:
                cur = trial
                changed = True
    return cur
