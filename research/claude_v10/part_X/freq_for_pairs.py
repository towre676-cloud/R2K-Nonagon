"""X.2 (Aronszajn-Krein): damped networks from O(n) intervention configurations. At each real s > 0, P(s) = K + sC + s^2 M is an
irreducible M-matrix. Baseline + two strengths per hidden node give, pointwise in s, the full first column g(s) = P(s)^{-1} e0
(Aronszajn-Krein / Sherman-Morrison), with no pair experiments. Then P(s) g(s) = e0 is LINEAR in (K, C, M). Exact rational test."""
import json, random, itertools
from fractions import Fraction as Fr
from math import isqrt
def solve(A, b):
    n = len(A); M = [list(A[i]) + [b[i]] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0); M[k], M[p] = M[p], M[k]
        for i in range(n):
            if i != k and M[i][k] != 0:
                f = M[i][k]/M[k][k]; M[i] = [x - f*y for x, y in zip(M[i], M[k])]
    return [M[i][n]/M[i][i] for i in range(n)]
def rank_and_solve(rows, rhs, nunk):
    """exact least-squares-free solve: Gaussian elimination on an overdetermined consistent system; returns rank and solution"""
    A = [r[:] + [b] for r, b in zip(rows, rhs)]; piv = []; r = 0
    for c in range(nunk):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; pv = A[r][c]; A[r] = [x/pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    consistent = all(A[i][nunk] == 0 for i in range(r, len(A)))
    sol = [Fr(0)]*nunk
    for i, c in enumerate(piv): sol[c] = A[i][nunk]
    return r, consistent, sol
def h(P, th):
    A = [row[:] for row in P]
    for j, t in th.items(): A[j][j] += t
    return solve(A, [Fr(1)] + [Fr(0)]*(len(P) - 1))[0]
def qsqrt(q):
    a, b = isqrt(q.numerator), isqrt(q.denominator); assert a*a == q.numerator and b*b == q.denominator; return Fr(a, b)
def net(n, E, g):
    K = [[Fr(0)]*n for _ in range(n)]
    for (i, j), w in E.items(): K[i][j] -= w; K[j][i] -= w; K[i][i] += w; K[j][j] += w
    for i in range(n): K[i][i] += g[i]
    return K
def rand_damped(n, rng):
    while True:
        E = {(i, j): Fr(rng.randint(1, 6), rng.randint(1, 3)) for i in range(n) for j in range(i+1, n) if rng.random() < 0.6}
        R = {0}; ch = True
        while ch:
            ch = False
            for i, j in E:
                if (i in R) != (j in R): R |= {i, j}; ch = True
        if len(R) == n: break
    gK = [Fr(rng.randint(0, 3)) for _ in range(n)]; gK[0] += 1
    Ec = {e: Fr(rng.randint(0, 5), rng.randint(1, 4)) for e in E}            # dampers on (a subset of) the spring edges, non-proportional
    gC = [Fr(rng.randint(0, 2), rng.randint(1, 3)) for _ in range(n)]
    Mv = [Fr(rng.randint(1, 5), rng.randint(1, 2)) for _ in range(n)]
    return net(n, E, gK), net(n, Ec, gC), [[Mv[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]
def unknown_index(n):
    idx = {}
    for i in range(n):
        for j in range(i, n): idx[('K', i, j)] = len(idx)
    for i in range(n):
        for j in range(i, n): idx[('C', i, j)] = len(idx)
    for i in range(n): idx[('M', i, i)] = len(idx)
    return idx
def first_column(Ps, n, t1, t2):
    """from 1 + 2(n-1) scalar port measurements at this s: g = P(s)^{-1} e0, exactly (Sherman-Morrison per node)"""
    h0 = h(Ps, {}); g = [h0]
    for j in range(1, n):
        d1, d2 = h0 - h(Ps, {j: t1}), h0 - h(Ps, {j: t2})
        inv_g2 = (1/d1 - 1/d2)/(1/t1 - 1/t2); g.append(qsqrt(1/inv_g2))   # G_0j(s) > 0 for real s > 0
    return g
def recover(K, C, M, svals, t1=Fr(1), t2=Fr(3)):
    n = len(K); idx = unknown_index(n); rows, rhs = [], []
    for s in svals:
        Ps = [[K[i][j] + s*C[i][j] + s*s*M[i][j] for j in range(n)] for i in range(n)]
        g = first_column(Ps, n, t1, t2)
        for i in range(n):          # row i of (K + sC + s^2 M) g = e0
            row = [Fr(0)]*len(idx)
            for j in range(n):
                a, b = min(i, j), max(i, j)
                row[idx[('K', a, b)]] += g[j]; row[idx[('C', a, b)]] += s*g[j]
            row[idx[('M', i, i)]] += s*s*g[i]
            rows.append(row); rhs.append(Fr(int(i == 0)))
    r, cons, sol = rank_and_solve(rows, rhs, len(idx))
    ok = r == len(idx) and cons and all(sol[idx[('K', min(i, j), max(i, j))]] == K[i][j] and sol[idx[('C', min(i, j), max(i, j))]] == C[i][j] for i in range(n) for j in range(n)) \
         and all(sol[idx[('M', i, i)]] == M[i][i] for i in range(n))
    return r, len(idx), ok
rng = random.Random(1010); rows = []
for t in range(40):
    n = rng.randint(2, 6); K, C, M = rand_damped(n, rng)
    # smallest number of real Laplace points that gives full rank (s = 0 included: static g at s = 0 too)
    for npts in range(1, 2*n + 3):
        sv = [Fr(0)] + [Fr(k, 2) for k in range(1, npts)]
        r, nu, ok = recover(K, C, M, sv)
        if ok: break
    configs = 1 + 2*(n - 1); rows.append({'n': n, 'unknowns': nu, 's_points_needed': npts, 'configurations': configs, 'pair_protocol_configurations': n*(n+1)//2,
                                          'scalar_measurements': configs*npts, 'exact': ok})
    print(rows[-1], flush=True)
# Part VII damped twins: path and triangle
c, s = Fr(12, 13), Fr(5, 13); Q = [[1, 0, 0], [0, c, s], [0, -s, c]]
def cong(A): return [[sum(Fr(Q[k][i])*A[k][l]*Fr(Q[l][j]) for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
Kp = net(3, {(0, 1): Fr(1), (1, 2): Fr(1)}, [Fr(3), Fr(1), Fr(1)]); Cp = net(3, {(0, 1): Fr(17, 55), (1, 2): Fr(161, 948)}, [Fr(101, 472), Fr(0), Fr(281, 915)])
I3 = [[Fr(int(i == j)) for j in range(3)] for i in range(3)]
tw = {}
for name, (KK, CC) in {'path': (Kp, Cp), 'triangle': (cong(Kp), cong(Cp))}.items():
    for npts in range(1, 8):
        r, nu, ok = recover(KK, CC, I3, [Fr(0)] + [Fr(k, 2) for k in range(1, npts)])
        if ok: break
    tw[name] = {'s_points_needed': npts, 'exact': ok}
out = {'random': rows, 'all_exact': all(r['exact'] for r in rows), 'twins': tw,
       'max_s_points_by_n': {n: max(r['s_points_needed'] for r in rows if r['n'] == n) for n in sorted({r['n'] for r in rows})}}
print(json.dumps({k: v for k, v in out.items() if k != 'random'}, indent=1)); json.dump(out, open('p10/freq_for_pairs.json', 'w'), indent=1, default=str)
