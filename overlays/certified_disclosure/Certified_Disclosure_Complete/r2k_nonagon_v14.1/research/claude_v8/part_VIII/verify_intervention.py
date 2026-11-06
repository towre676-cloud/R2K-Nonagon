#!/usr/bin/env python3
"""Independent verifier for the one-port intervention theorem (shares no code with intervene.py).
Measurements are simulated by Cramer's rule with fraction-free Bareiss determinants; the inverse is re-derived from
Sherman-Morrison (single node) and the 2x2 Woodbury identity (pairs). usage: verify_intervention.py [trials] [seed]"""
import sys, random, itertools, json
from fractions import Fraction as F
from math import isqrt
def det(A):
    A = [r[:] for r in A]; n = len(A); sign = 1; prev = F(1)
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None: return F(0)
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n): A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j])/prev
        prev = A[k][k]
    return sign*A[-1][-1]
def response(K, th):   # e0^T (K+Theta)^{-1} e0 = det(minor_00)/det
    A = [[K[i][j] + (th.get(i, 0) if i == j else 0) for j in range(len(K))] for i in range(len(K))]
    return det([r[1:] for r in A[1:]])/det(A) if len(K) > 1 else 1/A[0][0]
def rsqrt(q):
    a, b = isqrt(q.numerator), isqrt(q.denominator); assert a*a == q.numerator and b*b == q.denominator; return F(a, b)
def recover(n, meas, t1, t2, tp):
    h0 = meas[()]; G = {(0, 0): h0}
    for j in range(1, n):
        # Sherman-Morrison: h0 - h(theta) = theta G0j^2 / (1 + theta Gjj); two strengths -> linear system in (G0j^2, Gjj)
        D1, D2 = h0 - meas[((j, t1),)], h0 - meas[((j, t2),)]
        # D(1 + t Gjj) = t g  ->  t g - D t Gjj = D
        a11, a12, a21, a22 = t1, -D1*t1, t2, -D2*t2
        dd = a11*a22 - a12*a21; g = (D1*a22 - a12*D2)/dd; Gjj = (a11*D2 - a21*D1)/dd
        G[(0, j)] = rsqrt(g); G[(j, j)] = Gjj
    for i, j in itertools.combinations(range(1, n), 2):
        D = h0 - meas[((i, tp), (j, tp))]; p, q = G[(0, i)], G[(0, j)]; a, b = 1/tp + G[(i, i)], 1/tp + G[(j, j)]
        # Woodbury: D (ab - x^2) = p^2 b + q^2 a - 2 p q x ; choose the root with p b > q x and q a > p x (path product)
        cands = [x for x in ((p*q + s*rsqrt((D*a - p*p)*(D*b - q*q)))/D for s in (1, -1)) if p*b > q*x and q*a > p*x]
        assert len(set(cands)) == 1, 'root selection not unique'
        G[(i, j)] = cands[0]
    M = [[G[(min(i, j), max(i, j))] for j in range(n)] for i in range(n)]
    # invert by cofactors (independent of Gaussian elimination)
    dM = det(M); return [[(-1)**(i + j)*det([[M[r][c] for c in range(n) if c != i] for r in range(n) if r != j])/dM if n > 1 else 1/dM for j in range(n)] for i in range(n)]
if __name__ == '__main__':
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 40; rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    t1, t2, tp = F(1, 2), F(5), F(3, 2); results = []
    for t in range(trials):
        n = rng.randint(2, 6)
        while True:
            W = {(i, j): F(rng.randint(1, 7), rng.randint(1, 3)) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.55}
            reach = {0}; ch = True
            while ch:
                ch = False
                for (i, j) in W:
                    if (i in reach) != (j in reach): reach |= {i, j}; ch = True
            if len(reach) == n: break
        gr = [F(rng.randint(0, 4), rng.randint(1, 2)) for _ in range(n)]; gr[0] += 1
        K = [[F(0)]*n for _ in range(n)]
        for (i, j), w in W.items(): K[i][j] -= w; K[j][i] -= w; K[i][i] += w; K[j][j] += w
        for i in range(n): K[i][i] += gr[i]
        meas = {(): response(K, {})}
        for j in range(1, n):
            for tt in (t1, t2): meas[((j, tt),)] = response(K, {j: tt})
        for i, j in itertools.combinations(range(1, n), 2): meas[((i, tp), (j, tp))] = response(K, {i: tp, j: tp})
        results.append({'n': n, 'measurements': len(meas), 'minimum': n*(n + 1)//2, 'recovered': recover(n, meas, t1, t2, tp) == K})
    ok = all(r['recovered'] and r['measurements'] == r['minimum'] for r in results)
    print(json.dumps({'trials': trials, 'all_recovered_exactly_with_n(n+1)/2_measurements': ok, 'sizes': sorted({r['n'] for r in results})}))
    sys.exit(0 if ok else 1)
