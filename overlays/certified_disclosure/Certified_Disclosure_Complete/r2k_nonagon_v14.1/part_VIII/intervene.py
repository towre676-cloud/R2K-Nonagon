from pathlib import Path
HERE = Path(__file__).resolve().parent
"""VIII.1-3: one-port recovery of a grounded spring (and damped) network from calibrated grounding interventions.
Exact rational arithmetic throughout. Port = node 0. Measurement: h(Theta) = e0^T (K + Theta)^{-1} e0, Theta diagonal >= 0 on nodes 1..n-1."""
import json, random, itertools
from fractions import Fraction as Fr
from math import isqrt, comb
# ---------- exact linear algebra
def solve(A, b):
    n = len(A); M = [list(map(Fr, A[i])) + [Fr(b[i])] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0); M[k], M[p] = M[p], M[k]
        for i in range(n):
            if i != k and M[i][k] != 0:
                f = M[i][k]/M[k][k]; M[i] = [x - f*y for x, y in zip(M[i], M[k])]
    return [M[i][n]/M[i][i] for i in range(n)]
def inv(A):
    n = len(A); cols = [solve(A, [Fr(int(i == j)) for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]
def h(K, theta):  # theta: dict node -> value (Fraction or 'inf' not used)
    A = [row[:] for row in K]
    for j, t in theta.items(): A[j][j] += t
    return solve(A, [Fr(1)] + [Fr(0)]*(len(K) - 1))[0]
def qsqrt(q):
    q = Fr(q); a, b = isqrt(q.numerator), isqrt(q.denominator)
    if a*a != q.numerator or b*b != q.denominator: raise ValueError('not a rational square')
    return Fr(a, b)
# ---------- networks
def network(n, edges, grounds):
    K = [[Fr(0)]*n for _ in range(n)]
    for (i, j), w in edges.items(): K[i][j] -= w; K[j][i] -= w; K[i][i] += w; K[j][j] += w
    for i, g in enumerate(grounds): K[i][i] += g
    return K
def random_network(n, rng, p=0.6):
    while True:
        edges = {(i, j): Fr(rng.randint(1, 9), rng.randint(1, 4)) for i in range(n) for j in range(i+1, n) if rng.random() < p}
        adj = {i: set() for i in range(n)}
        for i, j in edges: adj[i].add(j); adj[j].add(i)
        seen, st = {0}, [0]
        while st:
            u = st.pop()
            for v in adj[u] - seen: seen.add(v); st.append(v)
        if len(seen) == n: break
    grounds = [Fr(rng.randint(0, 5), rng.randint(1, 3)) for _ in range(n)]
    if all(g == 0 for g in grounds): grounds[rng.randrange(n)] = Fr(1)
    return network(n, edges, grounds), edges, grounds
# ---------- the protocol
def protocol(n, th1, th2, thp):
    """list of experiments: baseline, two strengths per hidden node, one pair experiment per hidden pair. Count n(n+1)/2."""
    E = [{}]
    for j in range(1, n): E += [{j: th1}, {j: th2}]
    for i, j in itertools.combinations(range(1, n), 2): E.append({i: thp, j: thp})
    return E
def reconstruct(n, data, th1, th2, thp):
    """data: list of h values in protocol order. Returns K (exact) via the explicit inverse."""
    it = iter(data); h0 = next(it); G = [[None]*n for _ in range(n)]; G[0][0] = h0
    g, c = {}, {}
    for j in range(1, n):
        d1, d2 = h0 - next(it), h0 - next(it)          # d = theta g/(1 + theta Gjj): 1/d = (1/g)(1/theta) + Gjj/g
        inv_g = (1/d1 - 1/d2)/(1/th1 - 1/th2); Gjj_over_g = 1/d1 - inv_g/th1
        g[j] = 1/inv_g; G[j][j] = Gjj_over_g*g[j]; G[0][j] = G[j][0] = qsqrt(g[j])   # G_0j > 0 (connected grounded network)
    for i, j in itertools.combinations(range(1, n), 2):
        dl = h0 - next(it); a = 1/thp + G[i][i]; b = 1/thp + G[j][j]
        disc = (dl*a - g[i])*(dl*b - g[j])                       # factorized discriminant, >= 0
        G[i][j] = G[j][i] = (G[0][i]*G[0][j] - qsqrt(disc))/dl    # minus root: path-product inequality
    return inv(G), G
# ---------- checks
rng = random.Random(8080); OUT = {'random_networks': []}
th1, th2, thp = Fr(1), Fr(3), Fr(2)
for trial in range(60):
    n = rng.randint(2, 7); K, edges, grounds = random_network(n, rng)
    E = protocol(n, th1, th2, thp); data = [h(K, e) for e in E]
    Kr, G = reconstruct(n, data, th1, th2, thp)
    # both roots of the pair quadratic for the record, and which one is physical
    OUT['random_networks'].append({'n': n, 'experiments': len(E), 'n(n+1)/2': n*(n+1)//2, 'exact_recovery': Kr == K})
ok = all(r['exact_recovery'] and r['experiments'] == r['n(n+1)/2'] for r in OUT['random_networks'])
OUT['random_summary'] = {'networks': len(OUT['random_networks']), 'all_exact': ok, 'sizes': sorted(set(r['n'] for r in OUT['random_networks']))}
# repeated modes: complete graph K4 with equal springs and equal grounds has a triple eigenvalue
Krep = network(4, {(i, j): Fr(1) for i in range(4) for j in range(i+1, 4)}, [Fr(1)]*4)
Kr, _ = reconstruct(4, [h(Krep, e) for e in protocol(4, th1, th2, thp)], th1, th2, thp)
OUT['repeated_modes_K4'] = {'eigenvalues': 'grounded K4: 1 (simple) and 5 (multiplicity 3)', 'exact_recovery': Kr == Krep}
# the frontier path/triangle pair (same compliance) and the Part VII damped twin are separated by the protocol
Kpath = network(3, {(0, 1): Fr(1), (1, 2): Fr(1)}, [Fr(3), Fr(1), Fr(1)])
c, s = Fr(12, 13), Fr(5, 13); Q = [[1, 0, 0], [0, c, s], [0, -s, c]]
Ktri = [[sum(Fr(Q[k][i])*Kpath[k][l]*Fr(Q[l][j]) for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
dp = [h(Kpath, e) for e in protocol(3, th1, th2, thp)]; dt = [h(Ktri, e) for e in protocol(3, th1, th2, thp)]
OUT['path_vs_triangle'] = {'baseline_equal': dp[0] == dt[0], 'experiments_that_differ': sum(1 for x, y in zip(dp, dt) if x != y), 'of': len(dp),
                           'path_recovered': reconstruct(3, dp, th1, th2, thp)[0] == Kpath, 'triangle_recovered': reconstruct(3, dt, th1, th2, thp)[0] == Ktri}
# root selection: the + root is never physical (record both roots on all pairs of the random networks)
plus_ok = 0; total = 0; eq_cases = 0; rng2 = random.Random(5)
for trial in range(60):
    n = rng2.randint(3, 6); K, _, _ = random_network(n, rng2); G = inv(K)
    for i, j in itertools.combinations(range(1, n), 2):
        dl = G[0][0] - h(K, {i: thp, j: thp}); a = 1/thp + G[i][i]; b = 1/thp + G[j][j]
        r = qsqrt((dl*a - G[0][i]**2)*(dl*b - G[0][j]**2)); xm = (G[0][i]*G[0][j] - r)/dl; xp = (G[0][i]*G[0][j] + r)/dl
        total += 1; plus_ok += (xp == G[i][j]) and (r != 0)
        if xm != G[i][j]: raise RuntimeError('wrong physical root in exact campaign')
        # path-product inequalities that select the root
        if not (G[0][i]*G[j][j] >= G[0][j]*G[i][j] and G[0][j]*G[i][i] >= G[0][i]*G[i][j]): raise RuntimeError('path-product violation')
        eq_cases += (G[0][i]*G[j][j] == G[0][j]*G[i][j]) or (G[0][j]*G[i][i] == G[0][i]*G[i][j])
OUT['root_selection'] = {'pairs_checked': total, 'minus_root_correct': total, 'plus_root_ever_correct': plus_ok, 'path_product_holds': True, 'path_product_equality_cases (cut vertex)': eq_cases}
# Kron reduction: an inaccessible node. Star with hidden centre 3 vs triangle on {0,1,2} (Y-Delta): identical h for every Theta on {1,2}
Kstar = network(4, {(0, 3): Fr(1), (1, 3): Fr(2), (2, 3): Fr(3)}, [Fr(1), Fr(1), Fr(1), Fr(0)])
def kron(K, keep):
    drop = [i for i in range(len(K)) if i not in keep]
    Kdd = [[K[i][j] for j in drop] for i in drop]; Kdi = inv(Kdd)
    return [[K[i][j] - sum(K[i][p]*Kdi[a][b]*K[q][j] for a, p in enumerate(drop) for b, q in enumerate(drop)) for j in keep] for i in keep]
Kd = kron(Kstar, [0, 1, 2])
th_samples = [{}, {1: Fr(1)}, {2: Fr(5, 2)}, {1: Fr(2), 2: Fr(1, 3)}, {1: Fr(7), 2: Fr(7)}]
OUT['kron_ambiguity'] = {'reduced_triangle': [[str(x) for x in r] for r in Kd],
                         'h_equal_on_samples': all(h(Kstar, t) == h(Kd, t) for t in th_samples),
                         'Y_Delta_weights': 'edge_ij = w_i w_j / (w_0 + w_1 + w_2): 1*2/6, 1*3/6, 2*3/6 = 1/3, 1/2, 1'}
# damped extension: P(s) = K + s C + s^2 M at s = 0 and two positive s recovers K, C, M (Part VII damped twin separated)
def P(K, C, M, sv): return [[K[i][j] + sv*C[i][j] + (sv*sv*M[i][j]) for j in range(3)] for i in range(3)]
C = network(3, {(0, 1): Fr(17, 55), (1, 2): Fr(161, 948)}, [Fr(101, 472), Fr(0), Fr(281, 915)]); M = [[Fr(int(i == j)) for j in range(3)] for i in range(3)]
Ct = [[sum(Fr(Q[k][i])*C[k][l]*Fr(Q[l][j]) for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
res = {}
for name, (KK, CC) in {'damped_path': (Kpath, C), 'damped_triangle': (Ktri, Ct)}.items():
    rec = {}
    for sv in (Fr(0), Fr(1), Fr(2)):
        Ps = P(KK, CC, M, sv)
        try: rec[sv] = reconstruct(3, [h(Ps, e) for e in protocol(3, th1, th2, thp)], th1, th2, thp)[0]
        except ValueError: rec[sv] = None
    if None in rec.values(): res[name] = 'irrational square roots: exact rational path not applicable at this s'; continue
    K0 = rec[Fr(0)]; D1 = [[rec[Fr(1)][i][j] - K0[i][j] for j in range(3)] for i in range(3)]; D2 = [[rec[Fr(2)][i][j] - K0[i][j] for j in range(3)] for i in range(3)]
    Mr = [[(D2[i][j] - 2*D1[i][j])/2 for j in range(3)] for i in range(3)]; Cr = [[D1[i][j] - Mr[i][j] for j in range(3)] for i in range(3)]
    res[name] = {'K': K0 == KK, 'C': Cr == CC, 'M': Mr == M}
OUT['damped_extension'] = res
print(json.dumps(OUT, indent=1, default=str)[:3000]); json.dump(OUT, open(HERE / 'intervene.json', 'w'), indent=1, default=str)
