"""XI.3: a-posteriori certificates. Each reading is known only to within a deterministic bound eps (ball). The explicit inverse is
evaluated in Arb ball arithmetic and K = G^{-1} enclosed by arb_mat; an edge is CERTIFIED PRESENT when the lower bound of -K_ij
exceeds 0, CERTIFIED ABSENT when |K_ij| <= eps_gap (declared minimum-edge gap w_min/2) holds for the whole ball."""
import sys, json, time, itertools
sys.path.insert(0, 'p11')
import numpy as np, gridlib
from flint import arb, arb_mat, ctx
ctx.prec = 128
case = sys.argv[1]
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
K = gridlib.K_float(n, edges, ground); G = np.linalg.inv(K)
dg = np.diag(G).copy(); ok = dg > 0
C = np.where(np.outer(ok, ok), G/np.sqrt(np.abs(np.outer(dg, dg))), 0); port = int(np.argmax(np.where(ok, C.min(axis=1), -1)))
P = [port] + [i for i in range(n) if i != port]; K = K[np.ix_(P, P)]; G = G[np.ix_(P, P)]
E = {}
for (a, b), w in edges.items():
    i, j = P.index(a), P.index(b); E[(min(i, j), max(i, j))] = float(w)
wmin = min(abs(w) for w in E.values()); t = 1/np.median(np.diag(G)[1:]); h0 = G[0, 0]
rows = []
for rel in ([float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else (1e-15, 1e-13, 1e-11, 1e-9)):
    t0 = time.time(); eps = rel*h0; rng = np.random.default_rng(5)
    def reading(x): return arb(float(x + rng.uniform(-eps, eps)), eps)   # measured value with its deterministic error bound
    T = arb(t); T2 = arb(2*t)
    H0 = reading(h0); g = [None]*n; R = [[None]*n for _ in range(n)]; R[0][0] = H0
    for j in range(1, n):
        d1 = H0 - reading(h0 - t*G[0, j]**2/(1 + t*G[j, j])); d2 = H0 - reading(h0 - 2*t*G[0, j]**2/(1 + 2*t*G[j, j]))
        inv_g = (1/d1 - 1/d2)/(1/T - 1/T2); g[j] = 1/inv_g; R[j][j] = (1/d1 - inv_g/T)*g[j]; R[0][j] = R[j][0] = g[j].sqrt()
    for i, j in itertools.combinations(range(1, n), 2):
        a = 1/t + G[i, i]; b = 1/t + G[j, j]; x = G[i, j]; p = G[0, i]; q = G[0, j]
        dl = H0 - reading(h0 - (p*p*b + q*q*a - 2*p*q*x)/(a*b - x*x))
        A = 1/T + R[i][i]; B = 1/T + R[j][j]; disc = (dl*A - g[i])*(dl*B - g[j])
        if disc < 0: disc = arb(0)
        elif not disc > 0: disc = arb(0).union(disc.upper())  # contains 0: clip to [0, upper]
        R[i][j] = R[j][i] = (R[0][i]*R[0][j] - disc.sqrt())/dl
    try:
        Kr = arb_mat(R).inv(); fin = True
    except Exception as e: fin = False
    if fin:
        present = sum(1 for (i, j) in E if -Kr[i, j] > 0); gap_absent = 0; nonedge_pairs = 0; maxrel = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in E: maxrel = max(maxrel, float(Kr[i, j].rad())/E[(i, j)]); continue
                nonedge_pairs += 1; gap_absent += int(Kr[i, j].abs_upper() < wmin/2)
        contains_truth = all(Kr[i, j].contains(arb(-w)) for (i, j), w in E.items())
        rows.append({'eps_rel': rel, 'edges_certified_present': present, 'edges': len(E), 'nonedges_certified_below_gap': gap_absent, 'nonedges': nonedge_pairs,
                     'max_edge_rel_radius': maxrel, 'all_enclosures_contain_truth': contains_truth, 'seconds': time.time() - t0})
    else: rows.append({'eps_rel': rel, 'inverse_failed': True, 'seconds': time.time() - t0})
    print(case, rows[-1], flush=True)
json.dump({'case': case, 'n': n, 'w_min': wmin, 'rows': rows}, open('p11/grid_cert_%s.json' % case, 'w'), indent=1)
