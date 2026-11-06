"""XI.1: exact one-port disclosure of an IEEE/MATPOWER grid from n(n+1)/2 grounding-intervention readings (FLINT exact rationals)."""
import sys, json, time, itertools, random
sys.path.insert(0, 'p11')
from fractions import Fraction as Fr
from math import isqrt
import numpy as np
from flint import fmpq_mat, fmpq, fmpz
import gridlib
case = sys.argv[1]; t_start = time.time()
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
Kf = gridlib.K_float(n, edges, ground); Gf = np.linalg.inv(Kf)
rho = Gf[0:1, :]**0  # placeholder
# port: maximize the minimum port correlation rho_j = G_pj / sqrt(G_pp G_jj) (computable from public data only in a pilot; here reported)
D = np.sqrt(np.diag(Gf)); C = Gf/np.outer(D, D); port = int(np.argmax(C.min(axis=1)))
perm = [port] + [i for i in range(n) if i != port]; inv_perm = {v: k for k, v in enumerate(perm)}
E = {(min(inv_perm[a], inv_perm[b]), max(inv_perm[a], inv_perm[b])): w for (a, b), w in edges.items()}; Gr = [ground[perm[k]] for k in range(n)]
def fq(x): return fmpq(x.numerator, x.denominator)
K = fmpq_mat(n, n)
for (i, j), w in E.items(): K[i, i] += fq(w); K[j, j] += fq(w); K[i, j] -= fq(w); K[j, i] -= fq(w)
for i, g in enumerate(Gr): K[i, i] += fq(g)
t0 = time.time(); G = K.inv(); t_inv = time.time() - t0
def F(q): return Fr(int(q.p), int(q.q))
G00 = G[0, 0]; th1, th2, thp = fmpq(1), fmpq(2), fmpq(1)
# forward model: exact Sherman-Morrison / 2x2 Woodbury identities on the exact inverse (spot-checked below against direct exact solves)
def h_single(j, t): return G00 - t*G[0, j]**2/(1 + t*G[j, j])
def h_pair(i, j, t):
    a = 1/t + G[i, i]; b = 1/t + G[j, j]; x = G[i, j]; det = a*b - x*x
    return G00 - (G[0, i]**2*b + G[0, j]**2*a - 2*G[0, i]*G[0, j]*x)/det
t0 = time.time(); data = [G00]
for j in range(1, n): data += [h_single(j, th1), h_single(j, th2)]
for i, j in itertools.combinations(range(1, n), 2): data.append(h_pair(i, j, thp))
t_fwd = time.time() - t0
# independent spot check: direct exact solves of (K + Theta) x = e0
rng = random.Random(118); spot = 0; e0 = fmpq_mat(n, 1); e0[0, 0] = 1
for _ in range(12):
    i, j = sorted(rng.sample(range(1, n), 2)); A = fmpq_mat(K); A[i, i] += thp; A[j, j] += thp
    spot += int(A.solve(e0)[0, 0] == h_pair(i, j, thp))
    A = fmpq_mat(K); A[j, j] += th1; spot += int(A.solve(e0)[0, 0] == h_single(j, th1))
# reconstruction (the Part VIII explicit inverse, in FLINT rationals)
def qsqrt(q):
    p_, q_ = int(q.p), int(q.q); a, b = isqrt(p_), isqrt(q_); assert a*a == p_ and b*b == q_; return fmpq(a, b)
t0 = time.time(); it = iter(data); h0 = next(it); R = fmpq_mat(n, n); R[0, 0] = h0; g = {}
for j in range(1, n):
    d1, d2 = h0 - next(it), h0 - next(it); inv_g = (1/d1 - 1/d2)/(1/th1 - 1/th2); g[j] = 1/inv_g
    R[j, j] = (1/d1 - inv_g/th1)*g[j]; R[0, j] = R[j, 0] = qsqrt(g[j])
for i, j in itertools.combinations(range(1, n), 2):
    dl = h0 - next(it); a = 1/thp + R[i, i]; b = 1/thp + R[j, j]
    R[i, j] = R[j, i] = (R[0, i]*R[0, j] - qsqrt((dl*a - g[i])*(dl*b - g[j])))/dl
Krec = R.inv(); t_rec = time.time() - t0
exact = (Krec == K)
rec_edges = {(i, j) for i in range(n) for j in range(i+1, n) if Krec[i, j] != 0}
out = {'case': case, 'meta': {k: v for k, v in meta.items() if k != 'types'}, 'hidden_nodes_n': n, 'port_bus': nodes[port],
       'edges_true': len(E), 'grounded_nodes': sum(1 for g_ in Gr if g_ != 0), 'unknown_structural_scalars': len(E) + sum(1 for g_ in Gr if g_ != 0),
       'readings': len(data), 'readings_formula_n(n+1)/2': n*(n+1)//2, 'configurations': len(data),
       'exact_recovery_K': bool(exact), 'topology_precision': len(rec_edges & set(E))/max(len(rec_edges), 1), 'topology_recall': len(rec_edges & set(E))/len(E),
       'spot_checks_direct_solve': '%d/24' % spot, 'max_digits_in_G': max(len(str(G[i, j].q)) for i in range(0, n, 7) for j in range(0, n, 7)),
       'seconds': {'inverse': t_inv, 'forward': t_fwd, 'reconstruct': t_rec, 'total': time.time() - t_start}}
print(json.dumps(out, indent=1)); json.dump(out, open('p11/grid_exact_%s.json' % case, 'w'), indent=1)
