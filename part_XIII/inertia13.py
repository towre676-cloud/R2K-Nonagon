"""Exact negative index of K for the indefinite grids.
Upper bound (proved): K = K_+ + K_- with K_+ PSD (positive branches, grounds) and K_- = sum over negative merged
branches of w_e b_e b_e^T, so ind_-(K) <= #negative merged branches.
Lower bound (exact certificate): for every negative binary64 LDL^T pivot j, x_j = L^{-T} e_j rounded to a dyadic
rational vector; the exact rational Gram matrix X^T K X is shown strictly diagonally dominant with negative diagonal,
hence negative definite, so ind_-(K) >= number of such columns."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = ['../p11']
import gridlib
case = sys.argv[1]; t0 = time.time()
rec = json.load(open('record_%s.json' % case))
nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case); n = len(nodes)
idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
r = [a for a, c in edges]; c_ = [c for a, c in edges]
_, lab = connected_components(sp.coo_matrix((np.ones(len(r)), (r, c_)), shape=(n, n)), directed=False)
keep = [i for i in range(n) if lab[i] == lab[port]]; new = {o: k for k, o in enumerate(keep)}
E = {(new[a], new[c]): w for (a, c), w in edges.items() if a in new}; g = [ground[i] for i in keep]; n = len(keep)
kneg = sum(1 for w in E.values() if w < 0)
rows, cols, vals = [], [], []
for (a, c), w in E.items():
    w = float(w); rows += [a, c, a, c]; cols += [a, c, c, a]; vals += [w, w, -w, -w]
for i, x in enumerate(g):
    if x: rows.append(i); cols.append(i); vals.append(float(x))
A = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
pm = spl.splu(A, permc_spec='MMD_AT_PLUS_A').perm_c
B = A[pm][:, pm].tocsc()
lu = spl.splu(B, permc_spec='NATURAL', diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
assert (lu.perm_r == np.arange(n)).all() and (lu.perm_c == np.arange(n)).all()
d = lu.U.diagonal(); LT = lu.L.T.tocsr(); negj = np.where(d < 0)[0]
SCALE = 2 ** 48
X = []
for j in negj:
    e = np.zeros(n); e[j] = 1.0
    x = spl.spsolve_triangular(LT, e, lower=False); v = np.zeros(n); v[pm] = x
    nz = np.nonzero(np.abs(v) * SCALE >= 0.5)[0]
    X.append({int(i): int(round(v[i] * SCALE)) for i in nz})
t1 = time.time()
# exact Gram over Q, only pairs that share an edge or a grounded node
touch = defaultdict(list)
inc = defaultdict(list)
for (a, c), w in E.items(): inc[a].append((a, c)); inc[c].append((a, c))
for k, x in enumerate(X):
    es = set()
    for i in x: es.update(inc[i])
    for e in es: touch[e].append(k)
    for i in x:
        if g[i]: touch[('g', i)].append(k)
G = defaultdict(Fr)
for key, ks in touch.items():
    if key[0] == 'g':
        i = key[1]; wv = g[i]
        for p_ in ks:
            for q in ks:
                if p_ <= q: G[(p_, q)] += wv * X[p_].get(i, 0) * X[q].get(i, 0)
    else:
        a, c = key; wv = E[key]
        dl = {k: X[k].get(a, 0) - X[k].get(c, 0) for k in ks}
        for p_ in ks:
            if not dl[p_]: continue
            for q in ks:
                if p_ <= q and dl[q]: G[(p_, q)] += wv * dl[p_] * dl[q]
diag = [G.get((k, k), Fr(0)) for k in range(len(X))]
off = [Fr(0)] * len(X)
for (p_, q), v in G.items():
    if p_ != q: off[p_] += abs(v); off[q] += abs(v)
dominant = all(diag[k] < 0 and -diag[k] > off[k] for k in range(len(X)))
worst = max((float(off[k] / -diag[k]) for k in range(len(X)) if diag[k] < 0), default=0.0)
out = {'case': case, 'negative_merged_branches_upper_bound': kneg, 'binary64_negative_pivots': int(len(negj)),
       'witness_columns': len(X), 'gram_negative_definite_by_exact_dominance': dominant,
       'negative_index_exact': kneg if (dominant and len(X) == kneg) else None,
       'max_offdiag_ratio': worst, 'seconds': time.time() - t0, 'seconds_gram': time.time() - t1}
print(json.dumps(out)); json.dump(out, open('inertia_%s.json' % case, 'w'), indent=1)
