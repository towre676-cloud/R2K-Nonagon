"""Exact passivity status of each grid's DC susceptance matrix K (port component).
Positive definite: certified by structure (connected, every merged branch weight > 0, grounds >= 0, at least one > 0).
Indefinite: certified by an exact rational witness v with v^T K v < 0, v taken from a binary64 LDL^T column
(x = L^{-T} e_j at the most negative pivot) and then evaluated exactly over Q."""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
from fractions import Fraction as Fr
sys.path[:0] = ['../p11']
import gridlib
out = {}
for case in sys.argv[1:]:
    rec = json.load(open('record_%s.json' % case))
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case); n = len(nodes)
    idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
    r = [a for a, c in edges]; c_ = [c for a, c in edges]
    _, lab = connected_components(sp.coo_matrix((np.ones(len(r)), (r, c_)), shape=(n, n)), directed=False)
    keep = [i for i in range(n) if lab[i] == lab[port]]; new = {o: k for k, o in enumerate(keep)}
    E = {(new[a], new[c]): w for (a, c), w in edges.items() if a in new}; g = [ground[i] for i in keep]; n = len(keep)
    negw = sum(1 for w in E.values() if w < 0)
    if negw == 0 and all(x >= 0 for x in g) and any(x > 0 for x in g):
        out[case] = {'status': 'positive definite', 'certificate': 'structure: connected, all %d merged weights > 0, grounds >= 0 with %d positive' % (len(E), sum(1 for x in g if x > 0))}
        print(case, out[case]); continue
    rows, cols, vals = [], [], []
    for (a, c), w in E.items():
        w = float(w); rows += [a, c, a, c]; cols += [a, c, c, a]; vals += [w, w, -w, -w]
    for i, x in enumerate(g):
        if x: rows.append(i); cols.append(i); vals.append(float(x))
    A = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    pm = spl.splu(A, permc_spec='MMD_AT_PLUS_A').perm_c          # fill-reducing order, used symmetrically
    B = A[pm][:, pm].tocsc()
    lu = spl.splu(B, permc_spec='NATURAL', diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
    assert (lu.perm_r == np.arange(n)).all() and (lu.perm_c == np.arange(n)).all(), 'pivoting occurred'
    d = lu.U.diagonal(); L = lu.L.tocsr()                          # B = L U, L unit lower, U = D L^T
    j = int(np.argmin(d / np.abs(d).max()))
    e = np.zeros(n); e[j] = 1.0
    x = spl.spsolve_triangular(L.T.tocsr(), e, lower=False)        # x^T B x = D_j in exact arithmetic
    v = np.zeros(n); v[pm] = x                                     # B = P A P^T with (P y)_k = y[pm[k]]
    # exact evaluation over Q of v^T K v = sum_e w_e (v_a - v_c)^2 + sum_i g_i v_i^2, v_i = exact binary64 rationals
    V = [Fr(float(t)) for t in v]
    q = sum((w * (V[a] - V[c]) ** 2 for (a, c), w in E.items()), Fr(0)) + sum((x_ * V[i] ** 2 for i, x_ in enumerate(g) if x_), Fr(0))
    out[case] = {'status': 'indefinite' if q < 0 else 'undecided', 'negative_merged_weights': negw,
                 'binary64_negative_pivots': int((d < 0).sum()), 'witness_pivot': float(d[j]),
                 'exact_vKv_sign': -1 if q < 0 else (1 if q > 0 else 0), 'exact_vKv_float': float(q),
                 'witness_support': int((v != 0).sum())}
    print(case, out[case])
json.dump(out, open('passivity13.json', 'w'), indent=1)
