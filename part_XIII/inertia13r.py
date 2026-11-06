"""Rigorous-binary64 version of inertia13.py for the largest grid.
Witness: Y = the binary64 columns x_j = L^{-T} e_j (negative pivots j), with entries below 1e-13 max|x_j| set to
exactly zero; Y is an exact dyadic matrix, so the claim concerns Y^T K Y exactly. That Gram matrix is evaluated in
binary64 together with an a priori bound on its rounding error (Higham: every term carries at most gamma_6 relative
error from the rounded weight, difference and products, a sum of m terms adds gamma_m), so
|Gram_exact - Gram_computed| <= gamma_{m+6} * (|Dm|^T |W| |Dm| + |Y|^T |Gm| |Y|) entrywise, with m the number of
terms. Strict diagonal dominance with negative diagonal is checked against these bounds; it implies Y^T K Y is
negative definite, so the negative index is at least the number of columns. The upper bound is the count of
negative merged branches, as in inertia13.py."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
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
d = lu.U.diagonal(); LT = lu.L.T.tocsr(); negj = np.where(d < 0)[0]; k = len(negj)
ri, ci, vi = [], [], []
for col, j in enumerate(negj):
    e = np.zeros(n); e[j] = 1.0
    x = spl.spsolve_triangular(LT, e, lower=False); v = np.zeros(n); v[pm] = x
    v[np.abs(v) < 1e-13 * np.abs(v).max()] = 0.0
    nz = np.nonzero(v)[0]; ri += nz.tolist(); ci += [col] * len(nz); vi += v[nz].tolist()
Y = sp.csr_matrix((vi, (ri, ci)), shape=(n, k))
el = list(E); m = len(el)
Bi = sp.csr_matrix((np.r_[np.ones(m), -np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[[a for a, c in el], [c for a, c in el]])), shape=(m, n))
w = np.array([float(E[e]) for e in el]); gv = np.array([float(x) for x in g])
Dm = (Bi @ Y).tocsc()                                  # exact differences up to one rounding each
Gram = (Dm.T @ sp.diags(w) @ Dm + Y.T @ sp.diags(gv) @ Y).toarray()
Abs = (abs(Dm).T @ sp.diags(np.abs(w)) @ abs(Dm) + abs(Y).T @ sp.diags(np.abs(gv)) @ abs(Y)).toarray()
u = 2.0 ** -53; terms = m + n + 6; gamma = terms * u / (1 - terms * u)
err = 1.01 * gamma * Abs                               # 1% inflation covers the rounding of the bound itself
diag = np.diag(Gram); derr = np.diag(err)
off = np.abs(Gram).sum(1) - np.abs(diag) + err.sum(1) - derr
dominant = bool(((diag + derr) < 0).all() and ((-(diag + derr)) > off).all())
out = {'case': case, 'negative_merged_branches_upper_bound': kneg, 'binary64_negative_pivots': int(k), 'witness_columns': int(k),
       'gram_negative_definite_by_rigorous_dominance': dominant, 'negative_index_exact': kneg if (dominant and k == kneg) else None,
       'max_offdiag_ratio': float((off / -(diag + derr)).max()), 'max_relative_error_bound': float((derr / np.abs(diag)).max()),
       'witness_nnz': int(Y.nnz), 'method': 'binary64 with a priori rounding bound', 'seconds': time.time() - t0}
print(json.dumps(out)); json.dump(out, open('inertia_%s_rigorous_binary64.json' % case, 'w'), indent=1)
