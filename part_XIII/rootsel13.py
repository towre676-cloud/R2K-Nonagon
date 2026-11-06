"""Root selection for indefinite K (Proposition XIII.1.1): with a positive port column, every pair strength
theta_p < theta* = 1 / (||G||_2 (1 + rho)), rho = max_i G_0i / min_i G_0i, makes the Part VIII minus root the
true G_ij for every pair, and keeps Theta^{-1} + G_S positive definite. Computes theta* in binary64 and
compares with G_00 (the scale of the port reading)."""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
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
    E = {(new[a], new[c]): w for (a, c), w in edges.items() if a in new}; g = [ground[i] for i in keep]; n = len(keep); port = new[port]
    rows, cols, vals = [], [], []
    for (a, c), w in E.items():
        w = float(w); rows += [a, c, a, c]; cols += [a, c, c, a]; vals += [w, w, -w, -w]
    for i, x in enumerate(g):
        if x: rows.append(i); cols.append(i); vals.append(float(x))
    A = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    lu = spl.splu(A, permc_spec='MMD_AT_PLUS_A'); e = np.zeros(n); e[port] = 1; g0 = lu.solve(e)
    res = np.abs(A @ g0 - e).max()
    op = spl.LinearOperator((n, n), matvec=lu.solve, dtype=float)
    lam = spl.eigsh(op, k=1, which='LM', return_eigenvectors=False, tol=1e-8)[0]   # largest |eig| of G
    normG = abs(lam); rho = g0.max() / g0.min()
    theta = 1.0 / (normG * (1 + rho))
    out[case] = {'port_column_min': float(g0.min()), 'port_column_positive': bool(g0.min() > 0), 'rho': float(rho),
                 'norm_G_2': float(normG), 'theta_star': float(theta), 'G00': float(g0[port]), 'theta_star_times_G00': float(theta * g0[port]),
                 'solve_residual': float(res)}
    print(case, out[case], flush=True)
json.dump(out, open('rootsel13.json', 'w'), indent=1)
