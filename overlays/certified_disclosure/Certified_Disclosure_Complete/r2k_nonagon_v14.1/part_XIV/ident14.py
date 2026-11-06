"""Part XIV.3: modular first-order access diagnostics, not rational nullspace certificates.
The Jacobian is exact modulo p at the declared rational model. Full column rank
modulo p proves full rational rank, hence local identifiability at that model.
Deficient modular rank and modular kernel support alone do NOT certify rational
rank deficiency or individual coordinate identifiability. See the counterexample
J=[1,p] and the exact replacement in rational_identifiability.py.

Schemes: kronK has transfer readings only; hwK intervenes at declared shunt or
generator buses; portsK intervenes at every bus. All use K farthest-point ports.
"""
import sys, json, time
import numpy as np, flint
from crb14 import load, farthest_ports
from fractions import Fraction as Fr
import gridlib
p = 2147483647

def exact_G(case, D):
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
    idx = {b: k for k, b in enumerate(D['nodes'])}; old = {k: b for k, b in enumerate(nodes)}
    n = D['n']; K = flint.nmod_mat(n, n, p)
    def res(x): x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p
    E = {}
    for (a, c), w in edges.items():
        if old[a] in idx and old[c] in idx: E[(idx[old[a]], idx[old[c]])] = w
    for (a, c), w in E.items():
        r = res(w); K[a, a] += r; K[c, c] += r; K[a, c] -= r; K[c, a] -= r
    for i, gg in enumerate(ground):
        if gg and old[i] in idx: K[idx[old[i]], idx[old[i]]] += res(gg)
    Gm = K.inv()
    return np.array([[int(x) for x in row] for row in Gm.tolist()], dtype=object).astype(np.int64)

def rows(D, Xp, Xq):
    el = D['el']; a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el])
    dp = (Xp[a] - Xp[c]) % p; dq = (Xq[a] - Xq[c]) % p
    Rb = (dp * dq) % p; Rg = (Xp[D['Gs']] * Xq[D['Gs']]) % p
    return np.vstack([Rb, Rg]).T % p

def jac(D, Gp, ports, where):
    n = D['n']; pq = [(a, b) for ia, a in enumerate(ports) for b in ports[ia:]]
    pa = np.array([a for a, b in pq]); pb = np.array([b for a, b in pq]); out = [rows(D, Gp[:, pa], Gp[:, pb])]
    for j in where:
        f = pow(int(1 + Gp[j, j]) % p, -1, p)          # theta = 1:  (K + e_j e_j^T)^{-1} = G - g_j g_j^T / (1 + G_jj)
        cj = (Gp[j, pa] * f) % p; dj = (Gp[j, pb] * f) % p
        Xa = (Gp[:, pa] - (Gp[:, j][:, None] * cj[None, :]) % p) % p
        Xb = (Gp[:, pb] - (Gp[:, j][:, None] * dj[None, :]) % p) % p
        out.append(rows(D, Xa, Xb))
    return np.vstack(out)

def analyse(D, J):
    U = J.shape[1]; M = flint.nmod_mat(J.tolist(), p); r = M.rank()
    if r == U: return r, U, U, len(D['el'])
    N, nul = M.nullspace(); Nn = np.array([[int(N[i, k]) for k in range(nul)] for i in range(U)])
    ident = ~(Nn != 0).any(1); nb = len(D['el'])
    return r, U, int(ident.sum()), int(ident[:nb].sum())

if __name__ == '__main__':
    case = sys.argv[1]; schemes = sys.argv[2].split(',')
    D = load(case); Gp = exact_G(case, D); out = {'case': case, 'n': D['n'], 'hw_buses': len(D['hw']), 'branches': len(D['el']), 'schemes': {}}
    for sc in schemes:
        t0 = time.time()
        kind = ''.join(ch for ch in sc if ch.isalpha()); k = int(''.join(ch for ch in sc if ch.isdigit()))
        ports = farthest_ports(D, k)
        where = [] if kind == 'kron' else (D['hw'] if kind == 'hw' else list(range(D['n'])))
        J = jac(D, Gp, ports, where)
        r, U, iu, ib = analyse(D, J)
        s = {'readings': int(J.shape[0]), 'rank': r, 'unknowns': U, 'kernel_zero_coordinates_mod_p': iu, 'kernel_zero_branch_coordinates_mod_p': ib,
             'rational_coordinate_classification_certified': r == U,
             'scope': 'modular first-order diagnostic; deficient counts need rational witnesses',
             'branches': len(D['el']), 'seconds': time.time() - t0}
        out['schemes'][sc] = s; print(case, sc, s, flush=True)
    json.dump(out, open('ident14_%s.json' % case, 'w'), indent=1)
