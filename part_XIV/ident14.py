"""Part XIV.3: exact local identifiability under realistic access, modulo a prime.
The reading map y(w) is rational in the unknowns; its Jacobian at the true (rational) values is computed exactly
modulo p from G = K^{-1} mod p. Full column rank of the Jacobian mod p implies full rank over Q, hence local
identifiability of every unknown. When the rank is short, the unknowns that are identifiable to first order are
exactly those on which every kernel vector vanishes (computed from the FLINT nullspace mod p; a kernel vector mod p
may in principle exist only mod p, so the identifiable count is a lower bound).
Schemes (ports by farthest-point sampling from the Part XIII port, interventions of strength 1):
  kronK  : K ports, transfer compliances only, no interventions
  hwK    : K ports, one intervention at every bus that has a shunt or a generator in the source file
  portsK : K ports, one intervention at every bus"""
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
        s = {'readings': int(J.shape[0]), 'rank': r, 'unknowns': U, 'identifiable_unknowns': iu, 'identifiable_branches': ib,
             'branches': len(D['el']), 'seconds': time.time() - t0}
        out['schemes'][sc] = s; print(case, sc, s, flush=True)
    json.dump(out, open('ident14_%s.json' % case, 'w'), indent=1)
