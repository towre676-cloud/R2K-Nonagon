"""Part XV.3: maximum-likelihood estimation from ambient snapshots, against the bound of crb15.py.
Simulates N snapshots x ~ N(0, G V G + s^2 I) at every bus (chunked; only the sample covariance is kept) and one
injection step at a generator read with relative accuracy eps_inj. Fisher scoring on the Gaussian log-likelihood
plus the injection term, in relative parameters, from a nominal model 10% off and flat variances. Reports the median
relative branch error, the share within 1%, and the median ratio of error to the Cramer-Rao standard deviation."""
import sys, json, time
sys.path[:0] = ['../p14']
import numpy as np
from crb14 import load
import gridlib

def build(D, w, gv):
    n = D['n']; el = D['el']; nb = len(el); a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el])
    K = np.zeros((n, n)); wb = w[:nb]
    np.add.at(K, (a, a), wb); np.add.at(K, (c, c), wb); np.add.at(K, (a, c), -wb); np.add.at(K, (c, a), -wb)
    K[np.arange(n), np.arange(n)] += 0.0; K[D['Gs'], D['Gs']] += w[nb:]
    return np.linalg.inv(K), a, c

def factors(D, G, v, a, c):
    C = G @ (v[:, None] * G); Gs = D['Gs']
    U_w = G[:, a] - G[:, c]; H_w = C[:, a] - C[:, c]; U_g = G[:, Gs]; H_g = C[:, Gs]
    P = np.concatenate([np.stack([U_w, H_w], 1), np.stack([U_g, H_g], 1), np.stack([G, 0 * G], 1)], 2)
    Q = np.concatenate([np.stack([-H_w, -U_w], 1), np.stack([-H_g, -U_g], 1), np.stack([G, 0 * G], 1)], 2)
    return C, P, Q, U_w, U_g

def run(case, N, rho, eps_inj=1e-3, seed=11):
    D = load(case); n = D['n']; nb = len(D['el']); Gs = D['Gs']; rng = np.random.default_rng(seed)
    v0 = np.random.default_rng(3).integers(1, 101, size=n).astype(float); v0 /= v0.mean()
    G0, a, c = build(D, D['truth'], None); C0 = G0 @ (v0[:, None] * G0); s2 = rho * np.median(np.diag(C0))
    L = np.linalg.cholesky(C0 + s2 * np.eye(n)); Chat = np.zeros((n, n)); left = int(N)
    while left > 0:
        m = min(left, 20000); X = rng.standard_normal((m, n)) @ L.T; Chat += X.T @ X; left -= m
    Chat /= N
    bus, _ = gridlib.parse('../p11/data/%s.m' % case); idx = {b: k for k, b in enumerate(D['nodes'])}
    r = sorted(idx[int(float(b[0]))] for b in bus if int(float(b[1])) == 2 and int(float(b[0])) in idx)[0]
    sd_inj = eps_inj * np.abs(G0[:, r]).max(); yinj = G0[:, r] + rng.standard_normal(n) * sd_inj
    truth = np.concatenate([D['truth'], v0])
    th = np.concatenate([D['truth'] * (1 + 0.1 * rng.standard_normal(len(D['truth']))), np.ones(n)])
    t0 = time.time()
    for it in range(60):
        w, v = th[:len(D['truth'])], th[len(D['truth']):]
        G, a, c = build(D, w, None); C, P, Q, U_w, U_g = factors(D, G, v, a, c)
        Sig = C + s2 * np.eye(n); S = np.linalg.inv(Sig); U = P.shape[2]
        Pm = P.transpose(0, 2, 1).reshape(n, 2 * U); Qm = Q.transpose(0, 2, 1).reshape(n, 2 * U)
        M = (Qm.T @ S @ Pm).reshape(U, 2, U, 2); J = 0.5 * N * np.einsum('aibj,bjai->ab', M, M)
        W = S @ (Chat - Sig) @ S; g = 0.5 * N * np.einsum('nka,nm,mka->a', Q, W, P)
        gr = G[:, r]; dr = gr[a] - gr[c]
        Fi = np.zeros((n, U)); Fi[:, :nb] = -(U_w * dr[None, :]); Fi[:, nb:nb + len(Gs)] = -(U_g * gr[Gs][None, :])
        J += Fi.T @ Fi / sd_inj ** 2; g += Fi.T @ (yinj - gr) / sd_inj ** 2
        Jr = J * np.outer(th, th); gr_ = g * th
        step = np.linalg.solve(Jr + 1e-9 * np.diag(np.diag(Jr)), gr_)
        step = np.clip(step, -0.5, 0.5); th = th * (1 + step)
        if np.abs(step).max() < 1e-9: break
    err = np.abs(th - truth) / np.abs(truth); eb = err[:nb]
    cov = np.linalg.inv(Jr); sd = np.sqrt(np.clip(np.diag(cov), 0, None))[:nb]
    out = {'case': case, 'N': N, 'rho': rho, 'median_rel_err': float(np.median(eb)), 'within_1pct': float(np.mean(eb < 0.01)),
           'median_err_over_crb': float(np.median(eb / sd)), 'iterations': it + 1, 'seconds': time.time() - t0}
    print(out, flush=True); return out

if __name__ == '__main__':
    case = sys.argv[1]; res = []
    for N, rho in [(1e5, 0.0), (1e6, 0.0), (1e6, 0.01)]:
        res.append(run(case, N, rho))
    json.dump(res, open('est15_%s.json' % case, 'w'), indent=1)
