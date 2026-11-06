"""Part XV.2: how much ambient data the no-actuation route needs (Cramer-Rao, binary64).
Model: N independent snapshots of bus angles x ~ N(0, Sigma), Sigma = G V G + s^2 I, sensors at every bus, V = diag of
load-fluctuation variances (unknown), s^2 = sensor noise variance (known), plus ONE calibrated injection step at a
generator read at every bus with relative accuracy eps_inj (fixes the global scale of Proposition XV.1.1).
Fisher information: (N/2) tr(S dSigma_a S dSigma_b) with S = Sigma^{-1}; each dSigma_a = P_a Q_a^T with n x 2 factors:
  w_e: P = [u, h], Q = [-h, -u], u = G b_e, h = G V G b_e;   g_i: same with b = e_i;   v_i: P = Q = [g_i, 0].
Reports k_e with relative CRB of branch e = k_e / sqrt(N), the snapshots needed for a 1% median, and the time at 30 and
60 snapshots per second (phasor reporting rates)."""
import sys, json
sys.path[:0] = ['../p14']
import numpy as np
from crb14 import load
import gridlib

def run(case, rho, eps_inj=1e-3):
    D = load(case); G = D['G']; n = D['n']; el = D['el']; Gs = D['Gs']; nb = len(el)
    rng = np.random.default_rng(3); v = rng.integers(1, 101, size=n).astype(float); v /= v.mean()
    C = G @ (v[:, None] * G); s2 = rho * np.median(np.diag(C)); Sig = C + s2 * np.eye(n); S = np.linalg.inv(Sig)
    a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el])
    U_w = G[:, a] - G[:, c]; H_w = C[:, a] - C[:, c]                         # u_e, h_e for branches
    U_g = G[:, Gs]; H_g = C[:, Gs]
    P = [np.stack([U_w, H_w], 1), np.stack([U_g, H_g], 1), np.stack([G, np.zeros_like(G)], 1)]
    Q = [np.stack([-H_w, -U_w], 1), np.stack([-H_g, -U_g], 1), np.stack([G, np.zeros_like(G)], 1)]
    P = np.concatenate(P, 2); Q = np.concatenate(Q, 2)                         # (n, 2, U)
    U = P.shape[2]; Pm = P.transpose(0, 2, 1).reshape(n, 2 * U); Qm = Q.transpose(0, 2, 1).reshape(n, 2 * U)
    M = (Qm.T @ S @ Pm).reshape(U, 2, U, 2)
    J1 = 0.5 * np.einsum('aibj,bjai->ab', M, M)                               # Fisher per snapshot
    # one injection step at the first generator, every bus read with relative accuracy eps_inj
    bus, _ = gridlib.parse('../p11/data/%s.m' % case); idx = {b: k for k, b in enumerate(D['nodes'])}
    r = sorted(idx[int(float(b[0]))] for b in bus if int(float(b[1])) == 2 and int(float(b[0])) in idx)[0]
    gr = G[:, r]; dr = gr[a] - gr[c]
    Fi = np.zeros((n, U)); Fi[:, :nb] = -(U_w * dr[None, :]); Fi[:, nb:nb + len(Gs)] = -(U_g * gr[Gs][None, :])
    Fi /= (eps_inj * np.abs(gr).max())
    Ji = Fi.T @ Fi
    truth = np.concatenate([D['truth'], v]); out = {'case': case, 'rho_sensor_noise': rho, 'eps_inj': eps_inj, 'rows': []}
    for N in (1e3, 1e4, 1e5, 1e6):
        Jt = N * J1 + Ji
        Jr = Jt * np.outer(truth, truth)                                       # relative units
        cov = np.linalg.inv(Jr); sd = np.sqrt(np.clip(np.diag(cov), 0, None))[:nb]
        out['rows'].append({'N': N, 'median_rel_sd': float(np.median(sd)), 'p90_rel_sd': float(np.percentile(sd, 90)),
                            'within_1pct': float(np.mean(sd < 0.01)), 'hours_at_30Hz': N / 30 / 3600})
        print(case, rho, out['rows'][-1], flush=True)
    return out

if __name__ == '__main__':
    case = sys.argv[1]; res = [run(case, rho) for rho in (0.0, 0.01, 0.1)]
    json.dump(res, open('crb15_%s.json' % case, 'w'), indent=1)
