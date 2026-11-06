import numpy as np, json, mpmath as mp, time
R = json.load(open('results.json'))
mp.mp.dps = 25
def xi_ld(t):
    s = mp.mpf(0.5) + 1j * mp.mpc(t)
    return complex(1j * (1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)))
def moments(ld, c, a, b, Nn, nm=4):
    tot = np.zeros(nm, complex)
    for j in range(Nn // 2 + 1):
        th = 2 * np.pi * j / Nn
        tt = c + a * np.cos(th) + 1j * b * np.sin(th); dt = -a * np.sin(th) + 1j * b * np.cos(th)
        u = (tt - c) / a; T = [1, u]
        for k in range(2, nm): T.append(2 * u * T[-1] - T[-2])
        term = np.array(T) * ld(tt) * dt
        tot += term if (j == 0 or j == Nn // 2) else 2j * term.imag
    return tot / (2j * np.pi) * (2 * np.pi / Nn)
coarse = sorted([z[0] for z in R['E']['B1']['zeros_raw']] + [z[0] for z in R['E']['B2']['zeros_raw']])
edges = [10.0] + [(coarse[i] + coarse[i + 1]) / 2 for i in range(len(coarse) - 1)] + [76.5]
ref = R['E']['zeta_ref']
tiles = []; t0 = time.time()
for i in range(len(coarse)):
    lo, hi = edges[i], edges[i + 1]; c = (lo + hi) / 2; a = (hi - lo) / 2
    mu = moments(xi_ld, c, a, 0.45, 384)
    tz = c + a * (mu[1] / mu[0]).real
    # Newton-free second estimate from mu_2: u^2 = (mu2/mu0 + 1)/2 consistency
    cons = abs((mu[2] / mu[0]).real - (2 * ((tz - c) / a) ** 2 - 1))
    tiles.append({'lo': lo, 'hi': hi, 'mu0_dev': float(abs(mu[0] - 1)), 'mu1_imag': float(abs(mu[1].imag)),
                  't': float(tz), 'err': float(abs(tz - ref[i])), 'cheb_consistency': float(cons)})
R['E']['tiles'] = tiles; R['E']['tiles_time_s'] = time.time() - t0
# convergence study zero errors for B1 at several N (Chebyshev pencil)
import scipy.linalg as sla
def box_zeros(c, a, b, Nn, nmax):
    mu = moments(xi_ld, c, a, b, Nn, 2 * nmax + 3).real
    n_ = int(round(mu[0]))
    M0 = np.array([[(mu[i + j] + mu[abs(i - j)]) / 2 for j in range(n_)] for i in range(n_)])
    M1 = np.array([[(mu[i + j + 1] + mu[abs(i + j - 1)] + mu[abs(i - j) + 1] + mu[abs(abs(i - j) - 1)]) / 4 for j in range(n_)] for i in range(n_)])
    return np.sort((c + a * sla.eigvals(M1, M0)).real), mu[0], np.linalg.eigvalsh(M0).min()
conv = []
for Nn in (256, 512, 1024, 2048):
    zz, m0, emin = box_zeros(63.75, 12.75, 0.45, Nn, 9)
    conv.append([Nn, float(abs(m0 - 9)), float(np.max(np.abs(zz - np.array(ref[10:19])))), float(emin)])
R['E']['convB2'] = conv
# ---- Part D redo with composite panels
zeta_d, tau = 0.1, 1.0
fD = lambda s, k: s ** 2 + 2 * zeta_d * s + 1 + k * np.exp(-s * tau)
dfD = lambda s, k: 2 * s + 2 * zeta_d - k * tau * np.exp(-s * tau)
gx, gw = np.polynomial.legendre.leggauss(32)
def count(k, x0=0.0, x1=6.0, Y=40.0, panels=400):
    cs = [x0 - 1j * Y, x1 - 1j * Y, x1 + 1j * Y, x0 + 1j * Y, x0 - 1j * Y]; tot = 0
    for A_, B_ in zip(cs[:-1], cs[1:]):
        for p in range(panels):
            a0 = A_ + (B_ - A_) * p / panels; a1 = A_ + (B_ - A_) * (p + 1) / panels
            s = (a0 + a1) / 2 + (a1 - a0) / 2 * gx
            tot += np.sum(dfD(s, k) / fD(s, k) * gw) * (a1 - a0) / 2
    return tot / (2j * np.pi)
kc = R['D']['kappa_crit_analytic']
grid = np.linspace(0, 1.2, 121); cnt = [count(k) for k in grid]
R['D']['kappa_grid'] = grid.tolist(); R['D']['rhp_counts'] = [float(c.real) for c in cnt]
R['D']['count_integrality_dev'] = float(max(abs(c - round(c.real)) for c, k in zip(cnt, grid) if abs(k - kc) > 0.004))
# second crossing (count 2 -> 4)? find all crossing kappas from analytic candidates
wg_ = np.linspace(0.05, 12, 400001)
Gc = (-(1 - wg_ ** 2) - 2j * zeta_d * wg_) * np.exp(1j * wg_ * tau)
from scipy.optimize import brentq
cand = []
for i in np.where(np.sign(Gc.imag[:-1]) != np.sign(Gc.imag[1:]))[0]:
    wr = brentq(lambda wv: ((-(1 - wv ** 2) - 2j * zeta_d * wv) * np.exp(1j * wv * tau)).imag, wg_[i], wg_[i + 1], xtol=1e-16)
    kv = ((-(1 - wr ** 2) - 2j * zeta_d * wr) * np.exp(1j * wr * tau)).real
    if kv > 0: cand.append((float(kv), float(wr)))
R['D']['crossings'] = sorted(cand)[:6]
json.dump(R, open('results.json', 'w'), indent=1, default=float)
print(json.dumps(tiles, indent=0)); print(conv); print(R['D']['count_integrality_dev'], R['D']['crossings'])
print([ (round(k,2), round(c,6)) for k,c in zip(grid, R['D']['rhp_counts'])][::6])
