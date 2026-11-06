"""Bogomolny-Keating / Conrey-Snaith pair correlation -> number variance, against zeta zeros near the 10^6-th."""
import numpy as np, mpmath as mp, json, time
mp.mp.dps = 20
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].astype(float)
P = sieve(2_000_000); lp = np.log(P)
def A(eta):   # Euler product, eta complex
    x = P**(-(1+eta)); t = (1 - x)*(1 - 2/P + x)/(1 - 1/P)**2
    return complex(np.exp(np.sum(np.log(t.astype(complex)))))
def B(eta):
    return complex(np.sum((lp/(P**(1+eta) - 1))**2))
def R2(eps, T):
    Lg = np.log(T/(2*np.pi)); r = 2*np.pi*eps/Lg
    if eps == 0: return 0.0
    s = mp.mpc(1, r)
    dlog = mp.diff(lambda z: mp.zeta(z, derivative=1)/mp.zeta(z), s)
    term = complex(dlog) + complex(mp.power(T/(2*mp.pi), -1j*r)*mp.zeta(s)*mp.zeta(mp.conj(s)))*A(1j*r) - B(1j*r)
    return 1 + 2/Lg**2*term.real
def R2_gue(eps): return 1 - (np.sinc(eps))**2
T = 605744.0
xs = np.linspace(1e-4, 16, 1601); t0 = time.time()
r2 = np.array([R2(x, T) for x in xs]); tb = time.time() - t0
def sigma2(L, x, r):
    sel = x <= L; xx = x[sel]; yy = 1 - r[sel]
    return float(L - 2*np.trapezoid((L - xx)*yy, xx))
Ls = [0.5, 1, 2, 4, 8, 16]
emp = json.load(open('m5c.json'))
out = {'T': T, 'time_s': tb, 'rows': []}
# empirical bootstrap error bars (block bootstrap over 10 blocks of 2000 zeros)
from flint import acb
Z = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(1000000, 20000)])
U = Z/(2*np.pi)*np.log(Z/(2*np.pi*np.e)) + 7/8; U -= U[0]
rng = np.random.default_rng(98)
def s2emp(pts, L, M=6000):
    lo = rng.uniform(pts[0], pts[-1]-L, M); c = np.searchsorted(pts, lo+L) - np.searchsorted(pts, lo); return np.var(c)
for L in Ls:
    blocks = [s2emp(U[i*2000:(i+1)*2000], L) for i in range(10)]
    row = {'L': L, 'zeta': emp[str(float(L)) if str(float(L)) in emp else str(L)]['zeta'][0], 'zeta_block_se': float(np.std(blocks)/np.sqrt(10)),
           'BK': sigma2(L, xs, r2), 'GUE': sigma2(L, xs, np.array([R2_gue(x) for x in xs])), 'cue_sim': emp[str(float(L)) if str(float(L)) in emp else str(L)]['cue']}
    out['rows'].append(row); print(row, flush=True)
out['R2_sample'] = {str(round(x, 3)): float(v) for x, v in zip(xs[::100], r2[::100])}
json.dump(out, open('bk.json', 'w'), indent=1)
