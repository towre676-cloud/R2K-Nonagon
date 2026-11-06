"""High-precision R2 near the origin (all terms in mpmath, 50 digits) to bound the remainder of the origin correction."""
import numpy as np, mpmath as mp, json
mp.mp.dps = 50
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]
P = sieve(200000)   # tail beyond 2e5 moves A, B by < 1e-6 relative and does not affect the r^-2 cancellation
def R2(eps, T):
    Lg = mp.log(T/(2*mp.pi)); r = 2*mp.pi*eps/Lg; s = mp.mpc(1, r); eta = mp.mpc(0, r)
    dlog = mp.diff(lambda z: mp.zeta(z, derivative=1)/mp.zeta(z), s)
    A = mp.mpf(1); Bs = mp.mpf(0)
    for p in P:
        x = mp.power(p, -(1+eta)); A *= (1 - x)*(1 - mp.mpf(2)/p + x)/(1 - mp.mpf(1)/p)**2; Bs += (mp.log(p)/(mp.power(p, 1+eta) - 1))**2
    term = dlog + mp.power(T/(2*mp.pi), -1j*r)*mp.zeta(s)*mp.zeta(mp.conj(s))*A - Bs
    return 1 + 2/Lg**2*mp.re(term)
out = []
for T in (74920.0, 600269.0, 4992381.0):
    for eps in (1e-6, 1e-5, 1e-4):
        v = R2(mp.mpf(eps), mp.mpf(T)); out.append({'T': T, 'eps': eps, 'R2': float(v), 'R2_over_eps2': float(v/eps**2)}); print(out[-1], flush=True)
json.dump(out, open('p8/origin_hp.json', 'w'), indent=1)
