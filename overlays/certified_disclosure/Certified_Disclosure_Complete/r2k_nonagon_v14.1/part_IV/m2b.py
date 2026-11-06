import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 40; random.seed(7070)
def isreal(r): return abs(mp.im(r)) < mp.mpf(10)**-25*abs(r) and mp.re(r) < 0
rows = []
for regime in ((0.02, 0.4), (0.3, 1.5), (0.5, 3.0), (1.5, 6), (3, 15), (0.1, 8)):
    for t in range(140):
        n = random.choice([3, 4, 5, 6])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(*regime)) for _ in range(n)]
        for lev in range(n-1):
            D, N = chain_polys(m[lev:], k[lev:], c[lev:]); Dn, Nn = chain_polys(m[lev+1:], k[lev+1:], c[lev+1:])
            poles = mp.polyroots(list(reversed(trim(Dn))), maxsteps=500, extraprec=600)
            if any(isreal(r) for r in poles): continue
            i1 = lev + 1; c2 = c[i1+1] if i1+1 < n else 0
            X = 2*k[i1]/c[i1] - (c[i1] + c2)/m[i1]
            if X < 0: continue
            s0 = -k[i1]/c[i1]
            q, rem = pdiv(D, N); W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
            cands = [mp.re(r) for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
            up = [p for p in poles if mp.im(p) > 0]
            # pole pairs with real part in (s0, 0): max damping ratio, and 'aspect' |Im|/|Re - s0|
            inside = [p for p in up if s0 < mp.re(p) < 0]
            zmax_in = max([float(-mp.re(p)/abs(p)) for p in inside], default=0.0)
            zmax_all = max(float(-mp.re(p)/abs(p)) for p in up)
            asp = min([float(mp.im(p)/(mp.re(p) - s0)) for p in inside], default=1e9)
            rows.append({'multi': len(cands) > 1, 'zeta_in': zmax_in, 'zeta_all': zmax_all, 'aspect': asp, 'X': float(X), 'n_inside': len(inside)})
mul = [r for r in rows if r['multi']]; uni = [r for r in rows if not r['multi']]
print('levels', len(rows), 'multi', len(mul))
for r in mul: print('MULTI', {k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()})
import numpy as np
for key in ('zeta_in', 'zeta_all', 'aspect'):
    u = np.array([r[key] for r in uni]); mm = np.array([r[key] for r in mul])
    print(key, 'unique: min %.3f med %.3f max %.3f' % (u.min(), np.median(u), u.max()), '| multi: min %.3f max %.3f' % (mm.min(), mm.max()))
json.dump(rows, open('m2b.json', 'w'))
