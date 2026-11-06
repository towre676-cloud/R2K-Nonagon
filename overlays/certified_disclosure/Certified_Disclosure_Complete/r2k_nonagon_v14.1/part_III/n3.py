import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 40; random.seed(5050)
tab = {}; detail = []
def key(Xneg, realp, multi): return f"X{'<0' if Xneg else '>0'}|{'realpoles' if realp else 'noreal'}|{'multi' if multi else 'unique'}"
for regime in ((0.02, 0.4), (0.3, 1.5), (0.5, 3.0), (1.5, 6), (3, 15)):
    for t in range(60):
        n = random.choice([3, 4, 5])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(*regime)) for _ in range(n)]
        for lev in range(n-1):
            D, N = chain_polys(m[lev:], k[lev:], c[lev:])
            rts = mp.polyroots(list(reversed(trim(N))), maxsteps=500, extraprec=600)
            realp = any(abs(mp.im(r)) < mp.mpf(10)**-25*abs(r) for r in rts)
            q, rem = pdiv(D, N); W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
            cands = [mp.re(r) for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
            i1 = lev + 1; c2 = c[i1+1] if i1+1 < n else 0
            X = 2*k[i1]/c[i1] - (c[i1] + c2)/m[i1]
            s0 = -k[i1]/c[i1]
            left = sum(1 for x in cands if x < s0 - mp.mpf(10)**-15*abs(s0)); right = sum(1 for x in cands if x > s0 + mp.mpf(10)**-15*abs(s0))
            kk = key(X < 0, realp, len(cands) > 1); tab[kk] = tab.get(kk, 0) + 1
            if (X >= 0) and (not realp) and len(cands) > 1: detail.append({'X': float(X), 'left': left, 'right': right, 'n_cand': len(cands)})
print(json.dumps(tab, indent=1)); print('unexplained (X>0, no real poles, multi):', len(detail), detail[:6])
json.dump({'table': tab, 'unexplained': detail}, open('n3.json', 'w'), indent=1)
