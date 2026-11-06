import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 50; random.seed(4041)
tab = {'realpoles&multi': 0, 'realpoles&unique': 0, 'noreal&multi': 0, 'noreal&unique': 0}
for regime in ((0.02, 0.4), (0.3, 1.5), (0.5, 3.0), (3, 15)):
    for t in range(50):
        n = random.choice([3, 4, 5])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(*regime)) for _ in range(n)]
        D, N = chain_polys(m, k, c)
        for lev in range(n-1):
            # current level: (D,N); real negative roots of N = real poles of next-level H2
            rts = mp.polyroots(list(reversed(trim(N))), maxsteps=500, extraprec=600)
            has_real = any(abs(mp.im(r)) < mp.mpf(10)**-25*abs(r) for r in rts)
            q, rem = pdiv(D, N); W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
            cands = [r for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
            key = ('realpoles' if has_real else 'noreal') + ('&multi' if len(cands) > 1 else '&unique'); tab[key] += 1
            # advance one level using true parameters
            Dn, Nn = chain_polys(m[lev+1:], k[lev+1:], c[lev+1:])
            D, N = Dn, Nn
print(tab); json.dump(tab, open('p4c.json', 'w'))
