import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 40; random.seed(7070)
tab = {}; anomalies = []
def isreal(r): return abs(mp.im(r)) < mp.mpf(10)**-25*abs(r) and mp.re(r) < 0
for regime in ((0.02, 0.4), (0.3, 1.5), (0.5, 3.0), (1.5, 6), (3, 15), (0.1, 8)):
    for t in range(140):
        n = random.choice([3, 4, 5, 6])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(*regime)) for _ in range(n)]
        for lev in range(n-1):
            D, N = chain_polys(m[lev:], k[lev:], c[lev:])
            Dn, Nn = chain_polys(m[lev+1:], k[lev+1:], c[lev+1:])
            realpole = any(isreal(r) for r in mp.polyroots(list(reversed(trim(Dn))), maxsteps=500, extraprec=600))
            realzero = (len(trim(Nn)) > 1) and any(isreal(r) for r in mp.polyroots(list(reversed(trim(Nn))), maxsteps=500, extraprec=600))
            i1 = lev + 1; c2 = c[i1+1] if i1+1 < n else 0
            X = 2*k[i1]/c[i1] - (c[i1] + c2)/m[i1]
            q, rem = pdiv(D, N); W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
            cands = [r for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
            forced = (X < 0) or realpole or realzero; multi = len(cands) > 1
            key = ('forced' if forced else 'unforced') + '|' + ('multi' if multi else 'unique'); tab[key] = tab.get(key, 0) + 1
            if multi != forced: anomalies.append({'X': float(X), 'realpole': realpole, 'realzero': realzero, 'n_cand': len(cands)})
            sub = tab.setdefault('by_cause', {}); 
            if multi: sub[f"X<0:{X<0}|pole:{realpole}|zero:{realzero}"] = sub.get(f"X<0:{X<0}|pole:{realpole}|zero:{realzero}", 0) + 1
print(json.dumps(tab, indent=1)); print('anomalies', len(anomalies), anomalies[:8])
json.dump({'table': tab, 'anomalies': anomalies}, open('m2.json', 'w'), indent=1)
