import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 40; random.seed(6060)
res = {'levels_with_real_poles': 0, 'caseA': 0, 'caseA_confirmed': 0, 'caseB': 0, 'caseB_confirmed': 0, 'resR_positive': 0, 'resL_negative': 0, 'odd_real_pole_count': 0}
fails = []
for regime in ((0.3, 1.5), (0.5, 3.0), (1.5, 6), (3, 15)):
    for t in range(80):
        n = random.choice([3, 4, 5])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(*regime)) for _ in range(n)]
        for lev in range(n-1):
            D, N = chain_polys(m[lev:], k[lev:], c[lev:])
            Dn, Nn = chain_polys(m[lev+1:], k[lev+1:], c[lev+1:])     # H_{i+1} = Nn/Dn ; note Dn == N
            rts = mp.polyroots(list(reversed(trim(Dn))), maxsteps=500, extraprec=600)
            rp = sorted([mp.re(r) for r in rts if abs(mp.im(r)) < mp.mpf(10)**-25*abs(r)])
            if not rp: continue
            res['levels_with_real_poles'] += 1
            if len(rp) % 2: res['odd_real_pole_count'] += 1
            resid = lambda p: pev(Nn, p)/pev(pder(Dn), p)
            if resid(rp[-1]) > 0: res['resR_positive'] += 1
            if resid(rp[0]) < 0: res['resL_negative'] += 1
            i1 = lev + 1; s0 = -k[i1]/c[i1]
            q, rem = pdiv(D, N); W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
            cands = [mp.re(r) for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
            extra = [x for x in cands if abs(x - s0) > mp.mpf(10)**-15*abs(s0)]
            pL, pR = rp[0], rp[-1]
            if not (pR < s0 < 0):
                res['caseA'] += 1
                if any(pR < x < 0 for x in extra): res['caseA_confirmed'] += 1
                else: fails.append(('A', float(s0), [float(p) for p in rp], [float(x) for x in extra]))
            else:
                res['caseB'] += 1
                if any(pL < x < pR for x in extra): res['caseB_confirmed'] += 1
                else: fails.append(('B', float(s0), [float(p) for p in rp], [float(x) for x in extra]))
print(json.dumps(res, indent=1)); print('fails', fails[:5])
json.dump({'summary': res, 'fails': fails}, open('m1.json', 'w'), indent=1)
