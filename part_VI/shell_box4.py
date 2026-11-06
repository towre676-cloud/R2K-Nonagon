"""Precursor shell: adversarial search for unforced multiplicity (X >= 0, no real poles of the sub-chain, >= 2 distinct
negative real Wronskian candidates) minimizing the sub-chain's largest damping ratio zeta_max. Float search, exact
(Fraction / Hermite) confirmation of every reported point."""
import numpy as np, json, sys
from fractions import Fraction as Fr
from numpy.polynomial import polynomial as P
from p3v6 import positive_root_count_v6
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
def chain_polys(m, k, c, mul, add, one, zero):
    n = len(m); q = [[k[i] + (k[i+1] if i+1 < n else zero), c[i] + (c[i+1] if i+1 < n else zero), m[i]] for i in range(n)]
    e = [None] + [[k[i], c[i]] for i in range(1, n)]
    d = [None]*(n+2); d[n] = [one]; d[n+1] = [zero]
    for i in range(n-1, -1, -1):
        d[i] = add(mul(q[i], d[i+1]), [-x for x in mul(mul(e[i+1], e[i+1]), d[i+2])]) if i+1 < n else mul(q[i], d[i+1])
    return d[0], d[1]
def fmul(a, b): return list(P.polymul(a, b))
def fadd(a, b): return list(P.polyadd(a, b))
def qmul(a, b):
    r = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def qadd(a, b): n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def qtrim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a
def qdivmod(a, b):
    a = list(a); q = [Fr(0)]*max(len(a)-len(b)+1, 1)
    for i in range(len(a)-len(b), -1, -1):
        cf = a[i+len(b)-1]/b[-1]; q[i] = cf
        for j in range(len(b)): a[i+j] -= cf*b[j]
    return q, qtrim(a[:len(b)-1] or [Fr(0)])
def qder(a): return [i*a[i] for i in range(1, len(a))] or [Fr(0)]
def evaluate(m, k, c):
    D, N = chain_polys(m, k, c, fmul, fadd, 1.0, 0.0); D, N = np.array(D), np.array(N)
    q, rem = P.polydiv(D, N); W = P.polysub(P.polymul(rem, P.polyder(N)), P.polymul(P.polyder(rem), N))
    X = 2*k[1]/c[1] - (c[1] + (c[2] if len(c) > 2 else 0))/m[1]
    pN = P.polyroots(N); realN = np.sum(np.abs(pN.imag) < 1e-9*np.abs(pN))
    rW = P.polyroots(W); cand = np.sum((np.abs(rW.imag) < 1e-7*np.abs(rW)) & (rW.real < 0))
    up = pN[pN.imag > 1e-9*np.abs(pN)]; zmax = float(np.max(-up.real/np.abs(up))) if len(up) else 1.0
    pD = P.polyroots(D); realD = int(np.sum(np.abs(pD.imag) < 1e-9*np.abs(pD)))
    return X, int(realN), int(cand), zmax, realD
def exact(m, k, c):
    m, k, c = [Fr(x) for x in m], [Fr(x) for x in k], [Fr(x) for x in c]
    D, N = chain_polys(m, k, c, qmul, qadd, Fr(1), Fr(0)); D, N = qtrim(D), qtrim(N)
    q, rem = qdivmod(D, N); W = qtrim(qadd(qmul(rem, qder(N)), [-x for x in qmul(qder(rem), N)]))
    neg = lambda p: [p[i]*(-1)**i for i in range(len(p))]
    X = 2*k[1]/c[1] - (c[1] + (c[2] if len(c) > 2 else 0))/m[1]
    return {'X_exact_positive': bool(X > 0), 'neg_real_roots_W': positive_root_count_v6(neg(W))['distinct_positive_roots'],
            'neg_real_roots_N': positive_root_count_v6(neg(N))['distinct_positive_roots'], 'neg_real_roots_D': positive_root_count_v6(neg(D))['distinct_positive_roots']}
def sample_wide(L):
    m = np.exp(rng.uniform(np.log(.2), np.log(5), L)); k = np.exp(rng.uniform(np.log(.2), np.log(5), L)); c = np.exp(rng.uniform(np.log(.01), np.log(10), L))
    return m, k, c
def sample(L):
    m = rng.uniform(0.5, 2, L); k = rng.uniform(1, 5, L); lo_, hi_ = [(0.02, 0.4), (0.3, 1.5), (0.5, 3.0), (1.5, 6), (3, 15), (0.1, 8)][rng.integers(6)]; c = rng.uniform(lo_, hi_, L)
    return m, k, c
def ok(e): X, rN, cand, z, rD = e; return X > 0 and rN == 0 and cand >= 2
results = {}
for nsub in [int(x) for x in sys.argv[1].split(',')]:
    L = nsub + 1; hits = []; nsamp = 0; nadm = 0
    for _ in range(30000):
        m, k, c = sample(L); e = evaluate(m, k, c); nsamp += 1
        if e[0] > 0 and e[1] == 0: nadm += 1
        if ok(e): hits.append((e[3], m, k, c))
    hits.sort(key=lambda h: h[0])
    best = []
    for z0, m, k, c in hits[:12]:
        x = np.log(np.r_[m, k, c]); zc = z0; step = 0.05
        for it in range(3000):
            y = x + step*rng.standard_normal(len(x)); y = np.r_[np.clip(y[:L], np.log(.5), np.log(2)), np.clip(y[L:2*L], 0, np.log(5)), np.clip(y[2*L:], np.log(.02), np.log(15))]; e = evaluate(np.exp(y[:L]), np.exp(y[L:2*L]), np.exp(y[2*L:]))
            if ok(e) and e[3] < zc: x, zc = y, e[3]
            if it % 500 == 499: step *= 0.6
        mm, kk, cc = np.exp(x[:L]), np.exp(x[L:2*L]), np.exp(x[2*L:])
        ex = exact(mm, kk, cc); e = evaluate(mm, kk, cc)
        best.append({'zeta_max': zc, 'm': mm.tolist(), 'k': kk.tolist(), 'c': cc.tolist(), 'exact': ex, 'realD_float': e[4],
                     'exact_unforced_multiple': bool(ex['X_exact_positive'] and ex['neg_real_roots_N'] == 0 and ex['neg_real_roots_W'] >= 2)})
    best.sort(key=lambda b: b['zeta_max'])
    zs = sorted(h[0] for h in hits)
    results[nsub] = {'samples': nsamp, 'admissible_X_pos_no_real_poles': nadm, 'unforced_hits': len(hits),
                     'hit_zeta_quantiles': [float(np.quantile(zs, q)) for q in (0, .05, .5, 1)] if zs else None, 'optimized': best}
    print(nsub, json.dumps({k_: v for k_, v in results[nsub].items() if k_ != 'optimized'}), [round(b['zeta_max'], 4) for b in best], [b['exact_unforced_multiple'] for b in best], flush=True)
json.dump(results, open('shell_box4_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
