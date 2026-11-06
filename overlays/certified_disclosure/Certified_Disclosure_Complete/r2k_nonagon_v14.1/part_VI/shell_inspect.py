import json, numpy as np, sys, io, contextlib
src = open('shell_search.py').read().split('results = {}')[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
from numpy.polynomial import polynomial as P
R = json.load(open(sys.argv[1]))
for nsub, r in R.items():
    b = r['optimized'][0]; m, k, c = map(np.array, (b['m'], b['k'], b['c']))
    D, N = chain_polys(list(m), list(k), list(c), fmul, fadd, 1.0, 0.0); D, N = np.array(D), np.array(N)
    q, rem = P.polydiv(D, N); W = P.polysub(P.polymul(rem, P.polyder(N)), P.polymul(P.polyder(rem), N))
    print('nsub', nsub, 'zeta', b['zeta_max'], 'm', m.round(4), 'k', k.round(4), 'c', c.round(5))
    print('  N roots', np.round(P.polyroots(N), 5)); print('  D roots', np.round(P.polyroots(D), 5))
    print('  W roots', np.round(P.polyroots(W), 6), ' s0=-k1/c1', -k[1]/c[1], 'X', 2*k[1]/c[1] - (c[1] + (c[2] if len(c) > 2 else 0))/m[1])
    print('  exact', b['exact'], 'realD_float', b['realD_float'])
