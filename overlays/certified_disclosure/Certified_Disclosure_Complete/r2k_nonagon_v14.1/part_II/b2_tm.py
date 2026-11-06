import time, json
from tmint import *
from flint import arb_mat
lo, hi, n = 51.0, 76.5, 9
t0 = time.time(); m, st = moments(list(range(2*n-1)), lo, hi, 0.45, R2max=0.05); dt = time.time()-t0
mu = [x.real for x in m]
M = arb_mat(n, n, [(mu[i+j] + mu[abs(i-j)])/2 for i in range(n) for j in range(n)])
minors = [arb_mat(k, k, [M[i, j] for i in range(k) for j in range(k)]).det() for k in range(1, n+1)]
out = {'time': dt, 'segments': st, 'mu0': mu[0].str(18, radius=True), 'minors': [x.str(15, radius=True) for x in minors], 'all_positive': all(x > 0 for x in minors)}
print(out); json.dump(out, open('b2_tm.json', 'w'), indent=1)
