import time, json
from tmint import *
from flint import arb_mat
lo, hi = 7004.95, 7005.22; c = (lo+hi)/2; a = (hi-lo)/2
t0 = time.time(); m, st = moments([0, 1, 2], lo, hi, 0.45, R2max=0.02)
mu = [x.real for x in m]
M = arb_mat(2, 2, [mu[0], mu[1], mu[1], (mu[2]+mu[0])/2]); det = M.det()
z1 = acb.zeta_zero(6709).imag; z2 = acb.zeta_zero(6710).imag
gap = arb(a)*det.sqrt()
out = {'mu0': m[0].real.str(18, radius=True), 'mu0_imag': m[0].imag.str(5, radius=True), 'minors': [mu[0].str(18, radius=True), det.str(18, radius=True)],
       'gap_from_det': gap.str(18, radius=True), 'gap_ref': (z2-z1).str(18, radius=True), 'overlap': bool(gap.overlaps(z2-z1)),
       'segments': st, 'time': time.time()-t0, 'prec': 80, 'method': 'taylor-model'}
print(out); json.dump(out, open('lehmer_tm.json', 'w'), indent=1)
