import time, json, sys, tmint
from tmint import *
out = []
for nidx, mode in ((1000000, 'dft'), (100000000, 'dft')):
    tmint.MODE = mode
    z = [x.imag for x in acb.zeta_zeros(nidx-1, 3)]; zf = [float(x.mid()) for x in z]
    lo, hi = (zf[0]+zf[1])/2, (zf[1]+zf[2])/2; c, a = (lo+hi)/2, (hi-lo)/2
    for b, R2 in ((0.45, 0.2), (0.6, 0.2)):
        t0 = time.time(); m, st = tmint.moments_sym([0, 1], lo, hi, b, R2max=R2); dt = time.time()-t0
        loc = arb(c) + arb(a)*(m[1]/m[0]).real
        row = {'n': nidx, 'b': b, 'R2max': R2, 'time': dt, 'segments': st, 'mu0': m[0].real.str(12, radius=True),
               'zero': loc.str(20, radius=True), 'overlap': bool(loc.overlaps(z[1]))}
        out.append(row); print(row, flush=True)
json.dump(out, open('m3bench.json', 'w'), indent=1)
