import time, json
from tmint import *
out = []
for nidx in (1, 100, 1000, 3000, 6709, 100000, 1000000):
    z = [x.imag for x in acb.zeta_zeros(max(nidx-1, 1), 3)]
    zf = [float(x.mid()) for x in z]
    if nidx == 1: lo, hi = 10.0, (zf[0]+zf[1])/2; ref = z[0]; tz = zf[0]
    else: lo, hi = (zf[0]+zf[1])/2, (zf[1]+zf[2])/2; ref = z[1]; tz = zf[1]
    c, a = (lo+hi)/2, (hi-lo)/2
    t0 = time.time(); m, st = moments([0, 1], lo, hi, 0.45, R2max=0.02); dt = time.time()-t0
    loc = arb(c) + arb(a)*(m[1]/m[0]).real
    row = {'n': nidx, 't': tz, 'time': dt, 'segments': sum(st), 'mu0': m[0].real.str(12, radius=True),
           'count_is_1': bool(abs(m[0].real - 1) < 0.5), 'zero': loc.str(20, radius=True), 'overlap': bool(loc.overlaps(ref))}
    out.append(row); print(row, flush=True)
    json.dump(out, open('tmcost.json', 'w'), indent=1)
