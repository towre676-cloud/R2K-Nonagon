import time, json, sys, tmint
from tmint import *
tmint.MODE = 'dft'
out = []
for nidx in [int(x) for x in sys.argv[1:]]:
    t0 = time.time(); z = [x.imag for x in acb.zeta_zeros(max(nidx-1, 1), 3)]; tz0 = time.time()-t0
    zf = [float(x.mid()) for x in z]
    if nidx == 1: lo, hi = 10.0, (zf[0]+zf[1])/2; ref = z[0]; tz = zf[0]
    else: lo, hi = (zf[0]+zf[1])/2, (zf[1]+zf[2])/2; ref = z[1]; tz = zf[1]
    c, a = (lo+hi)/2, (hi-lo)/2
    t0 = time.time(); m, st = tmint.moments([0, 1], lo, hi, 0.45, R2max=0.02); dt = time.time()-t0
    loc = arb(c) + arb(a)*(m[1]/m[0]).real
    row = {'n': nidx, 't': tz, 'time': dt, 'zero_fetch_time': tz0, 'segments': sum(st), 'mu0': m[0].real.str(12, radius=True),
           'count_is_1': bool(abs(m[0].real - 1) < 0.5), 'zero': loc.str(22, radius=True), 'overlap': bool(loc.overlaps(ref)), 'mode': 'dft'}
    out.append(row); print(row, flush=True)
    json.dump(out, open('dftcost_%s.json' % sys.argv[1], 'w'), indent=1)
