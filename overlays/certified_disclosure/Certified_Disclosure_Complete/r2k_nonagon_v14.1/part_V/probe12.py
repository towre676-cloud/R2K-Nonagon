import time, tmint
from tmint import *
tmint.MODE = 'dft'; tmint.K = 60; tmint.HFRAC = 4; tmint.STEPDIV = 4
calls = {'n': 0, 'ok': 0, 'fail_crude': 0, 'fail_rouche': 0}
orig = tmint.seg_moments
def wrapped(A, D, x0, x1, ks, cc, a, R2):
    calls['n'] += 1; t0 = time.time(); r = orig(A, D, x0, x1, ks, cc, a, R2)
    print('R2=%.4g ok=%s %.1fs' % (R2, r is not None, time.time()-t0), flush=True)
    if r is not None: calls['ok'] += 1
    return r
tmint.seg_moments = wrapped
lo, hi = 1e12, 1e12 + 0.5; cc = acb((lo+hi)/2); a = acb((hi-lo)/2)
A = acb(lo, 0.45); B = acb(hi, 0.45)    # top side, Re s = 0.05
t0 = time.time()
try:
    D = B - A; L = 0.5; x0 = 0.0; R2 = 0.1; segs = 0
    while x0 < 1.0 and time.time() - t0 < 240:
        for trial in range(60):
            x1 = min(1.0, x0 + (R2/4)/L*0.999)
            r = wrapped(A, D, x0, x1, [0, 1], cc, a, R2) if False else tmint.seg_moments(A, D, x0, x1, [0, 1], cc, a, R2)
            if r is not None: break
            R2 /= 2
        x0 = x1; segs += 1; R2 = min(R2*1.5, 0.1)
    print('segments', segs, 'fraction of side', x0, 'calls', calls, 'time', time.time()-t0)
except Exception as e: print('err', e)
