import tmint
from tmint import *
from flint import acb, arb
LOGPI_ = arb.pi().log(); I_ = acb(0, 1)
def xild(s):
    z = acb_series([s, 1], prec=2).zeta()
    return 1/s + 1/(s-1) - LOGPI_/2 + (s/2).digamma()/2 + z[1]/z[0]
import flint
tmint.MODE = "series"; tmint.K = 60; tmint.HFRAC = 4; tmint.STEPDIV = 4
rec = []; orig = tmint.seg_moments
def wrapped(A, D, x0, x1, ks, cc, a, R2):
    r = orig(A, D, x0, x1, ks, cc, a, R2)
    if r is not None: rec.append((A, D, x0, x1, R2, r[0]))
    return r
tmint.seg_moments = wrapped
lo, hi = 7004.95, 7005.22
m, st = tmint.moments_sym([0], lo, hi, 0.6, R2max=0.2)
bad = 0
for (A, D, x0, x1, R2, r) in rec:
    def g(x, an, A=A, D=D):
        t = A + D*x; return I_*xild(acb(0.5) + I_*t)*D
    ex = acb.integral(g, x0, x1, rel_tol=2.0**-70)
    diff = (r - ex)
    if not r.overlaps(ex):
        bad += 1
        if bad <= 6: print('BAD side A=%s x=[%g,%g] R2=%.3g  rad=%.1e  |diff|=%.2e' % (A.str(6), x0, x1, R2, float(r.rad()), float(abs(diff).mid())))
print('segments', len(rec), 'bad', bad)
