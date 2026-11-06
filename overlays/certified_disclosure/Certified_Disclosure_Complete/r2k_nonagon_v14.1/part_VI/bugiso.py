import tmint, flint; flint.ctx.cap = 128
from tmint import *
flint.ctx.cap = 128
lo, hi = 7004.95, 7005.22
for sym in (True, False):
    for (K_, HF, SD) in ((40, 8, 16), (60, 4, 4)):
        for b in (0.45, 0.6):
            tmint.MODE = 'series'; tmint.K = K_; tmint.HFRAC = HF; tmint.STEPDIV = SD
            if sym: m, st = tmint.moments_sym([0], lo, hi, b, R2max=0.2)
            else: m, st = tmint.moments([0], lo, hi, b, R2max=0.2)
            v = m[0].real
            print('sym', sym, 'K', K_, 'HF', HF, 'b', b, v.str(18, radius=True), 'contains2', v.contains(2), 'segs', sum(st), flush=True)
