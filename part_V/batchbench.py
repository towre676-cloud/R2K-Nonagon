import time, json, tmint
from tmint import *
tmint.MODE = 'dft'
z = [x.imag for x in acb.zeta_zeros(999999, 3)]; zf = [float(x.mid()) for x in z]
lo, hi = (zf[0]+zf[1])/2, (zf[1]+zf[2])/2; c, a = (lo+hi)/2, (hi-lo)/2
out = []
for K_, HF, SD in ((40, 8, 16), (60, 4, 4), (60, 3, 3)):
    tmint.K = K_; tmint.HFRAC = HF; tmint.STEPDIV = SD
    t0 = time.time(); m, st = tmint.moments_sym([0, 1], lo, hi, 0.6, R2max=0.2); dt = time.time() - t0
    loc = arb(c) + arb(a)*(m[1]/m[0]).real
    r = {'K': K_, 'h_le_R1_over': HF, 'seg_len_R2_over': SD, 'segments': sum(st), 'evals': sum(st)*tmint.NDFT, 'time': dt,
         'mu0': m[0].real.str(12, radius=True), 'zero': loc.str(20, radius=True), 'overlap': bool(loc.overlaps(z[1]))}
    out.append(r); print(r, flush=True)
json.dump(out, open('batchbench.json', 'w'), indent=1)
