import time, json, tmint
from tmint import *
res = {}
for mode in ('series', 'dft'):
    tmint.MODE = mode
    t0 = time.time(); m, st = tmint.moments([0, 1], 10.0, 17.578382390253125, 0.45, R2max=0.02)
    res[mode] = [m[0].real.str(18, radius=True), m[1].real.str(18, radius=True), time.time()-t0, sum(st)]
    print(mode, res[mode], flush=True)
