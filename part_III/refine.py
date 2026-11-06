import time, json, sys, tmint
from tmint import *
from flint import arb_mat
tmint.MODE = sys.argv[3] if len(sys.argv) > 3 else 'series'
def cert(lo, hi, kmax=7):
    m, st = tmint.moments(list(range(kmax)), lo, hi, 0.45, R2max=0.02)
    mu = [x.real for x in m]; n = int(round(float(mu[0].mid()))); return mu, n, sum(st)
log = []; t0 = time.time()
stack = [(float(sys.argv[1]), float(sys.argv[2]))]; pairs = []
while stack:
    lo, hi = stack.pop()
    try: mu, n, seg = cert(lo, hi)
    except RuntimeError:
        stack.append((lo, hi + 0.01)); continue
    log.append({'lo': lo, 'hi': hi, 'n': n, 'segments': seg}); print(log[-1], flush=True)
    if n <= 1: continue
    if n == 2:
        M = arb_mat(2, 2, [mu[0], mu[1], mu[1], (mu[2] + mu[0])/2]); det = M.det(); a = (hi - lo)/2
        gap = arb(a)*det.sqrt(); pairs.append({'lo': lo, 'hi': hi, 'minor1': mu[0].str(10, radius=True), 'det': det.str(12, radius=True),
                                                'positive': bool(det > 0), 'gap': gap.str(16, radius=True), 'gap_mid': float(gap.mid())}); print(pairs[-1], flush=True)
        if float(gap.mid()) > 0.3*(hi - lo): continue  # well separated pair: stop
    mid = (lo + hi)/2
    # split at centroid-free midpoint; nudge if a zero sits on it
    stack.append((mid, hi)); stack.append((lo, mid))
z = [float(acb.zeta_zero(k).imag.mid()) for k in (6709, 6710)]
out = {'log': log, 'pairs': pairs, 'time': time.time() - t0, 'ref_gap': z[1] - z[0]}
json.dump(out, open('refine.json', 'w'), indent=1); print('time', out['time'])
