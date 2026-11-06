"""N6: root-free Lehmer scanning with rigorous tiles on a fixed grid (no zero positions used)."""
import time, json, sys, tmint
from tmint import *
from flint import arb_mat
tmint.MODE = sys.argv[3] if len(sys.argv) > 3 else 'series'
T0, T1, W = float(sys.argv[1]), float(sys.argv[2]), 2.0
KMAX = 11  # supports up to 6 zeros per tile
def disc_ref(n):
    u = [-1 + (2*j - 1)/n for j in range(1, n+1)]; d = arb(2)**((n-1)*(n-2))
    for i in range(n):
        for j in range(i+1, n): d *= arb(u[i] - u[j])**2
    return d
tiles = []; lo = T0; t_start = time.time()
while lo < T1 - 1e-9:
    hi = min(lo + W, T1); shift = 0
    while True:
        try:
            m, st = tmint.moments(list(range(KMAX)), lo, hi, 0.45, R2max=0.02); break
        except RuntimeError as e:
            hi += 0.05; shift += 1
    mu = [x.real for x in m]; n = int(round(float(mu[0].mid())))
    rec = {'lo': lo, 'hi': hi, 'edge_shifts': shift, 'mu0': mu[0].str(12, radius=True), 'n': n}
    if n >= 1:
        M = arb_mat(n, n, [(mu[i+j] + mu[abs(i-j)])/2 for i in range(n) for j in range(n)])
        minors = [arb_mat(k, k, [M[i, j] for i in range(k) for j in range(k)]).det() for k in range(1, n+1)]
        rec['all_real_simple'] = all(x > 0 for x in minors)
        det = minors[-1]; rec['det'] = det.str(10, radius=True)
        Q = det/disc_ref(n); rec['Q'] = float(Q.mid()); rec['Q_str'] = Q.str(8, radius=True)
        rec['centroid'] = float(((lo+hi)/2 + (hi-lo)/2*mu[1]/mu[0]).mid())
    rec['segments'] = sum(st); tiles.append(rec); print(rec, flush=True)
    lo = hi
out = {'window': [T0, T1], 'tile_width': W, 'mode': tmint.MODE, 'tiles': tiles, 'time': time.time()-t_start}
json.dump(out, open('scan_%d_%d.json' % (T0, T1), 'w'), indent=1)
