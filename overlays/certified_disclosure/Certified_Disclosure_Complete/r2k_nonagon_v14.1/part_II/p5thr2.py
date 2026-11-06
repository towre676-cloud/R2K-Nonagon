import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
rng = np.random.default_rng(92)
rec = []
for eta in (1e-9, 1e-7, 1e-5):
    for t in range(110):
        m, k, lam, w = chain(10)
        try: p, ww = pipeline(lam, w, eta); e = float(np.max(np.abs(p-lam)/lam))
        except Exception: e = np.inf
        rec.append([eta, float(w.min()/w.max()), e])
old = json.load(open('p5thr.json'))
rec = np.array(rec); x = rec[:, 1]/rec[:, 0]; fail = rec[:, 2] > 1e-3
bins = [(0, 0.1), (0.1, 1), (1, 3), (3, 10), (10, 30), (30, 100), (100, 1e3), (1e3, np.inf)]
out = {}
for lo, hi in bins:
    sel = (x >= lo) & (x < hi); out[f'[{lo},{hi})'] = [int(sel.sum()), int((fail & sel).sum())]
    print(f'x in [{lo:g},{hi:g}): trials {sel.sum():3d}  fail {(fail & sel).sum():3d}')
json.dump({'bins': out}, open('p5thr2.json', 'w'))
