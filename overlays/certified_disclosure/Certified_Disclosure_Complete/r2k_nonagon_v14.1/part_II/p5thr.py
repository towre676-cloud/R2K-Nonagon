import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
rng = np.random.default_rng(91)
rec = []
for eta in (1e-9, 1e-6):
    for t in range(120):
        m, k, lam, w = chain(10)
        try: p, ww = pipeline(lam, w, eta); e = float(np.max(np.abs(p-lam)/lam))
        except Exception: e = np.inf
        rec.append([eta, float(w.min()/w.max()), e])
rec = np.array(rec); x = rec[:, 1]/rec[:, 0]; fail = rec[:, 2] > 1e-3
out = {}
for lo, hi, nm in ((0, 1, 'x<1'), (1, 100, '1<=x<100'), (100, np.inf, 'x>=100')):
    sel = (x >= lo) & (x < hi); out[nm] = [int(sel.sum()), int((fail & sel).sum()), float(np.median(rec[sel & ~fail, 2]/rec[sel & ~fail, 0])) if (sel & ~fail).any() else None]
    print(nm, 'trials', sel.sum(), 'fail', (fail & sel).sum(), 'median pole/eta on success', out[nm][2])
json.dump(out, open('p5thr.json', 'w'))
