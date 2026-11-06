import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
rng = np.random.default_rng(93)
def crb_snr(lam, w, eta, N=100):
    om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, N); z = np.r_[(om+0.05j)**2, ((om+0.05j)**2).conj()]
    C = 1/(lam[None, :] - z[:, None]); H = C@w
    J = np.c_[C, w*C**2]                       # d/dw, d/dlam of sum w/(lam - z)
    sig = eta*np.abs(H)/np.sqrt(2); Jr = np.r_[J.real/sig[:, None], J.imag/sig[:, None]]
    cov = np.linalg.pinv(Jr.T@Jr); sw = np.sqrt(np.abs(np.diag(cov)))[:len(w)]; sl = np.sqrt(np.abs(np.diag(cov)))[len(w):]
    snr = w/sw; gaps = np.diff(lam); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
    return float(snr.min()), float(np.max(sl/mg))
rec = []
for eta in (1e-9, 1e-7, 1e-5, 1e-6):
    for t in range(90):
        m, k, lam, w = chain(10)
        snr, pole_frac = crb_snr(lam, w, eta)
        try: p, ww = pipeline(lam, w, eta); e = float(np.max(np.abs(p-lam)/lam))
        except Exception: e = np.inf
        rec.append([eta, float(w.min()/w.max()), snr, pole_frac, e])
rec = np.array(rec); fail = rec[:, 4] > 1e-3; x = rec[:, 1]/rec[:, 0]; snr = rec[:, 2]
def rate(v, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (v >= lo) & (v < hi); out.append([lo, hi, int(sel.sum()), int((fail & sel).sum())])
    return out
E = [0, 0.3, 1, 2, 3, 5, 10, 30, 1e9]
R = {'by_snr': rate(snr, E), 'by_x': rate(x, [0, 0.1, 1, 3, 10, 30, 100, 1e9])}
# best single threshold classifier
best = None
for thr in np.logspace(-1, 2, 200):
    acc = np.mean((snr < thr) == fail)
    if best is None or acc > best[1]: best = (float(thr), float(acc))
bestx = None
for thr in np.logspace(-2, 3, 300):
    acc = np.mean((x < thr) == fail)
    if bestx is None or acc > bestx[1]: bestx = (float(thr), float(acc))
R['best_snr_threshold'] = best; R['best_x_threshold'] = bestx; R['n'] = int(len(rec)); R['fails'] = int(fail.sum())
R['snr_over_x_median'] = float(np.median(snr/x))
print(json.dumps(R, indent=1)); json.dump(R, open('n5.json', 'w'), indent=1)
