import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
rng = np.random.default_rng(94)
rec = []
for eta in (1e-9, 1e-7, 1e-6, 1e-5):
    for t in range(90):
        m, k, lam, w = chain(10); n = 10; N = 100
        om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, N); mu = (om+0.05j)**2; nu = mu.conj()
        Hm, Hn = Hs(mu, lam, w), Hs(nu, lam, w)
        nm = eta*(rng.standard_normal(N)+1j*rng.standard_normal(N))/np.sqrt(2); nn = eta*(rng.standard_normal(N)+1j*rng.standard_normal(N))/np.sqrt(2)
        Hm2, Hn2 = Hm*(1+nm), Hn*(1+nn)
        L = lambda a, b: np.c_[(a[:, None]-b[None, :])/(mu[:, None]-nu[None, :]), (mu[:, None]*a[:, None]-nu[None, :]*b[None, :])/(mu[:, None]-nu[None, :])]
        Lt, Ln = L(Hm, Hn), L(Hm2, Hn2)
        sn = np.linalg.svd(Lt, compute_uv=False)[n-1]; E = np.linalg.norm(Ln - Lt, 2)
        # run pipeline with the same noisy pass-1 data replaced: reuse pipeline (fresh noise) for failure label
        try: p, ww = pipeline(lam, w, eta); e = float(np.max(np.abs(p-lam)/lam))
        except Exception: e = np.inf
        rec.append([eta, float(w.min()/w.max()), float(sn/E), e])
rec = np.array(rec); fail = rec[:, 3] > 1e-3; rho = rec[:, 2]; x = rec[:, 1]/rec[:, 0]
best = max(((float(th), float(np.mean((rho < th) == fail))) for th in np.logspace(-2, 2, 400)), key=lambda z: z[1])
bestx = max(((float(th), float(np.mean((x < th) == fail))) for th in np.logspace(-2, 3, 400)), key=lambda z: z[1])
bins = []
for lo, hi in ((0, 0.3), (0.3, 1), (1, 2), (2, 4), (4, 10), (10, 1e12)):
    sel = (rho >= lo) & (rho < hi); bins.append([lo, hi, int(sel.sum()), int((fail & sel).sum())])
out = {'n': int(len(rec)), 'fails': int(fail.sum()), 'best_rho_threshold': best, 'best_x_threshold': bestx, 'rho_bins': bins,
       'acc_rho_at_1': float(np.mean((rho < 1) == fail)), 'acc_rho_at_2': float(np.mean((rho < 2) == fail))}
print(json.dumps(out, indent=1)); json.dump(out, open('m4.json', 'w'), indent=1)
