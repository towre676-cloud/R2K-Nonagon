import numpy as np, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
rng = np.random.default_rng(78)
for eta in (1e-9, 1e-6):
    rows = []
    for t in range(60):
        m, k, lam, w = chain(10)
        # stage-1 only diagnostics
        n = 10; N = 100
        def sample(z): h = Hs(z, lam, w); return h*(1 + eta*(rng.standard_normal(len(z))+1j*rng.standard_normal(len(z)))/np.sqrt(2))
        om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, N); mu = (om+0.05j)**2; nu = mu.conj()
        p1 = projected(mu, nu, sample(mu), sample(nu), n)
        e1 = np.max(np.abs(p1-lam)/lam)
        p, ww = pipeline(lam, w, eta); e = np.max(np.abs(p-lam)/lam)
        # weakest mode visibility: residue relative to |H| at nearest sweep point
        vis = min(w[j]/0.05/np.max(np.abs(Hs(mu, lam, w))) for j in range(n))
        rows.append((e, e1, w.min()/w.max(), vis, np.min(np.diff(lam)/lam[1:])))
    rows = np.array(rows); bad = rows[:, 0] > 1e-3
    print('eta', eta, 'fail', bad.sum(), '/60')
    print('  stage1 err   fail med %.1e  ok med %.1e' % (np.median(rows[bad, 1]) if bad.any() else np.nan, np.median(rows[~bad, 1])))
    print('  wmin/wmax    fail med %.1e  ok med %.1e' % (np.median(rows[bad, 2]) if bad.any() else np.nan, np.median(rows[~bad, 2])))
    print('  rel min gap  fail med %.1e  ok med %.1e' % (np.median(rows[bad, 4]) if bad.any() else np.nan, np.median(rows[~bad, 4])))
