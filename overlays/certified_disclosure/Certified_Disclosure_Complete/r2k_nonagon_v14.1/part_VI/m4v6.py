"""M4 repaired: paired data (rho computed from the exact pass-1 samples the pipeline consumes), operational
(non-oracle) predictors, threshold fitted on a training half and scored on a held-out half with Wilson intervals."""
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5noise2.py').read().split('out = {}')[0])
def pipeline_paired(lam, w, eta, Hm, Hn, mu, nu, rng, theta=0.1, passes=4):
    n = len(lam)
    def sample(z): h = Hs(z, lam, w); return h*(1 + eta*(rng.standard_normal(len(z))+1j*rng.standard_normal(len(z)))/np.sqrt(2))
    p = projected(mu, nu, Hm, Hn, n); zs = np.r_[mu, nu]; H = np.r_[Hm, Hn]
    C = 1/(p[None, :]-zs[:, None]); wt = 1/np.abs(H)
    w0 = np.linalg.lstsq(np.r_[(C*wt[:, None]).real, (C*wt[:, None]).imag], np.r_[(H*wt).real, (H*wt).imag], rcond=None)[0]
    p, ww = gn(p, w0, zs, H, wt)
    for _ in range(passes):
        C = 1/(p[None, :]-zs[:, None]); J = np.c_[C, -ww*C**2]*(wt/eta)[:, None]; Jr = np.r_[J.real, J.imag]
        sg = np.sqrt(np.abs(np.diag(np.linalg.pinv(Jr.T@Jr))))[n:]
        gaps = np.diff(p); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
        sw = np.sqrt(np.abs(ww)); g = np.min(mg/sw); dl = np.minimum(np.maximum(theta*g*sw, 5*sg), theta*mg)
        znew = np.r_[p+1j*dl, p-1j*dl]; Hnew = sample(znew); zs = np.r_[zs, znew]; H = np.r_[H, Hnew]; wt = 1/np.abs(H)
        p, ww = gn(p, ww, zs, H, wt)
    return p
rng = np.random.default_rng(2026); rec = []
n = 10; N = 100
for eta in (1e-9, 1e-7, 1e-6, 1e-5):
    for t in range(100):
        m, k, lam, w = chain(n)
        om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, N); mu = (om+0.05j)**2; nu = mu.conj()
        Hm, Hn = Hs(mu, lam, w), Hs(nu, lam, w)
        cn = lambda: eta*(rng.standard_normal(N)+1j*rng.standard_normal(N))/np.sqrt(2)
        Hm2, Hn2 = Hm*(1+cn()), Hn*(1+cn())
        L = lambda a, b: np.c_[(a[:, None]-b[None, :])/(mu[:, None]-nu[None, :]), (mu[:, None]*a[:, None]-nu[None, :]*b[None, :])/(mu[:, None]-nu[None, :])]
        Lt, Ln = L(Hm, Hn), L(Hm2, Hn2)
        st = np.linalg.svd(Lt, compute_uv=False); sn_ = np.linalg.svd(Ln, compute_uv=False)
        E = np.linalg.norm(Ln - Lt, 2); rho_oracle = st[n-1]/E
        # operational: parametric bootstrap of the noise norm from the noisy data and the known eta (no access to truth)
        Eb = np.median([np.linalg.norm(L(Hm2*(1+cn()), Hn2*(1+cn())) - Ln, 2) for _ in range(8)])
        rho_op = sn_[n-1]/Eb; gapr = sn_[n-1]/sn_[n]
        try: p = pipeline_paired(lam, w, eta, Hm2, Hn2, mu, nu, rng); e = float(np.max(np.abs(p-lam)/lam))
        except Exception: e = np.inf
        rec.append([eta, rho_oracle, rho_op, gapr, e, st[n-1]/(2*E)])
rec = np.array(rec); fail = rec[:, 4] > 1e-3
idx = rng.permutation(len(rec)); tr, te = idx[:len(rec)//2], idx[len(rec)//2:]
def wilson(k, m, z=1.96):
    p = k/m; d = 1 + z*z/m; c = (p + z*z/(2*m))/d; h = z*np.sqrt(p*(1-p)/m + z*z/(4*m*m))/d; return [float(c-h), float(c+h)]
out = {'n': int(len(rec)), 'fails': int(fail.sum()), 'fail_def': 'max relative pole error > 1e-3', 'predictors': {}}
for j, name in ((1, 'rho_oracle'), (2, 'rho_operational_bootstrap'), (3, 'singular_gap_ratio')):
    x = rec[:, j]; grid = np.logspace(-2, 4, 600)
    th = max(grid, key=lambda g: np.mean((x[tr] < g) == fail[tr]))
    k = int(np.sum((x[te] < th) == fail[te])); m_ = len(te)
    fn = int(np.sum((x[te] >= th) & fail[te])); fp = int(np.sum((x[te] < th) & ~fail[te]))
    out['predictors'][name] = {'threshold_train': float(th), 'test_accuracy': k/m_, 'test_accuracy_wilson95': wilson(k, m_),
                               'test_missed_failures': fn, 'test_false_alarms': fp, 'test_n': m_}
# Weyl-Wedin one-sided sufficient condition: ||E|| < sigma_n/2 guarantees the numerical rank-n subspace is separated; check it never coexists with failure
weyl = rec[:, 5] > 1
out['weyl_sufficient'] = {'cases_condition_holds': int(weyl.sum()), 'failures_among_them': int((weyl & fail).sum()),
                          'failures_where_condition_fails': int((~weyl & fail).sum()), 'successes_where_condition_fails': int((~weyl & ~fail).sum())}
print(json.dumps(out, indent=1)); json.dump(out, open('m4v6.json', 'w'), indent=1)
