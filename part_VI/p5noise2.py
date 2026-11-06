import numpy as np, scipy.linalg as sla, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5.py').read().split('designs = [')[0])
rng = np.random.default_rng(78)
def gn(p, w, zs, H, wt, it=10):
    n = len(p)
    for _ in range(it):
        C = 1/(p[None, :]-zs[:, None]); r = (C@w - H)*wt; J = np.c_[C, -w*C**2]*wt[:, None]
        d = np.linalg.lstsq(np.r_[J.real, J.imag], -np.r_[r.real, r.imag], rcond=None)[0]; w = w + d[:n]; p = p + d[n:]
    o = np.argsort(p); return p[o], w[o]
def projected(mu, nu, Hm, Hn, n):
    L = (Hm[:, None]-Hn[None, :])/(mu[:, None]-nu[None, :]); Ls = (mu[:, None]*Hm[:, None]-nu[None, :]*Hn[None, :])/(mu[:, None]-nu[None, :])
    U, S, Vh = np.linalg.svd(np.c_[L, Ls]); Y = U[:, :n]
    U2, S2, Vh2 = np.linalg.svd(np.r_[L, Ls]); X = Vh2[:n].conj().T
    return np.sort(sla.eigvals(Y.conj().T@Ls@X, Y.conj().T@L@X).real)
def pipeline(lam, w, eta, N=100, theta=0.1, passes=4):
    n = len(lam)
    def sample(z): h = Hs(z, lam, w); return h*(1 + eta*(rng.standard_normal(len(z))+1j*rng.standard_normal(len(z)))/np.sqrt(2))
    om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, N); mu = (om+0.05j)**2; nu = mu.conj()
    Hm, Hn = sample(mu), sample(nu); p = projected(mu, nu, Hm, Hn, n)
    zs = np.r_[mu, nu]; H = np.r_[Hm, Hn]
    C = 1/(p[None, :]-zs[:, None]); wt = 1/np.abs(H)
    w0 = np.linalg.lstsq(np.r_[(C*wt[:, None]).real, (C*wt[:, None]).imag], np.r_[(H*wt).real, (H*wt).imag], rcond=None)[0]
    p, ww = gn(p, w0, zs, H, wt)
    for _ in range(passes):
        C = 1/(p[None, :]-zs[:, None]); J = np.c_[C, -ww*C**2]*(wt/eta)[:, None]; Jr = np.r_[J.real, J.imag]
        sg = np.sqrt(np.abs(np.diag(np.linalg.pinv(Jr.T@Jr))))[n:]
        gaps = np.diff(p); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
        sw = np.sqrt(np.abs(ww)); g = np.min(mg/sw); dl = np.minimum(np.maximum(theta*g*sw, 5*sg), theta*mg)
        znew = np.r_[p+1j*dl, p-1j*dl]; Hnew = sample(znew)
        zs = np.r_[zs, znew]; H = np.r_[H, Hnew]; wt = 1/np.abs(H)
        p, ww = gn(p, ww, zs, H, wt)
    return p, ww
out = {}
for eta in (1e-12, 1e-9, 1e-6, 1e-4, 1e-3):
    errs = []; rres = []
    for t in range(60):
        m, k, lam, w = chain(10)
        try: p, ww = pipeline(lam, w, eta); errs.append(np.max(np.abs(p-lam)/lam)); rres.append(np.max(np.abs(ww-w)/w))
        except Exception as e: errs.append(np.inf); rres.append(np.inf)
    out[eta] = [float(np.median(errs)), float(np.median(rres)), float(np.percentile(errs, 90)), float(np.max(errs))]
    print('eta %.0e  pole med %.1e  res med %.1e  p90 %.1e  max %.1e  (pole med/eta %.2g)' % (eta, *out[eta], out[eta][0]/eta), flush=True)
json.dump({str(k): v for k, v in out.items()}, open('p5noise2.json', 'w'))
