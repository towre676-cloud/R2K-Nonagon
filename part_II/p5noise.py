import numpy as np, json, warnings; warnings.filterwarnings('ignore')
exec(open('p5.py').read().split('designs = [')[0])
rng = np.random.default_rng(77)
def sigma_hat(p, w, zs, Hv, eta):
    C = 1/(p[None, :]-zs[:, None]); J = np.c_[C, -w*C**2]; Jr = np.r_[J.real, J.imag]
    scale = np.r_[np.abs(Hv), np.abs(Hv)]*eta/np.sqrt(2)
    Jw = Jr/scale[:, None]; cov = np.linalg.pinv(Jw.T@Jw); return np.sqrt(np.abs(np.diag(cov)))[len(p):]
def floored(lam, w, eta, theta=0.1, c=5.0, passes=5, tag=''):
    n = len(lam)
    def sample(z): h = Hs(z, lam, w); return h*(1 + eta*(rng.standard_normal(len(z))+1j*rng.standard_normal(len(z)))/np.sqrt(2))
    om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, n); mu = (om+0.05j)**2; nu = mu.conj()
    Hm, Hn = sample(mu), sample(nu); p, ww, _, _ = identify(mu, nu, Hm, Hn)
    zs = np.r_[mu, nu]; Hv = np.r_[Hm, Hn]
    for _ in range(passes):
        sg = sigma_hat(p, ww, zs, Hv, eta)
        gaps = np.diff(p); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
        sw = np.sqrt(np.abs(ww)); g = np.min(mg/sw); dl = np.minimum(np.maximum(theta*g*sw, c*sg), theta*mg)
        mu = p+1j*dl; nu = p-1j*dl; Hm, Hn = sample(mu), sample(nu); zs = np.r_[mu, nu]; Hv = np.r_[Hm, Hn]
        p, ww, _, c1 = identify(mu, nu, Hm, Hn)
    return p, ww
out = {}
for eta in (1e-12, 1e-9, 1e-6, 1e-4):
    errs = []; rres = []
    for t in range(60):
        m, k, lam, w = chain(10)
        try:
            p, ww = floored(lam, w, eta); errs.append(np.max(np.abs(p-lam)/lam)); rres.append(np.max(np.abs(ww-w)/w))
        except Exception: errs.append(np.inf); rres.append(np.inf)
    out[eta] = [float(np.median(errs)), float(np.median(rres)), float(np.percentile(errs, 90))]
    print('eta %.0e  pole med %.1e  res med %.1e  pole p90 %.1e  (pole/eta = %.1f)' % (eta, *out[eta], out[eta][0]/eta))
json.dump({str(k): v for k, v in out.items()}, open('p5noise.json', 'w'))
