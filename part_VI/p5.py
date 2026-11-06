"""P5: conditioning-optimal sampling morphisms for Löwner identification of chains."""
import numpy as np, scipy.linalg as sla, json
rng = np.random.default_rng(55)
def chain(n):
    m = rng.uniform(0.5, 2.0, n); k = rng.uniform(1.0, 5.0, n)
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] = k[i] + (k[i+1] if i+1 < n else 0)
        if i+1 < n: K[i, i+1] = K[i+1, i] = -k[i+1]
    M = np.diag(m); lam, Q = sla.eigh(K, M)
    w = Q[0]**2  # eigh(K,M) gives M-orthonormal Q; H = sum Q0k^2/(lam-z)
    return m, k, lam, w
def Hs(z, lam, w): return np.sum(w/(lam[None, :] - np.asarray(z)[:, None]), 1)
def lanczos(lv, wv):
    n = len(lv); Q = np.zeros((n, n)); a = np.zeros(n); b = np.zeros(n-1)
    Q[:, 0] = np.sqrt(wv/wv.sum()); prev = np.zeros(n); be = 0
    for j in range(n):
        v = lv*Q[:, j] - be*prev; a[j] = Q[:, j] @ v; v -= a[j]*Q[:, j]
        v -= Q[:, :j+1] @ (Q[:, :j+1].T @ v); v -= Q[:, :j+1] @ (Q[:, :j+1].T @ v)
        if j < n-1: be = np.linalg.norm(v); b[j] = be; prev = Q[:, j]; Q[:, j+1] = v/be
    return a, b
def recover(lam, w):
    a, b = lanczos(lam, w); m = [1/w.sum()]; k = [1/np.sum(w/lam)]
    for i in range(len(lam)-1):
        kn = a[i]*m[i] - k[i]; k.append(kn); m.append(kn**2/(b[i]**2*m[i]))
    return np.array(m), np.array(k)
def identify(mu, nu, Hm, Hn, polish=True, equil=True):
    n = len(mu)
    L = (Hm[:, None]-Hn[None, :])/(mu[:, None]-nu[None, :]); Ls = (mu[:, None]*Hm[:, None]-nu[None, :]*Hn[None, :])/(mu[:, None]-nu[None, :])
    c0 = np.linalg.cond(L)
    if equil:
        for _ in range(30):
            dr = 1/np.sqrt(np.max(np.abs(L), 1)); L = dr[:, None]*L; Ls = dr[:, None]*Ls
            dc = 1/np.sqrt(np.max(np.abs(L), 0)); L = L*dc; Ls = Ls*dc
    c1 = np.linalg.cond(L)
    p = np.sort(sla.eigvals(Ls, L).real)
    zs = np.r_[mu, nu]; H = np.r_[Hm, Hn]
    C = 1/(p[None, :]-zs[:, None]); w = np.linalg.lstsq(np.r_[C.real, C.imag], np.r_[H.real, H.imag], rcond=None)[0]
    if polish:
        for _ in range(8):
            C = 1/(p[None, :]-zs[:, None]); r = C@w - H
            J = np.c_[C, -w*C**2]; d = np.linalg.lstsq(np.r_[J.real, J.imag], -np.r_[r.real, r.imag], rcond=None)[0]
            w += d[:n]; p += d[n:]
    o = np.argsort(p); return p[o], w[o], c0, c1
designs = ['ring', 'sweep', 'cheb_sweep', 'resonance', 'balanced']
def run_design(dname, lam, w, noise=0.0, delta_frac=0.1):
    n = len(lam)
    def sample(z):
        h = Hs(z, lam, w); return h*(1 + noise*(rng.standard_normal(len(z)) + 1j*rng.standard_normal(len(z)))/np.sqrt(2))
    if dname == 'ring':
        th = np.linspace(0, 2*np.pi, 2*n, endpoint=False)+0.1; c, r = (lam[0]+lam[-1])/2, (lam[-1]-lam[0])/2*1.35
        z = c + r*np.exp(1j*th); mu, nu = z[0::2], z[1::2]
    elif dname in ('sweep', 'resonance', 'balanced'):
        om = np.linspace(np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1, n); mu = (om+0.05j)**2; nu = mu.conj()
    elif dname == 'cheb_sweep':
        a_, b_ = np.sqrt(lam[0])*0.8, np.sqrt(lam[-1])*1.1
        om = (a_+b_)/2 + (b_-a_)/2*np.cos(np.pi*(np.arange(n)+0.5)/n)[::-1]; mu = (om+0.05j)**2; nu = mu.conj()
    p, ww, c0, c1 = identify(mu, nu, sample(mu), sample(nu))
    if dname == 'resonance':  # second pass at resonance-centred samples
        gaps = np.diff(p); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
        dl = delta_frac*mg; mu = p + 1j*dl; nu = p - 1j*dl
        p, ww, c0, c1 = identify(mu, nu, sample(mu), sample(nu))
    if dname == 'balanced':  # residue-balanced resonance sampling, two refinement passes
        for _ in range(2):
            gaps = np.diff(p); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps])
            sw = np.sqrt(np.abs(ww)); g = np.min(mg/sw); dl = delta_frac*g*sw
            mu = p + 1j*dl; nu = p - 1j*dl
            p, ww, c0, c1 = identify(mu, nu, sample(mu), sample(nu))
    return p, ww, c0, c1
res = {d: {'pole': [], 'mass': [], 'cond0': [], 'cond1': []} for d in designs}
for trial in range(200):
    m, k, lam, w = chain(10)
    for d in designs:
        try:
            p, ww, c0, c1 = run_design(d, lam, w)
            mr, kr = recover(p, ww)
            res[d]['pole'].append(float(np.max(np.abs(p-lam)/lam))); res[d]['mass'].append(float(np.max(np.abs(mr-m)/m)))
            res[d]['cond0'].append(float(c0)); res[d]['cond1'].append(float(c1))
        except Exception as e:
            res[d]['pole'].append(np.inf); res[d]['mass'].append(np.inf); res[d]['cond0'].append(np.inf); res[d]['cond1'].append(np.inf)
summ = {}
for d in designs:
    summ[d] = {q: [float(np.nanmedian(np.array(res[d][q]))), float(np.max(np.array(res[d][q])))] for q in res[d]}
    print(d, {q: ['%.2e' % x for x in summ[d][q]] for q in summ[d]})
# noise robustness: pole and residue errors, recursion failures
noise_tab = {}
for nz in (1e-12, 1e-9, 1e-6):
    rows = {d: [] for d in designs}; rw = {d: [] for d in designs}
    for trial in range(60):
        m, k, lam, w = chain(10)
        for d in designs:
            try:
                p, ww, _, _ = run_design(d, lam, w, noise=nz)
                rows[d].append(float(np.max(np.abs(p-lam)/lam))); rw[d].append(float(np.max(np.abs(ww-w)/w)))
            except Exception:
                rows[d].append(np.inf); rw[d].append(np.inf)
    noise_tab[nz] = {d: [float(np.median(rows[d])), float(np.median(rw[d])), float(np.percentile(rows[d], 90))] for d in designs}
    print('noise', nz, {d: 'pole %.1e res %.1e p90 %.1e' % tuple(noise_tab[nz][d]) for d in designs}, flush=True)
# delta sweep for resonance design: Gershgorin quantity rho and achieved cond
dsw = []
m, k, lam, w = chain(10)
for df in (0.01, 0.03, 0.1, 0.3, 1.0):
    p, ww, c0, c1 = run_design('resonance', lam, w, delta_frac=df)
    gaps = np.diff(lam); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps]); dl = df*mg
    rho = max(sum(dl[i]/abs(lam[i]-lam[j]) for j in range(10) if j != i) for i in range(10))
    mu = lam + 1j*dl; S = (1j*dl)[:, None]/(mu[:, None]-lam[None, :])
    dsw.append([df, float(rho), float(np.linalg.cond(S)), float((1+rho)/(1-rho)) if rho < 1 else None, float(c1), float(np.max(np.abs(p-lam)/lam))])
    print('delta_frac %.2f rho %.3f condS %.3f gersh %s condL %.2e poleerr %.1e' % (df, rho, dsw[-1][2], dsw[-1][3], c1, dsw[-1][5]))
bsw = []
for df in (0.01, 0.03, 0.1, 0.3):
    p, ww, c0, c1 = run_design('balanced', lam, w, delta_frac=df)
    bsw.append([df, float(c1), float(np.max(np.abs(p-lam)/lam))]); print('balanced theta %.2f condL %.2e poleerr %.1e' % tuple(bsw[-1]))
json.dump({'balanced_sweep': bsw, 'summary': summ, 'noise': {str(k): v for k, v in noise_tab.items()}, 'delta_sweep': dsw}, open('p5.json', 'w'), indent=1)
