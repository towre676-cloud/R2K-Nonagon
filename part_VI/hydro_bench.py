"""Protected benchmark on the published Hydrophone synthetic 5-state plant (transcribed in the review bundle).
Train / selection / confirmation splits are frozen; noise is relative complex Gaussian per component.
Methods: (L) projected Loewner [reviewer baseline], (VF) relaxed vector fitting [Gustavsen], (R2K) Loewner seed +
weighted Gauss-Newton in pole-residue form + stability enforcement + exact Routh-Hurwitz certificate."""
import json, numpy as np
from fractions import Fraction as Fr
from scipy.linalg import svd, eigvals
from scipy.optimize import linear_sum_assignment
M = np.diag([0.012, 0.008])
K0 = np.diag(np.diag(M)*(2*np.pi*np.array([70., 110.]))**2) + 300*np.array([[1, -1], [-1, 1]])
C0 = np.diag([0.03, 0.025]) + .008*np.array([[1, -1], [-1, 1]])
a = np.array([1., -.6]); B = np.array([1., .35]); ki, km, tau, gain, cs = 800., 5000., 2., 1.3, 1.
def H(s):
    s = np.asarray(s, dtype=complex); kstar = ki + km*(s*tau)/(1 + s*tau)
    D = K0[None] + s[:, None, None]*cs*C0[None] + (s*s)[:, None, None]*M[None] + kstar[:, None, None]*np.outer(a, a)[None]
    q = np.linalg.solve(D, np.tile(B, (len(s), 1))[:, :, None])[:, :, 0]; return gain*q[:, 0]
A5 = np.zeros((5, 5)); A5[:2, 2:4] = np.eye(2); A5[2:4, :2] = -np.linalg.solve(M, K0 + (ki + km)*np.outer(a, a))
A5[2:4, 2:4] = -np.linalg.solve(M, cs*C0); A5[2:4, 4] = np.linalg.solve(M, km*a); A5[4, :2] = a/tau; A5[4, 4] = -1/tau
TRUE_POLES = np.linalg.eigvals(A5)
freq = np.r_[np.geomspace(.015, 2, 24), np.geomspace(65, 240, 52)]
idx_conf = np.arange(0, 76, 5); idx_sel = np.arange(2, 76, 5); idx_tr = np.array([i for i in range(76) if i not in set(idx_conf) | set(idx_sel)])
s_all = 2j*np.pi*freq; Y_all = H(s_all)
# ---------- exact Routh-Hurwitz over Q
def routh_hurwitz_exact(coeffs_desc):
    c = [Fr(float(x)) for x in coeffs_desc]
    if c[0] < 0: c = [-x for x in c]
    if any(x <= 0 for x in c): return False
    n = len(c) - 1; r0 = c[0::2]; r1 = c[1::2]
    rows = [r0, r1 + [Fr(0)]*(len(r0) - len(r1))]
    for _ in range(n - 1):
        a0, a1 = rows[-2], rows[-1]
        if a1[0] == 0: return False
        new = [(a1[0]*a0[j+1] - a0[0]*a1[j+1])/a1[0] for j in range(len(a0) - 1)] + [Fr(0)]
        rows.append(new)
    return all(r[0] > 0 for r in rows[:n+1])
def certify_stable(poles):
    Dc = np.real(np.poly(poles))   # real coefficients of prod (s - p)
    return routh_hurwitz_exact(Dc)
# ---------- methods: each returns (predict(s), poles)
def loewner_proj(s, y, k):
    z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; mu, nu = z[::2], z[1::2]; yl, yr = yy[::2], yy[1::2]
    L = (yl[:, None] - yr[None, :])/(mu[:, None] - nu[None, :]); Ls = (mu[:, None]*yl[:, None] - nu[None, :]*yr[None, :])/(mu[:, None] - nu[None, :])
    U = svd(np.c_[L, Ls], full_matrices=False)[0][:, :k]; V = svd(np.r_[L, Ls], full_matrices=False)[2].conj().T[:, :k]
    E = -U.conj().T@L@V; A = -U.conj().T@Ls@V; Bv = U.conj().T@yl; Cv = yr@V
    return (lambda ss: np.array([Cv@np.linalg.solve(x*E - A, Bv) for x in np.atleast_1d(ss)])), eigvals(A, E)
def pr_fit(s, y, poles, w):
    """residues for fixed poles (strictly proper, real-coefficient model), weighted LS."""
    z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; ww = np.r_[w, w]
    C = 1/(z[:, None] - poles[None, :])
    r = np.linalg.lstsq(C*ww[:, None], yy*ww, rcond=None)[0]; return r
def vector_fit(s, y, k, iters=12):
    w = 1/np.abs(y); z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; ww = np.r_[w, w]
    wmax = np.max(np.abs(s)); wmin = np.min(np.abs(s))
    beta = np.geomspace(wmin, wmax, (k + 1)//2); poles = []
    for b in beta: poles += [-b/100 + 1j*b, -b/100 - 1j*b]
    poles = np.array(poles[:k]) if k % 2 == 0 else np.array(poles[:k-1] + [-wmin])
    for _ in range(iters):
        C = 1/(z[:, None] - poles[None, :])
        Amat = np.c_[C, -C*yy[:, None]]*ww[:, None]
        x = np.linalg.lstsq(Amat, yy*ww, rcond=None)[0]; ct = x[k:]
        newp = np.linalg.eigvals(np.diag(poles) - np.outer(np.ones(k), ct))
        newp = np.where(newp.real > 0, -newp.real + 1j*newp.imag, newp)
        # enforce conjugate symmetry
        newp = np.sort_complex(newp); poles = newp
    r = pr_fit(s, y, poles, w)
    return (lambda ss: (1/(np.atleast_1d(ss)[:, None] - poles[None, :]))@r), poles
def r2k_fit(s, y, k, iters=15, seed='loewner'):
    w = 1/np.abs(y)
    _, p0 = loewner_proj(s, y, k) if seed == 'loewner' else vector_fit(s, y, k)
    p = np.where(p0.real > 0, -p0.real + 1j*p0.imag, p0)
    r = pr_fit(s, y, p, w)
    z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; ww = np.r_[w, w]
    lam = 1e-3
    def res(p, r): return ((1/(z[:, None] - p[None, :]))@r - yy)*ww
    cur = np.linalg.norm(res(p, r))
    for _ in range(iters):
        C = 1/(z[:, None] - p[None, :]); J = np.c_[C, r*C**2]*ww[:, None]
        rv = res(p, r); rr = np.r_[rv.real, rv.imag]
        JJ = np.r_[np.c_[J.real, -J.imag], np.c_[J.imag, J.real]]          # params [Re theta, Im theta], theta = (r, p)
        g = np.linalg.solve(JJ.T@JJ + lam*np.eye(4*k), -JJ.T@rr)
        dth = g[:2*k] + 1j*g[2*k:]; dr = dth[:k]; dp = dth[k:]
        pn = p + dp; pn = np.where(pn.real > 0, -pn.real + 1j*pn.imag, pn); rn = pr_fit(s, y, pn, w)
        new = np.linalg.norm(res(pn, rn))
        if new < cur: p, r, cur, lam = pn, rn, new, lam/3
        else: lam *= 10
    return (lambda ss: (1/(np.atleast_1d(ss)[:, None] - p[None, :]))@r), p
METHODS = {'loewner_proj': lambda s, y, k: loewner_proj(s, y, k), 'vector_fit': vector_fit, 'r2k': r2k_fit, 'vf_gn': lambda s, y, k: r2k_fit(s, y, k, seed='vf')}
def run(noise, trials=120, seed=0):
    rng = np.random.default_rng(seed); out = {m: [] for m in METHODS}
    for t in range(trials):
        n = noise*np.abs(Y_all)*(rng.normal(size=76) + 1j*rng.normal(size=76)); Yn = Y_all + n   # one frozen realization per trial
        for name, fn in METHODS.items():
            best = None
            for k in range(3, 9):
                try:
                    pred, poles = fn(s_all[idx_tr], Yn[idx_tr], k)
                    sel = np.median(np.abs(pred(s_all[idx_sel]) - Yn[idx_sel])/np.abs(Yn[idx_sel]))
                    if np.isfinite(sel) and (best is None or sel < best[0]): best = (sel, k, pred, poles)
                except Exception: pass
            if best is None: out[name].append(None); continue
            sel, k, pred, poles = best
            rc = np.abs(pred(s_all[idx_conf]) - Y_all[idx_conf])/np.abs(Y_all[idx_conf])
            stable_claim = bool(np.all(poles.real < 0)); cert = certify_stable(poles) if stable_claim else False
            if k == 5:
                ii, jj = linear_sum_assignment(np.abs(poles[:, None] - TRUE_POLES[None, :])); perr = float(np.max(np.abs(poles[ii] - TRUE_POLES[jj])/np.abs(TRUE_POLES[jj])))
            else: perr = None
            out[name].append({'k': k, 'conf_median': float(np.median(rc)), 'conf_max': float(np.max(rc)), 'stable_float': stable_claim, 'stable_certified': cert, 'pole_relerr': perr})
    summ = {}
    for name, rows in out.items():
        ok = [r for r in rows if r is not None]
        med = np.array([r['conf_median'] for r in ok]); mx = np.array([r['conf_max'] for r in ok])
        pe = [r['pole_relerr'] for r in ok if r['pole_relerr'] is not None]
        summ[name] = {'fits': len(ok), 'median_of_median_rel_err': float(np.median(med)), 'p90_of_median_rel_err': float(np.percentile(med, 90)),
                      'median_of_max_rel_err': float(np.median(mx)), 'frac_order5_selected': float(np.mean([r['k'] == 5 for r in ok])),
                      'median_pole_relerr_when_k5': float(np.median(pe)) if pe else None, 'frac_certified_stable': float(np.mean([r['stable_certified'] for r in ok])),
                      'false_stability_claims': int(sum(1 for r in ok if r['stable_float'] and not r['stable_certified']))}
    return summ
if __name__ == '__main__':
    res = {'splits': {'train': len(idx_tr), 'selection': len(idx_sel), 'confirmation': len(idx_conf)}, 'true_poles': [[float(p.real), float(p.imag)] for p in TRUE_POLES]}
    for noise in (0.0, 0.001, 0.01, 0.05):
        res[str(noise)] = run(noise, trials=1 if noise == 0 else 120, seed=int(noise*1e5) + 7)
        print(noise, json.dumps(res[str(noise)], indent=0), flush=True)
    json.dump(res, open('hydro_bench.json', 'w'), indent=1)
