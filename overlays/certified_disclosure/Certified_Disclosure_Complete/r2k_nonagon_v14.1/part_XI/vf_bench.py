"""XI.5: head-to-head with Vector Fitting (scikit-rf 2.2.0 implementation of Gustavsen-Semlyen with relaxation) on the measured
one-port and two-port responses shipped with scikit-rf. Frozen splits: every 5th frequency confirmation, 2::5 order selection, rest training.
R2K = projected Loewner seed + variable-projection refinement of poles (residues and constant eliminated) + pole reflection +
exact Routh-Hurwitz certificate over Q. Passivity of reflection models: max |S(j w)| on a dense grid (bounded-real necessary check)."""
import os, sys, json, time, warnings
from fractions import Fraction as Fr
from pathlib import Path as _P; _R = _P(__file__).resolve().parents[1]; sys.path[:0] = [str(_R/'part_X'), str(_R/'part_VI')]
import lossless_passivity as LP
import numpy as np, skrf
from skrf.vectorFitting import VectorFitting
from scipy.optimize import least_squares
from scipy.linalg import svd, eigvals
sys.path.insert(0, '.')
from hydro_bench import routh_hurwitz_exact
warnings.filterwarnings('ignore')
D = os.path.dirname(skrf.data.__file__)
sets = [('ring slot measured.s1p', 0, 0), ('ro,1.s1p', 0, 0), ('ro,2.s1p', 0, 0), ('ro,3.s1p', 0, 0), ('ntwk1.s2p', 0, 0), ('ntwk1.s2p', 1, 0)]
def split(N):
    conf = np.arange(0, N, 5); sel = np.arange(2, N, 5); tr = np.array([i for i in range(N) if i not in set(conf) | set(sel)]); return tr, sel, conf
def basis(s, poles):
    return np.c_[1/(s[:, None] - poles[None, :]), np.ones(len(s))]
def fit_residues(s, y, poles):
    A = basis(s, poles); z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; Az = basis(z, poles)
    c = np.linalg.lstsq(Az, yy, rcond=None)[0]; return c
def model(s, poles, c): return basis(s, poles)@c
def loewner_poles(s, y, k):
    z = np.r_[s, s.conj()]; yy = np.r_[y, y.conj()]; mu, nu = z[::2], z[1::2]; yl, yr = yy[::2], yy[1::2]
    L = (yl[:, None] - yr[None, :])/(mu[:, None] - nu[None, :]); Ls = (mu[:, None]*yl[:, None] - nu[None, :]*yr[None, :])/(mu[:, None] - nu[None, :])
    U = svd(np.c_[L, Ls], full_matrices=False)[0][:, :k]; V = svd(np.r_[L, Ls], full_matrices=False)[2].conj().T[:, :k]
    p = eigvals(-U.conj().T@Ls@V, -U.conj().T@L@V); p = p[np.isfinite(p)]
    p = np.where(p.real > 0, -p.real + 1j*p.imag, p); return p
def conj_closed(p):
    """split into real poles and upper-half complex poles; return parameter vector"""
    re = sorted([x.real for x in p if abs(x.imag) < 1e-9*max(1, abs(x))]); cp = [x for x in p if x.imag > 1e-9*max(1, abs(x))]
    return np.array(re), np.array(cp)
def unpack(theta, nr):
    r = theta[:nr]; cp = theta[nr:].reshape(-1, 2); cpx = cp[:, 0] + 1j*cp[:, 1]
    return np.r_[r, cpx, cpx.conj()]
def r2k_fit(s, y, k):
    p0 = loewner_poles(s, y, k); re, cp = conj_closed(p0)
    theta0 = np.r_[re, np.c_[cp.real, cp.imag].ravel()]; nr = len(re)
    w = 1/np.maximum(np.abs(y), 1e-3*np.max(np.abs(y)))
    def res(th):
        p = unpack(th, nr); c = fit_residues(s, y, p); e = (model(s, p, c) - y)*w; return np.r_[e.real, e.imag]
    sol = least_squares(res, theta0, method='lm', max_nfev=400) if len(theta0) else None
    p = unpack(sol.x, nr) if sol is not None else p0
    p = np.where(p.real > 0, -p.real + 1j*p.imag, p); c = fit_residues(s, y, p); return p, c
def exact_bounded_real(poles, resid, d):
    """full-axis passivity of a reflection model: S bounded real  <=>  Z = (1+S)/(1-S) positive real (exact Cayley/Brune, v10.1 checker),
    coefficients taken exactly from their binary64 values"""
    den = np.real(np.poly(poles))[::-1]; num = d*den.astype(complex)
    for k, r in enumerate(resid):
        num = num + r*np.r_[np.poly(np.delete(poles, k))[::-1], 0][:len(den)]
    num = np.real(num)
    N = [Fr(float(a)) + Fr(float(b)) for a, b in zip(den, num)]; Dd = [Fr(float(a)) - Fr(float(b)) for a, b in zip(den, num)]
    while N and N[-1] == 0: N.pop()
    while Dd and Dd[-1] == 0: Dd.pop()
    try: return bool(LP.certificate(N, Dd, reduce_common=True)['positive_real'])
    except Exception as e: return 'error: ' + str(e)[:60]
def certify(p):
    if len(p) == 0: return True
    return bool(routh_hurwitz_exact(np.real(np.poly(p))))
rows = []
for fn, i, j in sets:
    ntw = skrf.Network(os.path.join(D, fn)); f = ntw.f; y = ntw.s[:, i, j]; N = len(f); tr, sel, conf = split(N)
    fs = f.max(); s = 2j*np.pi*f/fs
    scale = np.sqrt(np.mean(np.abs(y[conf])**2))
    res = {'data': fn + ' S%d%d' % (i + 1, j + 1), 'points': N, 'band_GHz': [f.min()/1e9, f.max()/1e9]}
    # --- R2K
    t0 = time.time(); best = None
    for k in range(2, 21):
        try:
            p, c = r2k_fit(s[tr], y[tr], k); e = np.median(np.abs(model(s[sel], p, c) - y[sel]))
            if np.isfinite(e) and (best is None or e < best[0]): best = (e, k, p, c)
        except Exception: pass
    e, k, p, c = best; yc = model(s[conf], p, c); dense = 2j*np.pi*np.linspace(f.min(), f.max(), 4000)/fs
    res['r2k'] = {'order': len(p), 'conf_nrmse': float(np.sqrt(np.mean(np.abs(yc - y[conf])**2))/scale), 'conf_max_err_rel': float(np.max(np.abs(yc - y[conf]))/scale),
                  'stable_certified_exact': certify(p), 'max_abs_S_dense': float(np.max(np.abs(model(dense, p, c)))), 'seconds': time.time() - t0}
    t1 = time.time(); res['r2k']['passive_exact_full_axis'] = exact_bounded_real(p, c[:-1], c[-1]) if fn.endswith('.s1p') or (i == j) else 'n/a (transmission)'; res['r2k']['passivity_check_seconds'] = time.time() - t1
    # --- scikit-rf Vector Fitting: order grid, selection on the same split
    t0 = time.time(); bestv = None
    sub = skrf.Network(frequency=skrf.Frequency.from_f(f[tr], unit='hz'), s=y[tr][:, None, None])
    for nreal in (0, 1, 2):
        for ncx in range(1, 11):
            try:
                vf = VectorFitting(sub); vf.vector_fit(n_poles_real=nreal, n_poles_cmplx=ncx)
                ys = vf.get_model_response(0, 0, f[sel]); e = np.median(np.abs(ys - y[sel]))
                if np.isfinite(e) and (bestv is None or e < bestv[0]): bestv = (e, vf, nreal + 2*ncx)
            except Exception: pass
    e, vf, order = bestv; yc = vf.get_model_response(0, 0, f[conf]); poles = vf.poles
    pfull = np.r_[poles, poles[np.abs(poles.imag) > 0].conj()]
    try: passive = bool(vf.is_passive())
    except Exception: passive = None
    res['vector_fitting'] = {'order': int(order), 'conf_nrmse': float(np.sqrt(np.mean(np.abs(yc - y[conf])**2))/scale), 'conf_max_err_rel': float(np.max(np.abs(yc - y[conf]))/scale),
                             'stable_certified_exact': certify(pfull/(2*np.pi*fs)), 'max_abs_S_dense': float(np.max(np.abs(vf.get_model_response(0, 0, np.linspace(f.min(), f.max(), 4000))))),
                             'skrf_is_passive': passive, 'seconds': time.time() - t0}
    rr = vf.residues[0]; rfull = np.r_[rr, rr[np.abs(poles.imag) > 0].conj()]
    res['vector_fitting']['passive_exact_full_axis'] = exact_bounded_real(pfull/(2*np.pi*fs), rfull/(2*np.pi*fs), float(np.real(vf.constant_coeff[0]))) if (i == j) else 'n/a (transmission)'
    rows.append(res); print(json.dumps(res), flush=True)
json.dump(rows, open('p11/vf_bench.json', 'w'), indent=1)
