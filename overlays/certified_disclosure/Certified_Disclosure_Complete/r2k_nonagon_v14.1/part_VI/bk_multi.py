"""BK/ratios number variance at three disjoint heights (zeros from index 1e5, 1e6, 1e7), height-averaged prediction,
prime-cutoff and quadrature sensitivity, block standard errors."""
import numpy as np, mpmath as mp, json, time, sys
from flint import acb
mp.mp.dps = 20
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].astype(float)
PALL = sieve(2_000_000)
def make_R2(pmax):
    Pp = PALL[PALL <= pmax]; lp = np.log(Pp)
    def A(eta): x = Pp**(-(1+eta)); t = (1 - x)*(1 - 2/Pp + x)/(1 - 1/Pp)**2; return complex(np.exp(np.sum(np.log(t.astype(complex)))))
    def B(eta): return complex(np.sum((lp/(Pp**(1+eta) - 1))**2))
    cache = {}
    def R2(eps, T):
        Lg = np.log(T/(2*np.pi)); r = 2*np.pi*eps/Lg; s = mp.mpc(1, r)
        if r not in cache: cache[r] = complex(mp.diff(lambda z: mp.zeta(z, derivative=1)/mp.zeta(z), s))
        term = cache[r] + complex(mp.power(T/(2*mp.pi), -1j*r)*mp.zeta(s)*mp.zeta(mp.conj(s)))*A(1j*r) - B(1j*r)
        return 1 + 2/Lg**2*term.real
    return R2
xs = np.linspace(1e-4, 16, 1601); Ls = [0.5, 1, 2, 4, 8, 16]
def sigma2(L, x, r):
    sel = x <= L; xx = x[sel]; yy = 1 - r[sel]; return float(L - 2*np.trapezoid((L - xx)*yy, xx))
def s2(pts, L, step=0.01):
    lo = np.arange(pts[0], pts[-1]-L, step); c = np.searchsorted(pts, lo+L) - np.searchsorted(pts, lo); return float(np.var(c))
R2main = make_R2(2_000_000); R2low = make_R2(200_000)
out = {'heights': []}
for n0 in [int(x) for x in sys.argv[1].split(',')]:
    t0 = time.time(); Z = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(n0, 20000)]); tz = time.time() - t0
    U = Z/(2*np.pi)*np.log(Z/(2*np.pi*np.e)) + 7/8; U -= U[0]
    Ts = [float(np.quantile(Z, q)) for q in (0.1, 0.5, 0.9)]
    curves = {}
    for lab, T, R in (('T10', Ts[0], R2main), ('T50', Ts[1], R2main), ('T90', Ts[2], R2main), ('T50_p2e5', Ts[1], R2low)):
        curves[lab] = np.array([R(x, T) for x in xs])
    rows = []
    for L in Ls:
        emp = s2(U, L); blocks = [s2(U[i*2000:(i+1)*2000], L) for i in range(10)]; se = float(np.std(blocks, ddof=1)/np.sqrt(10))
        bk = sigma2(L, xs, curves['T50']); bkavg = float(np.mean([sigma2(L, xs, curves[k]) for k in ('T10', 'T50', 'T90')]))
        bkp = sigma2(L, xs, curves['T50_p2e5']); bkq = sigma2(L, xs[::2], curves['T50'][::2])
        gue = sigma2(L, xs, 1 - np.sinc(xs)**2)
        rows.append({'L': L, 'zeta': emp, 'se': se, 'BK_T50': bk, 'BK_height_avg': bkavg, 'prime_cutoff_shift': bkp - bk, 'quadrature_shift': bkq - bk,
                     'GUE': gue, 'z_BK_avg': (emp - bkavg)/se, 'frac_GUE_departure_explained': (gue - bkavg)/(gue - emp)})
    out['heights'].append({'n0': n0, 'T_quantiles': Ts, 'zeros_time_s': tz, 'rows': rows})
    print(n0, json.dumps(rows), flush=True)
    json.dump(out, open('bk_multi_%s.json' % sys.argv[1].replace(',', '_'), 'w'), indent=1)
