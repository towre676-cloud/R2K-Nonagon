import numpy as np, json, time, mpmath as mp
from flint import acb
from numpy.polynomial import chebyshev as Ch
def Qof(u):
    n = len(u)
    if n <= 1: return 1.0
    ub = -1 + (2*np.arange(1, n+1) - 1)/n
    ld = sum(np.log((u[i]-u[j])**2) for i in range(n) for j in range(i+1, n)); lr = sum(np.log((ub[i]-ub[j])**2) for i in range(n) for j in range(i+1, n))
    return float(np.exp(ld - lr))
# (1) consistency: rigorous scan Q vs Q from Arb zeros
S = json.load(open('scan_6990_7020.json'))['tiles']
zs = [float(z.imag.mid()) for z in acb.zeta_zeros(6690, 60)]
cons = []
for t in S:
    c, a = (t['lo']+t['hi'])/2, (t['hi']-t['lo'])/2
    u = np.array([(z - c)/a for z in zs if t['lo'] <= z < t['hi']])
    cons.append([t['lo'], t['n'], len(u), t['Q'], Qof(u)])
cons_err = max(abs(r[3]-r[4])/r[4] for r in cons)
print('scan-vs-zeros max rel diff in Q:', cons_err, ' counts equal:', all(r[1] == r[2] for r in cons))
# (2) statistics: zeta zeros near n=1e6, unfolded
t0 = time.time(); Z = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(1000000, 20000)]); tf = time.time() - t0
def Nsmooth(t): return t/(2*np.pi)*np.log(t/(2*np.pi*np.e)) + 7/8
U = Nsmooth(Z); U = U - U[0]
def tile_stats(pts, W=2.0):
    pts = np.sort(pts); L = pts[-1]; out = []
    for lo in np.arange(0.5, L - W, W):
        u = (pts[(pts >= lo) & (pts < lo+W)] - (lo + W/2))/(W/2); out.append((len(u), Qof(u)))
    return out
def cue(Nm, reps, rng):
    res = []
    for _ in range(reps):
        A = (rng.standard_normal((Nm, Nm)) + 1j*rng.standard_normal((Nm, Nm)))/np.sqrt(2)
        Qm, R = np.linalg.qr(A); Qm = Qm*(np.diag(R)/np.abs(np.diag(R)))
        th = np.sort(np.angle(np.linalg.eigvals(Qm))) % (2*np.pi); res += tile_stats(np.sort(th)*Nm/(2*np.pi))
    return res
rng = np.random.default_rng(95)
zeta_t = tile_stats(U); cue_t = cue(400, 25, rng); poi_t = tile_stats(np.cumsum(rng.exponential(1.0, 20000)))
def summarize(T, name):
    n = np.array([a for a, b in T]); q = np.array([b for a, b in T])
    d = {'tiles': len(T), 'count_dist': {int(k): float(np.mean(n == k)) for k in range(0, 6)}, 'count_var': float(np.var(n))}
    for k in (2, 3):
        sel = q[n == k]; d[f'logQ_n{k}_quantiles'] = [float(np.quantile(np.log(sel), p)) for p in (0.05, 0.25, 0.5, 0.75, 0.95)] if len(sel) else None
    d['P(Q<0.01 | n>=2)'] = float(np.mean(q[n >= 2] < 0.01)); return d, n, q
out = {'consistency': {'max_rel_diff': cons_err, 'rows': cons}, 'zeta_fetch_s': tf, 'window': [float(Z[0]), float(Z[-1])]}
from scipy.stats import ks_2samp
sz, nz, qz = summarize(zeta_t, 'zeta'); sc, nc, qc = summarize(cue_t, 'cue'); sp, np_, qp = summarize(poi_t, 'poisson')
out['zeta'] = sz; out['cue'] = sc; out['poisson'] = sp
for k in (2, 3):
    out[f'KS_n{k}_zeta_vs_cue'] = list(map(float, ks_2samp(np.log(qz[nz == k]), np.log(qc[nc == k]))))
    out[f'KS_n{k}_zeta_vs_poisson'] = list(map(float, ks_2samp(np.log(qz[nz == k]), np.log(qp[np_ == k]))))
print(json.dumps({k: v for k, v in out.items() if k != 'consistency'}, indent=1))
json.dump(out, open('m5.json', 'w'), indent=1)
