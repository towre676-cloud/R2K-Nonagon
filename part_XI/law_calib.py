"""Calibrate per-edge error against candidate resolution indices on IEEE 118 at sigma = 1e-9 (several noise draws)."""
import sys, json
sys.path.insert(0, 'p11'); import numpy as np, gridlib
src = open('p11/grid_noise.py').read()
case = sys.argv[1]; sys.argv = ['x', case]
exec(src.split('rows = []')[0])
sig = 1e-9; errs = []
for rep in range(20):
    s = sig*h0; nz = lambda k: rng.standard_normal(k)*s
    hh0 = h0 + nz(1)[0]; d1 = d1c - nz(n-1) + (hh0-h0); d2 = d2c - nz(n-1) + (hh0-h0); dp = dpc - nz(len(iu[0])) + (hh0-h0)
    inv_g = (1/d1 - 1/d2)/(1/t - 1/(2*t)); gg = 1/inv_g; Gjj = (1/d1 - inv_g/t)*gg
    R = np.zeros((n, n)); R[0, 0] = hh0; g0r = np.sqrt(np.abs(gg)); R[0, 1:] = R[1:, 0] = g0r; R[1:, 1:][np.diag_indices(n-1)] = Gjj
    aa = 1/t + Gjj; Ar = aa[iu[0]]; Br = aa[iu[1]]; disc = np.clip((dp*Ar - gg[iu[0]])*(dp*Br - gg[iu[1]]), 0, None)
    Xr = (g0r[iu[0]]*g0r[iu[1]] - np.sqrt(disc))/dp; S = R[1:, 1:]; S[iu] = Xr; S[(iu[1], iu[0])] = Xr
    Kr = np.linalg.inv(R); errs.append([(-Kr[e] - w)/w for e, w in E.items()])
sd = np.std(np.array(errs), axis=0)
idx = {}
for k, ((i, j), w) in enumerate(E.items()):
    gi, gj = G[0, i], G[0, j]; v = i if gi < gj else j; u = j if v == i else i
    idx.setdefault('far', []).append(h0*G[v, v]/(gi if v == i else gj)**2/w)
    idx.setdefault('far_noGvv', []).append(h0/(min(gi, gj))**2/w)
    idx.setdefault('both', []).append(h0*(1/gi**2 + 1/gj**2)/w*max(G[i, i], G[j, j]))
    idx.setdefault('deg_far', []).append(h0*G[v, v]/min(gi, gj)**2/w*abs(K[v, v]))
out = {}
for k, q in idx.items():
    q = np.array(q)*sig; r = sd/q; c = np.corrcoef(np.log(sd), np.log(q))[0, 1]
    out[k] = {'log_corr': float(c), 'ratio_median': float(np.median(r)), 'ratio_p95': float(np.quantile(r, .95)), 'ratio_max': float(r.max())}
print(case, json.dumps(out, indent=1)); json.dump(out, open('p11/law_calib_%s.json' % case, 'w'), indent=1)
