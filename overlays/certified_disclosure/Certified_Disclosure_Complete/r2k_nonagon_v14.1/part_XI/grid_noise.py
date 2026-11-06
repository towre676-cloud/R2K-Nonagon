"""XI.2: noisy one-port disclosure of MATPOWER grids (binary64), with the resolution law as an a-priori per-edge certificate.
Production schema of v9: common strength t (and 2t) per node, t at each pair; t = 1/median(G_jj) from a pilot."""
import sys, json, time
sys.path.insert(0, 'p11')
import numpy as np, gridlib
case = sys.argv[1]; levels = [0.0, 1e-12, 1e-9, 1e-6]; rng = np.random.default_rng(11)
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
K = gridlib.K_float(n, edges, ground); G = np.linalg.inv(K)
dg = np.diag(G).copy(); ok = dg > 0
C = np.where(np.outer(ok, ok), G/np.sqrt(np.abs(np.outer(dg, dg))), 0); port = int(np.argmax(np.where(ok, C.min(axis=1), -1)))
P = [port] + [i for i in range(n) if i != port]; K = K[np.ix_(P, P)]; G = G[np.ix_(P, P)]
E = {}
for (a, b), w in edges.items():
    i, j = P.index(a), P.index(b); E[(min(i, j), max(i, j))] = float(w)
t = 1/np.median(np.diag(G)[1:])
g0 = G[0, 1:]; gd = np.diag(G)[1:]; h0 = G[0, 0]
d1c = t*g0**2/(1 + t*gd); d2c = 2*t*g0**2/(1 + 2*t*gd)
iu = np.triu_indices(n - 1, 1)
a = 1/t + gd; A = a[iu[0]]; Bv = a[iu[1]]; X = G[1:, 1:][iu]; p = g0[iu[0]]; q = g0[iu[1]]
dpc = (p**2*Bv + q**2*A - 2*p*q*X)/(A*Bv - X**2)
m = 1 + 2*(n - 1) + len(iu[0])
wmin = min(abs(w) for w in E.values())
# resolution-law prediction per edge (far endpoint v = smaller port compliance): sd ~ c sigma h0 G_vv / G_0v^2, c = 63 (Table X.2.2)
CAL = json.load(open('p11/law_calib_case118.json'))['deg_far']; C_LAW = 2*CAL['ratio_max']   # calibrated on IEEE 118 only; held out elsewhere
def predicted_ok(sig):
    se = sig if sig > 0 else 1.1e-16; res = {}
    for (i, j), w in E.items():
        gi, gj = G[0, i], G[0, j]; v = i if gi < gj else j
        q = h0*G[v, v]*abs(K[v, v])/(min(gi, gj)**2*abs(w))
        res[(i, j)] = 3*C_LAW*se*q < 0.01
    return res
rows = []
for sig in levels:
    t0 = time.time(); s = sig*h0
    nz = lambda k: rng.standard_normal(k)*s
    hh0 = h0 + nz(1)[0]; d1 = d1c + nz(n - 1)*0 - nz(n - 1) + (hh0 - h0); d2 = d2c - nz(n - 1) + (hh0 - h0); dp = dpc - nz(len(iu[0])) + (hh0 - h0)
    with np.errstate(all='ignore'):
        inv_g = (1/d1 - 1/d2)/(1/t - 1/(2*t)); gg = 1/inv_g; Gjj = (1/d1 - inv_g/t)*gg
        R = np.zeros((n, n)); R[0, 0] = hh0; g0r = np.sqrt(np.abs(gg)); R[0, 1:] = R[1:, 0] = g0r; R[1:, 1:][np.diag_indices(n - 1)] = Gjj
        aa = 1/t + Gjj; Ar = aa[iu[0]]; Br = aa[iu[1]]
        disc = np.clip((dp*Ar - gg[iu[0]])*(dp*Br - gg[iu[1]]), 0, None)
        Xr = (g0r[iu[0]]*g0r[iu[1]] - np.sqrt(disc))/dp
        S = R[1:, 1:]; S[iu] = Xr; S[(iu[1], iu[0])] = Xr
        try: Kr = np.linalg.inv(R); finite = bool(np.all(np.isfinite(Kr)))
        except np.linalg.LinAlgError: Kr = None; finite = False
    if finite:
        thr = wmin/2; rec = {(i, j) for i in range(n) for j in range(i + 1, n) if -Kr[i, j] > thr}
        tp = len(rec & set(E)); prec = tp/max(len(rec), 1); recl = tp/len(E)
        errs = np.array([abs(-Kr[i, j] - w)/abs(w) for (i, j), w in E.items()])
        correct = {e: abs(-Kr[e] - w)/abs(w) < 0.01 for e, w in E.items()}
        pred = predicted_ok(sig); agree = np.mean([pred[e] == correct[e] for e in E]); safe = sum(1 for e in E if pred[e] and not correct[e])
        rows.append({'sigma_rel': sig, 'finite': True, 'precision': prec, 'recall': recl, 'median_edge_rel_err': float(np.median(errs)), 'frac_edges_within_1pct': float(np.mean([correct[e] for e in E])),
                     'law_predicted_resolvable': int(sum(pred.values())), 'law_agreement': float(agree), 'law_false_promises': int(safe), 'seconds': time.time() - t0})
    else: rows.append({'sigma_rel': sig, 'finite': False, 'seconds': time.time() - t0})
    print(case, rows[-1], flush=True)
rhos = np.sort(g0**2/(h0*gd))
out = {'law': 'relative sd(w_e) <= C sigma h0 G_vv K_vv / (G_0v^2 w_e), v = far endpoint, C = 2 x max ratio on IEEE 118', 'C_law': C_LAW, 'case': case, 'n_hidden': n, 'edges': len(E), 'readings': m, 'negative_weight_edges': sum(1 for w in E.values() if w < 0), 'common_strength_t': t,
       'port_correlation_sq_quantiles': {q: float(np.quantile(rhos, q)) for q in (0, 0.01, 0.1, 0.5)}, 'rows': rows}
json.dump(out, open('p11/grid_noise_%s.json' % case, 'w'), indent=1)
