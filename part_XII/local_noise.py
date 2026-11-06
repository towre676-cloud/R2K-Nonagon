"""XII.5: random reading noise on the known-support scheme (binary64), same noise model and levels as Part XI Table XI.4."""
import sys, json
sys.path.insert(0, 'p11')
import numpy as np, gridlib
case = sys.argv[1]; F = json.load(open('p12/fastbig_%s.json' % case))
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
K = gridlib.K_float(n, edges, ground); G = np.linalg.inv(K); h0 = G[0, 0]
support = sorted(edges); gsup = [i for i in range(n) if ground[i] != 0]; unk = [('w', e) for e in support] + [('g', i) for i in gsup]; U = len(unk)
col = {}
for k, (kind, x) in enumerate(unk):
    if kind == 'w': col.setdefault(x[0], []).append((k, x[1], 1)); col.setdefault(x[1], []).append((k, x[0], 1))
    else: col.setdefault(x, []).append((k, None, 0))
truth = np.array([float(edges[e]) for e in support] + [float(ground[i]) for i in gsup]); wmin = min(abs(float(w)) for w in edges.values())
t = 1/np.median(np.diag(G)[1:]); rng = np.random.default_rng(21); rows_out = []
pairs = [tuple(p) for p in F['pairs']]; plan = [tuple(p) for p in F['plan']]
for sig in (0.0, 1e-12, 1e-9, 1e-6):
    s = sig*h0; nz = lambda: rng.standard_normal()*s
    H0 = h0 + nz(); g0 = np.zeros(n); gd = np.zeros(n); g0[0] = H0
    for j in range(1, n):
        d1 = H0 - (h0 - t*G[0, j]**2/(1 + t*G[j, j]) + nz()); d2 = H0 - (h0 - 2*t*G[0, j]**2/(1 + 2*t*G[j, j]) + nz())
        ig = (1/d1 - 1/d2)/(1/t - 1/(2*t)); gg = 1/ig; gd[j] = (1/d1 - ig/t)*gg; g0[j] = np.sqrt(abs(gg))
    gd[0] = H0; Gx = {}
    for (i, j) in pairs:
        a = 1/t + G[i, i]; b = 1/t + G[j, j]; x = G[i, j]; p_ = G[0, i]; q = G[0, j]
        dl = H0 - (h0 - (p_*p_*b + q*q*a - 2*p_*q*x)/(a*b - x*x) + nz())
        A = 1/t + gd[i]; B = 1/t + gd[j]; disc = max((dl*A - g0[i]**2)*(dl*B - g0[j]**2), 0.0)
        Gx[(i, j)] = (g0[i]*g0[j] - np.sqrt(disc))/dl
    def Gv(a, b):
        if a == b: return gd[a]
        if a == 0: return g0[b]
        if b == 0: return g0[a]
        return Gx.get((min(a, b), max(a, b)))
    Arows = []; brh = []
    for (s_, i) in plan:
        r = np.zeros(U)
        for k, other, isw in col.get(i, []): r[k] = (Gv(i, s_) - Gv(other, s_)) if isw else Gv(i, s_)
        Arows.append(r); brh.append(float(i == s_))
    xhat = np.linalg.lstsq(np.array(Arows), np.array(brh), rcond=None)[0]
    rel = np.abs(xhat - truth)/np.abs(truth); w = xhat[:len(support)]
    rows_out.append({'sigma_rel': sig, 'median_rel_err': float(np.median(rel[:len(support)])), 'within_1pct': float(np.mean(rel[:len(support)] < 0.01)),
                     'branches_above_gap': int(np.sum(w > wmin/2)), 'branches': len(support)})
    print(case, rows_out[-1], flush=True)
json.dump({'case': case, 'readings': F['readings'], 'rows': rows_out}, open('p12/local_noise_%s.json' % case, 'w'), indent=1)
