"""XII.3: interval certificate for the 399-reading known-support scheme on IEEE 118 (same reading-error model as Part XI §3).
Readings are balls of radius eps*G00. From them, in Arb: G_0j and G_jj (two strengths per node), G_sj for each chosen pair
(Woodbury pair formula), then the selected square linear system for (branch weights, grounds) solved by arb_mat.solve."""
import sys, json, time
sys.path.insert(0, 'p11')
import numpy as np, gridlib
from flint import arb, arb_mat, ctx
ctx.prec = 128
case = 'case118'; L = json.load(open('p12/local_case118_known.json')); chosen = L['nominal']['chosen']
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
K = gridlib.K_float(n, edges, ground); G = np.linalg.inv(K); h0 = G[0, 0]
support = sorted(edges); gsup = [i for i in range(n) if ground[i] != 0]; unk = [('w', e) for e in support] + [('g', i) for i in gsup]; U = len(unk)
adj = {i: set() for i in range(n)}
for a, c in support: adj[a].add(c); adj[c].add(a)
def ball(s, R):
    B = {s}; fr = {s}
    for _ in range(R): fr = {v for u in fr for v in adj[u]} - B; B |= fr
    return B
t1, t2, tp = 1.0, 2.0, 1.0
rows_plan = [(0, i) for i in range(n)] + [(nodes.index(c[1]), i) for c in chosen[1:] for i in sorted(ball(nodes.index(c[1]), c[2]))]
needed_pairs = sorted({(min(s, j), max(s, j)) for c in chosen[1:] for s in [nodes.index(c[1])] for j in ball(s, c[2] + 1) - {0, s}})
# choose an independent square subsystem once (float pivoted QR on the exact-noise-free rows)
def row(Gf, s, i):
    r = np.zeros(U)
    for k, (kind, x) in enumerate(unk):
        if kind == 'w':
            a, c = x
            if i == a: r[k] = Gf(a, s) - Gf(c, s)
            elif i == c: r[k] = Gf(c, s) - Gf(a, s)
        elif x == i: r[k] = Gf(i, s)
    return r
A0 = np.array([row(lambda a, b: G[a, b], s, i) for s, i in rows_plan])
import scipy.linalg as sla
_, _, piv = sla.qr(A0.T, pivoting=True); sel = sorted(piv[:U])
truth = np.array([float(edges[e]) for e in support] + [float(ground[i]) for i in gsup])
readings = 1 + 2*(n - 1) + len(needed_pairs)
out = {'readings': readings, 'rows_available': len(rows_plan), 'unknowns': U, 'rows': []}
for rel in (1e-15, 1e-13, 1e-11, 1e-9):
    t0 = time.time(); eps = rel*h0; rng = np.random.default_rng(3)
    rd = lambda x: arb(float(x + rng.uniform(-eps, eps)), eps)
    H0 = rd(h0); Gb = {(0, 0): H0}
    for j in range(1, n):
        d1 = H0 - rd(h0 - t1*G[0, j]**2/(1 + t1*G[j, j])); d2 = H0 - rd(h0 - t2*G[0, j]**2/(1 + t2*G[j, j]))
        ig = (1/d1 - 1/d2)/(1/arb(t1) - 1/arb(t2)); g = 1/ig; Gb[(j, j)] = (1/d1 - ig/t1)*g; Gb[(0, j)] = Gb[(j, 0)] = g.sqrt()
    for (i, j) in needed_pairs:
        a = 1/tp + G[i, i]; b = 1/tp + G[j, j]; x = G[i, j]; p = G[0, i]; q = G[0, j]
        dl = H0 - rd(h0 - (p*p*b + q*q*a - 2*p*q*x)/(a*b - x*x))
        A = 1/arb(tp) + Gb[(i, i)]; B = 1/arb(tp) + Gb[(j, j)]; disc = (dl*A - Gb[(0, i)]**2)*(dl*B - Gb[(0, j)]**2)
        if disc < 0: disc = arb(0)
        elif not disc > 0: disc = arb(0).union(disc.upper())
        Gb[(i, j)] = Gb[(j, i)] = (Gb[(0, i)]*Gb[(0, j)] - disc.sqrt())/dl
    def Gf(a, b):
        return Gb.get((a, b), Gb.get((b, a)))
    try:
        M = arb_mat([[v for v in row_arb] for row_arb in ([[(Gf(rows_plan[r][1] if False else 0, 0))]] if False else [])]) if False else None
        Ar = []; br = []
        for r_ in sel:
            s, i = rows_plan[r_]; rr = [arb(0)]*U
            for k, (kind, x) in enumerate(unk):
                if kind == 'w':
                    a, c = x
                    if i == a: rr[k] = Gf(a, s) - Gf(c, s)
                    elif i == c: rr[k] = Gf(c, s) - Gf(a, s)
                elif x == i: rr[k] = Gf(i, s)
            Ar.append(rr); br.append([arb(int(i == s))])
        X = arb_mat(Ar).solve(arb_mat(br)); ok = True
    except Exception as e: ok = False; err = str(e)[:80]
    if ok:
        w = [X[k, 0] for k in range(U)]
        present = sum(1 for k in range(len(support)) if w[k] > 0); contains = all(w[k].contains(arb(truth[k])) for k in range(U))
        maxrel = max(float(w[k].rad())/abs(truth[k]) for k in range(U))
        out['rows'].append({'eps_rel': rel, 'branches_certified_present': present, 'branches': len(support), 'max_rel_radius': maxrel, 'truth_contained': contains, 'seconds': time.time() - t0})
    else: out['rows'].append({'eps_rel': rel, 'solve_failed': err})
    print(out['rows'][-1], flush=True)
json.dump(out, open('p12/local_cert_case118.json', 'w'), indent=1)
