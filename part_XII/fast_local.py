"""XII.4: scaled known-support scheme. Column 0 + single-row patches (s, R=0) accepted in order of cost (degree) when they raise the
rank mod p; then R=1 patches if needed. Full rank mod p proves full rank over Q, hence exact identifiability; noiseless binary64
recovery error reported."""
import sys, json, time
sys.path.insert(0, 'p11')
import numpy as np, gridlib
from flint import nmod_mat
PR = 2**61 - 1
from fractions import Fraction as Fr
case = sys.argv[1]
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
t0 = time.time()
def mp(x): x = Fr(x); return x.numerator % PR * pow(x.denominator % PR, -1, PR) % PR
Kp = [[0]*n for _ in range(n)]
for (i, j), w in edges.items():
    v = mp(w); Kp[i][i] = (Kp[i][i] + v) % PR; Kp[j][j] = (Kp[j][j] + v) % PR; Kp[i][j] = (Kp[i][j] - v) % PR; Kp[j][i] = (Kp[j][i] - v) % PR
for i, g in enumerate(ground): Kp[i][i] = (Kp[i][i] + mp(g)) % PR
Gp = nmod_mat(Kp, PR).inv(); t_inv = time.time() - t0
support = sorted(edges); gsup = [i for i in range(n) if ground[i] != 0]; unk = [('w', e) for e in support] + [('g', i) for i in gsup]; U = len(unk)
col = {}
for k, (kind, x) in enumerate(unk):
    if kind == 'w': col.setdefault(x[0], []).append((k, x[1], 1)); col.setdefault(x[1], []).append((k, x[0], 1))
    else: col.setdefault(x, []).append((k, None, 0))
adj = {i: set() for i in range(n)}
for a, c in support: adj[a].add(c); adj[c].add(a)
def row(s, i):
    r = [0]*U
    for k, other, isw in col.get(i, []):
        r[k] = (int(Gp[i, s]) - int(Gp[other, s])) % PR if isw else int(Gp[i, s])
    return r
rows = [row(0, i) for i in range(n)]; M = nmod_mat(rows, PR); r = M.rank(); readings = 1 + 2*(n - 1); srcs = 0
for s in sorted(range(1, n), key=lambda v: len(adj[v])):
    if r == U: break
    cand = rows + [row(s, s)]; r2 = nmod_mat(cand, PR).rank()
    if r2 > r: rows, r = cand, r2; readings += len(adj[s] - {0}); srcs += 1
stage1 = r
for s in range(1, n):
    if r == U: break
    B = {s} | adj[s]; cand = rows + [row(s, i) for i in sorted(B)]; r2 = nmod_mat(cand, PR).rank()
    if r2 > r:
        need = set().union(*[adj[i] | {i} for i in B]) - {0, s}; rows, r = cand, r2; readings += len(need); srcs += 1
out = {'case': case, 'n': n, 'unknowns': U, 'rank_mod_p': r, 'full_rank_proves_identifiability': r == U, 'rank_after_single_rows': stage1, 'sources': srcs,
       'readings_upper_bound': readings, 'full_protocol_readings': n*(n + 1)//2, 'saving_factor': n*(n + 1)/2/readings, 'readings_per_unknown': readings/U,
       'seconds': time.time() - t0}
print(json.dumps(out), flush=True); json.dump(out, open('p12/fast_%s.json' % case, 'w'), indent=1)
