"""XII.4 at scale: incremental row echelon mod p = 2^31-1 (numpy int64), G = K^{-1} mod p by FLINT nmod_mat."""
import sys, json, time
sys.path.insert(0, 'p11')
import numpy as np, gridlib
from flint import nmod_mat
from fractions import Fraction as Fr
P = 2147483647
case = sys.argv[1]; nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes); t0 = time.time()
def mp(x): x = Fr(x); return x.numerator % P * pow(x.denominator % P, -1, P) % P
Kp = [[0]*n for _ in range(n)]
for (i, j), w in edges.items():
    v = mp(w); Kp[i][i] = (Kp[i][i] + v) % P; Kp[j][j] = (Kp[j][j] + v) % P; Kp[i][j] = (Kp[i][j] - v) % P; Kp[j][i] = (Kp[j][i] - v) % P
for i, g in enumerate(ground): Kp[i][i] = (Kp[i][i] + mp(g)) % P
Gm = nmod_mat(Kp, P).inv(); G = np.array([[int(Gm[i, j]) for j in range(n)] for i in range(n)], dtype=np.int64); t_inv = time.time() - t0
support = sorted(edges); gsup = [i for i in range(n) if ground[i] != 0]; unk = [('w', e) for e in support] + [('g', i) for i in gsup]; U = len(unk)
col = {}
for k, (kind, x) in enumerate(unk):
    if kind == 'w': col.setdefault(x[0], []).append((k, x[1], 1)); col.setdefault(x[1], []).append((k, x[0], 1))
    else: col.setdefault(x, []).append((k, None, 0))
adj = {i: set() for i in range(n)}
for a, c in support: adj[a].add(c); adj[c].add(a)
def row(s, i):
    r = np.zeros(U + 1, dtype=np.int64); r[U] = int(i == s)
    for k, other, isw in col.get(i, []): r[k] = (G[i, s] - G[other, s]) % P if isw else G[i, s]
    return r
basis = {}  # pivot col -> normalized row
plan = []
def add(r):
    r = r.copy()
    for c, b in basis.items():
        if r[c]: r = (r - r[c]*b) % P
    nz = np.nonzero(r[:U])[0]
    if len(nz) == 0:
        assert r[U] == 0, 'inconsistent equation'
        return False
    c = int(nz[0]); inv = pow(int(r[c]), -1, P); basis[c] = (r*inv) % P; return True
for i in range(n): add(row(0, i)); plan.append((0, i))
rank0 = len(basis); readings = 1 + 2*(n - 1); pairs = set(); srcs = 0
for s in sorted(range(1, n), key=lambda v: len(adj[v])):
    if len(basis) == U: break
    if add(row(s, s)): plan.append((s, s)); srcs += 1; pairs |= {(min(s, j), max(s, j)) for j in adj[s] - {0}}
r1 = len(basis)
for s in range(1, n):
    if len(basis) == U: break
    gained = False
    for i in sorted({s} | adj[s]):
        if add(row(s, i)): gained = True; plan.append((s, i))
    if gained:
        srcs += 1; pairs |= {(min(s, j), max(s, j)) for i in {s} | adj[s] for j in adj[i] | {i} if j not in (0, s)}
readings += len(pairs)
# back-substitution mod p and comparison with the true weights reduced mod p
x = np.zeros(U, dtype=np.int64); piv = sorted(basis, reverse=True)
for c in piv:
    b = basis[c]; v = (b[U] - int(np.dot(b[:U] % P, x) % P)) % P  # b[c] = 1
    v = (int(b[U]) - sum(int(b[k])*int(x[k]) for k in np.nonzero(b[:U])[0] if k != c)) % P
    x[c] = v
truth = np.array([mp(edges[e]) for e in support] + [mp(ground[i]) for i in gsup], dtype=np.int64)
match = bool(len(basis) == U and np.all(x == truth))
out = {'case': case, 'n': n, 'unknowns': U, 'rank_column0': rank0, 'rank_after_single_rows': r1, 'rank_final_mod_p': len(basis),
       'identifiable_exact': len(basis) == U, 'recovered_equals_truth_mod_p': match, 'sources': srcs, 'readings': readings, 'full_protocol_readings': n*(n + 1)//2,
       'saving_factor': n*(n + 1)/2/readings, 'readings_per_unknown': readings/U, 'seconds_inverse_mod_p': t_inv, 'seconds_total': time.time() - t0}
out['plan'] = [[int(a), int(b)] for a, b in plan]; out['pairs'] = sorted([[int(a), int(b)] for a, b in pairs])
print(json.dumps({k: v for k, v in out.items() if k not in ('plan', 'pairs')}), flush=True); json.dump(out, open('p12/fastbig_%s.json' % case, 'w'), indent=1)
