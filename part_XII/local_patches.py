"""XII.2: known-support disclosure with LOCAL patches. Column 0 of G comes from the 2(n-1)+1 single-node readings. An extra source s
contributes only the rows i in a ball B_R(s); row i of K g_s = e_s needs G_js for j in N[i], so the pair readings needed are
(s, j) for j in B_{R+1}(s) minus {port, s}. Greedy selection by rank gain per reading (ranks mod p = 2^61-1), then exact
verification over Q of the recovered network. Also: the bridge-visibility count (ungrounded pendant blocks)."""
import sys, json, time, random
sys.path.insert(0, 'p11')
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpq, nmod_mat
import gridlib
PR = 2**61 - 1
case = sys.argv[1]; GM = sys.argv[2] if len(sys.argv) > 2 else 'known'
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
def fq(x): return fmpq(x.numerator, x.denominator)
def build(E, Gr):
    K = fmpq_mat(n, n)
    for (i, j), w in E.items(): K[i, i] += fq(w); K[j, j] += fq(w); K[i, j] -= fq(w); K[j, i] -= fq(w)
    for i, g in enumerate(Gr): K[i, i] += fq(g)
    return K
def modp(M): return nmod_mat([[int(M[i, j].p)*pow(int(M[i, j].q), -1, PR) % PR for j in range(M.ncols())] for i in range(M.nrows())], PR)
support = sorted(edges); gsup = [i for i in range(n) if ground[i] != 0] if GM == 'known' else list(range(n))
unk = [('w', e) for e in support] + [('g', i) for i in gsup]; U = len(unk)
adj = {i: set() for i in range(n)}
for a, c in support: adj[a].add(c); adj[c].add(a)
def ball(s, R):
    B = {s}; fr = {s}
    for _ in range(R):
        fr = {v for u in fr for v in adj[u]} - B; B |= fr
    return B
def rows_for(Gv, s, rows):
    out = []
    for i in rows:
        row = [0]*U
        for k, (kind, x) in enumerate(unk):
            if kind == 'w':
                a, c = x
                if i == a: row[k] = (Gv(a, s) - Gv(c, s)) % PR
                elif i == c: row[k] = (Gv(c, s) - Gv(a, s)) % PR
            elif x == i: row[k] = Gv(i, s)
        out.append(row)
    return out
def run(E, Gr, label):
    t0 = time.time(); K = build(E, Gr); Gp = modp(K).inv(); Gv = lambda i, j: int(Gp[i, j])
    rows = rows_for(Gv, 0, range(n)); r = nmod_mat(rows, PR).rank(); readings = 1 + 2*(n - 1); chosen = [('port column', 0, None, readings)]
    used = set()
    while r < U:
        best = None
        for s in range(1, n):
            for R in (0, 1, 2, 3):
                if (s, R) in used: continue
                B = ball(s, R); need = ball(s, R + 1) - {0, s}; cost = len(need)
                cand = rows + rows_for(Gv, s, sorted(B)); r2 = nmod_mat(cand, PR).rank()
                gain = r2 - r
                if gain > 0 and (best is None or gain/cost > best[0]): best = (gain/cost, s, R, cand, r2, cost)
        if best is None: break
        _, s, R, rows, r, cost = best; used.add((s, R)); readings += cost; chosen.append(('source', nodes[s], R, cost))
    # exact verification over Q on an independent subset of rows
    exact = None; zeros = []
    if r == U:
        A = nmod_mat(rows, PR); piv = []; Rr = A.transpose().rref()[0]
        for i in range(Rr.nrows()):
            for j in range(Rr.ncols()):
                if int(Rr[i, j]) != 0: piv.append(j); break
        G = K.inv(); srcs = [0] + [nodes.index(c[1]) for c in chosen[1:]]
        # rebuild the same rows exactly, in the same order
        qrows = []; qb = []
        def qrow(s, i):
            row = [fmpq(0)]*U
            for k, (kind, x) in enumerate(unk):
                if kind == 'w':
                    a, c = x
                    if i == a: row[k] = G[a, s] - G[c, s]
                    elif i == c: row[k] = G[c, s] - G[a, s]
                elif x == i: row[k] = G[i, s]
            return row
        for i in range(n): qrows.append(qrow(0, i)); qb.append(fmpq(int(i == 0)))
        for c in chosen[1:]:
            s = nodes.index(c[1])
            for i in sorted(ball(s, c[2])): qrows.append(qrow(s, i)); qb.append(fmpq(int(i == s)))
        Aq = fmpq_mat([qrows[i] for i in piv]); bq = fmpq_mat([[qb[i]] for i in piv]); x = Aq.solve(bq)
        Erec = {e: x[k, 0] for k, (kind, e) in enumerate(unk) if kind == 'w'}
        Grec = {i: x[len(support) + k, 0] for k, i in enumerate(gsup)}
        full = fmpq_mat(qrows) * x == fmpq_mat([[v] for v in qb])
        exact = full and all(Erec[e] == fq(E.get(e, Fr(0))) for e in support) and all(Grec.get(i, fmpq(0)) == fq(Gr[i]) for i in range(n))
        zeros = [[nodes[a], nodes[c]] for (a, c) in support if Erec[(a, c)] == 0]
    res = {'label': label, 'unknowns': U, 'rank': r, 'sources': len(chosen) - 1, 'readings': readings, 'full_protocol_readings': n*(n + 1)//2,
           'saving_factor': n*(n + 1)/2/readings, 'readings_per_unknown': readings/U, 'exact_recovery_and_all_equations': exact,
           'outaged_lines_detected': zeros, 'chosen': chosen, 'seconds': time.time() - t0}
    print(json.dumps({k: v for k, v in res.items() if k != 'chosen'}), flush=True); return res
# bridge-visibility count: bridges whose far side (away from port) holds no grounding
def bridges():
    disc = {}; low = {}; t = [0]; out = []
    sys.setrecursionlimit(10000)
    def dfs(u, pe):
        disc[u] = low[u] = t[0]; t[0] += 1
        for v in adj[u]:
            if v not in disc:
                dfs(v, u); low[u] = min(low[u], low[v])
                if low[v] > disc[u]: out.append((u, v))
            elif v != pe: low[u] = min(low[u], disc[v])
    dfs(0, None); return out
def side(u, v):  # nodes reachable from v without crossing (u, v)
    S = {v}; st = [v]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in S and not (x == v and y == u): S.add(y); st.append(y)
    return S
br = bridges(); blind = [(u, v) for u, v in br if all(ground[x] == 0 for x in side(u, v))]
out = {'case': case, 'ground_mode': GM, 'n': n, 'bridges': len(br), 'bridges_with_ungrounded_far_side': len(blind)}
out['nominal'] = run(edges, ground, 'nominal')
rng = random.Random(7)
for trial in range(2):
    drop = rng.sample([e for e in support if e not in {(min(u, v), max(u, v)) for u, v in br}], 3)
    E2 = {e: w for e, w in edges.items() if e not in drop}
    r_ = run(E2, ground, 'outage %s' % [[nodes[a], nodes[c]] for a, c in drop]); r_['true'] = [[nodes[a], nodes[c]] for a, c in drop]
    out.setdefault('outages', []).append(r_)
json.dump(out, open('p12/local_%s_%s.json' % (case, GM), 'w'), indent=1, default=str)
print(json.dumps({k: v for k, v in out.items() if k in ('bridges', 'bridges_with_ungrounded_far_side')}))
