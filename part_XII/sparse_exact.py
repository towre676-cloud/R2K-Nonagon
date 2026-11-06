"""XII.1: target-relative disclosure with known support (PerfectJet ch.7: exact target recovery iff ker O ⊆ ker T).
Unknowns: every branch weight on the nominal support and a ground at EVERY bus. Single-node interventions give the full first
column g_0 of G = K^{-1}; pair interventions (s, j) for one extra node s give column g_s. Each known column gives n linear
equations K g_s = e_s in the unknowns. Exact rank and exact solution over Q (FLINT)."""
import sys, json, time, random, itertools
sys.path.insert(0, 'p11')
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpq, nmod_mat, nmod
PR = 2**61 - 1
import gridlib
case = sys.argv[1] if len(sys.argv) > 1 else 'case118'
nodes, edges, ground, meta = gridlib.network('p11/data/%s.m' % case); n = len(nodes)
def fq(x): return fmpq(x.numerator, x.denominator)
def build(E, Gr):
    K = fmpq_mat(n, n)
    for (i, j), w in E.items(): K[i, i] += fq(w); K[j, j] += fq(w); K[i, j] -= fq(w); K[j, i] -= fq(w)
    for i, g in enumerate(Gr): K[i, i] += fq(g)
    return K
support = sorted(edges)                       # nominal one-line diagram
GROUND_MODE = sys.argv[2] if len(sys.argv) > 2 else 'all'
gsup = [i for i in range(n) if ground[i] != 0] if GROUND_MODE == 'known' else list(range(n))
unk = [('w', e) for e in support] + [('g', i) for i in gsup]
def system(G, cols):
    """rows: (K g_s)_i = delta_is for each known column s, linear in unknowns"""
    A = []; b = []
    for s in cols:
        for i in range(n):
            row = [fmpq(0)]*len(unk)
            for k, (kind, x) in enumerate(unk):
                if kind == 'w':
                    a, c = x
                    if i == a: row[k] = G[a, s] - G[c, s]
                    elif i == c: row[k] = G[c, s] - G[a, s]
                elif x == i: row[k] = G[i, s]
            A.append(row); b.append(fmpq(int(i == s)))
    return fmpq_mat(A), fmpq_mat([[v] for v in b])
def readings_for(cols):
    # baseline + 2 strengths at each of n-1 hidden nodes (column 0) + pairs (s, j), j not in {port, s}, for each extra column
    return 1 + 2*(n - 1) + sum(n - 2 for s in cols if s != 0)
out = {'case': case, 'n': n, 'unknowns': len(unk), 'branches': len(support), 'full_protocol_readings': n*(n+1)//2}
def modp(M):
    return nmod_mat([[int(M[i, j].p) * pow(int(M[i, j].q), -1, PR) % PR for j in range(M.ncols())] for i in range(M.nrows())], PR)
def run(E, Gr, label, rng):
    t0 = time.time(); K = build(E, Gr); Kp = modp(K); Gp = Kp.inv()
    class Wrap:  # index access mod p
        def __init__(s_, M): s_.M = M
        def __getitem__(s_, ij): return int(s_.M[ij[0], ij[1]])
    def system_p(cols):
        A = []
        for s_ in cols:
            for i in range(n):
                row = [0]*len(unk)
                for k, (kind, x) in enumerate(unk):
                    if kind == 'w':
                        a, c = x
                        if i == a: row[k] = (int(Gp[a, s_]) - int(Gp[c, s_])) % PR
                        elif i == c: row[k] = (int(Gp[c, s_]) - int(Gp[a, s_])) % PR
                    elif x == i: row[k] = int(Gp[i, s_])
                A.append(row)
        return nmod_mat(A, PR)
    cols = [0]; A = system_p(cols); r0 = r = A.rank()
    N = A.nullspace()[0]; nb = A.ncols() - r0
    rec0 = [k for k in range(len(unk)) if all(int(N[k, j]) == 0 for j in range(nb))]
    rec0_lines = sum(1 for k in rec0 if unk[k][0] == 'w'); rec0_grounds = sum(1 for k in rec0 if unk[k][0] == 'g')
    trace = [[1, r]]
    while r < len(unk):
        best = None
        for s_ in range(1, n):
            if s_ in cols: continue
            A2 = system_p(cols + [s_]); r2 = A2.rank()
            if best is None or r2 > best[1]: best = (s_, r2, A2)
            if r2 == len(unk): break
        if best[1] == r: break
        cols.append(best[0]); A, r = best[2], best[1]; trace.append([len(cols), r])
    exact = None; zeros_found = []
    if r == len(unk):
        # exact rational verification: independent rows chosen mod p, square system solved over Q
        piv = []; R_ = A.transpose().rref()[0]
        for i in range(R_.nrows()):
            for j in range(R_.ncols()):
                if int(R_[i, j]) != 0: piv.append(j); break
        G = K.inv(); Aq, bq = system(G, cols)
        Aq2 = fmpq_mat([[Aq[i, j] for j in range(Aq.ncols())] for i in piv]); bq2 = fmpq_mat([[bq[i, 0]] for i in piv])
        x = Aq2.solve(bq2)
        Erec = {e: x[k, 0] for k, (kind, e) in enumerate(unk) if kind == 'w'}; Grec = {i: x[len(support) + k, 0] for k, i in enumerate(gsup)}
        exact = all(Erec[e] == fq(E.get(e, Fr(0))) for e in support) and all(Grec.get(i, fmpq(0)) == fq(Gr[i]) for i in range(n))
        resid_ok = (Aq*x == bq)
        zeros_found = [e for e in support if Erec[e] == 0]; exact = exact and resid_ok
    res = {'label': label, 'rank_column0': r0, 'unknowns_recoverable_from_column0': len(rec0), 'lines_recoverable_from_column0': rec0_lines,
           'grounds_recoverable_from_column0': rec0_grounds, 'columns_needed': len(cols), 'rank_trace': trace,
           'readings': readings_for(cols), 'full_rank': r == len(unk), 'exact_recovery_and_all_equations': exact,
           'outaged_lines_detected': [[nodes[a], nodes[c]] for a, c in zeros_found], 'seconds': time.time() - t0}
    print(json.dumps(res), flush=True); return res
rng = random.Random(12)
out['nominal'] = run(edges, ground, 'nominal network', rng)
# outage scenario: 3 random in-service lines removed (not disconnecting the grid), recovered against the nominal support
for trial in range(3):
    while True:
        drop = rng.sample(support, 3); E2 = {e: w for e, w in edges.items() if e not in drop}
        adj = {i: set() for i in range(n)}
        for a, c in E2: adj[a].add(c); adj[c].add(a)
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in adj[u] - seen: seen.add(v); st.append(v)
        if len(seen) == n: break
    r_ = run(E2, ground, 'outage of %s' % [[nodes[a], nodes[c]] for a, c in drop], rng)
    r_['true_outages'] = [[nodes[a], nodes[c]] for a, c in drop]; out.setdefault('outages', []).append(r_)
out['ground_mode'] = GROUND_MODE
json.dump(out, open('p12/sparse_exact_%s_%s.json' % (case, GROUND_MODE), 'w'), indent=1, default=str)
