"""Synthetic test of the harmonic-pocket bound. A random grounded host network; a pocket B = I u T planted on it.
Readings: every column s outside I in full (all rows), the port column, and the diagonal row (s,s) of every interior s.
Prediction: deficiency >= e(B) - (|T|-1)(|B|-1) - |I|, with equality expected for generic weights."""
import random, itertools, json, flint
p = 2**31 - 1; random.seed(7)
def run(shape, nhost=30, trial=0):
    rnd = random.Random(1000 * trial + hash(shape) % 997)
    T, I, Bedges = SHAPES[shape]
    nB = len(T) + len(I); n = nhost + len(I)
    # host: random connected graph on 0..nhost-1, terminals are host nodes 0..|T|-1; interior nodes nhost..
    E = {}
    for v in range(1, nhost): E[(rnd.randrange(v), v)] = rnd.randrange(1, p)
    for _ in range(nhost): a, b = sorted(rnd.sample(range(nhost), 2)); E[(a, b)] = rnd.randrange(1, p)
    lab = {('t', k): k for k in range(len(T))}; lab.update({('i', k): nhost + k for k in range(len(I))})
    pocket = set()
    for a, b in Bedges:
        x, y = sorted((lab[a], lab[b])); E[(x, y)] = rnd.randrange(1, p); pocket.add((x, y))
    gr = {v: rnd.randrange(1, p) for v in rnd.sample(range(len(T), nhost), 4)}
    K = flint.nmod_mat(n, n, p)
    for (a, b), w in E.items(): K[a, a] += w; K[b, b] += w; K[a, b] -= w; K[b, a] -= w
    for v, g in gr.items(): K[v, v] += g
    G = K.inv(); edges = sorted(E); col = {e: k for k, e in enumerate(edges)}; gv = sorted(gr); U = len(edges) + len(gv)
    inc = {v: [] for v in range(n)}
    for e in edges: inc[e[0]].append(e); inc[e[1]].append(e)
    rows = []
    def row(i, s):
        r = [0] * U
        for e in inc[i]:
            o = e[1] if e[0] == i else e[0]; r[col[e]] = (int(G[i, s]) - int(G[o, s])) % p
        if i in gr: r[len(edges) + gv.index(i)] = int(G[i, s])
        rows.append(r)
    interior = set(range(nhost, n))
    for s in range(n):
        if s in interior: row(s, s)
        else:
            for i in range(n): row(i, s)
    M = flint.nmod_mat(rows, p); d = U - M.rank()
    # local characterisation: N(B) = {v on pocket edges : L_v h = 0 for h in H_B, e_s^T L_v phi_s = 0 for s in I}
    Bn = sorted(set(range(len(T))) | interior); bi = {v: k for k, v in enumerate(Bn)}; pe = sorted(e for e in E if e[0] in bi and e[1] in bi)   # induced edges of B
    Hcols = [[int(G[b, s]) for b in Bn] for s in range(n) if s not in interior]   # restrictions g_s|B, s outside I
    II = sorted(interior); KII = flint.nmod_mat(len(II), len(II), p)
    for a in range(len(II)):
        for b in range(len(II)): KII[a, b] = int(K[II[a], II[b]])
    KIIi = KII.inv()
    cond = []
    for h in Hcols:
        for b in Bn:
            r = []
            for (x, y) in pe:   # (L_e h)_b for edge e=(x,y): weight-1 Laplacian
                v = 0
                if b == x: v = (h[bi[x]] - h[bi[y]]) % p
                elif b == y: v = (h[bi[y]] - h[bi[x]]) % p
                r.append(v)
            cond.append(r)
    for a, s_ in enumerate(II):
        phi = {II[c]: int(KIIi[c, a]) for c in range(len(II))}
        r = []
        for (x, y) in pe:
            px, py = phi.get(x, 0), phi.get(y, 0); v = 0
            if s_ == x: v = (px - py) % p
            elif s_ == y: v = (py - px) % p
            r.append(v)
        cond.append(r)
    dimN = len(pe) - flint.nmod_mat(cond, p).rank()
    bound = len(pe) - (len(T) - 1) * (nB - 1) - len(I)
    return d, max(bound, 0), dimN
def K(nodes): return list(itertools.combinations(nodes, 2))
t = lambda k: ('t', k); i_ = lambda k: ('i', k)
SHAPES = {
 'K4, |T|=2': ([0, 1], [0, 1], K([t(0), t(1), i_(0), i_(1)])),
 'triangle off a cut vertex': ([0], [0, 1], K([t(0), i_(0), i_(1)])),
 'C4 off a cut vertex': ([0], [0, 1, 2], [(t(0), i_(0)), (i_(0), i_(1)), (i_(1), i_(2)), (i_(2), t(0))]),
 'diamond, |T|=2 at the tips': ([0, 1], [0, 1], [(t(0), i_(0)), (t(0), i_(1)), (i_(0), i_(1)), (i_(0), t(1)), (i_(1), t(1))]),
 'series chain, |T|=2': ([0, 1], [0, 1], [(t(0), i_(0)), (i_(0), i_(1)), (i_(1), t(1))]),
 'K5, |T|=2': ([0, 1], [0, 1, 2], K([t(0), t(1), i_(0), i_(1), i_(2)])),
 'K5, |T|=3': ([0, 1, 2], [0, 1], K([t(0), t(1), t(2), i_(0), i_(1)])),
 'K4, |T|=3': ([0, 1, 2], [0], K([t(0), t(1), t(2), i_(0)])),
 'leaf': ([0], [0], [(t(0), i_(0))]),
}
res = {}
for sh in SHAPES:
    out = [run(sh, trial=k) for k in range(3)]
    res[sh] = {'observed_deficiency': [o[0] for o in out], 'bound': [o[1] for o in out], 'dim_N_B': [o[2] for o in out]}
    print('%-30s observed %s  dim N(B) %s  bound %s' % (sh, [o[0] for o in out], [o[2] for o in out], [o[1] for o in out]))
json.dump(res, open('pocket_synth.json', 'w'), indent=1)
