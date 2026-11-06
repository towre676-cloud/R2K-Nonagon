"""Second independent check of the Part XIII identifiability claim, for grids too large for a dense row matrix.
Shares no code with tomo.c or record.py:
  * readings: dense FLINT inverse of K modulo a prime unused by record.py (default 998244353);
  * row system: rebuilt from the grid adjacency and the source list saved in record_CASE.json;
  * rank: sparse elimination in pure Python with a different pivot rule from tomo.c (rows taken shortest first,
    pivot on the LARGEST column index present, rows kept as dicts), so the elimination path is unrelated.
Usage: python3 verify13s.py CASE [prime]"""
import sys, json, time, gc
sys.path[:0] = ['../p11']
from fractions import Fraction as Fr
import flint, gridlib
case = sys.argv[1]; p = int(sys.argv[2]) if len(sys.argv) > 2 else 998244353
t0 = time.time()
rec = json.load(open('record_%s.json' % case))
nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
_i = {b: k for k, b in enumerate(nodes)}; _adj = {k: [] for k in range(len(nodes))}
for a, c in edges: _adj[a].append(c); _adj[c].append(a)
_seen = {_i[rec['port_bus']]}; _st = list(_seen)
while _st:
    u = _st.pop()
    for o in _adj[u]:
        if o not in _seen: _seen.add(o); _st.append(o)
keep = sorted(_seen); _new = {o: k for k, o in enumerate(keep)}
nodes = [nodes[i] for i in keep]; edges = {(_new[a], _new[c]): w for (a, c), w in edges.items() if a in _new}; ground = [ground[i] for i in keep]
n = len(nodes); idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
extra = [idx[b] for b in rec.get('extra_source_buses', [])]
wid = {idx[int(b)]: r for b, r in rec.get('widened_sources', {}).items()}
def res(x): x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p
K = flint.nmod_mat(n, n, p)
for (a, c), w in edges.items():
    r = res(w); K[a, a] += r; K[c, c] += r; K[a, c] -= r; K[c, a] -= r
for i, g in enumerate(ground):
    if g: K[i, i] += res(g)
G = K.inv(); del K; gc.collect(); t_inv = time.time() - t0
E = list(edges); Gs = [i for i in range(n) if ground[i] != 0]
col = {e: k for k, e in enumerate(E)}; gcol = {i: len(E) + k for k, i in enumerate(Gs)}; U = len(E) + len(Gs)
inc = {i: [] for i in range(n)}
for e in E: inc[e[0]].append(e); inc[e[1]].append(e)
need = []                                   # (i, s) pairs of rows
need += [(i, port) for i in range(n)]
need += [(s, s) for s in range(n) if s != port]
for s in extra:
    dist = {s: 0}; front = [s]
    for d in range(1, wid.get(s, 1) + 1):
        nxt = []
        for u in front:
            for e in inc[u]:
                o = e[1] if e[0] == u else e[0]
                if o not in dist: dist[o] = d; nxt.append(o)
        front = nxt
    need += [(i, s) for i in dist if i != s]
rows = []
for i, s in need:
    gi = int(G[i, s]); r = {}
    for e in inc[i]:
        o = e[1] if e[0] == i else e[0]; v = (gi - int(G[o, s])) % p
        if v: r[col[e]] = v
    if i in gcol and gi: r[gcol[i]] = gi
    if r: rows.append(r)
del G; gc.collect(); t_rows = time.time() - t0
rows.sort(key=len)
piv = {}                                    # column -> normalized row (dict), pivot = max column of the row
for r in rows:
    r = dict(r)
    while r:
        c = max(r)
        if c not in piv:
            inv = pow(r[c], -1, p); piv[c] = {k: v * inv % p for k, v in r.items()}; break
        f = r[c]; pr = piv[c]
        for k, v in pr.items():
            nv = (r.get(k, 0) - f * v) % p
            if nv: r[k] = nv
            else: r.pop(k, None)
rank = len(piv)
out = {'case': case, 'prime': p, 'rows': len(rows), 'unknowns': U, 'rank': rank, 'identifiable': rank == U,
       'agrees_with_record': (rank == U) == rec['identifiable'], 'pivot_rule': 'shortest row first, largest column',
       'seconds_inverse': t_inv, 'seconds_rows': t_rows - t_inv, 'seconds': time.time() - t0}
print(json.dumps(out)); json.dump(out, open('verify2_%s.json' % case, 'w'), indent=1)
