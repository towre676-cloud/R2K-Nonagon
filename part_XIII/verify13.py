"""Independent re-derivation of the Part XIII identifiability claim, sharing no code with tomo.c.
Dense FLINT nmod_mat inverse of K modulo a prime NOT used by record.py, the same reading set
(port column; diagonal (s,s) for every source; neighbour rows for the R=1 extra sources), rank by FLINT.
Usage: python3 verify13.py CASE [prime]"""
import sys, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'part_XI'))
from fractions import Fraction as Fr
import flint, gridlib
case = sys.argv[1]; p = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000007
t0 = time.time()
rec = json.load(open(HERE / ('record_%s.json' % case)))
nodes, edges, ground, meta = gridlib.network(ROOT / 'part_XIII/data' / (case + '.m'))
# same visibility restriction as record.py, recomputed here independently by breadth-first search from the port
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
def res(x): x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p
K = flint.nmod_mat(n, n, p)
for (a, c), w in edges.items():
    r = res(w); K[a, a] += r; K[c, c] += r; K[a, c] -= r; K[c, a] -= r
for i, g in enumerate(ground):
    if g: K[i, i] += res(g)
G = K.inv()
E = list(edges); Gs = [i for i in range(n) if ground[i] != 0]
col = {e: k for k, e in enumerate(E)}; gcol = {i: len(E) + k for k, i in enumerate(Gs)}; U = len(E) + len(Gs)
inc = {i: [] for i in range(n)}
for e in E: inc[e[0]].append(e); inc[e[1]].append(e)
rows = []
def row(i, s):
    r = {}
    gi = int(G[i, s])
    for e in inc[i]:
        o = e[1] if e[0] == i else e[0]
        r[col[e]] = (gi - int(G[o, s])) % p
    if i in gcol: r[gcol[i]] = gi
    rows.append(r)
for i in range(n): row(i, port)                     # port column, from single-node readings
for s in range(n):
    if s != port: row(s, s)                         # diagonal rows (equation i = s of column s)
wid = {idx[int(b)]: r for b, r in rec.get('widened_sources', {}).items()}
for s in extra:
    dist = {s: 0}; front = [s]
    for d in range(1, wid.get(s, 1) + 1):
        nxt = []
        for u in front:
            for e in inc[u]:
                o = e[1] if e[0] == u else e[0]
                if o not in dist: dist[o] = d; nxt.append(o)
        front = nxt
    for i in dist:
        if i != s: row(i, s)
M = flint.nmod_mat(len(rows), U, p)
for k, r in enumerate(rows):
    for c, v in r.items(): M[k, c] = v
rk = M.rank()
out = {'case': case, 'prime': p, 'rows': len(rows), 'unknowns': U, 'rank': rk, 'identifiable': rk == U,
       'agrees_with_record': (rk == U) == rec['identifiable'], 'seconds': time.time() - t0}
print(json.dumps(out)); json.dump(out, open(HERE / ('verify_%s.json' % case), 'w'), indent=1)
