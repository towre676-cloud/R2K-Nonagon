"""MATPOWER case -> grounded DC susceptance network (exact rationals). Reference bus = ground (Dirichlet)."""
import re
from fractions import Fraction as Fr
def parse(path):
    txt = open(path).read()
    def block(name):
        m = re.search(r'mpc\.%s\s*=\s*\[(.*?)\];' % name, txt, re.S)
        rows = []
        for line in m.group(1).split('\n'):
            line = line.split('%')[0].strip().rstrip(';').strip()
            if line: rows.append(line.split())
        return rows
    return block('bus'), block('branch')
def network(path):
    bus, br = parse(path)
    ids = [int(float(b[0])) for b in bus]; typ = {int(float(b[0])): int(float(b[1])) for b in bus}
    ref = [i for i in ids if typ[i] == 3][0]
    W = {}; neg = 0; zero = 0; off = 0
    for b in br:
        f, t, x, st = int(float(b[0])), int(float(b[1])), Fr(b[3]), int(float(b[10]))
        if st == 0: off += 1; continue
        if x == 0: zero += 1; continue
        if x < 0: neg += 1
        key = (min(f, t), max(f, t))
        if key[0] == key[1]: continue
        W[key] = W.get(key, Fr(0)) + 1/x
    nodes = [i for i in ids if i != ref]; idx = {b: k for k, b in enumerate(nodes)}
    n = len(nodes); edges = {}; ground = [Fr(0)]*n
    for (a, b), w in W.items():
        if a == ref: ground[idx[b]] += w
        elif b == ref: ground[idx[a]] += w
        else: edges[(idx[a], idx[b])] = w
    meta = {'buses': len(ids), 'reference_bus': ref, 'branches_in_service': len(br) - off, 'parallel_merged_edges': len(W),
            'negative_reactance_branches': neg, 'zero_reactance_branches': zero, 'out_of_service': off, 'types': typ}
    return nodes, edges, ground, meta
def K_float(n, edges, ground):
    import numpy as np
    K = np.zeros((n, n))
    for (i, j), w in edges.items(): w = float(w); K[i, i] += w; K[j, j] += w; K[i, j] -= w; K[j, i] -= w
    for i, g in enumerate(ground): K[i, i] += float(g)
    return K
