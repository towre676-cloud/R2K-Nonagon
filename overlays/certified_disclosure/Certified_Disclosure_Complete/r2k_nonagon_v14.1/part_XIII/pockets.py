"""Census of harmonic pockets and a necessary-condition test of Lemma XIII.2.2.
Pocket: node set B = I u T, interior I nonempty, ungrounded, every neighbour of I inside B; T = N(I) minus I.
For each record run, report pockets by (|I|,|T|), and whether each K4 pocket (|I|=|T|=2, B complete) has a source inside I."""
import sys, json, itertools
sys.path[:0] = ['../p11']
import gridlib
from collections import Counter
out = {}
for case in sys.argv[1:]:
    rec = json.load(open('record_%s.json' % case))
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
    idx = {b: k for k, b in enumerate(nodes)}; n = len(nodes)
    adj = {i: set() for i in range(n)}
    for a, c in edges: adj[a].add(c); adj[c].add(a)
    srcs = {idx[b] for b in rec['extra_source_buses']} | {idx[rec['port_bus']]}
    # interior components: connected components of the subgraph induced by ungrounded nodes whose
    # closed neighbourhood is small; we take I = connected sets of ungrounded nodes, grown greedily,
    # restricted to |I| <= 3 and |T| <= 3 (small pockets)
    k4 = []; small = Counter()
    for i in range(n):
        if ground[i]: continue
        for j in adj[i]:
            if j <= i or ground[j]: continue
            I = {i, j}; T = (adj[i] | adj[j]) - I
            if len(T) == 2:
                B = I | T
                complete = all(b in adj[a] for a, b in itertools.combinations(B, 2))
                small[('|I|=2,|T|=2', 'K4' if complete else 'other')] += 1
                if complete: k4.append((sorted(nodes[x] for x in I), sorted(nodes[x] for x in T), bool(I & srcs)))
    out[case] = {'K4_pockets': len(k4), 'K4_pockets_with_interior_source': sum(p[2] for p in k4),
                 'K4_pockets_without_interior_source': [p for p in k4 if not p[2]], 'identifiable': rec['identifiable'],
                 'other_two_interior_two_terminal_pockets': small[('|I|=2,|T|=2', 'other')]}
    print(case, {k: (v if not isinstance(v, list) else len(v)) for k, v in out[case].items()})
json.dump(out, open('pockets_census.json', 'w'), indent=1)
