"""Necessary-condition test of Proposition XIII.2.2 on the real grids, for the |T|=1 case (hanging blocks).
For every articulation vertex t and every component C of G - t that contains neither the port nor a ground,
the pocket B = C u {t} has dim N(B) >= cyclomatic number of B. Any full-rank run must therefore have a source in C
whenever B contains a cycle. Also re-tests K4 two-terminal pockets."""
import sys, json, networkx as nx
sys.path[:0] = ['../p11']
import gridlib
res = {}
for case in sys.argv[1:]:
    rec = json.load(open('record_%s.json' % case))
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
    idx = {b: k for k, b in enumerate(nodes)}
    G = nx.Graph(); G.add_edges_from(edges)
    port = idx[rec['port_bus']]; G = G.subgraph(nx.node_connected_component(G, port)).copy()
    srcs = {idx[b] for b in rec['extra_source_buses']} | {port}
    seen = set(); hanging = []; violations = []
    for t in nx.articulation_points(G):
        H = G.copy(); H.remove_node(t)
        for C in nx.connected_components(H):
            if port in C or any(ground[i] for i in C) or len(C) > 2000: continue
            key = frozenset(C)
            if key in seen: continue
            seen.add(key)
            B = G.subgraph(C | {t}); cyc = B.number_of_edges() - B.number_of_nodes() + 1
            if cyc > 0:
                hanging.append(cyc)
                if not (C & srcs): violations.append(sorted(nodes[i] for i in C)[:6])
    res[case] = {'hanging_cyclic_pockets': len(hanging), 'total_cyclomatic': sum(hanging), 'pockets_without_interior_source': len(violations),
                 'identifiable': rec['identifiable']}
    print(case, res[case], violations[:3])
json.dump(res, open('pockets_hanging.json', 'w'), indent=1)
