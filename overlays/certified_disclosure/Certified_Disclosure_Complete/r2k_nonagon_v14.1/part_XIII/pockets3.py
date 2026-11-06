"""Hanging-block census via the block-cut tree (linear time); same test as pockets2.py.
For each articulation vertex t and each child block of t (tree rooted at the port), the component C of G - t
below it is a pocket with T = {t} when C has no ground and does not contain the port. Its cycle rank is
|E(C u t)| - |C|. Records whether C contains a source."""
import sys, json, networkx as nx
sys.path[:0] = ['../p11']
import gridlib
sys.setrecursionlimit(1000000)
res = {}
for case in sys.argv[1:]:
    rec = json.load(open('record_%s.json' % case))
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
    idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
    G = nx.Graph(); G.add_edges_from(edges); G = G.subgraph(nx.node_connected_component(G, port)).copy()
    srcs = {idx[b] for b in rec['extra_source_buses']} | {port}
    blocks = [set(b) for b in nx.biconnected_components(G)]
    arts = set(nx.articulation_points(G))
    T = nx.Graph()
    for k, b in enumerate(blocks):
        for v in b & arts: T.add_edge(('b', k), ('a', v))
        T.add_node(('b', k))
    root = ('b', next(k for k, b in enumerate(blocks) if port in b)) if port not in arts else ('a', port)
    # aggregate below each tree node: set of non-articulation-owned nodes is handled by counting per block
    order = list(nx.dfs_preorder_nodes(T, root)); parent = {root: None}
    for u, v in nx.dfs_edges(T, root): parent[v] = u
    agg = {}   # tree node -> (nodes below excluding its parent articulation, edges below, grounded, sourced)
    for x in reversed(order):
        if x[0] == 'b':
            b = blocks[x[1]]; pa = parent[x][1] if parent[x] else None
            own = b - ({pa} if pa is not None else set())
            nn = len(own); ee = G.subgraph(b).number_of_edges(); gr = any(ground[i] for i in own); so = bool(own & srcs)
            for ch in T.neighbors(x):
                if ch != parent[x]:
                    a = agg[ch]; nn += a[0] - 1; ee += a[1]; gr = gr or a[2]; so = so or a[3]   # child art. already in own
            agg[x] = (nn, ee, gr, so)
        else:
            nn = 1; ee = 0; v = x[1]; gr = bool(ground[v]); so = v in srcs
            for ch in T.neighbors(x):
                if ch != parent[x]:
                    a = agg[ch]; nn += a[0]; ee += a[1]; gr = gr or a[2]; so = so or a[3]
            agg[x] = (nn, ee, gr, so)
    hanging = []; viol = 0
    for x in order:
        if x[0] == 'b' and parent[x] is not None:
            nn, ee, gr, so = agg[x]          # nn = |C| (excludes the parent articulation t), ee = |E(C u t)|
            if gr: continue
            cyc = ee - nn
            if cyc > 0:
                hanging.append(cyc); viol += int(not so)
    res[case] = {'hanging_cyclic_pockets': len(hanging), 'total_cyclomatic': sum(hanging), 'pockets_without_interior_source': viol}
    print(case, res[case], flush=True)
json.dump(res, open('pockets_hanging.json', 'w'), indent=1)
