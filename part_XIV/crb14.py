"""Part XIV.1: Cramer-Rao floor for known-support tomography (binary64).
Reading model: y = e_p^T (K + Theta)^{-1} e_q + noise, noise iid N(0, (sigma * G00)^2).
d y / d w_e = -(b_e^T x_p)(b_e^T x_q),  d y / d g_i = -x_p[i] x_q[i],  x_p = (K + Theta)^{-1} e_p.
Schemes:
  lean  : Part XIII one-port reading set (baseline, two single-node strengths per node, the pair list of record.py)
  full  : Part VIII n(n+1)/2 protocol (baseline, two strengths per node, every pair)
  ports k: k ports chosen by farthest-point sampling; baseline transfer compliances among ports, plus one single-node
           intervention at every node, each read at all port pairs.
Output per scheme: number of readings, and the relative CRB per unit sigma, k_e = G00 sqrt((F^T F)^{-1})_ee / |w_e|,
summarized by median and 90th percentile, and sigma_1pct = 0.01 / median k_e (the noise level at which an efficient
estimator reaches a 1% median error with one pass of the readings)."""
import sys, json, time
sys.path[:0] = ['../p11']
import numpy as np, gridlib
from scipy.sparse.csgraph import connected_components, shortest_path
import scipy.sparse as sp

def load(case):
    rec = json.load(open('../p13/record_%s.json' % case))
    nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case)
    n = len(nodes); idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
    r = [a for a, c in edges]; c_ = [c for a, c in edges]
    A = sp.coo_matrix((np.ones(len(r)), (r, c_)), shape=(n, n))
    _, lab = connected_components(A, directed=False)
    keep = [i for i in range(n) if lab[i] == lab[port]]; new = {o: k for k, o in enumerate(keep)}
    E = {(new[a], new[c]): w for (a, c), w in edges.items() if a in new}; g = [ground[i] for i in keep]
    nodes = [nodes[i] for i in keep]; n = len(nodes); idx = {b: k for k, b in enumerate(nodes)}; port = idx[rec['port_bus']]
    el = sorted(E); Gs = [i for i in range(n) if g[i] != 0]
    K = gridlib.K_float(n, E, g); G = np.linalg.inv(K)
    truth = np.array([float(E[e]) for e in el] + [float(g[i]) for i in Gs])
    adj = {i: set() for i in range(n)}
    for a, c in el: adj[a].add(c); adj[c].add(a)
    extra = [idx[b] for b in rec['extra_source_buses']]; wid = {idx[int(b)]: v for b, v in rec.get('widened_sources', {}).items()}
    bus, _ = gridlib.parse('../p11/data/%s.m' % case)
    hw = sorted(idx[int(float(b[0]))] for b in bus if int(float(b[0])) in idx and (float(b[5]) != 0 or int(float(b[1])) == 2))
    return dict(n=n, port=port, el=el, Gs=Gs, K=K, G=G, truth=truth, adj=adj, extra=extra, wid=wid, rec=rec, nodes=nodes, hw=hw)

def lean_pairs(D):
    n, port, adj = D['n'], D['port'], D['adj']; pairs = set()
    for s in range(n):
        if s == port: continue
        for j in adj[s]:
            if j != port: pairs.add((min(s, j), max(s, j)))
    for s in D['extra']:
        seen = {s}; front = [s]
        for _ in range(D['wid'].get(s, 1) + 1):
            front = [o for u in front for o in adj[u] if o not in seen and not seen.add(o)]
        for j in seen:
            if j not in (port, s): pairs.add((min(s, j), max(s, j)))
    return sorted(pairs)

def rows_for(D, x_p, x_q):
    """Jacobian row(s) for readings x_p^T (.) x_q: x_p, x_q are (n, r) arrays (r readings at once)."""
    el = D['el']; a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el])
    dp = x_p[a] - x_p[c]; dq = x_q[a] - x_q[c]
    R = np.vstack([-(dp * dq), -(x_p[D['Gs']] * x_q[D['Gs']])])   # (U, r)
    return (R * D['truth'][:, None]).T                            # columns scaled: unknowns in relative units

def fisher_lean_or_full(D, pairs, t):
    G, port, n = D['G'], D['port'], D['n']; U = len(D['el']) + len(D['Gs'])
    J = np.zeros((U, U)); m = 0
    g0 = G[:, port]
    R = rows_for(D, g0[:, None], g0[:, None]); J += R.T @ R; m += 1          # baseline
    for th in (t, 2 * t):                                                      # single-node, two strengths
        js = np.array([j for j in range(n) if j != port])
        X = g0[:, None] - G[:, js] * (th * G[port, js] / (1 + th * G[js, js]))[None, :]
        R = rows_for(D, X, X); J += R.T @ R; m += len(js)
    B = 4096
    for k0 in range(0, len(pairs), B):                                          # pair readings at strength t
        P = pairs[k0:k0 + B]; Xs = np.empty((n, len(P)))
        for k, (i, j) in enumerate(P):
            S = [i, j]; M = np.eye(2) / t + G[np.ix_(S, S)]
            Xs[:, k] = g0 - G[:, S] @ np.linalg.solve(M, G[S, port])
        R = rows_for(D, Xs, Xs); J += R.T @ R; m += len(P)
    return J, m

def fisher_ports(D, ports, t, where=None):
    G, n = D['G'], D['n']; U = len(D['el']) + len(D['Gs'])
    J = np.zeros((U, U)); m = 0
    pq = [(a, b) for ia, a in enumerate(ports) for b in ports[ia:]]
    pa = np.array([a for a, b in pq]); pb = np.array([b for a, b in pq])
    R = rows_for(D, G[:, pa], G[:, pb]); J += R.T @ R; m += len(pq)              # baseline transfer compliances
    for j in (range(n) if where is None else where):
        f = t / (1 + t * G[j, j])
        Xa = G[:, pa] - np.outer(G[:, j], f * G[j, pa]); Xb = G[:, pb] - np.outer(G[:, j], f * G[j, pb])
        R = rows_for(D, Xa, Xb); J += R.T @ R; m += len(pq)
    return J, m

def summarize(D, J, m):
    w, V = np.linalg.eigh(J); G00 = D['G'][D['port'], D['port']]
    lam_min = w[0]; tiny = w < w[-1] * 1e-15                       # directions beyond binary64 resolution of J
    wi = np.where(tiny, np.inf, w)
    var = (V ** 2 / wi[None, :]).sum(1); var[(V[:, tiny] ** 2).sum(1) > 1e-6] = np.inf
    ke = G00 * np.sqrt(var); nb = len(D['el']); kb = ke[:nb]
    fin = np.isfinite(kb)
    return {'readings': m, 'unknowns': len(D['truth']), 'median_k': float(np.median(kb)), 'p90_k': float(np.percentile(kb, 90)),
            'unresolved_branches': int((~fin).sum()), 'sigma_for_1pct_median': float(0.01 / np.median(kb)),
            'frac_1pct_at_1e-3': float(np.mean(kb * 1e-3 < 0.01)), 'frac_1pct_at_1e-4': float(np.mean(kb * 1e-4 < 0.01)),
            'frac_1pct_at_1e-6': float(np.mean(kb * 1e-6 < 0.01)), 'eig_ratio': float(w[0] / w[-1])}

def farthest_ports(D, k):
    n = D['n']; el = D['el']
    A = sp.coo_matrix((np.ones(len(el)), ([a for a, c in el], [c for a, c in el])), shape=(n, n))
    P = [D['port']]; dmin = shortest_path(A, directed=False, unweighted=True, indices=[D['port']])[0]
    while len(P) < k:
        nxt = int(np.argmax(dmin)); P.append(nxt)
        dmin = np.minimum(dmin, shortest_path(A, directed=False, unweighted=True, indices=[nxt])[0])
    return P

if __name__ == '__main__':
    case = sys.argv[1]; schemes = sys.argv[2].split(',') if len(sys.argv) > 2 else ['lean']
    D = load(case); t = 1 / np.median(np.diag(D['G'])); out = {'case': case, 'theta': float(t), 'hw_buses': len(D['hw']), 'schemes': {}}
    for sc in schemes:
        t0 = time.time()
        if sc == 'lean': J, m = fisher_lean_or_full(D, lean_pairs(D), t)
        elif sc == 'full':
            prs = [(i, j) for i in range(D['n']) for j in range(i + 1, D['n']) if D['port'] not in (i, j)]
            J, m = fisher_lean_or_full(D, prs, t)
        elif sc.startswith('hw'):            # hwK: K ports, interventions only at buses with a shunt or a generator
            k = int(sc[2:]); J, m = fisher_ports(D, farthest_ports(D, k), t, where=D['hw'])
        else:
            k = int(sc[5:]); J, m = fisher_ports(D, farthest_ports(D, k), t)
        s = summarize(D, J, m); s['seconds'] = time.time() - t0; out['schemes'][sc] = s
        print(case, sc, s, flush=True)
    json.dump(out, open('crb14_%s.json' % case, 'w'), indent=1)
