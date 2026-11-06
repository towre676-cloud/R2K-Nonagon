"""Part XIV.5: how many ports make every branch identifiable (exact, modulo 2^31-1).
Access: one intervention of strength 1 at every bus; each port reads its own driving-point compliance after every
intervention and at baseline (k(n+1) readings), and every pair of ports reads its transfer compliance at baseline.
Loop (PerfectJet ch.7 for ports): start from 5 farthest-point ports; while the Jacobian is rank-deficient, take the
FLINT kernel basis, and for each basis vector add as a port the endpoint of largest degree among the branches it
moves. Report the ports needed and the rank history."""
import sys, json, time
import numpy as np, flint
from ident14 import load, exact_G, rows, analyse, p
from crb14 import farthest_ports

def jac_dp(D, Gp, ports):
    n = D['n']; pq = [(a, b) for ia, a in enumerate(ports) for b in ports[ia:]]
    pa = np.array([a for a, b in pq]); pb = np.array([b for a, b in pq]); out = [rows(D, Gp[:, pa], Gp[:, pb])]
    P = np.array(ports)
    for j in range(n):
        f = pow(int(1 + Gp[j, j]) % p, -1, p); cj = (Gp[j, P] * f) % p
        X = (Gp[:, P] - (Gp[:, j][:, None] * cj[None, :]) % p) % p
        out.append(rows(D, X, X))
    return np.vstack(out)

case = sys.argv[1]; D = load(case); Gp = exact_G(case, D); deg = {i: len(D['adj'][i]) for i in range(D['n'])}
ports = farthest_ports(D, 5); hist = []; t0 = time.time()
for it in range(60):
    J = jac_dp(D, Gp, ports); U = J.shape[1]
    if J.shape[0] * U <= 2e7: M = flint.nmod_mat(J.tolist(), p)
    else:   # random compression R J (R uniform mod p, U+20 rows), chunked; rank(RJ) = U certifies rank(J) = U
        rng = np.random.default_rng(len(ports)); M = flint.nmod_mat(U + 20, U, p)
        for k0 in range(0, J.shape[0], 1500):
            Jc = flint.nmod_mat(J[k0:k0 + 1500].tolist(), p)
            Rc = flint.nmod_mat(rng.integers(0, p, size=(U + 20, Jc.nrows())).tolist(), p)
            M = M + Rc * Jc
    r = M.rank()
    hist.append({'ports': len(ports), 'readings': int(J.shape[0]), 'rank': r, 'unknowns': U, 'compressed': bool(J.shape[0] * U > 2e7)}); print(case, hist[-1], flush=True)
    if r == U: break
    N, nul = M.nullspace(); nb = len(D['el']); new = set()
    for k in range(nul):
        sup = [i for i in range(U) if int(N[i, k])]
        nodes = set()
        for i in sup:
            if i < nb: nodes |= set(D['el'][i])
            else: nodes.add(D['Gs'][i - nb])
        nodes -= set(ports)
        if nodes: new.add(max(nodes, key=lambda v: deg[v]))
    if not new: break
    ports = ports + sorted(new)
out = {'case': case, 'n': D['n'], 'ports_needed': len(ports) if hist[-1]['rank'] == hist[-1]['unknowns'] else None,
       'port_buses': [D['nodes'][q] for q in ports], 'history': hist, 'seconds': time.time() - t0}
json.dump(out, open('portloop14_%s.json' % case, 'w'), indent=1); print(json.dumps({k: out[k] for k in ('case', 'n', 'ports_needed', 'seconds')}))
