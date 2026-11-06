"""XIII: exact one-port known-support tomography at grid scale. C kernel (tomo.c) modulo three word primes, CRT, and
PerfectJet's rational reconstruction (Wang) -> exact rational susceptances, compared with the MATPOWER source exactly."""
import sys, json, time, subprocess, os
sys.path[:0] = ['../p11', '/tmp/claude-0/repos/PerfectJet/src']
from fractions import Fraction as Fr
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import gridlib
from perfectjet.modular import rational_reconstruction
import sympy
ALLP = []; _c = int(os.environ.get('PRIME_TOP', '2147483647'))
while len(ALLP) < 24:
    if sympy.isprime(_c): ALLP.append(_c)
    _c -= 2
PRIMES = ALLP[:1]
case = sys.argv[1]; t0 = time.time()
nodes, edges, ground, meta = gridlib.network('../p11/data/%s.m' % case); n = len(nodes)
deg = np.zeros(n, int)
for a, c in edges: deg[a] += 1; deg[c] += 1
port = int(np.argmax(deg))
# Islands joined to the rest only through the reference bus have G_{0i} = 0: no port reading sees them. Restrict to the port's component.
from scipy.sparse.csgraph import connected_components
_r = [a for a, c in edges]; _c = [c for a, c in edges]
ncomp, lab = connected_components(sp.coo_matrix((np.ones(len(_r)), (_r, _c)), shape=(n, n)), directed=False)
islands = {'components': int(ncomp), 'invisible_nodes': int((lab != lab[port]).sum()),
           'invisible_unknowns': sum(1 for a, c in edges if lab[a] != lab[port]) + sum(1 for i in range(n) if ground[i] and lab[i] != lab[port])}
if ncomp > 1:
    keep = [i for i in range(n) if lab[i] == lab[port]]; new = {o: k for k, o in enumerate(keep)}
    nodes = [nodes[i] for i in keep]; edges = {(new[a], new[c]): w for (a, c), w in edges.items() if a in new}
    ground = [ground[i] for i in keep]; port = new[port]; n = len(nodes)
# (dense K removed: the ordering below uses the sparse matrix only)
rows, cols, vals = [], [], []
for (a, c), w in edges.items():
    w = float(w); rows += [a, c, a, c]; cols += [a, c, c, a]; vals += [w, w, -w, -w]
for i, g in enumerate(ground):
    if g: rows.append(i); cols.append(i); vals.append(float(g))
A = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
ORDER = os.environ.get('ORDER', 'MMD_AT_PLUS_A'); perm = spl.splu(A, permc_spec=ORDER).perm_c.tolist(); iperm = {v: k for k, v in enumerate(perm)}
# physical-model checks in binary64 (reported, not used by the exact pipeline): definiteness of K, and the port column G_{0i}
_lu = spl.splu(A, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
_d = _lu.U.diagonal(); _e = np.zeros(n); _e[port] = 1.0; _x = _lu.solve(_e)
physical = {'K_positive_definite_binary64': bool((_d > 0).all()), 'negative_pivots': int((_d < 0).sum()),
            'port_column_min_over_max': float(np.abs(_x).min() / np.abs(_x).max()), 'port_column_all_positive': bool((_x > 0).all())}
E = sorted(edges, key=lambda e: min(iperm[e[0]], iperm[e[1]])); Gs = [i for i in range(n) if ground[i] != 0]
U = len(E) + len(Gs); t_prep = time.time() - t0
def res(x, p): x = Fr(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p
def run(p, extra):
    lines = ['%d %d %d %d %d' % (p, n, port, len(E), len(Gs)), ' '.join(map(str, perm))]
    lines += ['%d %d %d' % (a, c, res(edges[(a, c)], p)) for a, c in E] + ['%d %d' % (i, res(ground[i], p)) for i in Gs]
    lines += [str(len(extra)), ' '.join(map(str, extra)), str(len(extra)), ' '.join(str(radius.get(x, 1)) for x in extra)]
    r = subprocess.run(['./tomo'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-500:])
    out = r.stdout.split('\n'); rank, UU, nrows = map(int, out[0].split())
    if rank == UU: return rank, [int(x) for x in out[1:1 + UU]], [], r.stderr
    # free columns plus the supports of the kernel basis vectors (PerfectJet ch.7: ker O names the confounded unknowns)
    fr = [int(x.split()[1]) for x in out[1:] if x.startswith('free')]
    for x in out[1:]:
        if x.startswith('kern'): fr += [int(t) for t in x.split('|')[0].split()[2:]]
    return rank, None, sorted(set(fr)), r.stderr
adj = {i: set() for i in range(n)}
for a, c in E: adj[a].add(c); adj[c].add(a)
extra = []; log = []; radius = {}; last_rank = -1
for it in range(12):
    t1 = time.time(); rank, sol, free, err = run(PRIMES[0], extra); log.append({'iteration': it, 'extra_sources': len(extra), 'rank': rank, 'seconds': time.time() - t1, 'kernel_log': err.strip().split('\n')})
    print(case, log[-1]['iteration'], rank, U, len(extra), '%.1fs' % log[-1]['seconds'], flush=True)
    if rank == U: break
    new = set()
    for f in free:
        if f < len(E): new |= set(E[f])
        else: new.add(Gs[f - len(E)])
    new -= {port}
    if rank == last_rank:      # the R=1 rows have stalled: widen the row ball around the stalled sources (PerfectJet ch.7 loop)
        for x in new: radius[x] = radius.get(x, 1) + 1
    last_rank = rank; extra = sorted(set(extra) | new)
exact = None; mismatches = None; recon_log = []; sol_hash = None
def recon(sols, ps):
    M = 1
    for p in ps: M *= p
    out = []
    for k in range(U):
        x = 0
        for p, s in zip(ps, sols):
            Mi = M//p; x = (x + s[k]*Mi*pow(Mi, -1, p)) % M
        out.append(rational_reconstruction(x, M))
    return out
if rank == U:
    # early-termination CRT: add word primes until the reconstruction is defined everywhere and unchanged by one more prime.
    # The truth is consulted only after termination.
    sols = [sol]; prev = None
    from concurrent.futures import ThreadPoolExecutor
    cache = {}; rest = ALLP[1:]
    for k, p in enumerate(rest):
        if p not in cache:   # two primes at a time on two cores; the stopping rule below is unchanged
            batch = rest[k:k + 2]
            with ThreadPoolExecutor(2) as ex:
                for q, rr in zip(batch, ex.map(lambda q: run(q, extra), batch)): cache[q] = rr
        r = cache[p]; assert r[0] == U, 'rank drop mod %d' % p
        sols.append(r[1]); PRIMES.append(p)
        cur = recon(sols, PRIMES); ok = all(q is not None for q in cur)
        recon_log.append({'primes': len(PRIMES), 'all_defined': ok, 'stable': ok and cur == prev})
        if ok and cur == prev: break
        prev = cur if ok else None
    truth = [Fr(edges[e]) for e in E] + [Fr(ground[i]) for i in Gs]
    import hashlib
    canon = sorted(['e %d %d %s' % (min(nodes[e[0]], nodes[e[1]]), max(nodes[e[0]], nodes[e[1]]), q) for e, q in zip(E, cur)] + ['g %d %s' % (nodes[i], q) for i, q in zip(Gs, cur[len(E):])])
    sol_hash = hashlib.sha256('\n'.join(canon).encode()).hexdigest()
    mismatches = sum(int(q != t) for q, t in zip(cur, truth)); exact = mismatches == 0
    hbits = max(max(abs(t.numerator), t.denominator).bit_length() for t in truth)
pairs = set()
for s in range(n):
    if s == port: continue
    for j in adj[s]:
        if j != port: pairs.add((min(s, j), max(s, j)))
def ball(s, r):
    seen = {s}; front = [s]
    for _ in range(r):
        front = [o for u in front for o in adj[u] if o not in seen and not seen.add(o)]
    return seen
for s in extra:
    for j in ball(s, radius.get(s, 1) + 1):
        if j not in (port, s): pairs.add((min(s, j), max(s, j)))
readings = 1 + 2*(n - 1) + len(pairs)
out = {'case': case, 'buses': meta['buses'], 'hidden_nodes': n, 'nodes_total': meta['buses'] - 1, 'branches': len(E), 'grounds': len(Gs), 'unknowns': U, 'port_bus': nodes[port],
       'negative_reactance_branches': meta['negative_reactance_branches'], 'readings': readings, 'full_protocol_readings': n*(n + 1)//2,
       'readings_per_unknown': readings/U, 'saving_factor': n*(n + 1)/2/readings, 'rank_mod_p': rank, 'identifiable': rank == U,
       'primes': PRIMES, 'exact_rational_recovery': exact, 'rational_mismatches': mismatches, 'crt_log': recon_log, 'log2_modulus': sum(p.bit_length() for p in PRIMES), 'max_height_bits_of_truth': hbits if rank == U else None, 'iterations': log, 'islands': islands, 'physical': physical, 'extra_source_buses': [nodes[i] for i in extra], 'widened_sources': {str(nodes[i]): r for i, r in radius.items()}, 'port_index': port, 'seconds_total': time.time() - t0}
print(json.dumps({k: v for k, v in out.items() if k != 'iterations'}), flush=True)
out['ordering'] = ORDER; out['solution_sha256'] = sol_hash
json.dump(out, open('record_%s%s.json' % (case, os.environ.get('TAG', '')), 'w'), indent=1)
