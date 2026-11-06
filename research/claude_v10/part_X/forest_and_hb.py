"""X.4 Kirchhoff all-minors: G_ij as 2-forest ratio (exact enumeration). X.5 Hermite-Biehler/Brune PR certificate: Z = N/D (coprime)
is positive real iff N + D is strictly Hurwitz and Re Z(i w) >= 0; covers lossless (imaginary-axis) poles that Theorem VI.4.1 excludes."""
import sys, json, random, itertools, io, contextlib
sys.path.insert(0, '.')
from fractions import Fraction as Fr
with contextlib.redirect_stdout(io.StringIO()):
    from p3v6 import positive_root_count_v6, pr_certificate
from hydro_bench import routh_hurwitz_exact
OUT = {}
# ---------- Kirchhoff
def inv(A):
    n = len(A); M = [list(map(Fr, A[i])) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0); M[k], M[p] = M[p], M[k]; pv = M[k][k]; M[k] = [x/pv for x in M[k]]
        for i in range(n):
            if i != k and M[i][k] != 0: f = M[i][k]; M[i] = [x - f*y for x, y in zip(M[i], M[k])]
    return [r[n:] for r in M]
def comps(nodes, edges):
    par = {v: v for v in nodes}
    def fd(v):
        while par[v] != v: par[v] = par[par[v]]; v = par[v]
        return v
    for a, b in edges:
        ra, rb = fd(a), fd(b)
        if ra == rb: return None
        par[ra] = rb
    return fd
rng = random.Random(4); checks = 0; bad = 0
for t in range(25):
    n = rng.randint(2, 4); g = 'g'
    E = {}
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.7: E[(i, j)] = Fr(rng.randint(1, 5))
        if rng.random() < 0.6 or i == 0: E[(i, g)] = Fr(rng.randint(1, 4))
    nodes = list(range(n)) + [g]; el = list(E.items())
    K = [[Fr(0)]*n for _ in range(n)]
    for (a, b), w in el:
        if b == g: K[a][a] += w
        else: K[a][b] -= w; K[b][a] -= w; K[a][a] += w; K[b][b] += w
    try: G = inv(K)
    except StopIteration: continue
    trees = Fr(0); forests = {}
    for k in (n, n - 1):
        for sub in itertools.combinations(el, k):
            fd = comps(nodes, [e for e, _ in sub])
            if fd is None: continue
            w = Fr(1)
            for _, x in sub: w *= x
            if k == n: trees += w
            else:
                for i in range(n):
                    for j in range(n):
                        if fd(i) == fd(j) and fd(i) != fd(g): forests[(i, j)] = forests.get((i, j), Fr(0)) + w
    for i in range(n):
        for j in range(n):
            checks += 1; bad += (G[i][j] != forests.get((i, j), Fr(0))/trees)
OUT['kirchhoff_all_minors'] = {'entries_checked': checks, 'mismatches': bad, 'statement': 'G_ij = (weight of spanning 2-forests of the grounded graph with i, j in the tree not containing ground) / (weight of spanning trees)'}
# ---------- Hermite-Biehler / Brune
def hb_certificate(N, D):
    S = [Fr(0)]*max(len(N), len(D))
    for i, x in enumerate(N): S[i] += Fr(x)
    for i, x in enumerate(D): S[i] += Fr(x)
    while S and S[-1] == 0: S.pop()
    hur = routh_hurwitz_exact([float(x) for x in S[::-1]]) if len(S) > 1 else S[0] > 0
    c = pr_certificate(N, D)
    return {'N_plus_D_strictly_hurwitz': bool(hur), 'axis_nonneg': c['real_part_nonneg_on_axis'], 'PR_HB': bool(hur and c['real_part_nonneg_on_axis']), 'PR_VI41': c['positive_real']}
cases = {'1/(1-s)': ([1], [1, -1]), 's/(s^2+1) (lossless Foster)': ([0, 1], [1, 0, 1]), '1/s (capacitor)': ([1], [0, 1]), 's (inductor)': ([0, 1], [1]),
         '(s^2+1)/(s(s^2+4)) (Foster)': ([1, 0, 1], [0, 4, 0, 1]), '(s+2)/(s+1)': ([2, 1], [1, 1]), '(s-1)/(s+1) (non-PR)': ([-1, 1], [1, 1]),
         's/(s^2+0.1s+1)': ([0, 1], [1, Fr(1, 10), 1]), '(s^2+1)/(s^2+2) (not PR: Re<0 band)': ([1, 0, 1], [2, 0, 1])}
OUT['hb_cases'] = {k: hb_certificate(*v) for k, v in cases.items()}
for k, v in OUT['hb_cases'].items(): print(k, v)
print(OUT['kirchhoff_all_minors'])
json.dump(OUT, open('p10/forest_and_hb.json', 'w'), indent=1, default=str)
