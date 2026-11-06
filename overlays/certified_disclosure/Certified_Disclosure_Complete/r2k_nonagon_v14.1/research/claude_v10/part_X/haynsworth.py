"""X.1 (Haynsworth): inertia by Schur-complement pivots versus leading-minor determinants, in ball arithmetic, on boxes with many zeros."""
import sys; sys.path.insert(0, ".")
import json, sys, time, certgen, r2kverify
from flint import arb, arb_mat, acb, ctx
def mus_from(cert):  # independent full replay moments (verifier's own enclosures)
    v = r2kverify.verify(cert, full=True); return v
def pivots(M, n):
    """LDL^T without pivoting: d_k = Schur complement pivots; Haynsworth: In(M) = sum of signs of pivots (Sylvester via quotients)"""
    A = [[M[i][j] for j in range(n)] for i in range(n)]; d = []
    for k in range(n):
        p = A[k][k]; d.append(p)
        if not (p > 0 or p < 0): return d, False
        for i in range(k+1, n):
            f = A[i][k]/p
            for j in range(k+1, n): A[i][j] = A[i][j] - f*A[k][j]
    return d, True
def minors(M, n):
    out = []
    for k in range(1, n+1): out.append(arb_mat(k, k, [M[i][j] for i in range(k) for j in range(k)]).det())
    return out
res = []
for (lo, hi, nk, prec) in ((10.0, 51.0, 22, 80), (10.0, 76.5, 38, 80), (10.0, 100.0, 60, 80), (10.0, 100.0, 60, 128)):
    t0 = time.time()
    c = certgen.generate('xi', lo, hi, 0.45, list(range(nk)), True, mode='series', prec=prec)
    v = r2kverify.verify(c, full=False)   # count from recorded (generator) balls; inertia attempt below uses the same balls
    g = [arb(x[0], x[1]) for x in c['generator_moments']]
    ctx.prec = prec
    n = int(round(float(g[0].mid())))
    M = [[(g[i+j] + g[abs(i-j)])/2 for j in range(n)] for i in range(n)]
    mn = minors(M, n); ok_m = all(x > 0 for x in mn)
    pv, okp = pivots(M, n); ok_p = okp and all(x > 0 for x in pv)
    rel = lambda x: float(x.rad()/abs(x.mid())) if x.mid() != 0 else float('inf')
    res.append({'box': [lo, hi], 'zeros': n, 'moments': nk, 'prec_bits': prec, 'segments': len(c['segments']), 'gen_time': time.time() - t0,
                'minors_decide': ok_m, 'pivots_decide': ok_p, 'worst_minor_rel_radius': max(rel(x) for x in mn), 'worst_pivot_rel_radius': max(rel(x) for x in pv),
                'last_minor': mn[-1].str(5, radius=True), 'smallest_pivot': min(pv, key=lambda x: float(x.mid())).str(5, radius=True), 'verifier_light': v.get('status'), 'verifier_reason': v.get('reason', '')})
    print(res[-1], flush=True)
json.dump(res, open('p10/haynsworth.json', 'w'), indent=1)
