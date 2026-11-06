"""X.1b (Desnanot-Jacobi / Dodgson): leading minors by condensation versus Arb determinants on the 10-zero box; and the exact
Chebyshev-Vandermonde form M0 = sum_k v(u_k) v(u_k)^T, whose smallest eigenvalue sets the locality of the certificate."""
import sys; sys.path.insert(0, '.')
import json, certgen, numpy as np, mpmath as mp
from flint import arb, arb_mat, ctx, acb
c = certgen.generate('xi', 10.0, 51.0, 0.45, list(range(22)), True, mode='series', prec=80)
g = [arb(x[0], x[1]) for x in c['generator_moments']]; ctx.prec = 80; n = 10
M = [[(g[i+j] + g[abs(i-j)])/2 for j in range(n)] for i in range(n)]
# Dodgson: D^{(k+1)}_{ij} = (D^k_{ij} D^k_{i+1,j+1} - D^k_{i,j+1} D^k_{i+1,j}) / D^{(k-1)}_{i+1,j+1}; leading minor of order k+1 = D^{(k+1)}_{00}
prev = [[arb(1)]*(n+1) for _ in range(n+1)]; cur = [row[:] for row in M]; lead = [cur[0][0]]; ok = True
for k in range(1, n):
    m = n - k; nxt = [[None]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            den = prev[i+1][j+1]
            if not (den > 0 or den < 0): ok = False; break
            nxt[i][j] = (cur[i][j]*cur[i+1][j+1] - cur[i][j+1]*cur[i+1][j])/den
        if not ok: break
    if not ok: break
    prev, cur = cur, nxt; lead.append(cur[0][0])
direct = [arb_mat(k, k, [M[i][j] for i in range(k) for j in range(k)]).det() for k in range(1, n+1)]
rows = [{'k': k+1, 'dodgson': lead[k].str(6, radius=True) if k < len(lead) else None, 'direct': direct[k].str(6, radius=True),
         'rad_ratio_dodgson_over_direct': float(lead[k].rad()/direct[k].rad()) if k < len(lead) else None} for k in range(n)]
# locality law: lambda_min of the exact M0 from known zeros, boxes [10, T] with growing counts
mp.mp.dps = 40; Z = [float(acb.zeta_zero(k).imag.mid()) for k in range(1, 40)]
law = []
for T in (25, 35, 45, 51, 60, 70, 76.5, 85, 100):
    zs = [z for z in Z if 10 < z < T]; a = (T - 10)/2; cc = (T + 10)/2; u = [(z - cc)/a for z in zs]; m = len(zs)
    V = np.array([[np.cos(i*np.arccos(x)) for x in u] for i in range(m)])
    lam = float(np.linalg.eigvalsh(V@V.T)[0]); law.append({'box': [10, T], 'zeros': m, 'lambda_min_M0': lam})
out = {'dodgson_vs_direct': rows, 'dodgson_completed': ok, 'locality_law': law}
print(json.dumps(out, indent=1)); json.dump(out, open('p10/dodgson.json', 'w'), indent=1)
