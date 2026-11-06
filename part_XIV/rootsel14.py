"""Part XIV.4: root selection by two pair strengths (binary64).
For a pair (i,j) and strength theta, the pair reading delta gives Q(x) = delta x^2 - 2 G0i G0j x + (G0i^2 b + G0j^2 a - delta a b),
a = 1/theta + G_ii, b = 1/theta + G_jj, whose roots are x_true and s = 2 G0i G0j / delta - x_true. Two strengths give
two spurious roots that coincide only when delta_1 = delta_2, so the common root is x_true. This script reads every
lean pair of a grid at strengths theta and 2 theta (theta = c / median G_jj for c in a list), and counts
  minus_fail : pairs where the Part VIII minus root is not the true G_ij,
  two_fail   : pairs where the two-strength rule (the root of Q_1 closest to a root of Q_2) is not the true G_ij,
  min_gap    : the smallest |delta_1 - delta_2| / |delta_1| over the pairs (the rule's margin).
All quantities are computed from the model's G in binary64; 'not the true G_ij' means a relative miss above 1e-6."""
import sys, json
import numpy as np
from crb14 import load, lean_pairs

case = sys.argv[1]; D = load(case); G = D['G']; port = D['port']; pairs = lean_pairs(D)
P = np.array(pairs); i, j = P[:, 0], P[:, 1]
g0i, g0j, gii, gjj, x = G[port, i], G[port, j], G[i, i], G[j, j], G[i, j]
med = np.median(np.diag(G)); out = {'case': case, 'pairs': len(pairs), 'negative_branches': int(sum(1 for e in D['el'] if D['truth'][D['el'].index(e)] < 0)) if len(D['el']) < 20000 else None, 'rows': []}
def roots(th):
    a = 1 / th + gii; b = 1 / th + gjj
    dl = (g0i ** 2 * b + g0j ** 2 * a - 2 * g0i * g0j * x) / (a * b - x ** 2)
    disc = np.maximum((dl * a - g0i ** 2) * (dl * b - g0j ** 2), 0)
    return dl, (g0i * g0j - np.sqrt(disc)) / dl, (g0i * g0j + np.sqrt(disc)) / dl
for c in (0.001, 0.01, 0.1, 1.0, 10.0, 100.0):
    th = c / med
    d1, m1, p1 = roots(th); d2, m2, p2 = roots(2 * th)
    tol = 1e-6 * np.maximum(np.abs(x), 1e-300)
    minus_fail = int((np.abs(m1 - x) > tol).sum())
    c1 = np.stack([m1, p1]); c2 = np.stack([m2, p2])
    dist = np.abs(c1[:, None, :] - c2[None, :, :]).min(1)          # each root of Q_1 vs nearest root of Q_2
    pick = np.where(dist[0] <= dist[1], m1, p1)
    two_fail = int((np.abs(pick - x) > tol).sum())
    gap = np.abs(d1 - d2) / np.abs(d1)
    r = {'theta_times_median_Gjj': c, 'minus_fail': minus_fail, 'two_fail': two_fail, 'min_gap': float(gap.min()),
         'theta_G00': float(th * G[port, port])}
    out['rows'].append(r); print(case, r, flush=True)
json.dump(out, open('rootsel14_%s.json' % case, 'w'), indent=1)
