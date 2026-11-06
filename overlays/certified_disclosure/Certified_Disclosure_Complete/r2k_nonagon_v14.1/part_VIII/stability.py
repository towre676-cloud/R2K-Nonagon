from pathlib import Path
HERE = Path(__file__).resolve().parent
"""VIII.2: noise amplification of the intervention protocol, achievable design, and two-point lower-bound scales on hard families."""
import numpy as np, json, itertools
from scipy.optimize import minimize_scalar
rng = np.random.default_rng(2027)
def h(K, theta):
    A = K.copy()
    for j, t in theta.items(): A[j, j] += t
    return np.linalg.solve(A, np.eye(len(K))[0])[0]
def protocol(n, th):  # th[j] = (theta1, theta2) per node; pair strength = geometric mean of node first strengths
    E = [{}]
    for j in range(1, n): E += [{j: th[j][0]}, {j: th[j][1]}]
    for i, j in itertools.combinations(range(1, n), 2): E.append({i: th[i][0], j: th[j][0]})
    return E
def reconstruct(n, data, th):
    it = iter(data); h0 = next(it); G = np.zeros((n, n)); G[0, 0] = h0; g = {}
    for j in range(1, n):
        d1, d2 = h0 - next(it), h0 - next(it); t1, t2 = th[j]
        inv_g = (1/d1 - 1/d2)/(1/t1 - 1/t2); g[j] = 1/inv_g; G[j, j] = (1/d1 - inv_g/t1)*g[j]; G[0, j] = G[j, 0] = np.sqrt(max(g[j], 0))
    for i, j in itertools.combinations(range(1, n), 2):
        dl = h0 - next(it); a = 1/th[i][0] + G[i, i]; b = 1/th[j][0] + G[j, j]
        G[i, j] = G[j, i] = (G[0, i]*G[0, j] - np.sqrt(max((dl*a - g[i])*(dl*b - g[j]), 0)))/dl
    return np.linalg.inv(G), G
OUT = {}
# (a) single hidden node: worst-case relative error factor of G_0j^2 for strengths (u, inf), u = theta*G_jj; closed form 1/u + 2 + 2u
f = lambda u: 1/u + 2 + 2*u; r = minimize_scalar(f, bounds=(1e-3, 10), method='bounded')
OUT['clamp_design'] = {'u_opt': r.x, 'u_opt_exact': 2**-0.5, 'factor': r.fun, 'factor_exact': 2 + 2*2**0.5}
# (b) Gaussian two-point tail threshold Delta/2, NOT a pointwise SD lower bound.
# The theorem implies max-over-two-models RMSE >= Delta/4.
def two_node(w, g0=1.0, g1=1.0): return np.array([[g0 + w, -w], [-w, g1 + w]])
sigma = 1e-9; rows = []
for w in np.logspace(-3, 0, 7):
    K = two_node(w); G = np.linalg.inv(K); Gjj = G[1, 1]
    th = {1: (2**-0.5/Gjj, 4/Gjj)}   # first strength at the clamp-optimal scale; second large
    E = protocol(2, th); clean = np.array([h(K, e) for e in E]); est = []
    for t in range(400):
        Kr, _ = reconstruct(2, clean + sigma*rng.standard_normal(len(E)), th); est.append(Kr[1, 1])
    sd = float(np.std(est)); lb = sigma/(2*np.sqrt(len(E))*G[0, 1]**2)
    rows.append({'w': float(w), 'G01': float(G[0, 1]), 'sd_K11_protocol': sd, 'two_point_tail_threshold': lb, 'two_point_minimax_rmse_lower_bound': lb/2, 'sd_over_tail_threshold': sd/lb})
OUT['weak_link'] = rows
# (c) hard family 2: grounded path, distance from port. G_0j decays geometrically; resolution of the far grounding degrades as G_0j^-2
rows = []
for n in (3, 5, 7, 9):
    edges = np.ones(n - 1); grounds = np.full(n, 0.5)
    K = np.diag(grounds); 
    for i in range(n - 1): K[i, i] += 1; K[i+1, i+1] += 1; K[i, i+1] -= 1; K[i+1, i] -= 1
    G = np.linalg.inv(K); th = {j: (2**-0.5/G[j, j], 4/G[j, j]) for j in range(1, n)}
    E = protocol(n, th); clean = np.array([h(K, e) for e in E]); est = []
    sig = 1e-12
    for t in range(300):
        Kr, _ = reconstruct(n, clean + sig*rng.standard_normal(len(E)), th); est.append(Kr[n-1, n-1])
    sd = float(np.std(est)); lb = sig/(2*np.sqrt(len(E))*G[0, n-1]**2)
    rows.append({'n': n, 'experiments': len(E), 'G0_far': float(G[0, n-1]), 'sd_K_far': sd, 'two_point_tail_threshold': lb, 'two_point_minimax_rmse_lower_bound': lb/2, 'sd_over_tail_threshold': sd/lb})
OUT['distance'] = rows
# (d) intervention strength: error of G_0j^2 versus u = theta*G_jj (second strength large), Monte Carlo vs closed form
K = two_node(0.3); G = np.linalg.inv(K); rows = []
for u in (0.05, 0.2, 2**-0.5, 2, 10):
    th = {1: (u/G[1, 1], 1e6/G[1, 1])}; E = protocol(2, th); clean = np.array([h(K, e) for e in E]); est = []
    for t in range(2000):
        _, Gr = reconstruct(2, clean + 1e-10*rng.standard_normal(len(E)), th); est.append(Gr[0, 1]**2)
    rows.append({'u': u, 'rel_sd_G01sq_over_sigma': float(np.std(est)/G[0, 1]**2/1e-10), 'Gjj_over_G01sq': float(G[1, 1]/G[0, 1]**2),
                 'rms_model': float(G[1, 1]/G[0, 1]**2*np.sqrt((1/u+2)**2 + ((1+u)**2/u)**2 + u**2))})
OUT['strength'] = rows
OUT['statistical_contract'] = {'tail_threshold_is_not_sd_floor': True, 'comparison_is_descriptive': True, 'raw_baseline_shared': True}
print(json.dumps({'weak_link_rows': len(OUT['weak_link']), 'distance_rows': len(OUT['distance']), 'strength_rows': len(OUT['strength'])})); json.dump(OUT, open(HERE / 'stability.json', 'w'), indent=1)
