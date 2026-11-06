"""X.3 (Rayleigh / Kirchhoff): resolution of edge weights deep in the network. Exact edge sensitivity dh/dw_ij = -(G_0i - G_0j)^2;
two-point floor for any protocol: sd(w_ij) >= sigma / (2 sqrt(m) sup_Theta |(G_Theta)_0i - (G_Theta)_0j| max(...)) <= ... we use the
admissible bound |(G_Theta)_0i - (G_Theta)_0j| <= max(G_0i, G_0j) (entrywise monotonicity), and compare with the protocol."""
import numpy as np, json, itertools, sys
sys.path.insert(0, '.')
exec(open('p8/stability.py').read().split('OUT = {}')[0])
rng = np.random.default_rng(31)
out = {'path': [], 'sensitivity_identity_check': None}
# identity check: d h / d w_ij = -(G_0i - G_0j)^2 by finite differences on a random network
K = np.array([[3., -1, -1, 0], [-1, 3, -1, -1], [-1, -1, 3.5, -1], [0, -1, -1, 2.5]]); G = np.linalg.inv(K); e = 1e-7
E = np.zeros((4, 4)); E[1, 1] = E[3, 3] = 1; E[1, 3] = E[3, 1] = -1
fd = (np.linalg.inv(K + e*E)[0, 0] - np.linalg.inv(K - e*E)[0, 0])/(2*e)
out['sensitivity_identity_check'] = {'finite_difference': fd, 'formula': -(G[0, 1] - G[0, 3])**2}
n = 9; Kp = np.diag(np.full(n, 0.5))
for i in range(n - 1): Kp[i, i] += 1; Kp[i+1, i+1] += 1; Kp[i, i+1] -= 1; Kp[i+1, i] -= 1
G = np.linalg.inv(Kp); th = {j: (2**-0.5/G[j, j], 4/G[j, j]) for j in range(1, n)}
Ep = protocol(n, th); clean = np.array([h(Kp, x) for x in Ep]); sig = 1e-12; ests = []
for t in range(300):
    Kr, _ = reconstruct(n, clean + sig*rng.standard_normal(len(Ep)), th); ests.append(Kr)
ests = np.array(ests); m = len(Ep)
for k in range(n - 1):
    sd = float(np.std(-ests[:, k, k+1])); dG = G[0, k] - G[0, k+1]; mx = max(G[0, k], G[0, k+1])
    out['path'].append({'edge': [k, k+1], 'G0k': float(G[0, k]), 'G0k_minus_G0k1': float(dG), 'sd_w': sd,
                        'floor_max': sig/(2*np.sqrt(m)*mx**2), 'sd_times_G0k1_sq': sd*G[0, k+1]**2/sig, 'sd_times_dG_sq': sd*dG**2/sig, 'sd_times_G0k_sq': sd*G[0, k]**2/sig})
for r in out['path']: print(r)
lg = np.log([r['sd_w'] for r in out['path']])
for key in ('G0k', 'G0k_minus_G0k1'):
    x = np.log([r[key] for r in out['path']]); out['slope_vs_' + key] = float(np.polyfit(x[1:], lg[1:], 1)[0])
x = np.log([G[0, k+1] for k in range(n - 1)]); out['slope_vs_G0_far_end'] = float(np.polyfit(x[1:], lg[1:], 1)[0])
print({k: v for k, v in out.items() if k.startswith('slope') or k.startswith('sens')})
json.dump(out, open('p10/edge_resolution.json', 'w'), indent=1)
