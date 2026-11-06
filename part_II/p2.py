"""P2: Lehmer-pair collision under the de Bruijn-Newman heat flow (dt H = -dzz H), z = 2*gamma."""
import numpy as np, json, time
from flint import acb
t0 = time.time(); K = 300
zs = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(6709 - K, 2*K + 2)]); tz = time.time() - t0
i1 = K; z = 2*zs; g0 = z[i1+1] - z[i1]
dens = np.log(z[i1]/(4*np.pi))/(4*np.pi)
def S_of(z):
    oth = np.r_[z[:i1], z[i1+2:]]; a, b = z[i1], z[i1+1]
    return np.sum(1/((a - oth)*(b - oth))) + dens*(1/(a - z[0]) + 1/(z[-1] - b))
S0 = S_of(z); dt_lin = g0**2/(8 - 4*g0**2*S0)
def rhs(z):
    d = z[:, None] - z[None, :]; np.fill_diagonal(d, np.inf); return 2*np.sum(1/d, 1)
t = 0.0; h = -dt_lin/2000
while z[i1+1] - z[i1] > 0.15*g0:
    k1 = rhs(z); k2 = rhs(z + h/2*k1); k3 = rhs(z + h/2*k2); k4 = rhs(z + h*k3)
    z = z + h/6*(k1 + 2*k2 + 2*k3 + k4); t += h
g = z[i1+1] - z[i1]; S = S_of(z); rate = 8 - 4*g**2*S
t_coll = t - g**2/rate
out = {'gamma_6709': float(zs[i1]), 'gamma_6710': float(zs[i1+1]), 'gap_gamma': float(zs[i1+1]-zs[i1]), 'g_z': float(g0), 'S0': float(S0),
       'collision_time_linear': float(-dt_lin), 'collision_time_nbody': float(t_coll), 'window_zeros': 2*K + 2, 'zero_fetch_s': tz,
       'normalized_gap': float((zs[i1+1]-zs[i1])*np.log(zs[i1]/(2*np.pi))/(2*np.pi))}
# window-size convergence of S0
for KK in (50, 100, 200):
    zz = 2*zs[i1-KK:i1+KK+2]; j = KK; oth = np.r_[zz[:j], zz[j+2:]]; a, b = zz[j], zz[j+1]
    out[f'S0_window_{KK}'] = float(np.sum(1/((a-oth)*(b-oth))) + dens*(1/(a-zz[0]) + 1/(zz[-1]-b)))
print(out); json.dump(out, open('p2.json', 'w'), indent=1)
