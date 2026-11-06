import numpy as np, json
from flint import acb
rng = np.random.default_rng(97)
Z = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(1000000, 20000)])
def Nsm(t): return t/(2*np.pi)*np.log(t/(2*np.pi*np.e)) + 7/8
U = Nsm(Z); U -= U[0]
print('mean spacing check', (U[-1]-U[0])/(len(U)-1))
def sigma2(pts, L, M=20000):
    lo = rng.uniform(pts[0], pts[-1] - L, M); c = np.searchsorted(pts, lo + L) - np.searchsorted(pts, lo); return float(np.var(c)), float(np.mean(c))
cue = []
for _ in range(25):
    A = (rng.standard_normal((400, 400)) + 1j*rng.standard_normal((400, 400)))/np.sqrt(2); Q, R = np.linalg.qr(A); Q = Q*(np.diag(R)/np.abs(np.diag(R)))
    cue.append(np.sort(np.angle(np.linalg.eigvals(Q)) % (2*np.pi))*400/(2*np.pi))
def sigma2_cue(L, M=800):
    v = []
    for th in cue:
        lo = rng.uniform(0, 400 - L, M); v += list(np.searchsorted(th, lo + L) - np.searchsorted(th, lo))
    return float(np.var(v))
out = {}
for L in (0.5, 1, 2, 4, 8, 16):
    out[L] = {'zeta': sigma2(U, L), 'cue': sigma2_cue(L), 'gue_asym': float((np.log(2*np.pi*L) + np.euler_gamma + 1)/np.pi**2)}
    print(L, out[L], flush=True)
json.dump(out, open('m5c.json', 'w'), indent=1)
