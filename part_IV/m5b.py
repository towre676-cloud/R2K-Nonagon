import numpy as np, json
rng = np.random.default_rng(96)
def cue_counts(Nm, reps, W=2.0):
    v = []
    for _ in range(reps):
        A = (rng.standard_normal((Nm, Nm)) + 1j*rng.standard_normal((Nm, Nm)))/np.sqrt(2)
        Q, R = np.linalg.qr(A); Q = Q*(np.diag(R)/np.abs(np.diag(R)))
        th = np.sort(np.angle(np.linalg.eigvals(Q)) % (2*np.pi))*Nm/(2*np.pi)
        off = rng.uniform(0, W)
        for lo in np.arange(off, Nm - W + 1e-9, W): v.append(np.sum((th >= lo) & (th < lo + W)))
    return float(np.var(v))
out = {}
for Nm in (4, 6, 8, 10, 12, 16, 24, 50, 400):
    out[Nm] = cue_counts(Nm, max(200, 40000//Nm)); print(Nm, round(out[Nm], 4), flush=True)
t = 600269.68 + 5475
out['log_t_over_2pi'] = float(np.log(t/(2*np.pi)))
out['gue_sigma2_L2'] = float((np.log(2*np.pi*2) + np.euler_gamma + 1)/np.pi**2)
json.dump(out, open('m5b.json', 'w'), indent=1); print(out)
