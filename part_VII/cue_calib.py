"""VII.3: calibrate the number-variance estimator of bk_multi.py on CUE samples of the same size (20,000 unfolded points),
where the exact answer is known (sine-kernel Sigma^2, finite-N corrections O(L^2/N^2) negligible at N = 2000)."""
import numpy as np, json, time
from scipy.linalg import eigvalsh_tridiagonal
from scipy.special import sici
rng = np.random.default_rng(777)
Ls = [0.5, 1, 2, 4, 8, 16]
def s2(pts, L, step=0.01):
    lo = np.arange(pts[0], pts[-1]-L, step); c = np.searchsorted(pts, lo+L) - np.searchsorted(pts, lo); return float(np.var(c))
def sigma2_gue(L):  # exact sine-kernel number variance
    x = 2*np.pi*L; Si, Ci = sici(x)
    return float((np.log(x) + np.euler_gamma + 1 - Ci - x/np.pi*0 + (L*(2/np.pi)*0)) / np.pi**2 + L*(1 - 2/np.pi*Si) - (1 - np.cos(x))/np.pi**2)  # replaced below by quadrature
def sigma2_quad(L, n=200001):
    x = np.linspace(0, L, n); y = (L - x)*np.sinc(x)**2; return float(L - 2*np.trapezoid(y, x))
EX = {L: sigma2_quad(L) for L in Ls}
N = 40000; reps = 32; est = {L: [] for L in Ls}; ses = {L: [] for L in Ls}; t0 = time.time()
def gue_unfolded():
    # Dumitriu-Edelman beta=2 tridiagonal; spectrum scaled to semicircle on [-1,1]; central 20,000 eigenvalues unfolded by the exact semicircle CDF
    d = rng.standard_normal(N)*np.sqrt(2); e = np.sqrt(rng.chisquare(2*np.arange(N-1, 0, -1)))
    lam = np.sort(eigvalsh_tridiagonal(d, e))[N//2-10000:N//2+10000]/(2*np.sqrt(2*N))   # central 20,000, semicircle radius 1
    x = np.clip(lam, -1, 1); F = N*(0.5 + (x*np.sqrt(1 - x*x) + np.arcsin(x))/np.pi)
    return np.sort(F)
for r in range(reps):
    pts = gue_unfolded()
    for L in Ls:
        est[L].append(s2(pts, L)); blocks = [s2(pts[i*2000:(i+1)*2000], L) for i in range(10)]; ses[L].append(float(np.std(blocks, ddof=1)/np.sqrt(10)))
    print(r, '%.0fs' % (time.time() - t0), flush=True)
out = {'ensemble': 'GUE (Dumitriu-Edelman tridiagonal), central 20,000 of N, semicircle-unfolded', 'N': N, 'points_per_sample': 20000, 'replicates': reps, 'rows': []}
for L in Ls:
    e = np.array(est[L]); out['rows'].append({'L': L, 'exact': EX[L], 'mean_estimate': float(e.mean()), 'bias': float(e.mean() - EX[L]), 'bias_se': float(e.std(ddof=1)/np.sqrt(reps)),
                                              'true_sd_of_estimate': float(e.std(ddof=1)), 'mean_block_se': float(np.mean(ses[L])), 'se_ratio_true_over_block': float(e.std(ddof=1)/np.mean(ses[L]))})
    print(out['rows'][-1])
json.dump(out, open('p7/cue_calib.json', 'w'), indent=1)
