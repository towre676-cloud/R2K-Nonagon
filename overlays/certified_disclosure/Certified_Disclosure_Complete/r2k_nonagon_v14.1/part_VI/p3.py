"""P3: sample-domain / exact passivity certification and enforcement."""
import numpy as np, json, mpmath as mp
from fractions import Fraction as Fr
from numpy.polynomial import polynomial as P
from scipy.optimize import minimize
rng = np.random.default_rng(303)
OUT = {}
# ---------- exact weighted-Hermite count of positive roots
def newton_sums(c, K):
    """c: Fractions, ascending coeffs of degree-m poly. power sums s_0..s_{K-1} of its roots (Newton identities)."""
    m = len(c) - 1; a = [c[m - j] / c[m] for j in range(m + 1)]
    s = [Fr(m)]
    for k in range(1, K):
        if k <= m: s.append(-(sum(a[j] * s[k - j] for j in range(1, k)) + k * a[k]))
        else: s.append(-sum(a[j] * s[k - j] for j in range(1, m + 1)))
    return s
def bareiss_minors(M):
    """leading principal minors of a symmetric Fraction matrix via fraction-free elimination (no pivoting)."""
    n = len(M); A = [row[:] for row in M]; prev = Fr(1); minors = []
    for k in range(n):
        if A[k][k] == 0: return None
        minors.append(A[k][k])
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) / prev
        prev = A[k][k]
    return minors
def charpoly_FL(M):
    """Faddeev-LeVerrier over Q: coefficients c_0..c_n (descending powers) of det(xI - M)."""
    n = len(M); I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    Mk = [[Fr(0)] * n for _ in range(n)]; c = [Fr(1)]
    for k in range(1, n + 1):
        Mk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) + c[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AM = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c.append(-sum(AM[i][i] for i in range(n)) / k)
    return c
def signature(M):
    c = charpoly_FL(M)  # real-rooted: Descartes is exact
    def changes(seq):
        sq = [x for x in seq if x != 0]; return sum(1 for i in range(len(sq) - 1) if (sq[i] > 0) != (sq[i + 1] > 0))
    n = len(c) - 1
    pos = changes(c); neg = changes([c[k] * (-1) ** (n - k) for k in range(n + 1)])
    return pos - neg
def positive_root_count(coeffs):
    c = [Fr(float(x)) for x in coeffs]
    while c[-1] == 0: c.pop()
    m = len(c) - 1; s = newton_sums(c, 2 * m + 1)
    H0 = [[s[i + j] for j in range(m)] for i in range(m)]; H1 = [[s[i + j + 1] for j in range(m)] for i in range(m)]
    s0, s1 = signature(H0), signature(H1)
    return (s0, s1, None if (s0 is None or s1 is None) else (s0 + s1) // 2)
# ---------- ladder impedance (PR by construction)
def ladder(n):
    # Z = R0 + 1/(s C1 + G1 + 1/(R2 + s L2 + 1/(s C2 + G2 + ...)))
    num, den = np.array([1.0]), np.array([0.0])  # start with open circuit: Y_end = 0 -> represent impedance chain from the far end
    Zn, Zd = np.array([rng.uniform(0.5, 2)]), np.array([1.0])  # terminating resistor
    for k in range(n):
        C, G = rng.uniform(0.2, 1.5), rng.uniform(0.0, 0.05)
        Yn = P.polyadd(P.polymul([G, C], Zn), Zd); Yd = Zn  # Y = sC+G + 1/Z
        Zn, Zd = Yd, Yn  # Z = 1/Y
        if k < n - 1:
            R, L = rng.uniform(0.0, 0.05), rng.uniform(0.2, 1.5)
            Zn = P.polyadd(Zn, P.polymul([R, L], Zd))  # Z += R + sL
    R0 = 0.02; Zn = P.polyadd(Zn, R0 * Zd)
    lead = Zd[-1]; return Zn / lead, Zd / lead
def to_pole_residue(N, D):
    p = P.polyroots(D); dD = P.polyder(D)
    q, rem = P.polydiv(N, D); d = q[0] if len(q) else 0.0
    r = P.polyval(p, rem) / P.polyval(p, dD); return float(np.real(d)), p, r
def from_pole_residue(d, p, r):
    D = np.real(P.polyfromroots(p)); N = d * D
    for k in range(len(p)): N = P.polyadd(N, r[k] * P.polyfromroots(np.delete(p, k)))
    return np.real(N), D
def Zeval(d, p, r, s): return d + np.sum(r / (np.asarray(s)[:, None] - p), 1)
def reZpoly(N, D):
    Nm = N * (-1.0) ** np.arange(len(N)); Dm = D * (-1.0) ** np.arange(len(D))
    E = P.polyadd(P.polymul(N, Dm), P.polymul(Nm, D))  # even polynomial in s
    Pc = np.array([E[2 * j] * (-1.0) ** j for j in range((len(E) + 1) // 2)])  # P(x), x = omega^2
    return Pc
def certify(d, p, r):
    N, D = from_pole_residue(d, p, r); Pc = reZpoly(N, D)
    s0, s1, npos = positive_root_count(Pc)
    # exact roots location in high precision
    mp.mp.dps = 50; rts = mp.polyroots([mp.mpf(Fr(float(x)).numerator) / mp.mpf(Fr(float(x)).denominator) for x in Pc[::-1]], maxsteps=400, extraprec=400)
    pos = sorted([float(mp.re(z)) for z in rts if abs(mp.im(z)) < 1e-30 and mp.re(z) > 0])
    return {'sigH0': s0, 'sigH1': s1, 'n_pos_roots': npos, 'P0': float(Pc[0]), 'band_edges_omega': [float(np.sqrt(x)) for x in pos], 'deg': len(Pc) - 1}
nlad = 5
N0, D0 = ladder(nlad); d0, p0, r0 = to_pole_residue(N0, D0)
OUT['true_model'] = certify(d0, p0, r0)
# perturb residues: macromodel fitting error
om = np.logspace(-2, 1.5, 20000)
for attempt in range(400):
    xi = rng.standard_normal(len(r0)) + 1j * rng.standard_normal(len(r0))
    r1 = r0.copy()
    for k in range(len(p0)):
        j = int(np.argmin(np.abs(p0 - np.conj(p0[k]))))
        if j >= k: r1[k] = r0[k] * (1 + 0.08 * xi[k]); r1[j] = np.conj(r1[k]) if j != k else r1[k].real
    d1 = d0
    if Zeval(d1, p0, r1, 1j * om).real.min() < 0: break
OUT['perturbation_attempts'] = attempt + 1
OUT['perturbed'] = certify(d1, p0, r1)
om = np.logspace(-2, 1.5, 20000); reZ = Zeval(d1, p0, r1, 1j * om).real
OUT['perturbed_grid_min_reZ'] = float(reZ.min()); OUT['perturbed_grid_violating_fraction'] = float(np.mean(reZ < 0))
print('true', OUT['true_model']); print('pert', OUT['perturbed'], OUT['perturbed_grid_min_reZ'])
# ---------- enforcement: minimal change in |Z| on a log grid subject to Re Z >= eps at band points
pairs = [k for k in range(len(p0)) if p0[k].imag > 1e-12] + [k for k in range(len(p0)) if abs(p0[k].imag) <= 1e-12]
def unpack(x):
    r = r1.copy().astype(complex); i = 0
    for k in pairs:
        if abs(p0[k].imag) > 1e-12:
            r[k] = r1[k] + x[i] + 1j * x[i + 1]; j = int(np.argmin(np.abs(p0 - np.conj(p0[k])))); r[j] = np.conj(r[k]); i += 2
        else: r[k] = r1[k] + x[i]; i += 1
    return d1 + x[i], r
def ZE(x, s):
    dd, rr = unpack(x); return Zeval(dd, p0, rr, s)
nv = sum(2 if abs(p0[k].imag) > 1e-12 else 1 for k in pairs) + 1
og = np.logspace(-2, 1.5, 400); Z1g = Zeval(d1, p0, r1, 1j * og)
cons_om = []; hist = []
cur = OUT['perturbed']; x = np.zeros(nv)
for it in range(8):
    if cur['n_pos_roots'] == 0 and cur['P0'] > 0: break
    be = cur['band_edges_omega']
    for a_, b_ in zip(be[0::2], be[1::2]): cons_om += list(np.linspace(a_, b_, 12))
    if len(be) % 2 == 1: cons_om += list(np.linspace(be[-1], be[-1] * 3, 12))
    co = np.array(cons_om); eps = 1e-4
    obj = lambda x: np.sum(np.abs(ZE(x, 1j * og) - Z1g) ** 2 / np.abs(Z1g) ** 2)
    cons = {'type': 'ineq', 'fun': lambda x, co=co: ZE(x, 1j * co).real - eps}
    sol = minimize(obj, x, constraints=[cons], method='SLSQP', options={'maxiter': 500, 'ftol': 1e-14})
    x = sol.x; dE, rE = unpack(x); cur = certify(dE, p0, rE); cur['grid_min_reZ'] = float(Zeval(dE, p0, rE, 1j * om).real.min())
    hist.append({'iter': it, 'n_constraints': len(co), 'obj': float(sol.fun), 'cert': cur}); print(hist[-1])
OUT['enforcement'] = hist
dE, rE = unpack(x)
OUT['final_rel_change_Z'] = float(np.max(np.abs(Zeval(dE, p0, rE, 1j * og) - Z1g) / np.abs(Z1g)))
OUT['final_grid_min_reZ'] = float(Zeval(dE, p0, rE, 1j * om).real.min())
OUT['true_vs_pert_rel'] = float(np.max(np.abs(Zeval(d0, p0, r0, 1j * og) - Z1g) / np.abs(Z1g)))
print('final change', OUT['final_rel_change_Z'], 'grid min', OUT['final_grid_min_reZ'])
# ---------- rank-two conjecture, Foster/Pick setting
def foster_test(sampling, n=6, trials=50):
    rows = []
    for t in range(trials):
        om_ = np.sort(rng.uniform(0.5, 5, n)); rr = rng.uniform(0.3, 2, n); kneg = rng.integers(n); rr[kneg] = -rng.uniform(0.01, 0.2)
        if sampling == 'generic': s = rng.uniform(0.1, 1.0, 3 * n) + 1j * rng.uniform(-6, 6, 3 * n)
        else:
            gaps = np.diff(om_); mg = np.minimum(np.r_[gaps, np.inf], np.r_[np.inf, gaps]); dl = 0.02 * mg
            s = np.r_[dl + 1j * om_, dl - 1j * om_, 0.05 * mg + 1j * om_]
        A = []
        for k in range(n):
            gp = 1 / (s - 1j * om_[k]); gm = 1 / (s + 1j * om_[k]); A.append(np.outer(gp, gp.conj()) + np.outer(gm, gm.conj()))
        Q = np.array([[np.real(np.trace(A[k] @ A[l])) for l in range(n)] for k in range(n)])
        # min (r'-r)^T Q (r'-r), r' >= 0
        res = minimize(lambda y: (y - rr) @ Q @ (y - rr), np.maximum(rr, 0), jac=lambda y: 2 * Q @ (y - rr), bounds=[(0, None)] * n, method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-13})
        dr = res.x - rr; changed = int(np.sum(np.abs(dr) > 1e-6 * np.max(np.abs(rr))))
        Pc = sum(res.x[k] * A[k] for k in range(n)); mev = np.linalg.eigvalsh((Pc + Pc.conj().T) / 2).min() / np.abs(np.linalg.eigvalsh(Pc)).max()
        offdiag = np.max(np.abs(Q - np.diag(np.diag(Q))) / np.sqrt(np.outer(np.diag(Q), np.diag(Q))))
        others = np.delete(np.abs(dr), kneg); spill = float(others.max() / abs(dr[kneg]))
        rows.append([changed, 2 * changed, float(mev), float(offdiag), float(abs(res.x[kneg])), spill])
    rows = np.array(rows)
    return {'rank2_fraction': float(np.mean(rows[:, 0] == 1)), 'median_changed': float(np.median(rows[:, 0])), 'max_changed': int(rows[:, 0].max()),
            'min_eig_rel_after': float(rows[:, 2].min()), 'median_gram_coherence': float(np.median(rows[:, 3])), 'median_spill': float(np.median(rows[:, 5])), 'max_spill': float(np.max(rows[:, 5])), 'neg_residue_set_to_zero_frac': float(np.mean(rows[:, 4] < 1e-9))}
OUT['foster_generic'] = foster_test('generic'); OUT['foster_resonance'] = foster_test('resonance')
print('foster generic', OUT['foster_generic']); print('foster resonance', OUT['foster_resonance'])
json.dump(OUT, open('p3.json', 'w'), indent=1, default=float)
# --- diagnostics of the enforcement change
Zt = Zeval(d0, p0, r0, 1j * og); ZEg = Zeval(dE, p0, rE, 1j * og)
OUT['enf_change_normwise'] = float(np.max(np.abs(ZEg - Z1g)) / np.max(np.abs(Z1g)))
OUT['pert_vs_true_normwise'] = float(np.max(np.abs(Z1g - Zt)) / np.max(np.abs(Zt)))
OUT['enf_vs_true_normwise'] = float(np.max(np.abs(ZEg - Zt)) / np.max(np.abs(Zt)))
OUT['enf_change_relative_argmax_omega'] = float(og[np.argmax(np.abs(ZEg - Z1g) / np.abs(Z1g))]); OUT['minabsZ_on_grid'] = float(np.abs(Z1g).min())
OUT['rel_change_median'] = float(np.median(np.abs(ZEg - Z1g) / np.abs(Z1g)))
print({k: OUT[k] for k in ('enf_change_normwise', 'pert_vs_true_normwise', 'enf_vs_true_normwise', 'enf_change_relative_argmax_omega', 'minabsZ_on_grid', 'rel_change_median', 'true_vs_pert_rel')})
json.dump(OUT, open('p3.json', 'w'), indent=1, default=float)
oc = np.linspace(1.3, 2.4, 600)
OUT['curves'] = {'omega': oc.tolist(), 'true': Zeval(d0, p0, r0, 1j*oc).real.tolist(), 'perturbed': Zeval(d1, p0, r1, 1j*oc).real.tolist(), 'enforced': Zeval(dE, p0, rE, 1j*oc).real.tolist()}
json.dump(OUT, open('p3.json', 'w'), indent=1, default=float)
