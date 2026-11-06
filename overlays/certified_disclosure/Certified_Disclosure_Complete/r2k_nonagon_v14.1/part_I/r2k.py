"""R2K: Rank-Two Kernel Calculus -- numerical program.
Parts A-E. Writes results.json."""
import numpy as np, scipy.linalg as sla, json, time, mpmath as mp
from numpy.polynomial import polynomial as P
rng = np.random.default_rng(20261029)
R = {}
def rel(a, b): return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))
def inertia(M, tol=None):
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    tol = tol if tol is not None else 1e-9 * np.max(np.abs(ev))
    return (int(np.sum(ev > tol)), int(np.sum(ev < -tol)), int(np.sum(np.abs(ev) <= tol)))
def disp_sv(D): s = np.linalg.svd(D, compute_uv=False); return [float(x) for x in s / s[0]]

# ---------------- PART A: the nonagon ----------------
A = {}
n = 6
lam = np.sort(rng.uniform(-3, 3, n)); w = rng.uniform(0.5, 2, n) * rng.choice([-1, 1], n)
r = lambda z: np.sum(w / (np.asarray(z)[..., None] - lam), axis=-1)
x = np.linspace(-4.1, 4.3, n) + 0.37j; y = np.linspace(-3.9, 4.0, n) - 0.41j
Cx = 1 / (x[:, None] - lam); Cy = 1 / (y[:, None] - lam)
L = (r(x)[:, None] - r(y)[None, :]) / (x[:, None] - y[None, :])
A['lowner_cauchy_factor'] = rel(L, -Cx @ np.diag(w) @ Cy.T)
# Cauchy double alternant
Cxy = 1 / (x[:, None] - y[None, :])
num = np.prod([(x[j] - x[i]) * (y[i] - y[j]) for i in range(n) for j in range(i + 1, n)])
den = np.prod(x[:, None] - y[None, :])
A['cauchy_det'] = float(abs(np.linalg.det(Cxy) - num / den) / abs(num / den))
# closed form det L
def cdet(u, v):
    nn = len(u); num = np.prod([(u[j] - u[i]) * (v[i] - v[j]) for i in range(nn) for j in range(i + 1, nn)])
    return num / np.prod(u[:, None] - v[None, :])
detL_closed = (-1) ** n * cdet(x, lam) * np.prod(w) * cdet(y, lam)
A['lowner_det_closed'] = float(abs(np.linalg.det(L) - detL_closed) / abs(detL_closed))
# displacement
A['disp_lowner'] = disp_sv(np.diag(x) @ L - L @ np.diag(y))
A['disp_cauchy'] = disp_sv(np.diag(x) @ Cxy - Cxy @ np.diag(y))
# Hankel / Prony
mom = np.array([np.sum(w * lam ** j) for j in range(2 * n)])
H0 = np.array([[mom[i + j] for j in range(n)] for i in range(n)])
H1 = np.array([[mom[i + j + 1] for j in range(n)] for i in range(n)])
V = lam[:, None] ** np.arange(n)[None, :]
A['hankel_vandermonde'] = rel(H0, V.T @ np.diag(w) @ V)
pr = np.sort(sla.eigvals(H1, H0).real)
A['prony_pole_err'] = float(np.max(np.abs(pr - lam))); A['cond_H0'] = float(np.linalg.cond(H0)); A['cond_L'] = float(np.linalg.cond(L))
Z = np.diag(np.ones(n - 1), -1)
A['disp_hankel'] = disp_sv(Z @ H0 - H0 @ Z.T)
# Bezoutian = V^{-1} D W D V^{-T}; Lander inverse
p = P.polyfromroots(lam)
pk = [P.polyfromroots(np.delete(lam, k)) for k in range(n)]
q = sum(w[k] * pk[k] for k in range(n))
def bez(pc, qc):
    m = len(pc) - 1; pc = np.r_[pc, 0]; qc = np.r_[qc, np.zeros(m + 1 - len(qc)), 0]
    # (p(x)q(y)-p(y)q(x))/(x-y) via coefficient recurrence
    B = np.zeros((m, m))
    C = np.outer(pc[:m + 1], qc[:m + 1]) - np.outer(qc[:m + 1], pc[:m + 1])  # coeff x^i y^j
    # divide by (x - y): B_{i,j} with sum B_{ij} x^i y^j (x-y) = C
    for i in range(m - 1, -1, -1):
        for j in range(m):
            # coefficient of x^{i+1} y^j: B_{i,j} - B_{i+1,j-1} = C_{i+1,j}
            B[i, j] = C[i + 1, j] + (B[i + 1, j - 1] if (i + 1 < m and j >= 1) else 0)
    return B
Bz = bez(p, q)
dp = P.polyval(lam, P.polyder(p))
Vi = np.linalg.inv(V)  # columns = coeffs of Lagrange... V.T[i,k]=lam_k^i
A['bezout_factor'] = rel(Bz, Vi @ np.diag(dp * w * dp) @ Vi.T)
qt = sum((1 / (w[k] * dp[k] ** 2)) * pk[k] for k in range(n))
A['lander_inverse'] = rel(np.linalg.inv(H0), bez(p, qt))
xs = rng.uniform(-2, 2, n); ys = rng.uniform(-2, 2, n)
Bxy = (P.polyval(xs, p)[:, None] * P.polyval(ys, q)[None, :] - P.polyval(ys, p)[None, :] * P.polyval(xs, q)[:, None]) / (xs[:, None] - ys[None, :])
A['disp_bezout'] = disp_sv(np.diag(xs) @ Bxy - Bxy @ np.diag(ys))
# barycentric from Lowner null vector (AAA identity)
supp = np.linspace(-3.7, 3.9, n + 1) + 0.2j; test = np.linspace(-3.5, 3.6, n) - 0.3j
Lb = (r(test)[:, None] - r(supp)[None, :]) / (test[:, None] - supp[None, :])
beta = np.linalg.svd(Lb)[2].conj()[-1]
zz = rng.uniform(-5, 5, 200) + 1j * rng.uniform(-2, 2, 200)
fs = r(supp)
bary = (np.sum(beta * fs / (zz[:, None] - supp), 1)) / np.sum(beta / (zz[:, None] - supp), 1)
A['barycentric_err'] = rel(bary, r(zz))
# Pick <-> Nevanlinna via Cayley
cpos = rng.uniform(0.3, 1.5, n); cpos[[1, 4]] *= -1
f = lambda z: 0.3 + 0.7 * z + np.sum(cpos / (lam - np.asarray(z)[..., None]), -1)
Nz = 12; zp = rng.uniform(-3, 3, Nz) + 1j * rng.uniform(0.2, 2, Nz)
fz = f(zp)
Nev = (fz[:, None] - fz.conj()[None, :]) / (zp[:, None] - zp.conj()[None, :])
zeta_ = (zp - 1j) / (zp + 1j); S = (fz - 1j) / (fz + 1j)
Pick = (1 - S[:, None] * S.conj()[None, :]) / (1 - zeta_[:, None] * zeta_.conj()[None, :])
Dm = np.diag((zp + 1j) / (fz + 1j))
A['pick_cayley'] = rel(Pick, Dm @ Nev @ Dm.conj().T)
A['pick_inertia'] = inertia(Pick); A['nev_inertia'] = inertia(Nev)
A['pick_expected_neg'] = int(np.sum(cpos < 0))
A['disp_pick'] = disp_sv(Pick - np.diag(zeta_) @ Pick @ np.diag(zeta_.conj()))
# Christoffel-Darboux & Golub-Welsch (Legendre)
nL = 12
ak = np.array([k / np.sqrt(4 * k * k - 1) for k in range(1, nL + 1)])
def onp(xv, m):
    out = [np.full_like(xv, 1 / np.sqrt(2))]
    out.append(xv * out[0] / ak[0])
    for k in range(1, m):
        out.append((xv * out[k] - ak[k - 1] * out[k - 1]) / ak[k])
    return np.array(out)
xa = rng.uniform(-1, 1, 9); ya = rng.uniform(-1, 1, 9)
Px = onp(xa, nL); Py = onp(ya, nL)
Kdir = Px[:nL].T @ Py[:nL]
Kcd = ak[nL - 1] * (np.outer(Px[nL], Py[nL - 1]) - np.outer(Px[nL - 1], Py[nL])) / (xa[:, None] - ya[None, :])
A['cd_formula'] = rel(Kcd, Kdir)
A['disp_cd'] = disp_sv(np.diag(xa) @ Kdir - Kdir @ np.diag(ya))
Jl = np.diag(ak[:nL - 1], 1) + np.diag(ak[:nL - 1], -1)
ev, U = np.linalg.eigh(Jl); wgw = 2 * U[0] ** 2
xg, wg = np.polynomial.legendre.leggauss(nL)
A['golub_welsch_nodes'] = float(np.max(np.abs(ev - xg))); A['golub_welsch_weights'] = float(np.max(np.abs(wgw - wg)))
Pn = onp(xg, nL)[:nL]
Kn = Pn.T @ Pn
A['cd_at_nodes'] = rel(Kn, np.diag(1 / wg))
R['A'] = A

# ---------------- PART B: conservation of inertia ----------------
def lanczos(lv, wv):
    m = len(lv); Q = np.zeros((m, m)); a = np.zeros(m); b = np.zeros(m - 1)
    qv = np.sqrt(wv / wv.sum()); Q[:, 0] = qv; prev = np.zeros(m); beta_ = 0
    for k in range(m):
        v = lv * Q[:, k] - beta_ * prev
        a[k] = Q[:, k] @ v; v -= a[k] * Q[:, k]
        v -= Q[:, :k + 1] @ (Q[:, :k + 1].T @ v)
        if k < m - 1:
            beta_ = np.linalg.norm(v); b[k] = beta_; prev = Q[:, k]; Q[:, k + 1] = v / beta_
    return a, b
B = {'trials': 0, 'agree_lowner': 0, 'agree_hankel': 0, 'agree_bezout': 0, 'agree_cauchy_index': 0, 'agree_pick': 0}
for tr in range(2000):
    nn = rng.integers(2, 8)
    lv = np.sort(rng.uniform(-3, 3, nn))
    if np.min(np.diff(lv)) < 0.08: continue
    wv = rng.uniform(0.3, 2, nn) * rng.choice([-1, 1], nn)
    npos, nneg = int(np.sum(wv > 0)), int(np.sum(wv < 0))
    rr = lambda z: np.sum(wv / (np.asarray(z)[..., None] - lv), -1)
    rd = lambda z: -np.sum(wv / (np.asarray(z)[..., None] - lv) ** 2, -1)
    N = nn + 4; xp = np.sort(rng.uniform(-4, 4, N))
    if np.min(np.abs(xp[:, None] - lv)) < 0.05 or np.min(np.diff(xp)) < 0.05: continue
    B['trials'] += 1
    Ls = (rr(xp)[:, None] - rr(xp)[None, :]) / (xp[:, None] - xp[None, :] + np.eye(N))
    Ls[np.diag_indices(N)] = rd(xp)
    if inertia(Ls, 1e-10)[:2] == (nneg, npos): B['agree_lowner'] += 1
    mo = np.array([np.sum(wv * lv ** j) for j in range(2 * nn)])
    Hh = np.array([[mo[i + j] for j in range(nn)] for i in range(nn)])
    if inertia(Hh, 1e-12)[:2] == (npos, nneg): B['agree_hankel'] += 1
    pp = P.polyfromroots(lv); qq = sum(wv[k] * P.polyfromroots(np.delete(lv, k)) for k in range(nn))
    if inertia(bez(pp, qq), 1e-12)[:2] == (npos, nneg): B['agree_bezout'] += 1
    # Cauchy index from samples only: count pole crossings on a fine grid
    g = np.linspace(-3.5, 3.5, 200001); vals = rr(g); idx = 0
    s = np.sign(vals); ch = np.where(s[:-1] != s[1:])[0]
    for c in ch:
        if abs(vals[c]) > 1 and abs(vals[c + 1]) > 1:  # jump through infinity
            idx += 1 if vals[c] < 0 else -1
    if idx == npos - nneg: B['agree_cauchy_index'] += 1
    # Foster/Pick in RHP
    om = rng.uniform(0.5, 4, nn); rv = np.abs(wv) * np.sign(wv)
    Zf = lambda s_: np.sum(rv * 2 * np.asarray(s_)[..., None] / (np.asarray(s_)[..., None] ** 2 + om ** 2), -1)
    sp = rng.uniform(0.1, 1.0, 3 * nn) + 1j * rng.uniform(-5, 5, 3 * nn)
    Zs = Zf(sp); Pk = (Zs[:, None] + Zs.conj()[None, :]) / (sp[:, None] + sp.conj()[None, :])
    if inertia(Pk, 1e-11)[:2] == (2 * npos, 2 * nneg): B['agree_pick'] += 1
R['B'] = B

# ---------------- PART C: frequency-sweep structural identification ----------------
C = {}
nC = 10
m_true = rng.uniform(0.5, 2.0, nC); k_true = rng.uniform(1.0, 5.0, nC)
K = np.zeros((nC, nC))
for i in range(nC):
    K[i, i] = k_true[i] + (k_true[i + 1] if i + 1 < nC else 0)
    if i + 1 < nC: K[i, i + 1] = K[i + 1, i] = -k_true[i + 1]
M = np.diag(m_true)
Hfun = lambda z: np.array([np.linalg.solve(K - zz_ * M, np.eye(nC)[0])[0] for zz_ in np.atleast_1d(z)])
lam_true = np.sort(sla.eigh(K, M, eigvals_only=True))
C['lam_true'] = lam_true.tolist()
def lowner_poles(mu, nu, Hm, Hn):
    L_ = (Hm[:, None] - Hn[None, :]) / (mu[:, None] - nu[None, :])
    Ls_ = (mu[:, None] * Hm[:, None] - nu[None, :] * Hn[None, :]) / (mu[:, None] - nu[None, :])
    return sla.eigvals(Ls_, L_), L_, Ls_
# v1 ring design (kept as the baseline)
th = np.linspace(0, 2 * np.pi, 2 * nC, endpoint=False) + 0.1
ctr, rad = (lam_true[0] + lam_true[-1]) / 2, (lam_true[-1] - lam_true[0]) / 2 * 1.35
zr_ = ctr + rad * np.exp(1j * th)
p_ring, Lr_, _ = lowner_poles(zr_[0::2], zr_[1::2], Hfun(zr_[0::2]), Hfun(zr_[1::2]))
C['v1_ring_cond'] = float(np.linalg.cond(Lr_)); C['v1_ring_pole_relerr'] = float(np.max(np.abs(np.sort(p_ring.real) - lam_true) / lam_true))
# v2 conjugate frequency sweep z=(omega + i eta)^2, left = sweep, right = conj
wsw = np.linspace(np.sqrt(lam_true[0]) * 0.8, np.sqrt(lam_true[-1]) * 1.1, nC)
mu = (wsw + 0.05j) ** 2; nu = mu.conj(); zs = np.r_[mu, nu]; Hs = Hfun(zs)
p_sw, L_, Ls_ = lowner_poles(mu, nu, Hs[:nC], Hs[nC:])
C['v2_sweep_cond'] = float(np.linalg.cond(L_)); C['v2_sweep_pole_relerr'] = float(np.max(np.abs(np.sort(p_sw.real) - lam_true) / lam_true))
# v3 displacement-equilibrated pencil
for _ in range(30):
    dr = 1 / np.sqrt(np.max(np.abs(L_), 1)); L_ = dr[:, None] * L_; Ls_ = dr[:, None] * Ls_
    dc = 1 / np.sqrt(np.max(np.abs(L_), 0)); L_ = L_ * dc; Ls_ = Ls_ * dc
polesC = np.sort(sla.eigvals(Ls_, L_).real)
C['v3_equil_cond'] = float(np.linalg.cond(L_)); C['v3_equil_pole_relerr'] = float(np.max(np.abs(polesC - lam_true) / lam_true))
# v4 Cauchy-Gauss-Newton polish
Cm = 1 / (polesC[None, :] - zs[:, None])
wfit = np.linalg.lstsq(np.r_[Cm.real, Cm.imag], np.r_[Hs.real, Hs.imag], rcond=None)[0]
for _ in range(6):
    Cm = 1 / (polesC[None, :] - zs[:, None]); res = Cm @ wfit - Hs
    Jg = np.c_[Cm, -wfit * Cm ** 2]; Jr = np.r_[Jg.real, Jg.imag]
    dd = np.linalg.lstsq(Jr, -np.r_[res.real, res.imag], rcond=None)[0]; wfit += dd[:nC]; polesC += dd[nC:]
o_ = np.argsort(polesC); polesC = polesC[o_]; wfit = wfit[o_]
C['v4_fit_residual'] = float(np.max(np.abs(1 / (polesC[None, :] - zs[:, None]) @ wfit - Hs)))
C['pole_relerr_double'] = float(np.max(np.abs(polesC - lam_true) / lam_true))
# intrinsic floor: first-order forward error for data perturbation of size eps*|H|
Jg = np.c_[1 / (polesC[None, :] - zs[:, None]), -wfit / (polesC[None, :] - zs[:, None]) ** 2]
Jr = np.r_[Jg.real, Jg.imag]; Jp = np.linalg.pinv(Jr)
C['intrinsic_pole_floor'] = float(np.max(np.sum(np.abs(Jp[nC:]), 1) * np.finfo(float).eps * np.max(np.abs(Hs)) / lam_true))
C['residues'] = wfit.tolist()
C['lowner_cond'] = C['v2_sweep_cond']
a_, b_ = lanczos(polesC, wfit)
m1 = 1 / wfit.sum(); k1 = 1 / np.sum(wfit / polesC)
mr = [m1]; kr = [k1]
for i in range(nC - 1):
    kn = a_[i] * mr[i] - kr[i]; mn = kn ** 2 / (b_[i] ** 2 * mr[i]); kr.append(kn); mr.append(mn)
C['mass_relerr_double'] = float(np.max(np.abs(np.array(mr) - m_true) / m_true))
C['spring_relerr_double'] = float(np.max(np.abs(np.array(kr) - k_true) / k_true))
C['end_consistency_double'] = float(abs(a_[-1] * mr[-1] - kr[-1]) / kr[-1])
C['m_true'] = m_true.tolist(); C['k_true'] = k_true.tolist(); C['m_rec'] = mr; C['k_rec'] = kr
# mpmath 60-digit version
mp.mp.dps = 60
Km = mp.matrix(nC, nC); Mm = mp.matrix(nC, nC)
for i in range(nC):
    Km[i, i] = mp.mpf(k_true[i]) + (mp.mpf(k_true[i + 1]) if i + 1 < nC else 0); Mm[i, i] = mp.mpf(m_true[i])
    if i + 1 < nC: Km[i, i + 1] = Km[i + 1, i] = -mp.mpf(k_true[i + 1])
def Hmp(z):
    A_ = Km - z * Mm; return mp.lu_solve(A_, mp.matrix([1] + [0] * (nC - 1)))[0]
zsm = [mp.mpf(ctr) + mp.mpf(rad) * mp.expj(mp.mpf(t_)) for t_ in th]
mum, num_ = zsm[0::2], zsm[1::2]
Hmu = [Hmp(z) for z in mum]; Hnu = [Hmp(z) for z in num_]
Lm = mp.matrix(nC, nC); Lsm = mp.matrix(nC, nC)
for i in range(nC):
    for j in range(nC):
        Lm[i, j] = (Hmu[i] - Hnu[j]) / (mum[i] - num_[j])
        Lsm[i, j] = (mum[i] * Hmu[i] - num_[j] * Hnu[j]) / (mum[i] - num_[j])
evm = mp.eig(mp.inverse(Lm) * Lsm, left=False, right=False)
pm = sorted([mp.re(e) for e in evm])
lam_ref = sorted(mp.eig(mp.inverse(Mm) * Km, left=False, right=False), key=lambda e: mp.re(e))
C['pole_relerr_mp60'] = float(max(abs(pm[i] - mp.re(lam_ref[i])) / abs(mp.re(lam_ref[i])) for i in range(nC)))
# residues mp
Cmm = mp.matrix(2 * nC, nC); rhs = mp.matrix(2 * nC, 1); zall = mum + num_; Hall = Hmu + Hnu
for i in range(2 * nC):
    for k in range(nC): Cmm[i, k] = 1 / (pm[k] - zall[i])
    rhs[i] = Hall[i]
wm = mp.lu_solve(Cmm[:nC, :] if False else (Cmm.H * Cmm), Cmm.H * rhs)
wm = [mp.re(wm[k]) for k in range(nC)]
# mp lanczos
q = [mp.sqrt(wm[k] / sum(wm)) for k in range(nC)]
Qs = [q]; am = []; bm = []; prev = [0] * nC; bet = 0
for kk in range(nC):
    v = [pm[i] * Qs[kk][i] - bet * prev[i] for i in range(nC)]
    al = mp.fsum(Qs[kk][i] * v[i] for i in range(nC)); am.append(al)
    v = [v[i] - al * Qs[kk][i] for i in range(nC)]
    for Qj in Qs:
        c_ = mp.fsum(Qj[i] * v[i] for i in range(nC)); v = [v[i] - c_ * Qj[i] for i in range(nC)]
    if kk < nC - 1:
        bet = mp.sqrt(mp.fsum(vi ** 2 for vi in v)); bm.append(bet); prev = Qs[kk]; Qs.append([vi / bet for vi in v])
m1m = 1 / sum(wm); k1m = 1 / mp.fsum(wm[k] / pm[k] for k in range(nC))
mrm = [m1m]; krm = [k1m]
for i in range(nC - 1):
    kn = am[i] * mrm[i] - krm[i]; mn = kn ** 2 / (bm[i] ** 2 * mrm[i]); krm.append(kn); mrm.append(mn)
C['mass_relerr_mp60'] = float(max(abs(mrm[i] - mp.mpf(m_true[i])) / m_true[i] for i in range(nC)))
C['spring_relerr_mp60'] = float(max(abs(krm[i] - mp.mpf(k_true[i])) / k_true[i] for i in range(nC)))
R['C'] = C; print('A',json.dumps(R['A'],default=float)); print('B',R['B']); print('C',{k:v for k,v in C.items() if 'err' in k or 'cond' in k or 'cons' in k}, flush=True)
mp.mp.dps = 15

# ---------------- PART D: delay-system stability by contour moments ----------------
D = {}
zeta_d, tau = 0.1, 1.0
fD = lambda s, kap: s ** 2 + 2 * zeta_d * s + 1 + kap * np.exp(-s * tau)
dfD = lambda s, kap: 2 * s + 2 * zeta_d - kap * tau * np.exp(-s * tau)
gx, gw = np.polynomial.legendre.leggauss(400)
def rect_moments(kap, x0, x1, y0, y1, nm=8):
    corners = [x0 + 1j * y0, x1 + 1j * y0, x1 + 1j * y1, x0 + 1j * y1, x0 + 1j * y0]
    mus = np.zeros(nm, complex)
    for a0, a1 in zip(corners[:-1], corners[1:]):
        s = (a0 + a1) / 2 + (a1 - a0) / 2 * gx; ds = (a1 - a0) / 2 * gw
        g = dfD(s, kap) / fD(s, kap) * ds
        for k in range(nm): mus[k] += np.sum(s ** k * g)
    return mus / (2j * np.pi)
def rhp_count(kap): return rect_moments(kap, 0.0, 6.0, -40.0, 40.0, 1)[0]
def locate(kap, x0, x1, y0, y1):
    mu_ = rect_moments(kap, x0, x1, y0, y1, 12); nz = int(round(mu_[0].real))
    if nz == 0: return np.array([])
    H0_ = np.array([[mu_[i + j] for j in range(nz)] for i in range(nz)])
    H1_ = np.array([[mu_[i + j + 1] for j in range(nz)] for i in range(nz)])
    rt = sla.eigvals(H1_, H0_)
    for _ in range(6): rt = rt - fD(rt, kap) / dfD(rt, kap)
    return rt
kaps = np.linspace(0.0, 1.0, 101); counts = [float(rhp_count(k_).real) for k_ in kaps]
D['kappa_grid'] = kaps.tolist(); D['rhp_counts'] = [round(c, 12) for c in counts]
D['count_integrality_dev'] = float(max(abs(c - round(c)) for c in counts))
# crossing: kappa with rightmost root on imaginary axis
from scipy.optimize import brentq
def rightmost_re(kap):
    rt = locate(kap, -0.6, 0.6, 0.2, 3.0); return float(np.max(rt.real))
kc_dl = brentq(rightmost_re, 0.05, 0.95, xtol=1e-17, rtol=1e-15)
# analytic: 2 zeta w = kap sin(w tau), 1 - w^2 = -kap cos(w tau)
wg_ = np.linspace(0.05, 6, 200001)
Gc = (-(1 - wg_ ** 2) - 2j * zeta_d * wg_) * np.exp(1j * wg_ * tau)
cand = []
for i in np.where(np.sign(Gc.imag[:-1]) != np.sign(Gc.imag[1:]))[0]:
    wr = brentq(lambda wv: ((-(1 - wv ** 2) - 2j * zeta_d * wv) * np.exp(1j * wv * tau)).imag, wg_[i], wg_[i + 1], xtol=1e-16)
    kv = ((-(1 - wr ** 2) - 2j * zeta_d * wr) * np.exp(1j * wr * tau)).real
    if kv > 0: cand.append((kv, wr))
kc_val, wc = min(cand)
kc_an = float(np.hypot(1 - wc ** 2, 2 * zeta_d * wc))
D['kappa_crit_contour'] = kc_dl; D['kappa_crit_analytic'] = kc_an; D['omega_crit'] = float(wc)
D['kappa_crit_relerr'] = abs(kc_dl - kc_an) / kc_an
rt_c = locate(kc_an, -0.6, 0.6, 0.2, 3.0); D['root_at_crit'] = [float(rt_c[np.argmax(rt_c.real)].real), float(rt_c[np.argmax(rt_c.real)].imag)]
R['D'] = D

# ---------------- PART E: Chebyshev-Hermite box certificate ----------------
def cheb_moments(logderiv, c, a, b, Nn, nm, sym=True):
    th_ = 2 * np.pi * np.arange(Nn) / Nn
    tot = np.zeros(nm, complex)
    for j, t_ in enumerate(th_):
        if sym and j > Nn // 2: break
        tt = c + a * np.cos(t_) + 1j * b * np.sin(t_); dt = -a * np.sin(t_) + 1j * b * np.cos(t_)
        g = complex(logderiv(tt))
        u = (tt - c) / a; T = np.zeros(nm, complex); T[0] = 1; T[1] = u
        for k in range(2, nm): T[k] = 2 * u * T[k - 1] - T[k - 2]
        term = T * g * dt
        if (not sym) or j == 0 or j == Nn // 2: tot += term
        else: tot += 2j * term.imag
    return tot / (2j * np.pi) * (2 * np.pi / Nn)
def box_analysis(logderiv, c, a, b, Nn, name, nmax=24, sym=True):
    mus = cheb_moments(logderiv, c, a, b, Nn, 2 * nmax + 3, sym)
    out = {'name': name, 'c': c, 'a': a, 'b': b, 'N': Nn}
    out['mu0'] = [mus[0].real, mus[0].imag]
    n_ = int(round(mus[0].real)); out['count'] = n_
    mu = mus.real if False else mus
    def M0f(m):
        return np.array([[(mu[i + j] + mu[abs(i - j)]) / 2 for j in range(m)] for i in range(m)])
    def M1f(m):
        return np.array([[(mu[i + j + 1] + mu[abs(i + j - 1)] + mu[abs(i - j) + 1] + mu[abs(abs(i - j) - 1)]) / 4 for j in range(m)] for i in range(m)])
    M0 = M0f(n_); M1 = M1f(n_)
    out['moment_imag_max'] = float(np.max(np.abs(mus[:2 * n_ + 2].imag)))
    M0r, M1r = M0.real, M1.real
    ev0 = np.linalg.eigvalsh(M0r); out['M0_eigs'] = ev0.tolist()
    out['M0_inertia'] = [int(np.sum(ev0 > 1e-8)), int(np.sum(ev0 < -1e-8))]
    Mp = M0f(n_ + 1).real; out['M0_plus1_min_abs_eig'] = float(np.min(np.abs(np.linalg.eigvalsh(Mp))))
    out['cond_cheb'] = float(np.linalg.cond(M0r))
    # monomial Hankel of power sums in u, from Chebyshev moments via exact basis change
    from numpy.polynomial import chebyshev as Ch
    smon = []
    for k in range(2 * n_):
        cc = Ch.poly2cheb(np.eye(2 * n_)[k]); smon.append(np.sum(cc * mu[:len(cc)]))
    Hm = np.array([[smon[i + j] for j in range(n_)] for i in range(n_)]).real
    out['cond_monomial'] = float(np.linalg.cond(Hm))
    us = sla.eigvals(M1r, M0r)
    ts = np.sort_complex(c + a * us)
    out['zeros_raw'] = [[float(z.real), float(z.imag)] for z in ts]
    # Golub-Welsch in Chebyshev basis
    if out['M0_inertia'][1] == 0:
        Rc = np.linalg.cholesky(M0r).T
        J = np.linalg.solve(Rc.T, np.linalg.solve(Rc.T, M1r.T).T)
        J = (J + J.T) / 2
        evJ, UJ = np.linalg.eigh(J)
        out['gw_nodes'] = (c + a * evJ).tolist()
        out['gw_weights'] = (UJ[0] ** 2 * M0r[0, 0]).tolist()
    return out, mus
E = {}
mp.mp.dps = 25
def xi_logderiv_t(t):
    s = mp.mpf(0.5) + 1j * mp.mpc(t)
    g = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)
    return 1j * g
t0 = time.time()
b1, _ = box_analysis(xi_logderiv_t, 30.5, 20.5, 0.45, 2048, 'zeta B1 [10,51]')
b2, _ = box_analysis(xi_logderiv_t, 63.75, 12.75, 0.45, 2048, 'zeta B2 [51,76.5]')
refz = [float(mp.zetazero(k).imag) for k in range(1, 21)]
E['zeta_ref'] = refz
for bb in (b1, b2):
    zr = np.array([z[0] for z in bb['zeros_raw']]); ref = np.array([z for z in refz if bb['c'] - bb['a'] < z < bb['c'] + bb['a']])
    bb['ref'] = ref.tolist(); bb['zero_abs_err'] = float(np.max(np.abs(zr - ref))) if len(ref) == len(zr) else None
    bb['zero_imag_max'] = float(max(abs(z[1]) for z in bb['zeros_raw']))
    if 'gw_nodes' in bb:
        bb['gw_node_err'] = float(np.max(np.abs(np.array(bb['gw_nodes']) - ref)))
        bb['gw_weight_dev'] = float(np.max(np.abs(np.array(bb['gw_weights']) - 1)))
E['B1'] = b1; E['B2'] = b2; E['zeta_time_s'] = time.time() - t0
# convergence study for B1: count/moment error vs N
conv = []
for Nn in (256, 512, 1024):
    o, m_ = box_analysis(xi_logderiv_t, 30.5, 20.5, 0.45, Nn, 'conv')
    conv.append([Nn, abs(m_[0] - 10), o.get('zero_abs_err')])
E['conv'] = conv
# Davenport-Heilbronn
sq5 = mp.sqrt(5); kap = (mp.sqrt(10 - 2 * sq5) - 2) / (sq5 - 1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def dh_logderiv_t(t):
    s = mp.mpf(0.5) + 1j * mp.mpc(t)
    f0 = (1 - 1j * kap) / 2 * mp.dirichlet(s, chi) + (1 + 1j * kap) / 2 * mp.dirichlet(s, chib)
    f1 = (1 - 1j * kap) / 2 * mp.dirichlet(s, chi, derivative=1) + (1 + 1j * kap) / 2 * mp.dirichlet(s, chib, derivative=1)
    return 1j * f1 / f0
t0 = time.time()
dh, _ = box_analysis(dh_logderiv_t, 85.7, 1.0, 0.45, 768, 'DH box [84.7,86.7]', sym=False)
E['DH'] = dh; E['dh_time_s'] = time.time() - t0
dh_rt = [complex(z[0], z[1]) for z in dh['zeros_raw']]
pol = []
for z in dh_rt:
    s = mp.mpf(0.5) + 1j * mp.mpc(z)
    F = lambda ss: (1 - 1j * kap) / 2 * mp.dirichlet(ss, chi) + (1 + 1j * kap) / 2 * mp.dirichlet(ss, chib)
    sr = mp.findroot(F, s)
    pol.append([float(mp.re(sr)), float(mp.im(sr)), float(abs(sr - s))])
E['DH_polished_s'] = pol
R['E'] = E
json.dump(R, open('results.json', 'w'), indent=1, default=float)
print(json.dumps({k: (v if k != 'E' else {kk: (vv if not isinstance(vv, dict) else {q: vv[q] for q in vv if q not in ('M0_eigs',)}) for kk, vv in v.items()}) for k, v in R.items() if k != 'D'}, indent=1, default=float)[:12000])
print({k: R['D'][k] for k in R['D'] if k not in ('kappa_grid', 'rhp_counts')})
