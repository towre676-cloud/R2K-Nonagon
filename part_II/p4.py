"""P4: damped chains from one driving-point transfer function, 60 digits."""
import mpmath as mp, json, random
mp.mp.dps = 60
random.seed(404)
# ---------- polynomial helpers (ascending coefficient lists of mpf/mpc)
def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a
def padd(a, b): n = max(len(a), len(b)); return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def psc(a, c): return [c*x for x in a]
def pmul(a, b):
    r = [mp.mpf(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def pdiv(a, b):
    a = list(a); q = [mp.mpf(0)]*max(len(a)-len(b)+1, 1)
    for i in range(len(a)-len(b), -1, -1):
        c = a[i+len(b)-1]/b[-1]; q[i] = c
        for j in range(len(b)): a[i+j] -= c*b[j]
    return q, trim(a[:len(b)-1] or [mp.mpf(0)])
def pev(a, x): return mp.polyval(list(reversed(a)), x)
def pder(a): return [i*a[i] for i in range(1, len(a))] or [mp.mpf(0)]
def real(a): return [mp.re(x) for x in a]
# ---------- chain model
def chain_polys(m, k, c):
    """D = det Q(s), N = det Q_{2:n}(s) for Q = M s^2 + C s + K tridiagonal (dampers parallel to springs)."""
    n = len(m)
    q = [[k[i] + (k[i+1] if i+1 < n else 0), c[i] + (c[i+1] if i+1 < n else 0), m[i]] for i in range(n)]
    e = [None] + [[k[i], c[i]] for i in range(1, n)]
    dets = [None]*(n+2); dets[n] = [mp.mpf(1)]; dets[n+1] = [mp.mpf(0)]
    for i in range(n-1, -1, -1):
        dets[i] = padd(pmul(q[i], dets[i+1]), psc(pmul(pmul(e[i+1], e[i+1]), dets[i+2]), -1) if i+1 < n else [0])
    return dets[0], dets[1]
def H_eval(N, D, s): return pev(N, s)/pev(D, s)
# ---------- (b) general damping: Cauer-type peeling
def peel(N, D, n):
    m, k, c = [], [], []; log = []
    ki = pev(D, 0)/pev(N, 0)
    # c_1 = d/ds (D/N) at 0
    ci = (pev(pder(D), 0)*pev(N, 0) - pev(D, 0)*pev(pder(N), 0))/pev(N, 0)**2
    for i in range(n):
        q, rem = pdiv(D, N)          # D = (A s^2 + B s + C) N + rem
        A = q[2] if len(q) > 2 else mp.mpf(0); B = q[1]; Cc = q[0]
        m.append(A); k.append(ki); c.append(ci)
        if i == n-1:
            log.append({'level': i+1, 'last': True, 'B_minus_c': float(abs(B - ci)), 'C_minus_k': float(abs(Cc - ki)), 'rem_norm': float(max(abs(x) for x in rem))})
            break
        cn = B - ci
        W = padd(pmul(rem, pder(N)), psc(pmul(pder(rem), N), -1))
        cands = [r for r in mp.polyroots(list(reversed(trim(W))), maxsteps=500, extraprec=600) if abs(mp.im(r)) < mp.mpf(10)**-20*abs(r) and mp.re(r) < 0]
        best = None; resids = []
        for s0 in cands:
            s0 = mp.re(s0); kap = pev(rem, s0)/pev(N, s0); kn = -s0*cn
            r_ = abs(kn - (Cc - ki + kap))/abs(kn); resids.append(float(r_))
            if best is None or r_ < best[0]: best = (r_, s0, kap, kn)
        r_, s0, kap, kn = best
        mn = cn**2/kap
        e = [kn, cn]; Nn, rr = pdiv(padd(psc(N, kap), psc(rem, -1)), pmul(e, e))
        resids_sorted = sorted(resids)
        log.append({'level': i+1, 'n_wronskian_neg_real_roots': len(cands), 'best_resid': float(r_),
                    'second_best_resid': float(resids_sorted[1]) if len(resids_sorted) > 1 else None,
                    'division_remainder': float(max(abs(x) for x in rr))})
        D, N = N, trim(Nn); ki, ci = kn, cn
        # independent low-frequency check of the new level
        log[-1]['lowfreq_k_check'] = float(abs(pev(D, 0)/pev(N, 0) - ki)/ki)
    return m, k, c, log
OUT = {}
n = 5
m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(0.02, 0.4)) for _ in range(n)]
D, N = chain_polys(m, k, c)
mr, kr, cr, log = peel(N, D, n)
err = lambda a, b: float(max(abs(x - y)/abs(y) for x, y in zip(a, b)))
OUT['exact_polys'] = {'mass_err': err(mr, m), 'spring_err': err(kr, k), 'damper_err': err(cr, c), 'log': log}
print('from exact polys', OUT['exact_polys']['mass_err'], OUT['exact_polys']['spring_err'], OUT['exact_polys']['damper_err'])
for l in log: print('  ', l)
# ---------- end to end: samples -> Löwner (2n poles) -> poles/residues -> N, D -> peel
def lowner_mp(Hf, mu, nu):
    nn = len(mu); L = mp.matrix(nn, nn); Ls = mp.matrix(nn, nn)
    Hm = [Hf(z) for z in mu]; Hn = [Hf(z) for z in nu]
    for i in range(nn):
        for j in range(nn):
            L[i, j] = (Hm[i] - Hn[j])/(mu[i] - nu[j]); Ls[i, j] = (mu[i]*Hm[i] - nu[j]*Hn[j])/(mu[i] - nu[j])
    ev = mp.eig(mp.inverse(L)*Ls, left=False, right=False)
    zs = mu + nu; Hv = Hm + Hn
    Cm = mp.matrix(len(zs), nn); rhs = mp.matrix(len(zs), 1)
    for i, z in enumerate(zs):
        for j in range(nn): Cm[i, j] = 1/(z - ev[j])
        rhs[i] = Hv[i]
    res = mp.lu_solve(Cm.H*Cm, Cm.H*rhs)
    return list(ev), [res[j] for j in range(nn)]
Hf = lambda s: H_eval(N, D, s)
poles_true = sorted(mp.polyroots(list(reversed(D)), maxsteps=500, extraprec=600), key=lambda z: (float(mp.im(z)), float(mp.re(z))))
wmax = max(abs(mp.im(p)) for p in poles_true)
mu = [mp.mpc(0.05, w) for w in mp.linspace(0.1, float(wmax)*1.2, 2*n)]; nu = [mp.conj(z) for z in mu]
ev, res = lowner_mp(Hf, mu, nu)
ev_s = sorted(ev, key=lambda z: (float(mp.im(z)), float(mp.re(z))))
OUT['e2e_pole_err'] = float(max(abs(a - b)/abs(b) for a, b in zip(ev_s, poles_true)))
# rebuild N, D from poles/residues
Dm = [mp.mpf(1)]
for p in ev: Dm = pmul(Dm, [-p, 1])
Nm = [mp.mpf(0)]
for j, p in enumerate(ev):
    t = [res[j]]
    for l, p2 in enumerate(ev):
        if l != j: t = pmul(t, [-p2, 1])
    Nm = padd(Nm, t)
Dm = real(Dm); Nm = real(Nm)
OUT['e2e_dropped_leading_numerator'] = float(max(abs(x) for x in Nm[2*n-1:]))
Nm = trim(Nm[:2*n-1])
print('e2e pole err', OUT['e2e_pole_err'], 'dropped', OUT['e2e_dropped_leading_numerator'])
mr2, kr2, cr2, log2 = peel(Nm, Dm, n)
OUT['e2e'] = {'mass_err': err(mr2, m), 'spring_err': err(kr2, k), 'damper_err': err(cr2, c), 'log': log2}
print('end-to-end', OUT['e2e_pole_err'], OUT['e2e']['mass_err'], OUT['e2e']['spring_err'], OUT['e2e']['damper_err'])
for l in log2: print('  ', l)
# ---------- (a) proportional damping: C = alpha M + beta K, closed form + proportionality certificate
alpha, beta = mp.mpf('0.03'), mp.mpf('0.011')
cp = [None]*n  # proportional damping is not of parallel-damper form in general; build Q directly
def prop_H(s):
    Q = mp.matrix(n, n)
    for i in range(n):
        Kii = k[i] + (k[i+1] if i+1 < n else 0)
        Q[i, i] = m[i]*s**2 + (alpha*m[i] + beta*Kii)*s + Kii
        if i+1 < n: Q[i, i+1] = Q[i+1, i] = -(beta*s + 1)*k[i+1]
    return mp.lu_solve(Q, mp.matrix([1] + [0]*(n-1)))[0]
evp, resp = lowner_mp(prop_H, mu, nu)
up = [p for p in evp if mp.im(p) > 0]; up = sorted(up, key=lambda z: float(mp.im(z)))
lam = [abs(p)**2 for p in up]; damp = [-2*mp.re(p) for p in up]
wres = [resp[[j for j in range(2*n) if abs(evp[j] - p) < 1e-30][0]]*(p - mp.conj(p)) for p in up]
# LS fit damp = alpha + beta*lam
S1 = sum(lam); S2 = sum(x*x for x in lam); Sy = sum(damp); Sxy = sum(x*y for x, y in zip(lam, damp))
bfit = (n*Sxy - S1*Sy)/(n*S2 - S1**2); afit = (Sy - bfit*S1)/n
prop_resid = max(abs(dd - afit - bfit*l) for dd, l in zip(damp, lam))
# undamped chain from (lam, w)
def lanczos_mp(lv, wv):
    nn = len(lv); tot = sum(wv); q = [mp.sqrt(x/tot) for x in wv]; Qs = [q]; a = []; b = []; prev = [0]*nn; be = 0
    for j in range(nn):
        v = [lv[i]*Qs[j][i] - be*prev[i] for i in range(nn)]; al = mp.fsum(Qs[j][i]*v[i] for i in range(nn)); a.append(al)
        v = [v[i] - al*Qs[j][i] for i in range(nn)]
        for Qq in Qs: cc = mp.fsum(Qq[i]*v[i] for i in range(nn)); v = [v[i] - cc*Qq[i] for i in range(nn)]
        if j < nn-1: be = mp.sqrt(mp.fsum(x**2 for x in v)); b.append(be); prev = Qs[j]; Qs.append([x/be for x in v])
    return a, b
wre = [mp.re(x) for x in wres]
a_, b_ = lanczos_mp(lam, wre)
mrp = [1/sum(wre)]; krp = [1/mp.fsum(w/l for w, l in zip(wre, lam))]
for i in range(n-1):
    kn = a_[i]*mrp[i] - krp[i]; krp.append(kn); mrp.append(kn**2/(b_[i]**2*mrp[i]))
OUT['proportional'] = {'alpha_err': float(abs(afit - alpha)/alpha), 'beta_err': float(abs(bfit - beta)/beta), 'proportionality_residual': float(prop_resid),
                       'mass_err': err(mrp, m), 'spring_err': err(krp, k), 'w_imag_max': float(max(abs(mp.im(x)) for x in wres))}
print('proportional', OUT['proportional'])
# proportionality test applied to the non-proportional chain
up2 = sorted([p for p in ev if mp.im(p) > 0], key=lambda z: float(mp.im(z)))
lam2 = [abs(p)**2 for p in up2]; damp2 = [-2*mp.re(p) for p in up2]
S1 = sum(lam2); S2 = sum(x*x for x in lam2); Sy = sum(damp2); Sxy = sum(x*y for x, y in zip(lam2, damp2))
b2 = (n*Sxy - S1*Sy)/(n*S2 - S1**2); a2 = (Sy - b2*S1)/n
OUT['nonprop_proportionality_residual'] = float(max(abs(dd - a2 - b2*l) for dd, l in zip(damp2, lam2)))
print('non-proportional chain, proportionality residual', OUT['nonprop_proportionality_residual'])
json.dump(OUT, open('p4.json', 'w'), indent=1, default=float)
