"""P1/P2 rigorous: Arb ball arithmetic contour moments on rectangles."""
from flint import acb, arb, acb_series, ctx, arb_mat
import json, time, sys
ctx.prec = 80
I = acb(0, 1); LOGPI = arb.pi().log()
def xild(s):
    z = acb_series([s, 1], prec=2).zeta()
    return 1/s + 1/(s-1) - LOGPI/2 + (s/2).digamma()/2 + z[1]/z[0]
# Davenport-Heilbronn via Hurwitz zeta: L(s,chi) = 5^{-s} sum_a chi(a) zeta(s, a/5)
kap = ((10 - 2*arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
chi = {1: acb(1), 2: I, 3: -I, 4: acb(-1)}
def dhld(s):
    S = acb_series([s, 1], prec=2)
    p5 = (-S * acb(5).log()).exp()
    L1 = acb_series([0, 0], prec=2); L2 = acb_series([0, 0], prec=2)
    for a in (1, 2, 3, 4):
        hz = S.zeta(acb(a)/5)
        L1 += chi[a] * hz; L2 += chi[a].conjugate() * hz
    F = ((1 - I*kap)/2) * L1 + ((1 + I*kap)/2) * L2
    F = F * p5
    return F[1]/F[0]
def moments(ks, ld, c, a, lo, hi, b, tol=2.0**-62):
    out = []
    corners = [acb(lo, -b), acb(hi, -b), acb(hi, b), acb(lo, b)]
    for k in ks:
        tot = acb(0)
        for i in range(4):
            A, B = corners[i], corners[(i+1) % 4]
            def g(x, an, A=A, B=B, k=k):
                t = A + (B - A)*x; s = acb(0.5) + I*t
                return ((t - c)/a).chebyshev_t(k) * I * ld(s) * (B - A)
            tot += acb.integral(g, 0, 1, rel_tol=tol)
        out.append(tot / (2*acb.pi()*I))
    return out
def s(x): return x.str(20, radius=True)
R = {}
# ---- (a) 19 rigorous single-zero tiles
res = json.load(open('results.json'))
edges = [10.0] + [(res['E']['tiles'][i]['hi']) for i in range(18)] + [76.5]
tiles = []; t0 = time.time()
for i in range(19):
    lo, hi = edges[i], edges[i+1]; c = (lo+hi)/2; a = (hi-lo)/2
    m = moments([0, 1], xild, c, a, lo, hi, 0.45)
    cnt_ok = bool((m[0].real - 1).abs_upper() < 0.5)
    loc = c + a*(m[1]/m[0]).real
    ref = acb.zeta_zero(i+1).imag
    tiles.append({'j': i+1, 'mu0': s(m[0].real), 'mu0_imag': s(m[0].imag), 'count_is_1': cnt_ok,
                  'zero': loc.str(25, radius=True), 'rad': float(loc.rad()), 'arb_zeta_zero': ref.str(25, radius=True),
                  'overlap': bool(loc.overlaps(ref))})
    print(tiles[-1], flush=True)
R['tiles'] = tiles; R['tiles_time'] = time.time() - t0
# ---- (b) rigorous box certificates via Sylvester leading minors
def box(name, ld, lo, hi, b, n):
    c = (lo+hi)/2; a = (hi-lo)/2; t0 = time.time()
    m = moments(list(range(2*n - 1)), ld, c, a, lo, hi, b)
    mu = [x.real for x in m]
    M = arb_mat(n, n, [ (mu[i+j] + mu[abs(i-j)])/2 for i in range(n) for j in range(n)])
    minors = []
    for k in range(1, n+1):
        Mk = arb_mat(k, k, [M[i, j] for i in range(k) for j in range(k)]); minors.append(Mk.det())
    out = {'name': name, 'mu0': s(mu[0]), 'max_imag_rad': max(float(x.imag.abs_upper()) for x in m),
           'minors': [s(x) for x in minors], 'minor_signs': [ (1 if x > 0 else (-1 if x < 0 else 0)) for x in minors],
           'time': time.time() - t0, 'max_entry_rad': max(float(M[i, j].rad()) for i in range(n) for j in range(n))}
    print(out, flush=True); return out, M
R['B1'], _ = box('zeta [10,51]', xild, 10.0, 51.0, 0.45, 10)
R['B2'], _ = box('zeta [51,76.5]', xild, 51.0, 76.5, 0.45, 9)
# ---- (c) Davenport-Heilbronn rigorous inertia (1,1)
dh, Mdh = box('DH [84.7,86.7]', dhld, 84.7, 86.7, 0.45, 2)
R['DH'] = dh
det = Mdh.det(); R['DH']['det'] = s(det)
off = arb(1.0) / 2 * (-det).sqrt()  # a = 1
R['DH']['offline_distance'] = off.str(25, radius=True)
# ---- (d) Lehmer pair near t = 7005: rigorous 2-zero tile
z1 = acb.zeta_zero(6709).imag; z2 = acb.zeta_zero(6710).imag
R['lehmer_ref'] = [z1.str(25), z2.str(25)]
lo, hi = 7004.95, 7005.22
lp, Mlp = box('Lehmer tile [7004.95,7005.22]', xild, lo, hi, 0.45, 2)
a = (hi - lo)/2
gap = arb(a) * Mlp.det().sqrt()
R['lehmer'] = lp; R['lehmer']['gap_from_det'] = gap.str(25, radius=True); R['lehmer']['gap_ref'] = (z2 - z1).str(25, radius=True)
R['lehmer']['overlap'] = bool(gap.overlaps(z2 - z1))
json.dump(R, open('rig_p1.json', 'w'), indent=1)
print('DONE', R['lehmer'])
