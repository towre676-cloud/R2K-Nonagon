"""Taylor-model rigorous contour moments for d/dt log Xi(t), Xi(t) = xi(1/2 + i t).
Every segment: exact-point Taylor series (tight), crude-ball Cauchy tails (finite), Rouche zero-exclusion."""
from flint import acb, arb, acb_series, ctx
import flint
import time, json, sys
ctx.prec = 80
I = acb(0, 1); LOGPI = arb.pi().log()
K = 40
HFRAC = 8   # segment half-length bound: h <= R1/HFRAC
STEPDIV = 16  # segment length = R2/STEPDIV
FUNC = 'xi'   # 'xi' (Riemann xi, Taylor of zeta + exact gamma factor) or 'dh' (Davenport-Heilbronn, entire, no gamma factor)
KAP = ((10 - 2*arb(5).sqrt()).sqrt() - 2)/(arb(5).sqrt() - 1)
CHI = {1: acb(1), 2: acb(0, 1), 3: acb(0, -1), 4: acb(-1)}
def fpoint(s):
    if FUNC == 'xi': return s.zeta()
    L1 = acb(0); L2 = acb(0)
    for a_ in (1, 2, 3, 4):
        h = s.zeta(acb(a_)/5); L1 += CHI[a_]*h; L2 += CHI[a_].conjugate()*h
    return ((1 - acb(0, 1)*KAP)/2*L1 + (1 + acb(0, 1)*KAP)/2*L2)*(-s*acb(5).log()).exp()
def fseries(S):
    if FUNC == 'xi': return S.zeta()
    L1 = 0*S; L2 = 0*S
    for a_ in (1, 2, 3, 4):
        h = S.zeta(acb(a_)/5); L1 += CHI[a_]*h; L2 += CHI[a_].conjugate()*h
    return ((1 - acb(0, 1)*KAP)/2*L1 + (1 + acb(0, 1)*KAP)/2*L2)*(-S*acb(5).log()).exp()
def ball(c, r):  # complex box containing the disk |z-c|<=r
    return acb(arb(c.real.mid(), r), arb(c.imag.mid(), r))
def seg_moments(A, D, x0, x1, ks, cc, a, R2):
    """segment t = A + D*tau, tau in [x0, x1] (exact dyadics). integral of T_k(u) d/dt log Xi dt."""
    if flint.ctx.cap < K + 8: flint.ctx.cap = K + 8   # FIX v6: python-flint truncates series at ctx.cap (default 10)
    R1 = R2/2
    tauc = (arb(x0) + arb(x1))/2; ht = (arb(x1) - arb(x0))/2      # exact
    tc = A + D*tauc
    hw = (abs(D)*ht).upper()                                       # |w| <= hw on the segment
    if not (hw <= arb(R1)/HFRAC): return None
    sc = acb(0.5) + I*tc
    S = acb_series([sc, I], prec=K+2)
    zb = fpoint(ball(sc, R2))
    if not zb.is_finite(): return None
    M2 = zb.abs_upper(); q = arb(1)/2
    if MODE == 'series':
        Z = fseries(S); z = Z.coeffs()
    else:
        z = dft_coeffs(sc, arb(R1), M2, arb(R2)); Z = acb_series(z, prec=K+2)
    tail = M2*q**(K+2)/(1 - q)
    absz = [x.abs_upper() for x in z]
    up = sum((absz[j]*arb(R1)**j for j in range(1, len(z))), arb(0))
    m1 = z[0].abs_lower() - up - tail
    if not (m1 > 0): return None
    M1 = absz[0] + up + tail
    Mz = 4*M1/(arb(R1)*m1)   # FIX v6: Cauchy |Z'(w)| <= M1*R1/(R1-|w|)^2 = 4*M1/R1 on |w| <= R1/2
    sb = ball(sc, R1/2)
    G = I*(1/sb + 1/(sb - 1) - LOGPI/2 + (sb/2).digamma()/2) if FUNC == 'xi' else acb(0)
    if not G.is_finite(): return None
    Mf = Mz + G.abs_upper()
    Lg = (S.log() + (S - 1).log() - S*(LOGPI/2) + (S/2).lgamma() + Z.log()) if FUNC == 'xi' else Z.log()
    fc = Lg.derivative().coeffs()[:K+1]
    out = []; tb = ball(tc, R1/2)
    for k in ks:
        U = acb_series([(tc - cc)/a, acb(1)/a], prec=K+1)
        T0, T1 = acb_series([1], prec=K+1), U
        Tk = T0 if k == 0 else T1
        for _ in range(2, k+1): T0, T1 = T1, 2*U*T1 - T0; Tk = T1
        P = (Tk*acb_series(fc, prec=K+1)).coeffs()
        MT = ((tb - cc)/a).chebyshev_t(k).abs_upper()
        val = acb(0)
        for j in range(0, len(P), 2):          # w = D*(tau - tauc), dt = D dtau
            val += P[j] * D**(j+1) * 2*ht**(j+1)/(j+1)
        rho = hw/(arb(R1)/2)
        err = 2*hw * Mf*MT * rho**(K+1)/(1 - rho)
        r_ = err.upper()
        out.append(val + acb(arb(0, r_), arb(0, r_)))
    return out
MODE = 'series'; NDFT = 64; R2MIN = 2e-6
def dft_coeffs(sc, rho, M, R):
    # Z(w) = zeta(sc + i w); coefficients a_j, j <= K+1, from NDFT points on |w| = rho.
    # aliasing: |a~_j - a_j| <= M R^-j q^N / (1 - q^N), q = rho/R  (analytic Z on |w| <= R)
    N = NDFT; vals = []; roots = [(acb(2*k)/N).exp_pi_i() for k in range(N)]
    for k in range(N): vals.append(fpoint(sc + I*rho*roots[k]))
    q = rho/R; alias = M * q**N/(1 - q**N)
    out = []
    for j in range(K+2):
        acc = acb(0)
        for k in range(N): acc += vals[k]*roots[(-j*k) % N]
        aj = acc/N/rho**j
        e = (alias/R**j).upper()
        out.append(aj + acb(arb(0, e), arb(0, e)))
    return out
def side(A, B, ks, cc, a, R2max=0.05, stats=None):
    D = B - A; L = float(abs(complex(D.real.mid(), D.imag.mid()) if False else abs(complex(float(D.real.mid()), float(D.imag.mid())))))
    x0 = 0.0; tot = [acb(0) for _ in ks]; nseg = 0; R2 = R2max
    while x0 < 1.0:
        for trial in range(60):
            step = (R2/STEPDIV)/L*0.999
            x1 = min(1.0, x0 + step)
            r = seg_moments(A, D, x0, x1, ks, cc, a, R2)
            if r is not None: break
            R2 /= 2
            if R2 < R2MIN: raise RuntimeError('zero near contour at %r' % (x0,))
        else: raise RuntimeError('cannot certify segment at %r' % (x0,))
        tot = [p + q_ for p, q_ in zip(tot, r)]; x0 = x1; nseg += 1
        R2 = min(R2*1.5, R2max)
    if stats is not None: stats.append(nseg)
    return tot
def moments(ks, lo, hi, b, R2max=0.05):
    cc = (lo + hi)/2; a = (hi - lo)/2
    corners = [acb(lo, -b), acb(hi, -b), acb(hi, b), acb(lo, b)]
    tot = [acb(0) for _ in ks]; stats = []
    for i in range(4):
        r = side(corners[i], corners[(i+1) % 4], ks, acb(cc), acb(a), R2max, stats)
        tot = [x + y for x, y in zip(tot, r)]
    return [x/(2*acb.pi()*I) for x in tot], stats
if __name__ == '__main__':
    lo, hi = float(sys.argv[1]), float(sys.argv[2]); nk = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    t0 = time.time(); m, st = moments(list(range(nk)), lo, hi, 0.45)
    print([x.str(18, radius=True) for x in m], 'segments', st, 'time %.1f' % (time.time() - t0))
def moments_sym(ks, lo, hi, b, R2max=0.2):
    """conjugation-symmetric rectangle: mu_k = Im(integral over upper half-path)/pi, valid since Xi(conj t) = conj Xi(t)."""
    cc = (lo + hi)/2; a = (hi - lo)/2
    path = [acb(hi, 0), acb(hi, b), acb(lo, b), acb(lo, 0)]
    tot = [acb(0) for _ in ks]; stats = []
    for i in range(3):
        r = side(path[i], path[i+1], ks, acb(cc), acb(a), R2max, stats)
        tot = [x + y for x, y in zip(tot, r)]
    return [acb(x.imag/arb.pi()) for x in tot], stats
