"""X.0: rigorous closure of the origin remainder (v9 critique of Prop. VIII.0.1). R2 is analytic in eps near 0 (the r^-2 poles cancel).
On the circle |eps| = rho we evaluate R2 in Arb (the same truncated Euler products, primes <= 2e6), recover Taylor coefficients by a
DFT with a rigorous aliasing bound from a ball cover of the circle |eps| = R, and enclose int_0^a (L - x) R2(x) dx."""
import json, time, numpy as np, flint
from flint import acb, arb, acb_series, ctx
ctx.prec = 96; flint.ctx.cap = max(flint.ctx.cap, 8)
def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]
P = sieve(2_000_000); LP = [arb(p).log() for p in P]; IP = [arb(1)/p for p in P]
def AB(eta):
    A = acb(1); B = acb(0)
    for p, lp, ip in zip(P, LP, IP):
        x = (-(1 + eta)*lp).exp()                   # p^{-(1+eta)}
        A *= (1 - x)*(1 - 2*ip + x)/((1 - ip)**2)
        B += (lp*x/(1 - x))**2                        # (log p / (p^{1+eta} - 1))^2
    return A, B
def f(r, T):
    s = acb(1) + acb(0, 1)*r; z = acb_series([s, 1], prec=3).zeta().coeffs(); z0, z1, z2 = z[0], z[1], z[2]
    dlog = 2*z2/z0 - (z1/z0)**2
    zm = (acb(1) - acb(0, 1)*r).zeta(); A, B = AB(acb(0, 1)*r)
    return dlog + (-(acb(0, 1)*r)*(T/(2*arb.pi())).log()).exp()*z0*zm*A - B
def R2(eps, T):
    Lg = (T/(2*arb.pi())).log(); r = 2*arb.pi()*eps/Lg
    return 1 + (f(r, T) + f(-r, T))/Lg**2
out = []; a = arb(1)/10000; N = 64
for T in (74920.0, 600269.0, 4992381.0):
    t0 = time.time(); T = arb(T); rho = arb(1)/40; Rr = arb(1)/20
    w = [(acb(2*k)/N).exp_pi_i() for k in range(N)]
    vals = [R2(rho*w[k], T) for k in range(N)]
    # crude bound of |R2| on |eps| <= R: cover the circle |eps| = R by 64 boxes (maximum principle)
    M = arb(0)
    for k in range(64):
        c0 = (acb(2*k)/64).exp_pi_i()*Rr; box = acb(arb(c0.real.mid(), (Rr*arb.pi()/32).upper()), arb(c0.imag.mid(), (Rr*arb.pi()/32).upper()))
        M = M.max(R2(box, T).abs_upper())
    q = rho/Rr; alias = M*q**N/(1 - q**N)
    coef = []
    for j in range(N):
        acc = acb(0)
        for k in range(N): acc += vals[k]*w[(-j*k) % N]
        e = (alias/Rr**j).upper(); coef.append(acc/N/rho**j + acb(arb(0, e), arb(0, e)))
    # enclosure of R2 on [0, a]: |R2(x) - sum_{j<N} c_j x^j| <= M (x/R)^N / (1 - x/R)
    tail = M*(a/Rr)**N/(1 - a/Rr)
    res = {'T': float(T.mid()), 'M_on_R_circle': float(M.upper()), 'c0': coef[0].str(8, radius=True), 'c1': coef[1].str(8, radius=True), 'c2': coef[2].str(8, radius=True),
           'c3': coef[3].str(8, radius=True), 'seconds': time.time() - t0}
    for L in (8, 16):
        I = acb(0)
        for j in range(N): I += coef[j]*(L*a**(j+1)/(j+1) - a**(j+2)/(j+2))     # int_0^a (L - x) x^j dx
        I += acb(arb(0, (tail*L*a).upper()), 0)
        omitted = 2*L*a - a*a - 2*I                                               # 2 int_0^a (L - x)(1 - R2) dx
        res['L%d' % L] = {'omitted': omitted.real.str(15, radius=True), 'remainder_2int_LR2': (2*I).real.str(6, radius=True)}
    out.append(res); print(res, flush=True)
json.dump(out, open('p10/origin_rigorous.json', 'w'), indent=1)
