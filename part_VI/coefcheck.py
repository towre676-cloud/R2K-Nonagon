from flint import acb, arb, acb_series, ctx
ctx.prec = 80
I = acb(0, 1); LOGPI = arb.pi().log(); K = 60
def F_point(t):  # d/dt log Xi(t) = i xi'/xi(1/2 + i t)
    s = acb(0.5) + I*t; z = acb_series([s, 1], prec=2).zeta()
    return I*(1/s + 1/(s-1) - LOGPI/2 + (s/2).digamma()/2 + z[1]/z[0])
tc = acb(7005.22, 0.00624)
sc = acb(0.5) + I*tc
S = acb_series([sc, I], prec=K+2)
for name, Lg in (('full', S.log() + (S - 1).log() - S*(LOGPI/2) + (S/2).lgamma() + S.zeta().log()),
                 ('zeta only', S.zeta().log()), ('gamma only', (S/2).lgamma()), ('logs', S.log() + (S - 1).log() - S*(LOGPI/2))):
    fc = Lg.derivative().coeffs()[:K+1]
    print(name, 'c0', fc[0].str(12), 'c1', fc[1].str(8))
fc = (S.log() + (S - 1).log() - S*(LOGPI/2) + (S/2).lgamma() + S.zeta().log()).derivative().coeffs()[:K+1]
for w in (acb(0, 0.001), acb(0, 0.00624), acb(0.00624, 0)):
    poly = sum((fc[j]*w**j for j in range(len(fc))), acb(0))
    ex = F_point(tc + w)
    print('w=', w.str(4), 'poly', poly.str(15), 'exact', ex.str(15), 'diff', float(abs(poly - ex).mid()))
# compare zeta-only part at w via point
z0 = (S.zeta()).coeffs()
w = acb(0, 0.00624)
print('Z poly', sum((z0[j]*w**j for j in range(len(z0))), acb(0)).str(15), 'Z point', (sc + I*w).zeta().str(15))
