"""Root-free rigorous tiles at large height, plus independent Hardy-Z sign-change counts."""
import time, json, sys, tmint
from tmint import *
from flint import arb_mat
tmint.MODE = 'dft'; tmint.K = 60; tmint.HFRAC = 4; tmint.STEPDIV = 4
T0, W, ntiles, b, tag = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
R2 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
KMAX = 9
def hardyZ(t):
    s = acb(0.5, t); th = (acb(0.25, t/2)).lgamma().imag - arb(t)/2*arb.pi().log()
    return (acb(0, 1)*th).exp()*s.zeta()
out = []; lo = T0
for i in range(ntiles):
    hi = lo + W; shifts = 0
    while True:
        try:
            t0 = time.time(); m, st = tmint.moments_sym(list(range(KMAX)), lo, hi, b, R2max=R2); dt = time.time() - t0; break
        except RuntimeError: hi += 0.02; shifts += 1
    mu = [x.real for x in m]; n = int(round(float(mu[0].mid())))
    rec = {'lo': lo, 'hi': hi, 'b': b, 'n': n, 'mu0': mu[0].str(10, radius=True), 'segments': sum(st), 'time': dt, 'edge_shifts': shifts}
    if n >= 1:
        M = arb_mat(n, n, [(mu[i_+j] + mu[abs(i_-j)])/2 for i_ in range(n) for j in range(n)])
        mins = [arb_mat(k, k, [M[i_, j] for i_ in range(k) for j in range(k)]).det() for k in range(1, n+1)]
        rec['all_real_simple'] = all(x > 0 for x in mins); rec['minors'] = [x.str(8, radius=True) for x in mins]
    # independent count: sign changes of Hardy Z on a fine grid
    t1 = time.time(); N = int(W*120); vals = [float(hardyZ(lo + (hi - lo)*k/N).real.mid()) for k in range(N+1)]
    rec['hardyZ_sign_changes'] = sum(1 for a_, b_ in zip(vals[:-1], vals[1:]) if (a_ > 0) != (b_ > 0)); rec['hardyZ_time'] = time.time() - t1
    out.append(rec); print(rec, flush=True); json.dump(out, open('hscan_%s.json' % tag, 'w'), indent=1)
    lo = hi
