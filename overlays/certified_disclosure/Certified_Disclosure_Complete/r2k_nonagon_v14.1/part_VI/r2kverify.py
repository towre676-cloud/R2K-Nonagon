#!/usr/bin/env python3
"""r2kverify: independent checker for r2k contour certificates. Imports nothing from the generator (tmint/certgen).

usage: r2kverify.py CERT.json [--full] [--evaluator auto|series|dft]

light mode : integrity hash, exact tiling and orientation of every side, recomputation of moments and Sylvester minors
             from the recorded per-segment balls.
full mode  : the record is treated as HINTS only (segmentation + disk radii). Every segment enclosure is re-derived
             here with an independent Taylor-model implementation; the verdict uses only the verifier's own enclosures,
             and recorded balls must overlap them.
Exit status 0 = CERTIFIED, 1 = REJECTED.  A JSON verdict is printed."""
import json, sys, hashlib
from flint import acb, arb, acb_series, arb_mat, ctx
import flint

class Reject(Exception): pass

def need(cond, msg):
    if not cond: raise Reject(msg)

def fx(h): return float.fromhex(h)

def integrity(cert):
    body = dict(cert); claimed = body.pop('sha256', None)
    need(claimed is not None, 'missing sha256')
    need(hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest() == claimed, 'sha256 mismatch: record altered')

def expected_path(cert):
    r = cert['rectangle']; lo, hi, b = fx(r['lo']), fx(r['hi']), fx(r['b'])
    need(lo < hi and 0 < b, 'degenerate rectangle')
    if cert['symmetric_half_path']:
        need(cert['function'] == 'xi', 'symmetric half-path is only valid for xi (real on the critical line)')
        return [(hi, 0.0), (hi, b), (lo, b), (lo, 0.0)], lo, hi, b
    return [(lo, -b), (hi, -b), (hi, b), (lo, b), (lo, -b)], lo, hi, b

def tiling(cert, path):
    sides = [(path[i], path[i+1]) for i in range(len(path) - 1)]
    groups = {i: [] for i in range(len(sides))}
    for s in cert['segments']:
        A = tuple(s['A']); D = tuple(s['D'])
        idx = [i for i, (P, Q) in enumerate(sides) if (P[0], P[1]) == A and (Q[0]-P[0], Q[1]-P[1]) == D]
        need(len(idx) == 1, 'segment not attached to exactly one side with correct orientation: A=%r D=%r' % (A, D))
        groups[idx[0]].append(s)
    for i, segs in groups.items():
        need(len(segs) > 0, 'side %d has no segments' % i)
        xs = [(fx(s['x0']), fx(s['x1'])) for s in segs]
        need(xs[0][0] == 0.0, 'side %d does not start at 0' % i)
        need(xs[-1][1] == 1.0, 'side %d does not end at 1' % i)
        for (a0, a1), (b0, b1) in zip(xs[:-1], xs[1:]):
            need(a1 == b0, 'side %d: gap or overlap between segments (%r vs %r)' % (i, a1, b0))
        for a0, a1 in xs: need(a0 < a1, 'side %d: empty or reversed segment' % i)
    return sides, groups

# ---------------- the functions (re-implemented here, independently of the generator) ----------------
FUNCTION = 'xi'
def _dh_kappa(): r5 = arb(5).sqrt(); return ((10 - 2*r5).sqrt() - 2)/(r5 - 1)
def F_point(s):
    """zeta(s) for xi (the gamma factor is handled exactly in the log-derivative), or the Davenport-Heilbronn function."""
    if FUNCTION == 'xi': return s.zeta()
    k = _dh_kappa(); I = acb(0, 1); w1 = (1 - I*k)/2; w2 = (1 + I*k)/2
    chi = [None, acb(1), I, -I, acb(-1)]
    tot = acb(0)
    for a_ in range(1, 5): tot += (w1*chi[a_] + w2*chi[a_].conjugate())*s.zeta(acb(a_)/5)
    return tot*acb(5)**(-s)
def F_series(S, K):
    if FUNCTION == 'xi': return S.zeta().coeffs()
    k = _dh_kappa(); I = acb(0, 1); w1 = (1 - I*k)/2; w2 = (1 + I*k)/2
    chi = [None, acb(1), I, -I, acb(-1)]
    tot = acb_series([0], prec=K+2)
    for a_ in range(1, 5): tot = tot + S.zeta(acb(a_)/5)*(w1*chi[a_] + w2*chi[a_].conjugate())
    return (tot*(S*(-acb(5).log())).exp()).coeffs()

def ball_from(rec):
    return acb(arb(rec[0], rec[1]), arb(rec[2], rec[3]))

# ---------------- independent segment enclosure ----------------
LOGPI = None
def zeta_coeffs(sc, R1, R2, K, evaluator, NDFT):
    I = acb(0, 1)
    if evaluator == 'series':
        S = acb_series([sc, I], prec=K+2); return F_series(S, K)
    # DFT from point values with rigorous aliasing bound (analytic on |w|<=R2)
    crude = F_point(acb(arb(sc.real.mid(), R2), arb(sc.imag.mid(), R2)))
    need(crude.is_finite(), 'crude bound not finite')
    M = crude.abs_upper(); N = NDFT; rho = arb(R1); q = rho/arb(R2)
    w = [(acb(2*k)/N).exp_pi_i() for k in range(N)]
    vals = [F_point(sc + I*rho*w[k]) for k in range(N)]
    out = []
    for j in range(K+2):
        acc = acb(0)
        for k in range(N): acc += vals[k]*w[(-j*k) % N]
        e = (M*q**N/(1 - q**N)/arb(R2)**j).upper()
        out.append(acc/N/rho**j + acb(arb(0, e), arb(0, e)))
    return out

def zeta_series_xi(sc, z, K):
    """d/dw log xi(sc + i w) as a power series of length K+1, given Taylor coeffs z of zeta."""
    global LOGPI
    if LOGPI is None: LOGPI = arb.pi().log()
    I = acb(0, 1); S = acb_series([sc, I], prec=K+2); Z = acb_series(z, prec=K+2)
    L = (S.log() + (S - 1).log() - S*(LOGPI/2) + (S/2).lgamma() + Z.log()) if FUNCTION == 'xi' else Z.log()
    return L.derivative().coeffs()[:K+1]

def segment_enclosure(A, D, x0, x1, R2, ks, cc, a, K, evaluator, NDFT):
    I = acb(0, 1); R1 = R2/2
    flint.ctx.cap = max(flint.ctx.cap, K + 8)
    tauc = (arb(x0) + arb(x1))/2; ht = (arb(x1) - arb(x0))/2
    tc = A + D*tauc; hw = (abs(D)*ht).upper()
    rho = hw/(arb(R1)/2); need(rho < 1, 'segment too long for its disk (rho >= 1)')
    sc = acb(0.5) + I*tc
    crude = F_point(acb(arb(sc.real.mid(), R2), arb(sc.imag.mid(), R2)))
    need(crude.is_finite(), 'crude bound not finite on disk')
    M2 = crude.abs_upper()
    z = zeta_coeffs(sc, R1, R2, K, evaluator, NDFT)
    half = arb(1)/2; tail = M2*half**(K+2)/(1 - half)
    absz = [x.abs_upper() for x in z]
    up = sum((absz[j]*arb(R1)**j for j in range(1, len(z))), arb(0))
    m1 = z[0].abs_lower() - up - tail
    need(m1 > 0, 'Rouche exclusion fails: zeta may vanish within R1 of the contour')
    M1 = absz[0] + up + tail; Mz = 4*M1/(arb(R1)*m1)
    sb = acb(arb(sc.real.mid(), R1/2), arb(sc.imag.mid(), R1/2))
    global LOGPI
    if LOGPI is None: LOGPI = arb.pi().log()
    G = I*(1/sb + 1/(sb - 1) - LOGPI/2 + (sb/2).digamma()/2) if FUNCTION == 'xi' else acb(0)
    need(G.is_finite(), 'gamma part not finite')
    Mf = Mz + G.abs_upper()
    fc = zeta_series_xi(sc, z, K)
    tb = acb(arb(tc.real.mid(), R1/2), arb(tc.imag.mid(), R1/2))
    res = []
    for k in ks:
        U = acb_series([(tc - cc)/a, acb(1)/a], prec=K+1)
        T0, T1 = acb_series([1], prec=K+1), U; Tk = T0 if k == 0 else T1
        for _ in range(2, k+1): T0, T1 = T1, 2*U*T1 - T0; Tk = T1
        P = (Tk*acb_series(fc, prec=K+1)).coeffs()
        MT = ((tb - cc)/a).chebyshev_t(k).abs_upper()
        v = acb(0)
        for j in range(0, len(P), 2): v += P[j]*D**(j+1)*2*ht**(j+1)/(j+1)
        e = (2*hw*Mf*MT*rho**(K+1)/(1 - rho)).upper()
        res.append(v + acb(arb(0, e), arb(0, e)))
    return res

def decide(cert, totals, ks):
    I = acb(0, 1)
    if cert['symmetric_half_path']: mus = [x.imag/arb.pi() for x in totals]; imag_ok = True
    else:
        m = [x/(2*acb.pi()*I) for x in totals]; mus = [x.real for x in m]
        imag_ok = all(x.imag.contains(0) for x in m)
    n_ball = mus[0]; n = int(round(float(n_ball.mid())))
    need((n_ball - n).abs_upper() < arb(1)/2, 'zero count not an unambiguous integer: %s' % n_ball.str(10, radius=True))
    need(n_ball.contains(n), 'UNSOUND: count enclosure %s excludes the integer %d (error bound failure)' % (n_ball.str(20, radius=True), n))
    verdict = {'count': n, 'mu0': n_ball.str(15, radius=True), 'imag_parts_contain_zero': imag_ok}
    if n >= 1:
        need(2*n - 2 <= max(ks), 'not enough moments recorded for a %d-zero box' % n)
        M = arb_mat(n, n, [(mus[i+j] + mus[abs(i-j)])/2 for i in range(n) for j in range(n)])
        mins = [arb_mat(k, k, [M[i, j] for i in range(k) for j in range(k)]).det() for k in range(1, n+1)]
        signs = []
        for x in mins:
            if x > 0: signs.append(1)
            elif x < 0: signs.append(-1)
            else: signs.append(0)
        need(0 not in signs, 'a Sylvester minor contains 0: inertia undecided')
        seq = [1] + signs; neg = sum(1 for i in range(n) if seq[i] != seq[i+1])
        verdict.update({'minors': [x.str(10, radius=True) for x in mins], 'inertia': [n - neg, neg],
                        'all_zeros_simple_on_line': neg == 0, 'conjugate_pairs_off_line': neg})
    return verdict

def verify(cert, full=False, evaluator='auto'):
    out = {'file_format': cert.get('format'), 'mode': 'full' if full else 'light'}
    try:
        need(cert.get('format') in ('r2k-contour-certificate/1', 'r2k-contour-certificate/2'), 'unknown format')
        integrity(cert)
        global FUNCTION
        need(cert['function'] in ('xi', 'dh'), 'unsupported function'); FUNCTION = cert['function']
        path, lo, hi, b = expected_path(cert)
        sides, groups = tiling(cert, path)
        ks = cert['moments_k']; need(ks == list(range(len(ks))), 'moments must be k = 0..K')
        prec = cert['method']['prec_bits']; ctx.prec = prec
        cc = acb(fx(cert['chebyshev_center'])); a = acb(fx(cert['chebyshev_halfwidth']))
        need(fx(cert['chebyshev_center']) == (lo + hi)/2 and fx(cert['chebyshev_halfwidth']) == (hi - lo)/2, 'Chebyshev map inconsistent with rectangle')
        totals = [acb(0) for _ in ks]; overlaps = 0; checked = 0
        K = cert['method']['K']; NDFT = cert['method']['NDFT']
        if evaluator == 'auto': evaluator = 'series' if abs(hi) < 2e5 else 'dft'
        out['evaluator'] = evaluator if full else None
        for i, segs in groups.items():
            P, Q = sides[i]; A = acb(P[0], P[1]); D = acb(Q[0] - P[0], Q[1] - P[1])
            for s in segs:
                rec = [ball_from(c) for c in s['contrib']]; need(len(rec) == len(ks), 'segment has wrong number of contributions')
                if full:
                    mine = segment_enclosure(A, D, fx(s['x0']), fx(s['x1']), fx(s['R2']), ks, cc, a, K, evaluator, NDFT)
                    for r_, m_ in zip(rec, mine): need(r_.overlaps(m_), 'recorded segment enclosure disagrees with recomputation')
                    overlaps += 1; use = mine
                else: use = rec
                totals = [t + u for t, u in zip(totals, use)]; checked += 1
        out['segments_checked'] = checked; out['segments_recomputed'] = overlaps if full else 0
        out.update(decide(cert, totals, ks)); out['status'] = 'CERTIFIED'
        r = cert['rectangle']
        out['scope'] = 'zeros of %s with %s <= Im s <= %s and %s <= Re s <= %s' % (cert['function'], repr(lo), repr(hi), repr(0.5 - b), repr(0.5 + b))
    except Reject as e:
        out['status'] = 'REJECTED'; out['reason'] = str(e)
    return out

if __name__ == '__main__':
    cert = json.load(open(sys.argv[1])); full = '--full' in sys.argv
    ev = 'auto'
    if '--evaluator' in sys.argv: ev = sys.argv[sys.argv.index('--evaluator') + 1]
    v = verify(cert, full, ev); print(json.dumps(v, indent=1)); sys.exit(0 if v['status'] == 'CERTIFIED' else 1)
