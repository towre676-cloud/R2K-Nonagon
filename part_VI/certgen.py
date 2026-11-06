"""certgen: run the Taylor-model contour engine (tmint) and emit a portable proof record.
The record is a HINT object: exact segmentation, disk radii, method parameters, and the generator's own
per-segment enclosures. The independent verifier (r2kverify.py) recomputes every segment itself."""
import json, hashlib, time, sys, platform
import tmint
from tmint import acb, arb
import flint

def _ball(x):  # acb -> [re_mid, re_rad, im_mid, im_rad] as exact-ish decimal strings
    return [x.real.mid().str(40, radius=False), float(x.real.rad()), x.imag.mid().str(40, radius=False), float(x.imag.rad())]

def generate(function, lo, hi, b, ks, symmetric, mode='dft', K=60, HFRAC=4, STEPDIV=4, R2max=0.2, prec=80, NDFT=64):
    assert function == 'xi' or not symmetric
    tmint.FUNC = function
    tmint.MODE = mode; tmint.K = K; tmint.HFRAC = HFRAC; tmint.STEPDIV = STEPDIV; tmint.NDFT = NDFT
    flint.ctx.prec = prec; tmint.ctx.prec = prec
    rec = []; orig = tmint.seg_moments
    def wrapped(A, D, x0, x1, ks_, cc, a, R2):
        r = orig(A, D, x0, x1, ks_, cc, a, R2)
        if r is not None:
            rec.append({'A': [float(A.real.mid()), float(A.imag.mid())], 'D': [float(D.real.mid()), float(D.imag.mid())],
                        'x0': float(x0).hex(), 'x1': float(x1).hex(), 'R2': float(R2).hex(), 'contrib': [_ball(v) for v in r]})
        return r
    tmint.seg_moments = wrapped
    t0 = time.time()
    try:
        if function == 'xi' and symmetric:
            m, st = tmint.moments_sym(ks, lo, hi, b, R2max=R2max)
            path = [[hi, 0.0], [hi, b], [lo, b], [lo, 0.0]]; closed = False
        else:
            m, st = tmint.moments(ks, lo, hi, b, R2max=R2max)
            path = [[lo, -b], [hi, -b], [hi, b], [lo, b], [lo, -b]]; closed = True
    finally:
        tmint.seg_moments = orig
    dt = time.time() - t0
    src = open(tmint.__file__, 'rb').read()
    cert = {
        'format': 'r2k-contour-certificate/2',
        'function': function, 'variable': 't, s = 1/2 + i t',
        'rectangle': {'lo': float(lo).hex(), 'hi': float(hi).hex(), 'b': float(b).hex()},
        'symmetric_half_path': symmetric, 'path': [[float(p[0]).hex(), float(p[1]).hex()] for p in path],
        'moments_k': ks, 'chebyshev_center': ((lo + hi)/2).hex(), 'chebyshev_halfwidth': ((hi - lo)/2).hex(),
        'method': {'mode': mode, 'K': K, 'HFRAC': HFRAC, 'NDFT': NDFT, 'prec_bits': prec},
        'segments': rec,
        'generator_moments': [_ball(x) for x in m],
        'generator_segments_per_side': st, 'generator_time_s': dt,
        'provenance': {'tmint_sha256': hashlib.sha256(src).hexdigest(), 'python_flint': flint.__version__, 'python': platform.python_version()},
        'function_definition': ('xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s)' if function == 'xi' else
                            'Davenport-Heilbronn f(s) = 5^(-s) sum_a [((1-i k)/2) chi(a) + ((1+i k)/2) conj chi(a)] zeta(s, a/5), chi mod 5 with chi(2) = i, k = (sqrt(10-2 sqrt5)-2)/(sqrt5-1)'),
    'scope': ('Claims follow from the moments: mu_0 counts zeros of xi(1/2+it) in the rectangle; if all Sylvester minors of the '
                  'Chebyshev-Hermite matrix are positive, every zero in the rectangle is simple and on the critical line. '
                  'The rectangle covers Re s in [1/2-b, 1/2+b] and Im s in [lo, hi].')}
    body = json.dumps(cert, sort_keys=True).encode()
    cert['sha256'] = hashlib.sha256(body).hexdigest()
    return cert

if __name__ == '__main__':
    lo, hi, b, nk, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    mode = sys.argv[6] if len(sys.argv) > 6 else 'dft'
    c = generate('xi', lo, hi, b, list(range(nk)), True, mode=mode)
    json.dump(c, open(out, 'w'), indent=1)
    print(out, 'segments', len(c['segments']), 'time %.1f' % c['generator_time_s'], 'mu0', c['generator_moments'][0][:2])
