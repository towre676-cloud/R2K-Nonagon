"""v6: regenerate every Taylor-model result of Parts II-V with the corrected engine, as portable certificates."""
import json, sys, time, os, flint
import certgen
from flint import acb
group = sys.argv[1]
os.makedirs('certs', exist_ok=True)
log = []
def tile(name, lo, hi, b, nk, mode):
    fn = 'certs/%s.json' % name
    if os.path.exists(fn): print('skip', fn); return
    t0 = time.time()
    try:
        c = certgen.generate('xi', lo, hi, b, list(range(nk)), True, mode=mode)
        json.dump(c, open(fn, 'w'), indent=1)
        rec = {'name': name, 'lo': lo, 'hi': hi, 'b': b, 'segments': len(c['segments']), 'time': time.time() - t0, 'mu0': c['generator_moments'][0][:2]}
    except RuntimeError as e:
        rec = {'name': name, 'lo': lo, 'hi': hi, 'b': b, 'error': str(e)}
    print(rec, flush=True); log.append(rec); json.dump(log, open('certs/rerun_%s_log.json' % group, 'w'), indent=1)
def single_zero(n, b, mode):
    z = [float(x.imag.mid()) for x in acb.zeta_zeros(max(n - 1, 1), 3)]
    if n == 1: lo, hi = 10.0, (z[0] + z[1])/2
    else: lo, hi = (z[0] + z[1])/2, (z[1] + z[2])/2
    tile('zero_%d_b%s' % (n, str(b).replace('.', 'p')), lo, hi, b, 2, mode)
if group == 'A':
    tile('II_lehmer_b0p45', 7004.95, 7005.22, 0.45, 3, 'series')
    tile('II_lehmer_b0p6', 7004.95, 7005.22, 0.6, 3, 'series')
    tile('II_B2_box', 51.0, 76.5, 0.45, 17, 'series')
    for n in (1, 100, 1000, 3000, 6709): single_zero(n, 0.45, 'series')
    for i, lo in enumerate(range(6990, 7020, 2)): tile('III_scan7005_%02d' % i, float(lo), float(lo + 2), 0.45, 11, 'series')
    for lo, hi in ((7005.0, 7006.0), (7005.0, 7005.5), (7005.0, 7005.25), (7005.0, 7005.125)): tile('III_refine_%g_%g' % (lo, hi), lo, hi, 0.45, 3, 'series')
    for n in (100000, 1000000): single_zero(n, 0.45, 'dft'); single_zero(n, 0.6, 'dft')
elif group == 'B':
    for n in (10000000, 100000000): single_zero(n, 0.45, 'dft')
    single_zero(100000000, 0.6, 'dft')
    for i, lo in enumerate(range(600262, 600272, 2)): tile('III_scan6e5_%02d' % i, float(lo), float(lo + 2), 0.45, 11, 'dft')
    for i, lo in enumerate(range(42653544, 42653550)): tile('IV_fsscan_%02d' % i, float(lo), float(lo + 1), 0.6, 11, 'dft')
    tile('V_e10', 1e10, 1e10 + 0.6, 0.45, 5, 'dft')
elif group == 'C':
    tile('V_e12a', 1e12, 1e12 + 0.5, 0.45, 5, 'dft')
elif group == 'D':
    tile('V_e12b', 1e12 + 0.5, 1e12 + 1.0, 0.45, 5, 'dft')
