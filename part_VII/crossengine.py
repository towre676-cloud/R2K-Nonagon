"""VII.1: our Taylor-model engine + independent verifier on the frontier bundle's 45 direct-Arb leaf tiles over [10,100], full strip."""
import json, glob, os, time, sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
import certgen, r2kverify
src = '/tmp/claude-0/frontier_in/NONAGON_FRONTIER/independent_certification/tiles'
os.makedirs('certs_p7', exist_ok=True); rows = []
for fn in sorted(glob.glob(src + '/tile_*.json')):
    t = json.load(open(fn)); g = t['geometry']
    lo, hi, b = Fr(g['lo']), Fr(g['hi']), Fr(g['b'])
    assert float(lo) == lo and float(hi) == hi, 'non-binary64 endpoint'
    their = {'status': t['status'], 'count': t.get('zero_count', t.get('count'))}
    out = 'certs_p7/x_%s.json' % os.path.basename(fn)[5:9]; t0 = time.time()
    try:
        c = certgen.generate('xi', float(lo), float(hi), float(b), [0, 1, 2], True, mode='series')
        json.dump(c, open(out, 'w'), indent=1)
        v = r2kverify.verify(c, full=True)
        rows.append({'tile': os.path.basename(fn), 'lo': str(lo), 'hi': str(hi), 'b': str(b), 'their_status': t['status'], 'their_raw': {k: t[k] for k in t if k in ('zero_count', 'count', 'inertia', 'all_real_simple')},
                     'ours_status': v['status'], 'ours_count': v.get('count'), 'ours_inertia': v.get('inertia'), 'ours_mu0': v.get('mu0'), 'segments': len(c['segments']), 'seconds': time.time() - t0})
    except Exception as e:
        rows.append({'tile': os.path.basename(fn), 'lo': str(lo), 'hi': str(hi), 'their_status': t['status'], 'ours_status': 'GENERATOR_FAILED', 'error': str(e)})
    print(rows[-1], flush=True)
json.dump(rows, open('p7/crossengine.json', 'w'), indent=1)
