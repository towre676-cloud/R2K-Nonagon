"""Lehmer-pair gap from v6 certificates vs an Arb reference ball (no float subtraction)."""
import json, glob
from flint import acb, arb, ctx
ctx.prec = 200
z1, z2 = acb.zeta_zero(6709), acb.zeta_zero(6710); ref = z2.imag - z1.imag
out = {'reference_gap_arb': ref.str(30, radius=True), 'tiles': []}
for fn in sorted(glob.glob('certs/III_refine_*.json')) + ['certs/II_lehmer_b0p45.json', 'certs/II_lehmer_b0p6.json']:
    c = json.load(open(fn)); g = c['generator_moments']
    mu = [arb(x[0], x[1]) for x in g]   # generator_moments store mu_k directly
    lo, hi = float.fromhex(c['rectangle']['lo']), float.fromhex(c['rectangle']['hi']); a = (arb(hi) - arb(lo))/2
    det = mu[0]*(mu[2] + mu[0])/2 - mu[1]**2; gap = a*det.sqrt()
    out['tiles'].append({'file': fn, 'lo': lo, 'hi': hi, 'det': det.str(15, radius=True), 'gap': gap.str(20, radius=True),
                         'contains_reference': bool(gap.overlaps(ref)), 'abs_diff_upper': float((gap - ref).abs_upper())})
print(json.dumps(out, indent=1)); json.dump(out, open('refinev6.json', 'w'), indent=1)
