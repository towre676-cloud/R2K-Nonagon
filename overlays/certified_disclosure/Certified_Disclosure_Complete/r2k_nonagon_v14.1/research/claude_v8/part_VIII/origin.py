"""VIII.0: the omitted origin of the number-variance quadrature. Sigma^2(L) = L - 2 int_0^L (L-x)(1-R2(x)) dx; the code integrated
from a = 1e-4. Omitted piece: 2 int_0^a (L-x)(1-R2(x)) dx = 2La - a^2 - 2 int_0^a (L-x) R2(x) dx, and 0 <= R2(x) <= c x^2 near 0."""
import json, numpy as np, mpmath as mp, sys
sys.argv = ['x', '0']
mp.mp.dps = 20
src = open('bk_multi.py').read().split("R2main = make_R2")[0]
exec(src)
R2 = make_R2(2_000_000); a = 1e-4; out = {'a': a, 'rows': []}
B = json.load(open('bk_multi_100000_1000000_10000000.json'))
for h in B['heights']:
    T = h['T_quantiles'][1]
    xs0 = np.linspace(1e-7, a, 41); r = np.array([R2(x, T) for x in xs0])
    c_est = float(np.max(np.abs(r[1:])/xs0[1:]**2))           # empirical constant in |R2| <= c x^2 on (0, a]
    for row in h['rows']:
        L = row['L']; lead = 2*L*a - a*a
        omitted = 2*np.trapezoid((L - xs0)*(1 - r), xs0) + 2*L*1e-7  # [0,1e-7] piece: (1-R2) = 1 to 1e-14
        rem_bound = 2*L*c_est*a**3/3
        corr = row['BK_height_avg'] - omitted
        out['rows'].append({'n0': h['n0'], 'L': L, 'leading_2La_minus_a2': lead, 'omitted_numeric': float(omitted), 'remainder_bound': rem_bound,
                            'BK_old': row['BK_height_avg'], 'BK_corrected': float(corr), 'zeta': row['zeta'], 'se': row['se'],
                            'residual_old': row['zeta'] - row['BK_height_avg'], 'residual_corrected': float(row['zeta'] - corr),
                            'GUE_old': row['GUE'], 'GUE_corrected': row['GUE'] - lead})
    out.setdefault('c_est', {})[str(h['n0'])] = c_est
# exact GUE check: the stored GUE column vs exact sine-kernel from 0
def s2exact(L, n=400001):
    x = np.linspace(0, L, n); return float(L - 2*np.trapezoid((L - x)*(1 - (1 - np.sinc(x)**2)), x))
out['gue_check'] = [{'L': r['L'], 'stored': r['GUE_old'], 'corrected': r['GUE_corrected'], 'exact_from_0': s2exact(r['L'])} for r in out['rows'][:6]]
for r in out['rows']: print(r['n0'], r['L'], 'lead %.8f omitted %.8f bound %.1e  resid old %+.5f new %+.5f (se %.4f)' % (r['leading_2La_minus_a2'], r['omitted_numeric'], r['remainder_bound'], r['residual_old'], r['residual_corrected'], r['se']))
for g in out['gue_check']: print('GUE', g)
json.dump(out, open('p8/origin.json', 'w'), indent=1)
