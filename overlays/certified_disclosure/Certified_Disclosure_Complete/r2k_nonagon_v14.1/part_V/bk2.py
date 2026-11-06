import numpy as np, json
from flint import acb
D = json.load(open('bk.json'))
Z = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(1000000, 20000)])
U = Z/(2*np.pi)*np.log(Z/(2*np.pi*np.e)) + 7/8; U -= U[0]
def s2(pts, L, step=0.01):
    lo = np.arange(pts[0], pts[-1]-L, step); c = np.searchsorted(pts, lo+L) - np.searchsorted(pts, lo); return float(np.var(c))
rows = []
for r in D['rows']:
    L = r['L']; full = s2(U, L); blocks = [s2(U[i*4000:(i+1)*4000], L) for i in range(5)]
    se = float(np.std(blocks, ddof=1)/np.sqrt(5))
    rows.append({'L': L, 'zeta_dense': full, 'se': se, 'BK': r['BK'], 'GUE': r['GUE'], 'z_BK': (full - r['BK'])/se if se > 0 else None,
                 'frac_of_GUE_departure_explained': (r['GUE'] - r['BK'])/(r['GUE'] - full) if r['GUE'] != full else None})
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()})
D['dense_rows'] = rows; json.dump(D, open('bk.json', 'w'), indent=1)
