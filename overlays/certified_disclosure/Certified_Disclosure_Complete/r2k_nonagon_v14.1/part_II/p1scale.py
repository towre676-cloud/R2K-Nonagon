import numpy as np, json
from flint import acb
from numpy.polynomial import chebyshev as Ch
out = {}
for start in (1, 1000, 6700, 100000):
    zs = [float(z.imag.mid()) for z in acb.zeta_zeros(start, 26)]
    rows = []
    for n in range(1, 21):
        lo = (zs[0] + (zs[0] - (zs[1]-zs[0]))) / 2 if True else 0
        # box: from half-gap below first to half-gap above last
        lo = zs[0] - (zs[1]-zs[0])/2; hi = zs[n-1] + (zs[n]-zs[n-1])/2
        c, a = (lo+hi)/2, (hi-lo)/2; u = (np.array(zs[:n]) - c)/a
        V = np.array([[Ch.chebval(x, np.eye(n)[k]) for k in range(n)] for x in u])
        M = V.T @ V; ev = np.linalg.eigvalsh(M)
        logdisc = (n-1)*(n-2)*np.log(2) + 2*sum(np.log(abs(u[i]-u[j])) for i in range(n) for j in range(i+1, n))
        tr = np.trace(M)
        lb = np.exp(logdisc - (n-1)*np.log(tr/(n-1))) if n > 1 else 1.0
        rows.append([n, float(ev.min()), float(lb), float(np.exp(logdisc)), float(np.sum(np.log(ev)) - logdisc)])
    out[start] = rows
    print('start', start, 'height', round(zs[0], 2))
    for r in rows[::3] + [rows[-1]]: print('  n=%2d  lmin=%.3e  bound=%.3e  disc=%.3e  logdet-logdisc=%.1e' % tuple(r))
# chebyshev-spaced reference: same n, u = chebyshev points
cheb = []
for n in range(1, 21):
    u = np.cos(np.pi*(np.arange(n)+0.5)/n); V = np.array([[np.cos(k*np.arccos(x)) for k in range(n)] for x in u]); cheb.append(float(np.linalg.eigvalsh(V.T@V).min()))
out['cheb_points'] = cheb
# fit log lmin ~ alpha + beta n + gamma n^2 for start=1000
r = np.array(out[1000]); nn = r[:, 0]; ll = np.log(r[:, 1])
coef = np.polyfit(nn[1:], ll[1:], 2); out['fit_1000_quadratic'] = coef.tolist()
coef1 = np.polyfit(nn[1:], ll[1:], 1); out['fit_1000_linear'] = coef1.tolist()
print('cheb points lmin', ['%.2f' % x for x in cheb[::4]]); print('quad fit', coef, 'lin', coef1)
json.dump(out, open('p1scale.json', 'w'), indent=1)
