"""Part XIV.2: an efficient estimator for noisy readings, and its distance from the Cramer-Rao floor.
Readings are simulated from the source model with iid Gaussian noise sigma * G00 (G00 = port compliance).
The estimator is Levenberg-Marquardt on the raw readings, in relative parameters, started from a nominal
model in which every branch and ground is off by an independent 10% (a utility's database values).
It reports the median relative branch error, the median of the per-branch ratio |error| / CRB, and the share of
branches within 1%. Schemes as in crb14.py: 'lean' (one port) or 'portsK'."""
import sys, json, time
import numpy as np
from crb14 import load, lean_pairs, farthest_ports, rows_for
import gridlib

def model(D, w, scheme, t, pairs=None, ports=None):
    n = D['n']; el = D['el']; nb = len(el)
    E = {e: w[k] for k, e in enumerate(el)}; g = np.zeros(n); g[D['Gs']] = w[nb:]
    K = np.zeros((n, n))
    a = np.array([e[0] for e in el]); c = np.array([e[1] for e in el]); wb = w[:nb]
    np.add.at(K, (a, a), wb); np.add.at(K, (c, c), wb); np.add.at(K, (a, c), -wb); np.add.at(K, (c, a), -wb)
    K[np.arange(n), np.arange(n)] += g
    G = np.linalg.inv(K); Dc = dict(D); Dc['truth'] = w      # rows_for scales columns by the current values
    ys, Fs = [], []
    if scheme == 'lean':
        port = D['port']; g0 = G[:, port]
        ys.append([G[port, port]]); Fs.append(rows_for(Dc, g0[:, None], g0[:, None]))
        js = np.array([j for j in range(n) if j != port])
        for th in (t, 2 * t):
            f = th / (1 + th * G[js, js])
            X = g0[:, None] - G[:, js] * (f * G[port, js])[None, :]
            ys.append(X[port]); Fs.append(rows_for(Dc, X, X))
        Xs = np.empty((n, len(pairs)))
        for k, (i, j) in enumerate(pairs):
            S = [i, j]; M = np.eye(2) / t + G[np.ix_(S, S)]
            Xs[:, k] = g0 - G[:, S] @ np.linalg.solve(M, G[S, port])
        ys.append(Xs[port]); Fs.append(rows_for(Dc, Xs, Xs))
    else:
        pq = [(p, q) for ip, p in enumerate(ports) for q in ports[ip:]]
        pa = np.array([p for p, q in pq]); pb = np.array([q for p, q in pq])
        ys.append(G[pa, pb]); Fs.append(rows_for(Dc, G[:, pa], G[:, pb]))
        for j in range(n):
            f = t / (1 + t * G[j, j])
            Xa = G[:, pa] - np.outer(G[:, j], f * G[j, pa]); Xb = G[:, pb] - np.outer(G[:, j], f * G[j, pb])
            ys.append(Xb[pa, np.arange(len(pq))]); Fs.append(rows_for(Dc, Xa, Xb))
    return np.concatenate([np.ravel(y) for y in ys]), np.vstack(Fs)


def init_linear(D, y, t, pairs):
    """Data-only starting point for the lean scheme: Part VIII closed forms on the noisy readings, then the Part XII
    linear system in least squares (Part XII estimator)."""
    n, port, el = D['n'], D['port'], D['el']; js = [j for j in range(n) if j != port]; m = len(js)
    H0 = y[0]; h1 = y[1:1 + m]; h2 = y[1 + m:1 + 2 * m]; hp = y[1 + 2 * m:]
    g0 = np.zeros(n); gd = np.zeros(n); g0[port] = H0; gd[port] = H0
    d1 = H0 - h1; d2 = H0 - h2
    ig = (1 / d1 - 1 / d2) / (1 / t - 1 / (2 * t)); gg = 1 / ig
    gd[js] = (1 / d1 - ig / t) * gg; g0[js] = np.sqrt(np.abs(gg))
    Gx = {}
    for k, (i, j) in enumerate(pairs):
        dl = H0 - hp[k]; A = 1 / t + gd[i]; B = 1 / t + gd[j]
        disc = max((dl * A - g0[i] ** 2) * (dl * B - g0[j] ** 2), 0.0); Gx[(i, j)] = (g0[i] * g0[j] - np.sqrt(disc)) / dl
    def Gv(a, b):
        if a == b: return gd[a]
        if a == port: return g0[b]
        if b == port: return g0[a]
        return Gx.get((min(a, b), max(a, b)))
    col = {}; nb = len(el)
    for k, (a, c) in enumerate(el): col.setdefault(a, []).append((k, c)); col.setdefault(c, []).append((k, a))
    gi = {i: nb + k for k, i in enumerate(D['Gs'])}
    rowsA, rhs = [], []
    plan = [(port, i) for i in range(n)] + [(s, s) for s in range(n) if s != port]
    for s in D['extra']:
        seen = {s}; front = [s]
        for _ in range(D['wid'].get(s, 1)):
            front = [o for u in front for o in D['adj'][u] if o not in seen and not seen.add(o)]
        plan += [(s, i) for i in seen if i != s]
    U = nb + len(D['Gs'])
    for s, i in plan:
        r = np.zeros(U); ok = True
        gis = Gv(i, s)
        for k, o in col.get(i, []):
            gos = Gv(o, s)
            if gis is None or gos is None: ok = False; break
            r[k] = gis - gos
        if not ok: continue
        if i in gi: r[gi[i]] = gis
        rowsA.append(r); rhs.append(float(i == s))
    return np.linalg.lstsq(np.array(rowsA), np.array(rhs), rcond=None)[0]

def run(case, scheme, sigmas, seed=5):
    D = load(case); t = 1 / np.median(np.diag(D['G'])); truth = D['truth']; nb = len(D['el'])
    pairs = lean_pairs(D) if scheme == 'lean' else None
    ports = farthest_ports(D, int(scheme[5:])) if scheme.startswith('ports') else None
    crb = json.load(open('crb14_%s.json' % case))['schemes'].get(scheme)
    y0, F0 = model(D, truth, scheme, t, pairs, ports)
    Jt = F0.T @ F0; wv, V = np.linalg.eigh(Jt); ok = wv > wv[-1] * 1e-15
    var = (V[:, ok] ** 2 / wv[ok]).sum(1); G00 = D['G'][D['port'], D['port']]; kcrb = G00 * np.sqrt(var)
    rng = np.random.default_rng(seed); out = []
    for sig in sigmas:
        y = y0 + rng.standard_normal(len(y0)) * sig * G00
        t0 = time.time(); best = None
        for start in range(7):                      # start 0: data-only linear estimate (lean); then 10%-off nominal models
            if start == 0 and scheme == 'lean':
                w = init_linear(D, y, t, pairs)
                if not np.isfinite(w).all(): w = truth * (1 + 0.1 * rng.standard_normal(len(truth)))
            else:
                w = truth * (1 + 0.1 * rng.standard_normal(len(truth)))
            lam = 1e-3
            yc, F = model(D, w, scheme, t, pairs, ports); cost = np.sum((y - yc) ** 2)
            for it in range(300):
                A = F.T @ F; gvec = F.T @ (y - yc)
                step = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-30), gvec)
                wn = w * (1 + step)
                try: yn, Fn = model(D, wn, scheme, t, pairs, ports); cn = np.sum((y - yn) ** 2)
                except np.linalg.LinAlgError: cn = np.inf
                if not np.isfinite(cn): cn = np.inf
                if cn < cost:
                    rel = (cost - cn) / cost; w, yc, F, cost = wn, yn, Fn, cn; lam = max(lam / 5, 1e-12)
                    if rel < 1e-10: break
                else:
                    lam *= 10
                    if lam > 1e12: break
            if best is None or cost < best[0]: best = (cost, w, it, start)
            if cost / (sig * G00) ** 2 / len(y) < 2: break     # a fit at the noise level: accept
        cost, w, it, start = best
        err = np.abs(w - truth) / np.abs(truth); eb = err[:nb]
        ratio = eb / (kcrb[:nb] * sig)
        r = {'sigma': sig, 'median_rel_err': float(np.median(eb)), 'p90_rel_err': float(np.percentile(eb, 90)),
             'within_1pct': float(np.mean(eb < 0.01)), 'median_err_over_crb': float(np.median(ratio)),
             'crb_median_pred': float(np.median(kcrb[:nb]) * sig), 'iterations': it + 1, 'seconds': time.time() - t0,
             'chi2_per_reading': float(cost / (sig * G00) ** 2 / len(y)), 'starts_used': start + 1}
        out.append(r); print(case, scheme, r, flush=True)
    return out

if __name__ == '__main__':
    case, scheme = sys.argv[1], sys.argv[2]; sig = [float(x) for x in sys.argv[3].split(',')]
    res = run(case, scheme, sig)
    try: allr = json.load(open('gn14_%s.json' % case))
    except Exception: allr = {}
    allr[scheme] = res; json.dump(allr, open('gn14_%s.json' % case, 'w'), indent=1)
