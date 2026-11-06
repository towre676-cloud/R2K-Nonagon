"""Measurement-design test at a fixed budget of 45 real-frequency measurements (hardware-realizable, s = i*omega)."""
import json, numpy as np
exec(open('hydro_bench.py').read().split("if __name__ == '__main__':")[0])
def adaptive_design(Yfun, rng, noise, budget=45, pilot=25):
    f_p = np.geomspace(.015, 240, pilot); s_p = 2j*np.pi*f_p
    y_p = Yfun(s_p) + noise*np.abs(Yfun(s_p))*(rng.normal(size=pilot) + 1j*rng.normal(size=pilot))
    _, poles = vector_fit(s_p, y_p, 5)
    up = [p for p in poles if p.imag > 1e-6]
    extra = []
    per = (budget - pilot)//max(len(up), 1)
    for p in up:
        w0, dmp = abs(p.imag), max(abs(p.real), 1e-3*abs(p.imag))
        for j in range(per): extra.append(w0 + dmp*(-2 + 4*j/(per - 1)))   # omega within +-2 half-bandwidths
    while len(extra) < budget - pilot: extra.append(abs(up[0].imag) if up else 1.0)
    s_x = 1j*np.array(extra[:budget - pilot])
    y_x = Yfun(s_x) + noise*np.abs(Yfun(s_x))*(rng.normal(size=len(s_x)) + 1j*rng.normal(size=len(s_x)))
    return np.r_[s_p, s_x], np.r_[y_p, y_x]
res = {}
for noise in (0.01, 0.05):
    rng = np.random.default_rng(int(noise*1e4) + 11); rows = {('sweep', m): [] for m in ('vector_fit', 'r2k')}
    rows.update({('adaptive', m): [] for m in ('vector_fit', 'r2k')})
    for t in range(80):
        n = noise*np.abs(Y_all)*(rng.normal(size=76) + 1j*rng.normal(size=76)); Yn = Y_all + n
        sA, yA = s_all[idx_tr], Yn[idx_tr]
        sB, yB = adaptive_design(H, rng, noise)
        for design, (ss, yy) in (('sweep', (sA, yA)), ('adaptive', (sB, yB))):
            for m in ('vector_fit', 'r2k'):
                best = None
                for k in range(3, 9):
                    try:
                        pred, poles = METHODS[m](ss, yy, k)
                        sel = np.median(np.abs(pred(s_all[idx_sel]) - Yn[idx_sel])/np.abs(Yn[idx_sel]))
                        if np.isfinite(sel) and (best is None or sel < best[0]): best = (sel, k, pred, poles)
                    except Exception: pass
                if best is None: rows[(design, m)].append([float('inf'), float('inf'), False]); continue
                sel, k, pred, poles = best
                rc = np.abs(pred(s_all[idx_conf]) - Y_all[idx_conf])/np.abs(Y_all[idx_conf])
                rows[(design, m)].append([float(np.median(rc)), float(np.max(rc)), bool(certify_stable(poles))])
    res[str(noise)] = {f'{d}/{m}': {'median_of_median': float(np.median([r[0] for r in v])), 'median_of_max': float(np.median([r[1] for r in v])),
                                    'certified_stable': float(np.mean([r[2] for r in v]))} for (d, m), v in rows.items()}
    print(noise, json.dumps(res[str(noise)], indent=0), flush=True)
json.dump(res, open('hydro_design.json', 'w'), indent=1)
