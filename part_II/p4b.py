import mpmath as mp, json, random
exec(open('p4.py').read().split('OUT = {}')[0])
mp.mp.dps = 50; random.seed(4040)
stats = {}
for regime, (clo, chi) in (('light', (0.02, 0.4)), ('heavy', (0.5, 3.0)), ('very_heavy', (3.0, 15.0))):
    multi = 0; levels = 0; maxerr = 0; gaps = []; fails = 0
    for t in range(70):
        n = random.choice([3, 4, 5, 6])
        m = [mp.mpf(random.uniform(0.5, 2)) for _ in range(n)]; k = [mp.mpf(random.uniform(1, 5)) for _ in range(n)]; c = [mp.mpf(random.uniform(clo, chi)) for _ in range(n)]
        D, N = chain_polys(m, k, c)
        try:
            mr, kr, cr, log = peel(N, D, n)
            e = max(float(max(abs(x-y)/abs(y) for x, y in zip(a, b))) for a, b in ((mr, m), (kr, k), (cr, c))); maxerr = max(maxerr, e)
            for l in log:
                if 'n_wronskian_neg_real_roots' in l:
                    levels += 1
                    if l['n_wronskian_neg_real_roots'] > 1: multi += 1; gaps.append(l['second_best_resid'])
        except Exception as ex: fails += 1
    stats[regime] = {'levels': levels, 'levels_with_multiple_candidates': multi, 'min_second_best_resid': (min(gaps) if gaps else None), 'max_param_err': maxerr, 'failures': fails}
    print(regime, stats[regime], flush=True)
json.dump(stats, open('p4b.json', 'w'), indent=1)
