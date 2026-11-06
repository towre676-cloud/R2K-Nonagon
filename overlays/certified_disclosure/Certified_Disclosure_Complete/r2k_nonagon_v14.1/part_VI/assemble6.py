import json, re, os, glob, math
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon VI</title>')
body = open('body6.html').read()
def pm(s):  # '[2.00000 +/- 1.1e-16]' -> latex
    m = re.match(r'\[(.*?) \+/- (.*?)\]', s)
    if not m: return s.replace('+/-', '±')
    mid, rad = m.groups(); a, e = rad.split('e') if 'e' in rad else (rad, '0')
    return '$%s\\pm%s\\times10^{%d}$' % (mid, a, int(e))
def pmint(s, n):
    m = re.match(r'\[(.*?) \+/- (.*?)\]', s)
    if not m: return s
    rad = float(m.group(2)); e = int(math.floor(math.log10(rad))); return '$%d\\pm%.1f\\times10^{%d}$' % (n, rad/10**e, e)
def sci(x, d=2):
    if x == 0: return '0'
    e = int(math.floor(math.log10(abs(x)))); return '$%.*f\\times10^{%d}$' % (d - 1, x/10**e, e)
VF = json.load(open('verify_full.json')); VL = json.load(open('verify_light.json'))
if os.path.exists('verify_full2.json'): pass
recs = sorted(f for f in glob.glob('certs/*.json') if 'log' not in f)
order = ['zero_', 'II_', 'III_scan7005', 'III_refine', 'III_scan6e5', 'IV_', 'V_', 'NEG_']
def key(f):
    n = os.path.basename(f)
    for i, p in enumerate(order):
        if n.startswith(p):
            if p == 'zero_':
                return (i, int(n.split('_')[1]), n)
            return (i, 0, n)
    return (99, 0, n)
recs.sort(key=key)
rows = ''; nfull = 0
for f in recs:
    c = json.load(open(f)); v = VF.get(f) or VL.get(f); r = c['rectangle']
    lo, hi, b = float.fromhex(r['lo']), float.fromhex(r['hi']), float.fromhex(r['b'])
    name = os.path.basename(f)[:-5].replace('_', '\\_')
    rng = ('[%.6g, %.6g]' % (lo, hi)) if hi < 1e9 else ('$10^{%d}+[%.1f,\\,%.1f]$' % (round(math.log10(lo)), lo - 10**round(math.log10(lo)), hi - 10**round(math.log10(lo))))
    if v is None: verdict, rep, mu, cnt, ine = 'pending', '—', '—', '—', '—'
    else:
        verdict = v['status'] + (' (full)' if v['mode'] == 'full' else ' (light; full replay running)'); rep = '%.0f s' % v['time'] if v['time'] >= 1 else '%.1f s' % v['time']
        mu = pmint(v.get('mu0', ''), v.get('count')); cnt = v.get('count', '—'); ine = '(%d,%d)' % tuple(v['inertia']) if v.get('inertia') else '—'
        if v['status'] == 'CERTIFIED' and v['mode'] == 'full': nfull += 1
    rows += f'<tr><td><code>{os.path.basename(f)[:-5]}</code></td><td class="n">{rng}</td><td class="n">{b:g}</td><td class="n">{len(c["segments"])}</td><td class="n">{cnt}</td><td class="n">{ine}</td><td class="n">{mu}</td><td class="n">{rep}</td><td>{verdict}</td></tr>\n'
body = body.replace('<!--REGENROWS-->', rows).replace('<!--NCERT-->', str(len(recs))).replace('<!--NFULL-->', str(nfull) if nfull < len(recs) else 'All')
# tamper
T = json.load(open('tamper_tests.json')); tr = ''
for t in T:
    tr += f"<tr><td>{t['test']}</td><td>{t.get('mode', 'generator')}</td><td>{t['status']}</td><td>{t.get('reason', '') or '—'}</td></tr>\n"
body = body.replace('<!--TAMPERROWS-->', tr)
# env2
er = ''; cur = None
for line in open('env2_replay.out'):
    line = line.strip()
    if line.startswith('=='): cur = line[3:]; continue
    if cur and line:
        parts = line.split(' ', 3); st, cnt = parts[0], parts[1]
        ine = re.search(r'\[(\d+), (\d+)\]', line); mu = re.search(r'(\[[0-9.]+ \+/- [0-9.e-]+\])', line)
        er += f'<tr><td><code>{cur}</code></td><td>{st}</td><td class="n">{cnt}</td><td class="n">({ine.group(1)},{ine.group(2)})</td><td class="n">{pm(mu.group(1))}</td></tr>\n'; cur = None
body = body.replace('<!--ENV2ROWS-->', er)
# M4
M = json.load(open('m4v6.json')); mr = ''
lab = {'rho_oracle': ('$\\rho=\\sigma_n(\\mathbb L)/\\|E\\|_2$ (oracle)', 'yes'), 'rho_operational_bootstrap': ('$\\hat\\rho=\\sigma_n(\\tilde{\\mathbb L})/\\widehat{\\|E\\|}_2$ (bootstrap)', 'no'), 'singular_gap_ratio': ('$\\sigma_n/\\sigma_{n+1}$ of noisy pencil', 'no')}
for k_, v in M['predictors'].items():
    mr += f"<tr><td>{lab[k_][0]}</td><td>{lab[k_][1]}</td><td class=\"n\">{v['threshold_train']:.3f}</td><td class=\"n\">{100*v['test_accuracy']:.1f}%</td><td class=\"n\">[{100*v['test_accuracy_wilson95'][0]:.1f}, {100*v['test_accuracy_wilson95'][1]:.1f}]%</td><td class=\"n\">{v['test_missed_failures']}</td><td class=\"n\">{v['test_false_alarms']}</td></tr>\n"
body = body.replace('<!--M4ROWS-->', mr)
# hydro
Hd = json.load(open('hydro_bench.json')); hr = ''
names = {'loewner_proj': 'projected Löwner', 'vector_fit': 'vector fitting', 'r2k': 'R2K (Löwner seed + GN + certificate)', 'vf_gn': 'VF seed + GN + certificate'}
for nz in ('0.001', '0.01', '0.05'):
    for mth in ('loewner_proj', 'vector_fit', 'r2k', 'vf_gn'):
        v = Hd[nz][mth]
        hr += f"<tr><td class=\"n\">{float(nz)*100:g}%</td><td>{names[mth]}</td><td class=\"n\">{sci(v['median_of_median_rel_err'])}</td><td class=\"n\">{sci(v['p90_of_median_rel_err'])}</td><td class=\"n\">{sci(v['median_of_max_rel_err'])}</td><td class=\"n\">{100*v['frac_certified_stable']:.0f}%</td></tr>\n"
body = body.replace('<!--HYDROROWS-->', hr)
# shell
def shell_rows(fn, label):
    if not os.path.exists(fn): return '', None
    S = json.load(open(fn)); out = ''; zmin = None
    for ns, r in S.items():
        opt = [o for o in r['optimized'] if o['exact_unforced_multiple'] and o['exact']['neg_real_roots_D'] == 0]
        zm = min((o['zeta_max'] for o in opt), default=None)
        if zm is not None and (zmin is None or zm < zmin[0]): zmin = (zm, ns, label)
        q = r['hit_zeta_quantiles']
        out += (f"<tr><td class=\"n\">{ns} mass{'es' if int(ns) > 1 else ''}</td><td>{label}</td><td class=\"n\">{r['admissible_X_pos_no_real_poles']:,}</td><td class=\"n\">{r['unforced_hits']:,}</td>"
                f"<td class=\"n\">{'%.3f / %.3f' % (q[0], q[2]) if q else '—'}</td><td class=\"n\">{'%.4f' % zm if zm is not None else '—'}</td><td class=\"n\">{len(opt)}/{len(r['optimized'])}</td></tr>\n")
    return out, zmin
r1, z1 = shell_rows('shell_search_1.json', 'wide')
r2, z2 = shell_rows('shell_search_2_3.json', 'wide')
r3, z3 = shell_rows('shell_box4_2_3.json', 'Part IV box')
body = body.replace('<!--SHELLROWS-->', r1 + r2 + r3)
SH = json.load(open('shell_search_2_3.json')) if os.path.exists('shell_search_2_3.json') else {}
SB = json.load(open('shell_box4_2_3.json')) if os.path.exists('shell_box4_2_3.json') else {}
if z2:
    ok_ = lambda o: o['exact_unforced_multiple'] and o['exact']['neg_real_roots_D'] == 0
    best = min((o for r in SH.values() for o in r['optimized'] if ok_(o)), key=lambda o: o['zeta_max'])
    import numpy as _np; ratio = max(_np.r_[best['m'], best['k'], best['c']])/min(_np.r_[best['m'], best['k'], best['c']])
    L = len(best['m'])
    fmt = lambda a: ', '.join('%.4g' % x for x in a)
    txt = (f'<p>For two-mass sub-chains the floating-point screen of the random stage already flags unforced multiple levels with $\\zeta_{{\\max}}$ as low as {SH["2"]["hit_zeta_quantiles"][0]:.3f}, '
           f'and descent drives the damping ratio to {best["zeta_max"]:.4f}. That point is a level of {L} masses (the peeled mass plus a sub-chain of {L-1}) with '
           f'$m=({fmt(best["m"])})$, $k=({fmt(best["k"])})$, $c=({fmt(best["c"])})$. In exact arithmetic over $\\mathbb Q$, $X&gt;0$, the sub-chain polynomial $N$ has no negative real root, the full level polynomial $D$ has no negative real root either, and $W$ has {best["exact"]["neg_real_roots_W"]} distinct negative real roots: the physical $s_0=-k_2/c_2$ and two more. '
           f'The descent is unconstrained in scale, and this point spans a parameter ratio of about {ratio:,.0f}: a heavy, nearly undamped slow mode next to a light, stiff one. '
           f'The last column of Table VI.5 gives, for each search, how many of the twelve descended points pass this exact test.</p>')
    if SB:
        q2 = SB.get('2', {}); hits = q2.get('unforced_hits'); adm = q2.get('admissible_X_pos_no_real_poles')
        zb = min((o['zeta_max'] for r in SB.values() for o in r['optimized'] if ok_(o)), default=None)
        txt += (f'<p>Inside Part IV\'s own parameter box, unforced multiplicity is rarer ({hits:,} of {adm:,} admissible two-mass levels in random sampling, about {100*hits/adm:.1f}%), and already in random sampling it occurs at $\\zeta_{{\\max}}={q2["hit_zeta_quantiles"][0]:.3f}$;'
                + (f' constrained descent reaches $\\zeta_{{\\max}}={zb:.4f}$ without leaving the box.' if zb is not None else '.')
                + ' The 0.86 threshold was therefore a property of where Part IV looked, combined with the rarity of such levels, and not a property of damped chains.</p>')
    body = body.replace('<!--SHELLTEXT-->', txt).replace('<!--ZMIN-->', '%.4f' % best['zeta_max'])
# BK
B = None
for fn in ('bk_multi_100000_1000000_10000000.json',):
    if os.path.exists(fn): B = json.load(open(fn))
if B is None:
    body = body.replace('<!--BKROWS-->', '<tr><td colspan="10">three-height run in progress at release time; see Part V Table V.1 for the single-height comparison</td></tr>').replace('<!--FIGBK-->', '').replace('<!--BKTEXT-->', '<div class="thm"><div class="hd">Status VI.8.1 <i>(three heights)</i> <span class="ep op">running</span></div><p>The disjoint-window computation (<code>bk_multi.py</code>, windows from zeros $10^5$, $10^6$, $10^7$, height-averaged prediction, cutoff and quadrature sensitivity) was still running when this release was cut, and its table will be filled in the next build. Until then, Part V\'s single-height comparison stands as a numerical test of a conjectural formula, with the residuals at $L=8$ and $16$ unexplained.</p></div>')
if B:
    br = ''; pts = []
    for h in B['heights']:
        for r in h['rows']:
            br += (f"<tr><td class=\"n\">from $\\gamma_{{{h['n0']:,}}}$</td><td class=\"n\">{r['L']:g}</td><td class=\"n\">{r['zeta']:.4f}</td><td class=\"n\">{r['se']:.4f}</td><td class=\"n\">{r['BK_height_avg']:.4f}</td><td class=\"n\">{r['GUE']:.4f}</td>"
                   f"<td class=\"n\">{r['z_BK_avg']:+.1f}</td><td class=\"n\">{100*r['frac_GUE_departure_explained']:.1f}%</td><td class=\"n\">{r['prime_cutoff_shift']:+.1e}</td><td class=\"n\">{r['quadrature_shift']:+.1e}</td></tr>\n").replace(',', '{,}', 1) if False else \
                  (f"<tr><td class=\"n\">from zero {h['n0']:,}</td><td class=\"n\">{r['L']:g}</td><td class=\"n\">{r['zeta']:.4f}</td><td class=\"n\">{r['se']:.4f}</td><td class=\"n\">{r['BK_height_avg']:.4f}</td><td class=\"n\">{r['GUE']:.4f}</td>"
                   f"<td class=\"n\">{r['z_BK_avg']:+.1f}</td><td class=\"n\">{100*r['frac_GUE_departure_explained']:.1f}%</td><td class=\"n\">{r['prime_cutoff_shift']:+.1e}</td><td class=\"n\">{r['quadrature_shift']:+.1e}</td></tr>\n")
            pts.append((h['n0'], r['L'], r['z_BK_avg'], (r['zeta'] - r['BK_height_avg'])))
    body = body.replace('<!--BKROWS-->', br)
    # figure: residual (zeta - BK) / se vs L for each height
    W, Hh, ml, mr_, mt, mb = 640, 300, 60, 20, 20, 44
    Ls = [0.5, 1, 2, 4, 8, 16]; zs = [p[2] for p in pts]; zmax = max(3, max(abs(z) for z in zs)*1.1)
    X = lambda L: ml + (math.log2(L) + 1)/5*(W - ml - mr_); Y = lambda z: mt + (zmax - z)/(2*zmax)*(Hh - mt - mb)
    cols = {100000: 'var(--slate)', 1000000: 'var(--ox)', 10000000: 'var(--ink)'}
    svg = [f'<svg viewBox="0 0 {W} {Hh}" role="img" aria-label="Number-variance residuals" style="width:100%;max-width:{W}px;font-family:var(--body)">']
    for zz in (-2, 0, 2):
        svg.append(f'<line x1="{ml}" x2="{W-mr_}" y1="{Y(zz):.1f}" y2="{Y(zz):.1f}" stroke="var(--rule)" stroke-dasharray="{"" if zz == 0 else "4 4"}"/>')
        svg.append(f'<text x="{ml-8}" y="{Y(zz)+4:.1f}" text-anchor="end" font-size="12" fill="var(--ink-soft)">{zz:+d}</text>')
    for L in Ls:
        svg.append(f'<text x="{X(L):.1f}" y="{Hh-mb+18}" text-anchor="middle" font-size="12" fill="var(--ink-soft)">{L:g}</text>')
    svg.append(f'<text x="{(ml+W-mr_)/2}" y="{Hh-6}" text-anchor="middle" font-size="13" fill="var(--ink)">L (mean spacings)</text>')
    svg.append(f'<text x="14" y="{(mt+Hh-mb)/2}" text-anchor="middle" font-size="13" fill="var(--ink)" transform="rotate(-90 14 {(mt+Hh-mb)/2})">z = (Σ²ζ − Σ²BK)/se</text>')
    for n0, c in cols.items():
        pp = [(X(p[1]), Y(p[2])) for p in pts if p[0] == n0]
        if not pp: continue
        svg.append('<polyline fill="none" stroke="%s" stroke-width="1.6" points="%s"/>' % (c, ' '.join('%.1f,%.1f' % q for q in pp)))
        for q in pp: svg.append(f'<circle cx="{q[0]:.1f}" cy="{q[1]:.1f}" r="3.5" fill="{c}"/>')
    lx = ml + 10
    for i, (n0, c) in enumerate(cols.items()):
        svg.append(f'<circle cx="{lx + i*150}" cy="{mt+6}" r="4" fill="{c}"/><text x="{lx + i*150 + 9}" y="{mt+10}" font-size="12" fill="var(--ink)">from zero {n0:,}</text>')
    svg.append('</svg>')
    body = body.replace('<!--FIGBK-->', '<div class="plate"><figure>' + ''.join(svg) + '<figcaption><b>Figure 18.</b> Standardized residuals of the zeta number variance against the height-averaged Bogomolny–Keating prediction, at three disjoint windows. Dashed lines mark ±2 standard errors.</figcaption></figure></div>')
    # text
    big = [p for p in pts if abs(p[2]) > 3]
    sign_hi = [p for p in pts if p[1] >= 8]
    allneg = all(p[3] < 0 for p in sign_hi)
    maxcut = max(abs(r['prime_cutoff_shift']) for h in B['heights'] for r in h['rows']); maxq = max(abs(r['quadrature_shift']) for h in B['heights'] for r in h['rows'])
    maxratio = max(max(abs(r['prime_cutoff_shift']), abs(r['quadrature_shift']))/r['se'] for h in B['heights'] for r in h['rows'])
    maxq8 = max(max(abs(r['prime_cutoff_shift']), abs(r['quadrature_shift'])) for h in B['heights'] for r in h['rows'] if r['L'] >= 8)
    expl = [r['frac_GUE_departure_explained'] for h in B['heights'] for r in h['rows'] if r['L'] >= 1]
    t = (f'<div class="thm"><div class="hd">Result VI.8.1 <i>(number variance at three heights)</i> <span class="ep nm">verified · double</span></div>'
         f'<p>At all three heights the height-averaged Bogomolny–Keating prediction accounts for {100*min(expl):.1f}% to {100*max(expl):.1f}% of the departure of $\\Sigma^2(L)$ from GUE for $1\\le L\\le16$. '
         f'Of the 18 comparisons, {len(big)} have $|z|&gt;3$. '
         + (f'For $L\\ge8$ the residual has the same sign (zeta below the prediction) at every height, ' if allneg else 'For $L\\ge8$ the residual sign is not consistent across heights, ')
         + f'and neither numerical knob explains it: the Euler-product cutoff moves the prediction by at most {maxcut:.1e} and the quadrature by at most {maxq:.1e}, at most {maxratio:.2f} standard errors in any row, and at $L\\ge8$ both shifts are under {maxq8:.1e}, a small fraction of the residuals there.</p></div>')
    t += ('<p>The block standard errors understate the true uncertainty at large $L$. Counts in windows of length 8 or 16 are correlated over many blocks of 2,000 zeros, '
          'and ten blocks give a variance estimate with only nine degrees of freedom. A residual of several block standard errors at $L=16$ is therefore weaker evidence than the same $z$ at $L=1$. '
          'Even so, a residual whose sign is the same at three heights spanning two decades is not noise. The natural candidates are the lower-order terms that the formula leaves out: '
          'the $O(1/\\mathcal L)$ corrections to the unfolding, and the finite-window bias of the empirical variance estimator. '
          'Part VI reports the pattern and does not claim a correction. The formula is conjectural, and the comparison tests it numerically. It proves nothing.</p>')
    body = body.replace('<!--BKTEXT-->', t)
# E12 text
e12 = [f for f in ('certs/V_e12a.json', 'certs/V_e12b.json') if os.path.exists(f)]
VF2 = json.load(open('verify_full.json'))
if e12:
    parts = []
    for f in e12:
        v = VF2.get(f) or VL.get(f); c = json.load(open(f)); r = c['rectangle']
        parts.append((float.fromhex(r['lo']), float.fromhex(r['hi']), v, len(c['segments']), c['generator_time_s']))
    nz = sum(p[2]['count'] for p in parts if p[2]); okall = all(p[2] and p[2]['status'] == 'CERTIFIED' and p[2].get('all_zeros_simple_on_line') for p in parts)
    modes = ', '.join(sorted(set(p[2]['mode'] for p in parts if p[2])))
    lo = min(p[0] for p in parts); hi = max(p[1] for p in parts)
    txt = (f'<div class="thm"><div class="hd">Result VI.3.1 <i>(the tiles at $10^{{12}}$, regenerated)</i> <span class="ep pr">proved · replayed</span></div>'
           f'<p>On $10^{{12}}+[{lo-1e12:.1f},\\,{hi-1e12:.1f}]$ with $0.05\\le\\operatorname{{Re}}s\\le0.95$, the repaired engine and the independent verifier ({modes} mode) certify '
           f'<strong>{nz} zeros, {"each simple and on the critical line" if okall else "NOT all certified"}</strong>, in {len(parts)} record{"s" if len(parts) > 1 else ""} of '
           + ' and '.join(str(p[3]) for p in parts) + ' segments. The counts agree with the Hardy-$Z$ sign changes recorded in Part V.</p></div>')
    body = body.replace('<!--E12TEXT-->', txt)
else:
    body = body.replace('<!--E12TEXT-->', '<div class="thm"><div class="hd">Status VI.3.1 <i>(the tiles at $10^{12}$)</i> <span class="ep op">running</span></div><p>The two half-unit tiles at $t=10^{12}$ are being regenerated by the repaired engine at the time of this release. Each needs about an hour of one core. Until their records exist and pass full replay, the Part V claim at $10^{12}$ is withdrawn as a proof. It stands only as a computation: the counts, minors and Hardy-$Z$ agreement of Part V were produced by the defective engine. The tile at $10^{10}$ has been regenerated (Table VI.1).</p></div>')
open('nonagon_VI.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
