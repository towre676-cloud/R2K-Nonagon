import json, re, os, glob
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon V</title>')
body = open('body5.html').read(); F = json.load(open('figs5.json'))
for k, v in F.items(): body = body.replace(f'<!--{k}-->', v)
D = json.load(open('bk.json'))
rows = ''.join(f"<tr><td class=\"n\">{r['L']:g}</td><td class=\"n\">{r['zeta_dense']:.4f}</td><td class=\"n\">{r['se']:.4f}</td><td class=\"n\">{r['BK']:.4f}</td><td class=\"n\">{r['GUE']:.4f}</td><td class=\"n\">{(100*r['frac_of_GUE_departure_explained']):.1f}%</td></tr>\n" for r in D['dense_rows'])
body = body.replace('<!--BKROWS-->', rows)
H = []
for f in ('hscan_e10.json', 'hscan_e12a.json', 'hscan_e12b.json'):
    if os.path.exists(f): H += json.load(open(f))
hr = ''
for r in H:
    mins = ', '.join(m.replace('+/-', '±') for m in r.get('minors', []))
    hr += f"<tr><td class=\"n\">[{r['lo']:.1f}, {r['hi']:.1f}]</td><td class=\"n\">{r['n']}</td><td class=\"n\">{r['mu0'].replace('+/-', '±')}</td><td class=\"n\">{mins}</td><td>{'yes' if r.get('all_real_simple') else 'NO'}</td><td class=\"n\">{r['hardyZ_sign_changes']}</td><td class=\"n\">{r['segments']}</td><td class=\"n\">{r['time']/60:.1f} min</td></tr>\n"
body = body.replace('<!--HROWS-->', hr)
e12 = [r for r in H if r['lo'] >= 1e12 - 1]
if e12:
    nz = sum(r['n'] for r in e12); ok = all(r.get('all_real_simple') for r in e12); agree = all(r['n'] == r['hardyZ_sign_changes'] for r in e12)
    lo = min(r['lo'] for r in e12); hi = max(r['hi'] for r in e12)
    txt = (f'<div class="thm"><div class="hd">Result V.3.1 <i>(certificates at $t=10^{{12}}$)</i> <span class="ep pr">proved · arb</span></div>'
           f'<p>On $t\\in[{lo:.1f},\\,{hi:.1f}]$, with $0.05\\le\\operatorname{{Re}}s\\le0.95$, the box criterion certifies <strong>{nz} zeros of $\\zeta$, {"every one simple and on the critical line" if ok else "NOT all real"}</strong>, by strictly positive Sylvester minors of ball matrices assembled from Taylor-model contour moments. '
           f'The independent Hardy-$Z$ counts {"agree tile by tile" if agree else "DISAGREE"}. Each tile ran on its own core, in {max(r["time"] for r in e12)/60:.0f} minutes for the slower one. '
           'The certificate uses no Gram points, no Turing method and no zero list.</p></div>')
    body = body.replace('<!--E12TEXT-->', txt)
comp = ('<p>The program is complete. It began with a structural claim: ten equations scattered across two and a half centuries are rank-two displacement kernels with one conserved quantity. It ends with that claim proved, made computable, and pushed to the edge of what is known about the zeros of $\\zeta$. ' 
        'The <strong>Chebyshev–Hermite box criterion</strong> started in Part I as a floating-point observation at $t\\le76$. It is now a rigorous certificate over the full critical strip, computed in seconds at the millionth zero and in minutes at $t=10^{12}$. Its discriminant identity measures off-line distances, as it did for Davenport–Heilbronn, and gaps of near-collisions, as it did for the classic Lehmer pair, both to sixteen or more digits without locating a zero. ' 
        'Its meter carries the full random-matrix and arithmetic statistics of the zeros: CUE locally, and a number-variance plateau now explained in level as well as shape. ' 
        'On the engineering side, the same rank-two structure gives single-sweep identification of damped mechanical chains at 48 digits, with a proved forcing law for its branch structure. It gives exact passivity certificates over $\\mathbb Q$, sampling with condition number near 1, and an identifiability threshold that comes from perturbation theory with no free constant.</p>'
        '<p>What is left is engineering, not mathematics. A batched multipoint evaluator for the Riemann–Siegel sum would replace the 64-fold repetition of a $\\sqrt t$ computation in every Taylor segment. With it, the box criterion becomes a production verification method whose certificates are inertia statements, independent of sign-change counting.</p>')
body = body.replace('<!--COMPLETION-->', comp)
open('nonagon_V.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
