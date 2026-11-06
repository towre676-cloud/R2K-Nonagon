import json, re, math
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon X</title>')
body = open('body10.html').read()
H = json.load(open('p10/haynsworth.json'))
body = body.replace('<!--HAYROWS-->', ''.join(
    f"<tr><td class=\"n\">[{r['box'][0]:g}, {r['box'][1]:g}]</td><td class=\"n\">{r['zeros']}</td><td class=\"n\">{r['moments']}</td><td class=\"n\">{r['prec_bits']}</td><td>{'yes' if r['minors_decide'] else 'no'}</td><td>{'yes' if r['pivots_decide'] else 'no'}</td>"
    f"<td class=\"n\">{('%.2e' % r['worst_minor_rel_radius']) if math.isfinite(r['worst_minor_rel_radius']) else '∞'}</td><td class=\"n\">{r['worst_pivot_rel_radius']:.2e}</td></tr>\n" for r in H))
D = json.load(open('p10/dodgson.json'))
body = body.replace('<!--LAWROWS-->', ''.join(f"<tr><td class=\"n\">[{r['box'][0]}, {r['box'][1]}]</td><td class=\"n\">{r['zeros']}</td><td class=\"n\">{r['lambda_min_M0']:.2e}</td></tr>\n" for r in D['locality_law']))
F = json.load(open('p10/freq_for_pairs.json')); byn = {}
for r in F['random']: byn.setdefault(r['n'], []).append(r)
body = body.replace('<!--FFPROWS-->', ''.join(f"<tr><td class=\"n\">{n}</td><td class=\"n\">{len(v)}</td><td class=\"n\">{n*(n+2)}</td><td class=\"n\">{max(x['s_points_needed'] for x in v)}</td><td class=\"n\">{2*n-1}</td><td class=\"n\">{n*(n+1)//2}</td><td>{'all' if all(x['exact'] for x in v) else 'NO'}</td></tr>\n" for n, v in sorted(byn.items())))
E = json.load(open('p10/edge_resolution.json'))
body = body.replace('<!--EDGEROWS-->', ''.join(f"<tr><td class=\"n\">({r['edge'][0]},{r['edge'][1]})</td><td class=\"n\">{r['G0k']:.4f}</td><td class=\"n\">{r['G0k_minus_G0k1']:.4f}</td><td class=\"n\">{r['sd_w']:.2e}</td><td class=\"n\">{r['sd_times_G0k_sq']:.1f}</td><td class=\"n\">{r['sd_times_G0k1_sq']:.1f}</td><td class=\"n\">{r['sd_times_dG_sq']:.1f}</td></tr>\n" for r in E['path']))
B = json.load(open('p10/forest_and_hb.json'))
yn = lambda x: 'yes' if x else 'no'
body = body.replace('<!--HBROWS-->', ''.join(f"<tr><td>{k}</td><td>{yn(v['N_plus_D_strictly_hurwitz'])}</td><td>{yn(v['axis_nonneg'])}</td><td>{'PR' if v['PR_HB'] else 'not PR'}</td><td>{'PR' if v['PR_VI41'] else 'rejected'}</td></tr>\n" for k, v in B['hb_cases'].items()))
try:
    O = json.load(open('p10/origin_rigorous.json'))
    rows = ''.join(f"<tr><td class=\"n\">{r['T']:.4g}</td><td class=\"n\">{r['c0'].split(' + ')[0]}</td><td class=\"n\">{r['c2'].split(' + ')[0]}</td><td class=\"n\">{r['L8']['omitted']}</td><td class=\"n\">{r['L16']['omitted']}</td></tr>\n" for r in O)
    txt = ('<div class="tw"><table><caption>Table X.0. <span>Rigorous enclosures (Arb, 96 bits, DFT of 64 points, aliasing bounded). $c_0$ and $c_2$ are the real parts of the Taylor coefficients of $R_2$ at $\\varepsilon=0$; the omitted piece is $2\\int_0^{a}(L-x)(1-R_2)\\,dx$ with $a=10^{-4}$.</span></caption>'
           '<thead><tr><th>$t$</th><th>$c_0=R_2(0)$</th><th>$c_2=R_2\'\'(0)/2$</th><th>omitted, $L=8$</th><th>omitted, $L=16$</th></tr></thead><tbody>' + rows + '</tbody></table></div>')
    txt += ('<div class="thm"><div class="hd">Theorem X.0.1 <i>(the origin correction, certified)</i> <span class="ep pr">proved · arb</span></div><p>For the truncated Bogomolny–Keating formula as evaluated in Parts V–IX, at each of the three heights, the omitted contribution lies in the balls of Table X.0. Each ball has radius below $4\\times10^{-18}$ and lies below $2La-a^2$ ($1.59999\\times10^{-3}$ at $L=8$, $3.19999\\times10^{-3}$ at $L=16$) by exactly the certified remainder $2\\int_0^a(L-x)R_2\\,dx\\in[1.3,3.2]\\times10^{-11}$. The enclosures give $R_2(0)=0$ and $R_2\'(0)=0$ to within $3\\times10^{-15}$, and $R_2\'\'(0)/2=2.5108249$, $2.7623003$ and $2.9139476$ to about $2\\times10^{-8}$. These differ in the fifth digit from the Part VIII samples (2.51078 and so on), which used primes below $2\\times10^5$ rather than $2\\times10^6$: the curvature at the origin is sensitive to the Euler-product cutoff at that level, while the omitted integral is not, beyond $10^{-14}$. This replaces Proposition VIII.0.1, whose constant was only sampled.</p></div>')
except FileNotFoundError:
    txt = '<p><em>The rigorous enclosure run had not completed at assembly time.</em></p>'
body = body.replace('<!--ORIGINRIG-->', txt)
open('nonagon_X.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
