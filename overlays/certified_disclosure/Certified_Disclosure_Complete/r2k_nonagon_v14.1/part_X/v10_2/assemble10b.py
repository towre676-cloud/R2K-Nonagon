import json, re
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon X·2</title>')
body = open('body10b.html').read()
P = json.load(open('p10b/passivity_crosscheck.json'))
yn = lambda b: 'PR' if b else 'reject'
names = {'(s^2+1)/(s^2+s+1)  [tangency at w=1]': '$(s^2+1)/(s^2+s+1)$, tangency at $\\omega=1$'}
rows = ''
for r in P['cases']:
    if 'tangency' in r['Z'] and 'w=1' not in r['Z']: continue
    rows += f"<tr><td>{names.get(r['Z'], '$' + r['Z'] + '$')}</td><td>{yn(r['claude_VI41'])}</td><td>{yn(r['claude_X51'])}</td><td>{yn(r['v10_1_cayley_brune'])}</td></tr>\n"
for k, v in P['part_VI_ladder'].items():
    rows += f"<tr><td>Part VI ladder, {k} model</td><td>—</td><td>{yn(v['claude_X51'])}</td><td>{yn(v['v10_1'])}</td></tr>\n"
body = body.replace('<!--PASSROWS-->', rows)
E = json.load(open('p10b/edge_vs_compensated.json'))
body = body.replace('<!--EDGEROWS-->', ''.join(f"<tr><td class=\"n\">({r['edge'][0]},{r['edge'][1]})</td><td class=\"n\">{r['sd_protocol']:.2e}</td><td class=\"n\">{r['compensated_efficient_sd']:.2e}</td><td class=\"n\">{r['ratio_to_efficient']:.1f}</td><td class=\"n\">{r['ratio_to_old_floor']:.0f}</td></tr>\n" for r in E))
body = body.replace('as in the v10.1 estimation example, attains the efficient value to within Monte Carlo error.', 'attains the efficient value within the one-parameter diagonal-congruence family: the v10.1 estimation example, 80,000 replicates, gives an empirical-to-efficient ratio of 1.00004.')
body = body.replace('For every interior edge, the protocol', 'For every interior edge from $(1,2)$ to $(6,7)$, the protocol')
C = json.load(open('p10b/composition_scaling.json')); seen = set(); cr = ''
for r in C:
    if r['tiles'] in seen: continue
    seen.add(r['tiles']); fl = r['flat_parent_inertia']
    cr += f"<tr><td class=\"n\">{r['tiles']}</td><td class=\"n\">{r['zeros']}</td><td>{'yes' if r['blocks_decided'] else 'no'}</td><td class=\"n\">{('(%d,%d)' % tuple(r['block_inertia'])) if r['block_inertia'] else '—'}</td><td class=\"n\">{('(%d,%d)' % tuple(fl)) if fl else 'undecided'}</td><td class=\"n\">{r['max_cross_abs_upper']:.1e}</td></tr>\n"
body = body.replace('<!--COMPROWS-->', cr)
txt = ('<div class="thm"><div class="hd">Result X·2.3.1 <i>(composition obeys the locality law)</i> <span class="ep nm">verified · arb</span></div>'
       '<p>At 192 bits, the flat parent form in the power basis is decided up to 5 zeros (2 tiles) and no further. Composition decides the parent inertia up to 10 zeros (4 tiles), correctly $(N,0)$ in every case. From 12 zeros on, the sign of $\\det V$ itself is undecided. The basis of child root polynomials is a resultant-type matrix, and its determinant shrinks as tiles are added. The cross blocks, exactly zero in exact arithmetic, grow in enclosure from $10^{-14}$ at two tiles to $10^{-4}$ at fifteen.</p></div>'
       '<p>Composition therefore does not escape the locality law. It moves the limit for the power-basis form from 5 to 10 zeros at 192 bits, and then meets the same Chebyshev–Vandermonde conditioning in the transformation matrix. What it contributes is a consistency theorem: the parent moments overlap the sum of the child moments, and for small unions the parent inertia is the sum of the child inertias. For certifying large ranges, the conjunction of child certificates is already the complete proof. The zeros of a union of disjoint closed boxes, each with a zero-free certified boundary, are the union of their zeros, and reality and simplicity are local properties. No parent form is needed, and at scale none can be decided. The two Part X branches agree on this once it is stated: v10.1 composes small unions algebraically, and the locality law explains why large unions are certified tile by tile.</p>')
body = body.replace('<!--COMPTEXT-->', txt)
body = body.replace('<!--COMPABSTRACT-->', 'Composing $k$ tiles of the Part III scan with the v10.1 root-polynomial basis at 192 bits decides the parent inertia up to 10 zeros, against 5 for the flat power-basis form, and fails from 12 zeros on because the basis determinant becomes undecidable. Composition obeys the locality law rather than escaping it. Large ranges are certified, completely and without any parent form, as the conjunction of tile certificates.')
open('nonagon_X2.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
