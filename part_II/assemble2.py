import json, re, ast, os
head = open('css_head.txt').read(); body = open('body2.html').read(); F = json.load(open('figs2.json'))
for k, v in F.items(): body = body.replace(f'<!--{k}-->', v)
L = open('rig_p1.log').read()
tiles = [ast.literal_eval(line) for line in L.splitlines() if line.startswith("{'j':")]
rows = ''.join(f"<tr><td>{t['j']}</td><td class=\"n\">{t['mu0']}</td><td class=\"n\">{t['zero']}</td><td class=\"n\">{t['arb_zeta_zero']}</td><td>{'yes' if t['overlap'] else 'NO'}</td></tr>\n" for t in tiles)
tab = ('<div class="tw"><table>\n<caption>Table II.3. <span>Rigorous single-zero tiles: Arb enclosures of $\\mu_0$ and of the zero $c+a\\mu_1/\\mu_0$, against Arb\'s independent <code>zeta_zero</code>.</span></caption>\n'
       '<thead><tr><th>$j$</th><th>$\\boldsymbol\\mu_0$</th><th>tile enclosure of $\\gamma_j$</th><th><code>zeta_zero</code></th><th>overlap</th></tr></thead>\n<tbody>\n' + rows + '</tbody></table></div>')
body = body.replace('<!--TILETABLE-->', tab)
from flint import arb, ctx
ctx.prec = 80
det = arb('-0.380731407483928332', '3.46e-19')
dist = (-det).sqrt()/2
body = body.replace('<!--DHDIST-->', '$|\\operatorname{Re}s-\\tfrac12|\\in' + dist.str(19, radius=True).replace('+/-', '\\pm') + '$')
lp = json.load(open('lehmer_tm.json')); f = lambda x: x.replace('+/-', '\\pm')
txt = ("<p>The tile $t\\in[7004.95,\\,7005.22]\\times[-0.45,0.45]$ around $\\gamma_{6709},\\gamma_{6710}$ was certified by the Taylor-model integrator of Theorem II.1.7, "
       f"in {lp['time']:.0f} s with {sum(lp['segments'])} segments at 80 bits. It gives $\\boldsymbol\\mu_0 = {f(lp['mu0'])}$ and Sylvester minors ${f(lp['minors'][0])}$ and ${f(lp['minors'][1])}$, both strictly positive, "
       "so both zeros are simple and on the critical line. The discriminant meter encloses their gap, without locating either zero, as "
       f"$$\\gamma_{{6710}}-\\gamma_{{6709}} = a\\sqrt{{\\det M_0}}\\ \\in\\ {f(lp['gap_from_det'])},$$ which overlaps the difference of Arb's independent <code>zeta_zero</code> enclosures, ${f(lp['gap_ref'])}$.</p>")
body = body.replace('<!--LEHMERTM-->', txt)
tc = json.load(open('tmcost.json'))
rows = ''.join(f"<tr><td class=\"n\">{r['n']}</td><td class=\"n\">{r['t']:.4f}</td><td class=\"n\">{r['segments']}</td><td class=\"n\">{r['time']:.1f}</td><td class=\"n\">{f(r['mu0'])}</td><td class=\"n\">{f(r['zero'])}</td><td>{'yes' if r['overlap'] else 'NO'}</td></tr>\n" for r in tc)
tab = ('<div class="tw"><table>\n<caption>Table II.4. <span>Taylor-model certification of single-zero tiles up the critical line: segments, wall time on one core, the count enclosure, the zero enclosure, and agreement with Arb\'s <code>zeta_zero</code>.</span></caption>\n'
       '<thead><tr><th>zero $n$</th><th>$\\gamma_n$</th><th>segments</th><th>s</th><th>$\\boldsymbol\\mu_0$</th><th>enclosure of $\\gamma_n$</th><th>overlap</th></tr></thead>\n<tbody>\n' + rows + '</tbody></table></div>')
body = body.replace('<!--TMCOST-->', tab)
open('nonagon_II.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body), len(head + body))
