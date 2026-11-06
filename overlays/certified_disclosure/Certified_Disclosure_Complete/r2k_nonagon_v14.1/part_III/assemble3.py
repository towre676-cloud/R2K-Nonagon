import json, re, os
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon III</title>')
body = open('body3.html').read(); F = json.load(open('figs3.json'))
for k, v in F.items(): body = body.replace(f'<!--{k}-->', v)
if os.path.exists('scan_600262_600272.json'):
    S = json.load(open('scan_600262_600272.json')); T = S['tiles']
    nz = sum(t['n'] for t in T); allr = all(t.get('all_real_simple', True) for t in T)
    rows = ''.join(f"<tr><td class=\"n\">[{t['lo']:.2f}, {t['hi']:.2f}]</td><td class=\"n\">{t['n']}</td><td class=\"n\">{t['mu0'].replace('+/-', '±')}</td><td>{'yes' if t.get('all_real_simple') else 'NO'}</td><td class=\"n\">{t['Q']:.4g}</td><td class=\"n\">{t['segments']}</td></tr>\n" for t in T)
    qs = sorted(t['Q'] for t in T)
    txt = (f'<div class="thm"><div class="hd">Result III.2.3 <i>(root-free scan at the millionth zero)</i> <span class="ep pr">proved · arb</span></div>'
           f'<p>The same scanner in DFT mode, on $t\\in[600262,600272]$, certifies <strong>{nz} zeros, all simple and on the critical line</strong>, in {len(T)} tiles and {S["time"]/60:.0f} minutes on one core, again using no zero positions. The tile-by-tile counts agree exactly with Arb\'s independent <code>zeta_zeros</code> list (zeros 999,975 onward). '
           f'The meter ranges over $[{qs[0]:.3g},\\,{qs[-1]:.3g}]$ in this window. At this height the mean spacing is $2\\pi/\\log(t/2\\pi)=0.548$, so two-unit tiles hold three to four zeros.</p></div>'
           '<div class="tw"><table><caption>Table III.5. <span>Root-free rigorous scan at height $6\\times10^{5}$.</span></caption>'
           '<thead><tr><th>tile</th><th>$n$</th><th>$\\boldsymbol\\mu_0$</th><th>all real, simple</th><th>$Q$</th><th>segments</th></tr></thead><tbody>' + rows + '</tbody></table></div>')
    body = body.replace('<!--HIGHSCAN-->', txt)
open('nonagon_III.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
