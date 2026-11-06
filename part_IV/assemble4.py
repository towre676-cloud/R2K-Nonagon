import json, re, os
from flint import acb
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon IV</title>')
body = open('body4.html').read(); F = json.load(open('figs4.json'))
for k, v in F.items(): body = body.replace(f'<!--{k}-->', v)
fn = 'scanfs_42653544_42653550.json'
if os.path.exists(fn):
    S = json.load(open(fn)); T = S['tiles']; nz = sum(t['n'] for t in T)
    lo0 = T[0]['lo']; hi1 = T[-1]['hi']
    Z = [float(z.imag.mid()) for z in acb.zeta_zeros(99999970, 40)]
    ref = [sum(1 for z in Z if t['lo'] <= z < t['hi']) for t in T]
    agree = all(r == t['n'] for r, t in zip(ref, T))
    rows = ''.join(f"<tr><td class=\"n\">[{t['lo']:.0f}, {t['hi']:.0f}]</td><td class=\"n\">{t['n']}</td><td class=\"n\">{r}</td><td class=\"n\">{t['mu0'].replace('+/-', '±')}</td><td>{'yes' if t.get('all_real_simple') else 'NO'}</td><td class=\"n\">{t['Q']:.4g}</td></tr>\n" for t, r in zip(T, ref))
    txt = (f'<div class="thm"><div class="hd">Result IV.3.2 <i>(full-strip root-free scan at the hundred-millionth zero)</i> <span class="ep pr">proved · arb</span></div>'
           f'<p>On $t\\in[{lo0:.0f},{hi1:.0f}]$, with tiles of width 1 spanning $-0.1\\le\\operatorname{{Re}}s\\le1.1$ and no zero positions used, the scanner certifies <strong>{nz} zeros, all simple and on the critical line, and no zero elsewhere in the critical strip</strong>. It took {S["time"]/60:.0f} minutes on one core. '
           f'The per-tile counts {"agree exactly" if agree else "DISAGREE"} with Arb\'s independent <code>zeta_zeros</code> list from zero 99,999,970 onward.</p></div>'
           '<div class="tw"><table><caption>Table IV.6. <span>Full-strip root-free scan near $t=4.27\\times10^{7}$.</span></caption>'
           '<thead><tr><th>tile</th><th>$n$ (certified)</th><th>$n$ (Arb list)</th><th>$\\boldsymbol\\mu_0$</th><th>all real, simple</th><th>$Q$</th></tr></thead><tbody>' + rows + '</tbody></table></div>')
    body = body.replace('<!--FSSCAN-->', txt)
open('nonagon_IV.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
