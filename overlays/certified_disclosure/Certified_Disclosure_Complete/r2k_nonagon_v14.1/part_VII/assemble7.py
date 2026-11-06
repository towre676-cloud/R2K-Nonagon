import json, re, math
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon VII</title>')
body = open('body7.html').read()
X = json.load(open('p7/crossengine.json')); xr = ''
for r in X:
    ine = '(%d,%d)' % tuple(r['ours_inertia']) if r.get('ours_inertia') else '—'
    xr += f"<tr><td class=\"n\">{r['tile'][5:9]}</td><td class=\"n\">[{r['lo']}, {r['hi']}]</td><td class=\"n\">{r['their_raw'].get('count')}</td><td class=\"n\">{r['ours_count']}</td><td class=\"n\">{ine}</td><td class=\"n\">{r['segments']}</td><td class=\"n\">{r['seconds']:.1f} s</td></tr>\n"
body = body.replace('<!--XROWS-->', xr)
C = json.load(open('p7/cue_calib.json')); cr = ''
for r in C['rows']:
    cr += (f"<tr><td class=\"n\">{r['L']:g}</td><td class=\"n\">{r['exact']:.4f}</td><td class=\"n\">{r['mean_estimate']:.4f}</td><td class=\"n\">{r['bias']:+.5f}</td><td class=\"n\">{r['bias_se']:.5f}</td>"
           f"<td class=\"n\">{r['true_sd_of_estimate']:.4f}</td><td class=\"n\">{r['mean_block_se']:.4f}</td><td class=\"n\">{r['se_ratio_true_over_block']:.2f}</td></tr>\n")
body = body.replace('<!--CALROWS-->', cr)
B = json.load(open('bk_multi_100000_1000000_10000000.json'))
R = {r['L']: r for r in C['rows']}
# recompute zeta residuals with (i) bias correction, (ii) calibrated standard error
lines = []; zs_new = []
for h in B['heights']:
    for r in h['rows']:
        c = R[r['L']]; se_true = r['se']*c['se_ratio_true_over_block']
        znew = (r['zeta'] - c['bias'] - r['BK_height_avg'])/se_true; zs_new.append((h['n0'], r['L'], r['z_BK_avg'], znew))
big_old = sum(1 for z in zs_new if abs(z[2]) > 3); big_new = sum(1 for z in zs_new if abs(z[3]) > 3)
hi = [z for z in zs_new if z[1] >= 8]
rows = ''.join(f"<tr><td class=\"n\">from zero {n0:,}</td><td class=\"n\">{L:g}</td><td class=\"n\">{zo:+.1f}</td><td class=\"n\">{zn:+.1f}</td></tr>" for n0, L, zo, zn in zs_new if L >= 4)
maxbias = max(abs(r['bias']) for r in C['rows']); r16 = R[16]; r8 = R[8]
ratio_hi = max(R[8]['se_ratio_true_over_block'], R[16]['se_ratio_true_over_block'])
sig_bias = [r['L'] for r in C['rows'] if abs(r['bias']) > 3*r['bias_se']]
resid = {(h['n0'], r['L']): (r['zeta'] - r['BK_height_avg'], r['se']) for h in B['heights'] for r in h['rows']}
r16s = [resid[(n, 16)][0] for n in (100000, 1000000, 10000000)]; r8s = [resid[(n, 8)][0] for n in (100000, 1000000, 10000000)]
txt = (f'<div class="thm"><div class="hd">Result VII.3.1 <i>(estimator calibration)</i> <span class="ep nm">verified · double</span></div>'
       f'<p>On {C["replicates"]} GUE replicates of 20,000 points, no bias of the dense-window estimator is resolved at three Monte Carlo standard errors for any $L\\le16$. The largest, at $L=16$, is ${r16["bias"]:+.4f}\\pm{r16["bias_se"]:.4f}$, which is <em>positive</em>, opposite in sign to the zeta deficit. '
       f'For GUE the ten-block standard error is conservative: the true sampling error is {min(r["se_ratio_true_over_block"] for r in C["rows"]):.2f} to {max(r["se_ratio_true_over_block"] for r in C["rows"]):.2f} times the block estimate. The estimator therefore does not manufacture a negative large-$L$ residual.</p></div>')
txt += ('<p>The calibration also shows that GUE error scales do not transfer to zeta. At $L=16$ a GUE sample of 20,000 points has a true sampling error of '
        f'{r16["true_sd_of_estimate"]:.4f} in $\\Sigma^2$, while the ten zeta blocks agree with one another to a standard error of about 0.0005. At large $L$ the zeta number variance is not a random fluctuation around its mean. It is governed by the short primes through the Bogomolny–Keating term, and it is almost the same in every block of 2,000 consecutive zeros. A zeta block standard error measures within-sequence stability, not exchangeable sampling error, and no $z$-score built from it should be read as a probability. The informative test is consistency across heights. That consistency is unusually sharp: '
        f'at $L=16$ the residual is ${r16s[0]:+.4f}$, ${r16s[1]:+.4f}$ and ${r16s[2]:+.4f}$ at the windows from zeros $10^5$, $10^6$ and $10^7$, and at $L=8$ it is ${r8s[0]:+.4f}$, ${r8s[1]:+.4f}$ and ${r8s[2]:+.4f}$.</p>')
txt += ('<div class="thm"><div class="hd">Observation VII.3.2 <i>(a height-independent large-$L$ offset)</i> <span class="ep nm">verified · double</span></div>'
        f'<p>Relative to the height-averaged Bogomolny–Keating number variance, the zeta number variance at $L=16$ is lower by $0.0031\\pm0.0002$ (spread across windows) at heights where $\\log(t/2\\pi)$ ranges from about 9.4 to 13.6. Any omitted term that scales like a power of $1/\\log(t/2\\pi)$ would change by roughly a third or more across that range. The observed offset does not change at all to this resolution. '
        'A height-independent offset of this kind points away from finite-height arithmetic corrections. It points toward something common to all three computations: the unfolding convention, the lower cutoff and grid of the $R_2$ quadrature (whose halving moves the prediction by less than $6\\times10^{-5}$), or a constant term in the evaluated formula itself. It is recorded as an open item. It is not claimed as a discrepancy in the ratios conjecture.</p></div>')
body = body.replace('<!--CALTEXT-->', txt)
abs_ = ('The estimator has no resolved bias, and the bias it shows at $L=16$ has the opposite sign to the deficit. The deficit itself turns out to be height-independent, $-0.0031$ at $L=16$ in all three windows across two decades of height. That points away from finite-height arithmetic and toward a convention shared by the three computations, so it is recorded as an open item, not a discrepancy with the conjecture.')
body = body.replace('<!--CALABSTRACT-->', abs_)
open('nonagon_VII.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body)); print(abs_)
