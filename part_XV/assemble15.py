import json, re, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon XV</title>')
body = open('body15.html').read()
lab = {'case118': 'IEEE 118', 'case300': 'IEEE 300', 'case1354pegase': 'PEGASE 1354'}
J = lambda f: json.load(open(f)) if os.path.exists(f) else None
pct = lambda x: f'{100 * x:.0f}%' if x >= 0.01 or x == 0 else f'{100 * x:.1f}%'
pc2 = lambda x: f'{100 * x:.2f}%'
def sci(x):
    e = int(f'{x:e}'.split('e')[1]); m = x / 10 ** e
    return f'{m:.0f}×10<sup>{e}</sup>' if abs(m - round(m)) < 1e-9 else f'{m:.1f}×10<sup>{e}</sup>'
names = {'inj1': 'one generator step', 'inj': 'steps at every generator', 'amb': 'ambient covariance', 'amb1': 'ambient + one calibrated step',
         'amb1l': 'ambient (distance ≤ 2 entries) + one step', 'amb1s': 'ambient (random entries) + one step', 'both': 'ambient + steps at every generator'}
rows = ''; A = {}
for c in lab:
    d = J('p15/ambient15_%s.json' % c)
    if not d: continue
    A[c] = d
    for sc, s in d['schemes'].items():
        kind, arg = sc.split(':')
        rows += (f"<tr><td>{lab[c]}</td><td>{names[kind]}</td><td class=\"n\">{arg}</td><td class=\"n\">{s['readings']:,}</td>"
                 f"<td class=\"n\">{s['rank']:,} / {s['unknowns']:,}</td><td class=\"n\">{s['identifiable_branches']:,} of {s['branches']:,}</td></tr>\n")
def ib(c, sc):
    s = A[c]['schemes'][sc]; return f"{s['identifiable_branches']:,} of {s['branches']:,}"
full = [c for c in A for sc in ('amb1:all', 'amb1l:all') if sc in A[c]['schemes'] and A[c]['schemes'][sc]['rank'] == A[c]['schemes'][sc]['unknowns']]
full = sorted(set(full), key=list(lab).index)
amb1list = '; '.join(f"{lab[c]} {A[c]['schemes'][sc]['rank']:,} of {A[c]['schemes'][sc]['unknowns']:,}" for c in full for sc in ('amb1:all', 'amb1l:all') if sc in A[c]['schemes'])
ambdef = '; '.join(f"{lab[c]} {A[c]['schemes']['amb:all']['rank']:,} of {A[c]['schemes']['amb:all']['unknowns']:,}" for c in A if 'amb:all' in A[c]['schemes'])
crbrows = ''; C = {}
for c in lab:
    d = J('p15/crb15_%s.json' % c)
    if not d: continue
    C[c] = d
    for blk in d:
        for r in blk['rows']:
            crbrows += (f"<tr><td>{lab[c]}</td><td class=\"n\">{blk['rho_sensor_noise']:g}</td><td class=\"n\">{sci(r['N'])}</td><td class=\"n\">{r['hours_at_30Hz']:.2g}</td>"
                        f"<td class=\"n\">{pc2(r['median_rel_sd'])}</td><td class=\"n\">{pc2(r['p90_rel_sd'])}</td><td class=\"n\">{pct(r['within_1pct'])}</td></tr>\n")
def crb(c, rho, N, key):
    for blk in C[c]:
        if blk['rho_sensor_noise'] == rho:
            for r in blk['rows']:
                if r['N'] == N: return r[key]
est = ''; ratios = []; e118 = None
for c in lab:
    d = J('p15/est15_%s.json' % c)
    if not d: continue
    for r in d:
        ratios.append(r['median_err_over_crb'])
        if c == 'case118' and r['N'] == 1e5 and r['rho'] == 0: e118 = r['median_rel_err']
        est += (f"<tr><td>{lab[c]}</td><td class=\"n\">{sci(r['N'])}</td><td class=\"n\">{r['rho']:g}</td><td class=\"n\">{pc2(r['median_rel_err'])}</td>"
                f"<td class=\"n\">{pct(r['within_1pct'])}</td><td class=\"n\">{r['median_err_over_crb']:.2f}</td></tr>\n")
rep = {'<!--IDROWS-->': rows, '<!--CRBROWS-->': crbrows, '<!--ESTROWS-->': est,
       '<!--INJALL118-->': ib('case118', 'inj:all'), '<!--INJALL300-->': ib('case300', 'inj:all'),
       '<!--AMB1GRIDS-->': (', '.join(lab[c] for c in full[:-1]) + ' and ' + lab[full[-1]]) if len(full) > 1 else lab[full[0]], '<!--AMB1LIST-->': amb1list, '<!--AMBDEF-->': ambdef,
       '<!--AMB1E5-->': ', '.join(f"{pc2(crb(c, 0.0, 1e5, 'median_rel_sd'))} ({lab[c]})" for c in C),
       '<!--AMBRHO-->': ', '.join(f"{pc2(crb(c, 0.01, 1e6, 'median_rel_sd'))} ({lab[c]})" for c in C),
       '<!--RHO01MED-->': ', '.join(f"{pc2(crb(c, 0.01, 1e6, 'median_rel_sd'))} ({lab[c]})" for c in C),
       '<!--FRAC1E5-->': ' to '.join([pct(min(crb(c, 0.0, 1e5, 'within_1pct') for c in C)), pct(max(crb(c, 0.0, 1e5, 'within_1pct') for c in C))]),
       '<!--ESTRATIO-->': f"{min(ratios):.2f}–{max(ratios):.2f}", '<!--E118-->': pc2(e118)}
for k, v in rep.items(): body = body.replace(k, str(v))
open('nonagon_XV.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
