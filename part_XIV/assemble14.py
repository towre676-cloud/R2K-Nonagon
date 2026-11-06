import json, re, os, glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon XIV</title>')
body = open('body14.html').read()
lab = {'case118': 'IEEE 118', 'case300': 'IEEE 300', 'case1354pegase': 'PEGASE 1354', 'case2869pegase': 'PEGASE 2869',
       'case6515rte': 'RTE 6515', 'case9241pegase': 'PEGASE 9241', 'case_ACTIVSg10k': 'ACTIVSg 10k'}
order = list(lab)
def J(f):
    return json.load(open(f)) if os.path.exists(f) else None
def sci(x):
    if x == 0: return '0'
    e = int(f'{x:e}'.split('e')[1]); m = x / 10 ** e
    return f'{m:.1f}×10<sup>{e}</sup>'
pct = lambda x: f'{100 * x:.0f}%' if x >= 0.01 or x == 0 else f'{100 * x:.1f}%'
crb = {c: J('p14/crb14_%s.json' % c) for c in order}
r1 = r2 = ''
lean_sig = []; lean_fr3 = []
for c in order:
    d = crb.get(c)
    if not d: continue
    for sc in ('lean', 'full'):
        s = d['schemes'].get(sc)
        if not s: continue
        if sc == 'lean': lean_sig.append(s['sigma_for_1pct_median']); lean_fr3.append(s['frac_1pct_at_1e-3'])
        r1 += (f"<tr><td>{lab[c]}</td><td>{'one port, Part XIII' if sc == 'lean' else 'one port, full protocol'}</td><td class=\"n\">{s['readings']:,}</td>"
               f"<td class=\"n\">{s['median_k']:,.0f}</td><td class=\"n\">{s['p90_k']:,.0f}</td><td class=\"n\">{sci(s['sigma_for_1pct_median'])}</td>"
               f"<td class=\"n\">{pct(s['frac_1pct_at_1e-4'])}</td><td class=\"n\">{pct(s['frac_1pct_at_1e-3'])}</td></tr>\n")
    for sc, s in d['schemes'].items():
        if not sc.startswith('ports'): continue
        r2 += (f"<tr><td>{lab[c]}</td><td class=\"n\">{sc[5:]}</td><td class=\"n\">{s['readings']:,}</td><td class=\"n\">{s['median_k']:,.1f}</td>"
               f"<td class=\"n\">{sci(s['sigma_for_1pct_median'])}</td><td class=\"n\">{pct(s['frac_1pct_at_1e-4'])}</td><td class=\"n\">{pct(s['frac_1pct_at_1e-3'])}</td>"
               f"<td class=\"n\">{s['unresolved_branches']}</td></tr>\n")
gn = ''; ratios = []
for c in order:
    d = J('p14/gn14_%s.json' % c)
    if not d: continue
    for sc in ('lean', 'ports10', 'ports20'):
        for r in d.get(sc, []):
            if r['chi2_per_reading'] < 3: ratios.append(r['median_err_over_crb'])
            gn += (f"<tr><td>{lab[c]}</td><td>{'one port' if sc == 'lean' else sc[5:] + ' ports'}</td><td class=\"n\">{sci(r['sigma'])}</td>"
                   f"<td class=\"n\">{sci(r['median_rel_err'])}</td><td class=\"n\">{sci(r['p90_rel_err'])}</td><td class=\"n\">{pct(r['within_1pct'])}</td>"
                   f"<td class=\"n\">{r['median_err_over_crb']:.2f}</td><td class=\"n\">{r['chi2_per_reading']:.2f}</td></tr>\n")
g118 = J('p14/gn14_case118.json')
p20 = [r for r in g118['ports20'] if r['sigma'] == 1e-3][0]; l6 = [r for r in g118['lean'] if r['sigma'] == 1e-6][0]
idr = ''; idd = {}
for c in order:
    d = J('p14/ident14_%s.json' % c)
    if not d: continue
    idd[c] = d
    for sc, s in d['schemes'].items():
        kind = ''.join(ch for ch in sc if ch.isalpha()); k = ''.join(ch for ch in sc if ch.isdigit())
        acc = {'kron': 'ports only, no interventions', 'hw': 'interventions at shunt and generator buses', 'ports': 'intervention at every bus'}[kind]
        idr += (f"<tr><td>{lab[c]}</td><td>{acc}</td><td class=\"n\">{k}</td><td class=\"n\">{s['readings']:,}</td>"
                f"<td class=\"n\">{s['rank']:,} / {s['unknowns']:,}</td><td class=\"n\">{s['identifiable_branches']:,} of {s['branches']:,}</td></tr>\n")
def ib(c, sc):
    s = idd[c]['schemes'][sc]; return f"{s['identifiable_branches']:,} of {s['branches']:,} ({100 * s['identifiable_branches'] / s['branches']:.0f}%)"
pl = [J('p14/portloop14_%s.json' % c) for c in order]; pl = [x for x in pl if x]
pltxt = '; '.join(f"{x['ports_needed']} ports on {lab[x['case']]} ({x['n']:,} buses)" if x['ports_needed'] else f"no completion within the loop budget on {lab[x['case']]}" for x in pl)
rs = ''
for c in order:
    d = J('p14/rootsel14_%s.json' % c)
    if not d: continue
    for r in d['rows']:
        if r['theta_times_median_Gjj'] < 0.1: continue
        rs += (f"<tr><td>{lab[c]}</td><td class=\"n\">{d['pairs']:,}</td><td class=\"n\">{r['theta_times_median_Gjj']:g}</td><td class=\"n\">{r['theta_G00']:.2g}</td>"
               f"<td class=\"n\">{r['minus_fail']}</td><td class=\"n\">{r['two_fail']}</td><td class=\"n\">{sci(r['min_gap'])}</td></tr>\n")
rep = {'<!--CRB1ROWS-->': r1, '<!--CRB2ROWS-->': r2, '<!--GNROWS-->': gn, '<!--IDROWS-->': idr, '<!--RSROWS-->': rs,
       '<!--LEANSIG-->': f"{sci(min(lean_sig))} to {sci(max(lean_sig))}", '<!--SIGRANGE-->': f"{sci(min(lean_sig))} and {sci(max(lean_sig))}",
       '<!--MAXFRAC1E3-->': pct(max(lean_fr3)), '<!--GNRATIO-->': f"{min(ratios):.2f}–{max(ratios):.2f}",
       '<!--P20ERR-->': f"{100 * p20['median_rel_err']:.2f}%", '<!--P20FRAC-->': pct(p20['within_1pct']), '<!--GN118E6-->': sci(l6['median_rel_err']),
       '<!--HW118-->': ib('case118', 'hw20'), '<!--HW300-->': ib('case300', 'hw20'), '<!--HW1354-->': ib('case1354pegase', 'hw5'),
       '<!--KRON118-->': idd['case118']['schemes']['kron60']['identifiable_branches'], '<!--KRON300-->': idd['case300']['schemes']['kron60']['identifiable_branches'],
       '<!--KRON1354-->': idd['case1354pegase']['schemes']['kron60']['identifiable_branches'], '<!--PORTLOOP-->': pltxt}
for k, v in rep.items(): body = body.replace(k, str(v))
open('nonagon_XIV.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
