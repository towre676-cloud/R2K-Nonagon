import json, re
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon XII</title>')
body = open('body12.html').read(); lab = {'case118': 'IEEE 118', 'case300': 'IEEE 300', 'case1354pegase': 'PEGASE 1354', 'case2869pegase': 'PEGASE 2869'}
rows = ''
G = json.load(open('p12/local_cert_case118.json'))
rows += f"<tr><td>IEEE 118 (greedy patches)</td><td class=\"n\">179</td><td class=\"n\">{G['readings']}</td><td class=\"n\">6,903</td><td class=\"n\">{6903/G['readings']:.1f}×</td><td class=\"n\">{G['readings']/179:.2f}</td><td>yes</td><td>yes (over $\\mathbb Q$)</td><td class=\"n\">—</td></tr>\n"
for c in lab:
    d = json.load(open('p12/fastbig_%s.json' % c))
    rows += (f"<tr><td>{lab[c]}</td><td class=\"n\">{d['unknowns']:,}</td><td class=\"n\">{d['readings']:,}</td><td class=\"n\">{d['full_protocol_readings']:,}</td><td class=\"n\">{d['saving_factor']:.1f}×</td><td class=\"n\">{d['readings_per_unknown']:.2f}</td>"
             f"<td>{'yes' if d['identifiable_exact'] else 'no'}</td><td>{'yes (mod $p$)' if d['recovered_equals_truth_mod_p'] else 'no'}</td><td class=\"n\">{d['seconds_total']:.1f} s</td></tr>\n")
body = body.replace('<!--SCALEROWS-->', rows)
F = json.load(open('p11/grid_cert_case118.json')); full = {r['eps_rel']: r for r in F['rows']}
cr = ''.join(f"<tr><td class=\"n\">{r['eps_rel']:.0e}</td><td class=\"n\">{r['branches_certified_present']}/173</td><td class=\"n\">{r['max_rel_radius']:.2e}</td><td class=\"n\">{full[r['eps_rel']]['edges_certified_present']}/173</td><td class=\"n\">{full[r['eps_rel']]['max_edge_rel_radius']:.2e}</td></tr>\n" for r in G['rows'])
body = body.replace('<!--CERTROWS-->', cr)
nr = ''
for c in lab:
    d = json.load(open('p12/local_noise_%s.json' % c))
    for r in d['rows']:
        nr += f"<tr><td>{lab[c]}</td><td class=\"n\">{d['readings']:,}</td><td class=\"n\">{('%.0e' % r['sigma_rel']) if r['sigma_rel'] else '0'}</td><td class=\"n\">{r['median_rel_err']:.1e}</td><td class=\"n\">{100*r['within_1pct']:.1f}%</td><td class=\"n\">{r['branches_above_gap']:,}/{r['branches']:,}</td></tr>\n"
body = body.replace('<!--NOISEROWS-->', nr)
open('nonagon_XII.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
