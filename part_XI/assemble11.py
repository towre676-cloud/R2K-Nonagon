import json, re, os
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon XI</title>')
body = open('body11.html').read()
cases = ['case118', 'case300', 'case1354pegase', 'case2869pegase']; lab = {'case118': 'IEEE 118', 'case300': 'IEEE 300', 'case1354pegase': 'PEGASE 1354', 'case2869pegase': 'PEGASE 2869'}
gr = ''
for c in cases:
    d = json.load(open('p11/grid_noise_%s.json' % c)); q = d['port_correlation_sq_quantiles']
    gr += f"<tr><td>{lab[c]}</td><td class=\"n\">{d['n_hidden']:,}</td><td class=\"n\">{d['edges']:,}</td><td class=\"n\">{d['negative_weight_edges']}</td><td class=\"n\">{d['readings']:,}</td><td class=\"n\">{('%.3g' % q['0']) if q['0'] > 0 else 'n/a (out of class)'}</td><td class=\"n\">{q['0.5']:.3g}</td></tr>\n"
body = body.replace('<!--GRIDROWS-->', gr)
er = ''; peg = ''
for c in cases:
    fn = 'p11/grid_exact_%s.json' % c
    if not os.path.exists(fn): continue
    d = json.load(open(fn))
    er += f"<tr><td>{lab[c]}</td><td class=\"n\">{d['readings']:,}</td><td>{'yes' if d['exact_recovery_K'] else 'NO'}</td><td class=\"n\">{d['topology_precision']:.3f}</td><td class=\"n\">{d['topology_recall']:.3f}</td><td class=\"n\">{d['spot_checks_direct_solve']}</td><td class=\"n\">{d['max_digits_in_G']}</td><td class=\"n\">{d['seconds']['total']:.1f} s</td></tr>\n"
    if c == 'case1354pegase' and d['exact_recovery_K']: peg = f"It recovers <strong>PEGASE 1354 exactly</strong> from {d['readings']:,} readings in {d['seconds']['total']/60:.0f} minutes. "
body = body.replace('<!--EXACTROWS-->', er)
body = body.replace('<!--PEGEXACT-->', peg or 'Exact rational disclosure of the PEGASE cases exceeded this environment within the time budget, so they are treated in binary64 and in ball arithmetic below. ')
cr = ''
for c in cases:
    fn = 'p11/grid_cert_%s.json' % c
    if not os.path.exists(fn): continue
    d = json.load(open(fn))
    for r in d['rows']:
        if r.get('inverse_failed'): cr += f"<tr><td>{lab[c]}</td><td class=\"n\">{r['eps_rel']:.0e}</td><td colspan=\"5\">interval inverse failed</td></tr>\n"; continue
        cr += f"<tr><td>{lab[c]}</td><td class=\"n\">{r['eps_rel']:.0e}</td><td class=\"n\">{r['edges_certified_present']}/{r['edges']}</td><td class=\"n\">{r['nonedges_certified_below_gap']:,}/{r['nonedges']:,}</td><td class=\"n\">{r['max_edge_rel_radius']:.2e}</td><td>{'yes' if r['all_enclosures_contain_truth'] else 'NO'}</td><td class=\"n\">{r['seconds']:.1f} s</td></tr>\n"
body = body.replace('<!--CERTROWS-->', cr)
nr = ''
for c in cases:
    d = json.load(open('p11/grid_noise_%s.json' % c))
    for r in d['rows']:
        nr += (f"<tr><td>{lab[c]}</td><td class=\"n\">{('%.0e' % r['sigma_rel']) if r['sigma_rel'] else '0'}</td><td class=\"n\">{r['precision']:.3f}</td><td class=\"n\">{r['recall']:.3f}</td><td class=\"n\">{r['median_edge_rel_err']:.1e}</td>"
               f"<td class=\"n\">{100*r['frac_edges_within_1pct']:.1f}%</td><td class=\"n\">{r['law_predicted_resolvable']}</td><td class=\"n\">{r['law_false_promises']}</td></tr>\n")
body = body.replace('<!--NOISEROWS-->', nr)
V = json.load(open('p11/vf_bench.json')); vr = ''
for r in V:
    for m, nm in (('r2k', 'R2K'), ('vector_fitting', 'Vector Fitting')):
        x = r[m]; pe = x.get('passive_exact_full_axis'); pe = ('yes' if pe is True else 'no' if pe is False else pe)
        sk = x.get('skrf_is_passive', '—'); sk = 'yes' if sk is True else 'no' if sk is False else sk
        vr += (f"<tr><td>{r['data'] if m == 'r2k' else ''}</td><td>{nm}</td><td class=\"n\">{x['order']}</td><td class=\"n\">{x['conf_nrmse']:.3g}</td><td class=\"n\">{x['conf_max_err_rel']:.3g}</td><td class=\"n\">{x['seconds']:.1f} s</td>"
               f"<td>{'yes' if x['stable_certified_exact'] else 'no'}</td><td>{pe}</td><td>{sk if m == 'vector_fitting' else '—'}</td></tr>\n")
body = body.replace('<!--VFROWS-->', vr)
open('nonagon_XI.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
