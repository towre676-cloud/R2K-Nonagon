import json, re, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon XIII</title>')
body = open('body13.html').read()
lab = {'case118': 'IEEE 118', 'case300': 'IEEE 300', 'case1354pegase': 'PEGASE 1354', 'case2869pegase': 'PEGASE 2869',
       'case6515rte': 'RTE 6515', 'case9241pegase': 'PEGASE 9241', 'case_ACTIVSg10k': 'ACTIVSg 10k',
       'case13659pegase': 'PEGASE 13659', 'case_ACTIVSg25k': 'ACTIVSg 25k', 'case_ACTIVSg70k': 'ACTIVSg 70k'}
R = {c: json.load(open('p13/record_%s.json' % c)) for c in lab if os.path.exists('p13/record_%s.json' % c)}
R = {c: d for c, d in R.items() if d.get('crt_log') is not None or d['case'] == 'case118'}
def nnz(d):
    for it in reversed(d['iterations']):
        for l in it['kernel_log']:
            m = re.search(r'nnz\(L\)\s*=?\s*(\d+)', l)
            if m: return int(m.group(1))
    return None
rows = crt = ver = alts = inert = ''
for c, d in R.items():
    ex = d['exact_rational_recovery']
    rows += (f"<tr><td>{lab[c]}</td><td class=\"n\">{d['buses']:,}</td><td class=\"n\">{d['hidden_nodes']:,}</td><td class=\"n\">{d['unknowns']:,}</td><td class=\"n\">{d['readings']:,}</td>"
             f"<td class=\"n\">{d['full_protocol_readings']:,}</td><td class=\"n\">{d['saving_factor']:,.0f}×</td><td class=\"n\">{d['readings_per_unknown']:.2f}</td>"
             f"<td class=\"n\">{len(d.get('extra_source_buses', [])):,}</td><td class=\"n\">{d['physical']['negative_pivots']:,}</td><td>{'yes' if d['identifiable'] else 'no'}</td><td>{'yes' if ex else ('no' if ex is False else '—')}</td><td class=\"n\">{d['seconds_total']:,.1f} s</td></tr>\n")
    if d['identifiable']:
        crt += (f"<tr><td>{lab[c]}</td><td class=\"n\">{d['unknowns']:,}</td><td class=\"n\">{d['max_height_bits_of_truth']}</td><td class=\"n\">{len(d['primes'])}</td>"
                f"<td class=\"n\">{d['log2_modulus']}</td><td class=\"n\">{d['rational_mismatches']}</td></tr>\n")
    v1 = json.load(open('p13/verify_%s.json' % c)) if os.path.exists('p13/verify_%s.json' % c) else None
    v2 = json.load(open('p13/verify2_%s.json' % c)) if os.path.exists('p13/verify2_%s.json' % c) else None
    if v1 or v2:
        agree = all(v['agrees_with_record'] for v in (v1, v2) if v)
        r1 = f"{v1['rank']:,}" if v1 else '—'; r2 = f"{v2['rank']:,}" if v2 else '—'
        s1 = f"{v1['seconds']:.0f} s" if v1 else '—'; s2 = f"{v2['seconds']:.0f} s" if v2 else '—'
        ver += (f"<tr><td>{lab[c]}</td><td class=\"n\">{d['unknowns']:,}</td><td class=\"n\">{r1}</td><td class=\"n\">{r2}</td>"
                f"<td>{'yes' if agree else 'NO'}</td><td class=\"n\">{s1} / {s2}</td></tr>\n")
    alt = 'p13/record_%s_alt.json' % c
    if os.path.exists(alt):
        b = json.load(open(alt))
        alts += (f"<tr><td>{lab[c]}</td><td class=\"n\">{nnz(d):,}</td><td class=\"n\">{nnz(b):,}</td><td class=\"n\">{len(d['extra_source_buses']):,} / {len(b['extra_source_buses']):,}</td>"
                 f"<td class=\"n\">{len(d['primes'])} / {len(b['primes'])}</td><td class=\"n\">{b['rank_mod_p']:,}</td><td>{'yes' if b['exact_rational_recovery'] else 'no'}</td>"
                 f"<td>{'yes' if (b['solution_sha256'] and b['solution_sha256'] == d.get('solution_sha256')) else 'NO'}</td></tr>\n")
    inr = 'p13/inertia_%s.json' % c
    if not os.path.exists(inr) and os.path.exists('p13/inertia_%s_rigorous_binary64.json' % c): inr = 'p13/inertia_%s_rigorous_binary64.json' % c
    if os.path.exists(inr):
        q = json.load(open(inr))
        if 'rigorous' in inr: q['negative_index_exact'] = (str(q['negative_index_exact']) + ' (rigorous binary64)') if q['negative_index_exact'] is not None else None
        inert += (f"<tr><td>{lab[c]}</td><td class=\"n\">{q['negative_merged_branches_upper_bound']:,}</td><td class=\"n\">{q['negative_index_exact'] if q['negative_index_exact'] is not None else '≥ ' + str(q['witness_columns'])}</td>"
                  f"<td class=\"n\">{q['max_offdiag_ratio']:.1e}</td><td class=\"n\">{q['seconds']:.0f} s</td></tr>\n")
    elif d['physical']['negative_pivots'] == 0:
        inert += f"<tr><td>{lab[c]}</td><td class=\"n\">0</td><td class=\"n\">0 (positive definite by structure)</td><td class=\"n\">—</td><td class=\"n\">—</td></tr>\n"
ok = [d for d in R.values() if d['exact_rational_recovery']]
big = max(ok, key=lambda d: d['buses'])
tot_u = sum(d['unknowns'] for d in ok)
cap = json.load(open('p13/caption13.json'))
body = (body.replace('<!--RECROWS-->', rows).replace('<!--CRTROWS-->', crt).replace('<!--VERROWS-->', ver)
        .replace('<!--MAXBUSES-->', f"{big['buses']:,}").replace('<!--MAXFULL-->', f"{big['full_protocol_readings']:,}")
        .replace('<!--MAXREAD-->', f"{big['readings']:,}").replace('<!--MAXNNZ-->', f"{nnz(big):,}" if nnz(big) else 'the reported number of')
        .replace('<!--SUBTITLE-->', cap['subtitle']).replace('<!--ABSRESULT-->', cap['absresult']).replace('<!--RECTEXT-->', cap['rectext']))
S = json.load(open('p13/pocket_synth.json'))
syn = ''.join(f"<tr><td>{k}</td><td class=\"n\">{', '.join(map(str, v['observed_deficiency']))}</td><td class=\"n\">{', '.join(map(str, v['dim_N_B']))}</td><td class=\"n\">{', '.join(map(str, v['bound']))}</td></tr>\n" for k, v in S.items())
vm = max([d['buses'] for c, d in R.items() if os.path.exists('p13/verify2_%s.json' % c) or os.path.exists('p13/verify_%s.json' % c)])
body = (body.replace('<!--SYNTHROWS-->', syn).replace('<!--ALTROWS-->', alts).replace('<!--INERTIAROWS-->', inert)
        .replace('<!--VERMAX-->', f'{vm:,}').replace('<!--SURETEXT-->', cap.get('suretext', '').replace('<!--VERMAX2-->', f'{vm:,}')).replace('<!--SUREABS-->', cap.get('sureabs', '')))
open('nonagon_XIII.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body), 'exact grids', len(ok), 'unknowns', tot_u, 'largest', big['case'])
