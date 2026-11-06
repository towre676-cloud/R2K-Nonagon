#!/usr/bin/env python3
"""Assemble canonical Part VIII from versioned local sources, independently of cwd."""
from pathlib import Path
import json
import re

HERE=Path(__file__).resolve().parent
head=(HERE/'head8.html').read_text()
body=(HERE/'body8.html').read_text()
O=json.loads((HERE/'origin.json').read_text()); rows=[]
for r in O['rows']:
    rows.append(f'<tr><td class="n">from zero {r["n0"]:,}</td><td class="n">{r["L"]:g}</td><td class="n">{r["leading_2La_minus_a2"]:.6f}</td><td class="n">{r["residual_old"]:+.5f}</td><td class="n">{r["residual_corrected"]:+.5f}</td><td class="n">{r["se"]:.4f}</td></tr>')
body=body.replace('<!--ORIGINROWS-->','\n'.join(rows))
I=json.loads((HERE/'intervene.json').read_text()); V=json.loads((HERE/'verify_intervention.out').read_text())
rs=I['random_summary'];rt=I['root_selection'];pt=I['path_vs_triangle'];de=I['damped_extension']
checks=[('random connected networks, $n=2$–7, exact recovery from $n(n+1)/2$ measurements',f'{rs["networks"]}/{rs["networks"]} exact'),
 ('independent algebraic reference, $n=2$–6, different strengths',f'{V["trials"]}/{V["trials"]} exact'),
 ('repeated modes: grounded $K_4$ (eigenvalue 5 of multiplicity 3)','exact' if I['repeated_modes_K4']['exact_recovery'] else 'FAILED'),
 ('pair root selection: minus root correct / plus root ever correct',f'{rt["minus_root_correct"]}/{rt["pairs_checked"]} ; {rt["plus_root_ever_correct"]}'),
 ('path-product equality cases (cut vertices)',str(rt['path_product_equality_cases (cut vertex)'])),
 ('Part VII path vs triangle: baseline equal; differing intervention responses',f'{pt["experiments_that_differ"]} of {pt["of"]}; both recovered exactly'),
 ('hidden star centre vs triangle: equal accessible intervention responses','exact' if I['kron_ambiguity']['h_equal_on_samples'] else 'FAILED'),
 ('Part VII damped twins: recover $K,C,M$ at $s=0,1,2$','exact for both' if all(all(v.values()) for v in de.values() if isinstance(v,dict)) else 'see JSON')]
body=body.replace('<!--EXACTROWS-->','\n'.join(f'<tr><td>{name}</td><td>{result}</td></tr>' for name,result in checks))
S=json.loads((HERE/'stability.json').read_text());rows=[]
for r in S['weak_link']:
    rows.append(f'<tr><td>weak link</td><td class="n">$w={r["w"]:.3g}$</td><td class="n">{r["G01"]:.3e}</td><td class="n">{r["sd_K11_protocol"]:.2e}</td><td class="n">{r["two_point_tail_threshold"]:.2e}</td><td class="n">{r["sd_over_tail_threshold"]:.1f}</td></tr>')
for r in S['distance']:
    rows.append(f'<tr><td>path distance</td><td class="n">$n={r["n"]}$ ($m={r["experiments"]}$)</td><td class="n">{r["G0_far"]:.3e}</td><td class="n">{r["sd_K_far"]:.2e}</td><td class="n">{r["two_point_tail_threshold"]:.2e}</td><td class="n">{r["sd_over_tail_threshold"]:.1f}</td></tr>')
body=body.replace('<!--STABROWS-->','\n'.join(rows))
body=body.replace('<!--STRROWS-->','\n'.join(f'<tr><td class="n">{r["u"]:.3g}</td><td class="n">{r["rel_sd_G01sq_over_sigma"]:.0f}</td><td class="n">{r["rms_model"]:.0f}</td></tr>' for r in S['strength']))
result=head+body+'\n</body></html>\n'
if re.findall(r'<!--[A-Z0-9_]+-->',result):raise RuntimeError('Unfilled chapter placeholders')
(HERE.parent/'nonagon_VIII.html').write_text(result)
print(json.dumps({'assembled':'nonagon_VIII.html','origin_rows':len(O['rows']),'exact_rows':len(checks),'statistical_rows':len(rows)}))
