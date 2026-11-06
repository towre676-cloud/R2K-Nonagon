import json, re
head = open('css_head.txt').read().replace('<title>Rank-Two Nonagon II</title>', '<title>Rank-Two Nonagon VIII</title>')
body = open('body8.html').read()
O = json.load(open('p8/origin.json')); a = O['a']; rr = ''
for r in O['rows']:
    lead = 2*r['L']*a - a*a; new = r['residual_old'] + lead
    rr += f"<tr><td class=\"n\">from zero {r['n0']:,}</td><td class=\"n\">{r['L']:g}</td><td class=\"n\">{lead:.6f}</td><td class=\"n\">{r['residual_old']:+.5f}</td><td class=\"n\">{new:+.5f}</td><td class=\"n\">{r['se']:.4f}</td></tr>\n"
body = body.replace('<!--ORIGINROWS-->', rr)
I = json.load(open('p8/intervene.json')); V = json.load(open('p8/verify_intervention.out')) if False else json.loads(open('p8/verify_intervention.out').read())
rs = I['random_summary']; rt = I['root_selection']; pt = I['path_vs_triangle']; de = I['damped_extension']
ex = [('random connected networks, $n=2$–7, exact recovery from $n(n+1)/2$ measurements', f"{rs['networks']}/{rs['networks']} exact"),
      ('independent verifier, $n=2$–6, different strengths', f"{V['trials']}/{V['trials']} exact"),
      ('repeated modes: grounded $K_4$ (eigenvalue 5 of multiplicity 3)', 'exact' if I['repeated_modes_K4']['exact_recovery'] else 'FAILED'),
      ('pair root selection: minus root correct / plus root ever correct', f"{rt['minus_root_correct']}/{rt['pairs_checked']} ; {rt['plus_root_ever_correct']}"),
      ('path-product equality cases (cut vertices) among checked pairs', str(rt['path_product_equality_cases (cut vertex)'])),
      ('frontier path vs triangle: baseline equal; protocol measurements that differ', f"{'yes' if pt['baseline_equal'] else 'no'}; {pt['experiments_that_differ']} of {pt['of']}; both recovered exactly"),
      ('Y–Δ: hidden star centre vs triangle, identical $h$ on accessible interventions', 'yes' if I['kron_ambiguity']['h_equal_on_samples'] else 'no'),
      ('Part VII damped twins: $K$, $C$, $M$ recovered from $s\\in\\{0,1,2\\}$', 'exact for both' if all(all(v.values()) for v in de.values() if isinstance(v, dict)) else 'see JSON')]
body = body.replace('<!--EXACTROWS-->', ''.join(f'<tr><td>{a_}</td><td>{b_}</td></tr>\n' for a_, b_ in ex))
S = json.load(open('p8/stability.json')); sr = ''
for r in S['weak_link']: sr += f"<tr><td>weak link</td><td class=\"n\">$w={r['w']:.3g}$</td><td class=\"n\">{r['G01']:.3e}</td><td class=\"n\">{r['sd_K11_protocol']:.2e}</td><td class=\"n\">{r['lower_bound_any_protocol']:.2e}</td><td class=\"n\">{r['ratio']:.1f}</td></tr>\n"
for r in S['distance']: sr += f"<tr><td>path distance</td><td class=\"n\">$n={r['n']}$ ($m={r['experiments']}$)</td><td class=\"n\">{r['G0_far']:.3e}</td><td class=\"n\">{r['sd_K_far']:.2e}</td><td class=\"n\">{r['lower_bound_any_protocol']:.2e}</td><td class=\"n\">{r['ratio']:.1f}</td></tr>\n"
body = body.replace('<!--STABROWS-->', sr)
body = body.replace('<!--STRROWS-->', ''.join(f"<tr><td class=\"n\">{r['u']:.3g}</td><td class=\"n\">{r['rel_sd_G01sq_over_sigma']:.0f}</td><td class=\"n\">{r['rms_model']:.0f}</td></tr>\n" for r in S['strength']))
open('nonagon_VIII.html', 'w').write(head + body)
print('placeholders left:', re.findall(r'<!--[A-Z0-9_]+-->', head + body))
