from pathlib import Path
from fractions import Fraction as Q
import json,math
from flint import arb
HERE=Path(__file__).resolve().parent

def main():
 s=json.loads((HERE/'finite_window_summary.json').read_text());p=json.loads((HERE/'pair_closure_results.json').read_text());m=json.loads((HERE/'fixed_continuous_height_results.json').read_text())
 assert s['status']=='PASS' and s['all_fixed_target_bounds_below_1e_minus_6']
 assert p['status']=='PASS' and p['selftest']['independent_exact_cases']==200
 assert m['status']=='complete' and m['all_targets_pass'] and len(m['results'])==30
 for name in ['fixed_preparation.json','moment_bridge_checks.json','sine_bias_results.json']:
  assert json.loads((HERE/name).read_text())['status']=='PASS'
 for w in s['windows']:
  assert len(w['rows'])==6
  for r in w['rows']:
   assert arb(r['combined_deterministic_radius'])<arb('1e-6')
  assert [t['batches'] for t in w['dependence_tiers']]==[5,10,20,40,80]
  if w['geometry']=='fixed_span':
   d=json.loads((HERE/f"fixed_data_{w['n0']}.json").read_text());assert d['span']==18000 and d['event_and_pair_moments_exactly_agree']
   assert arb(d['lower_buffer'])>0 and arb(d['upper_buffer'])>0
   assert all(0<float.fromhex(x)<18000 for x in d['points_hex']) and len(d['points_hex'])==d['point_count']
   for r in d['rows']:
    L=Q(r['L']);v=Q(r['variance_exact']);delta=Q(r['mean_exact'])-L
    reconstructed=arb(r['translation_variance_ball'])+arb(Q(r['mean_correction_exact']).numerator)/Q(r['mean_correction_exact']).denominator+arb(r['edge_correction_ball'])
    assert reconstructed.contains(arb(v.numerator)/v.denominator)
 target=next(w for w in s['windows'] if w['geometry']=='fixed_span' and w['n0']==20000000)
 b=json.loads((HERE/'sine_bias_results.json').read_text());cb=json.loads((HERE/'moment_bridge_checks.json').read_text())
 out={'status':'PASS','exact_original_window_cases':30,'exact_fixed_coordinate_cases':30,'fixed_coordinate_model_certificates':30,
 'points_in_fixed_coordinate_spans':sum(w['point_count'] for w in s['windows'] if w['geometry']=='fixed_span'),
 'max_fixed_combined_deterministic_radius':max(float(arb(r['combined_deterministic_radius']).upper()) for w in s['windows'] if w['geometry']=='fixed_span' for r in w['rows']),
 'fixed_index_20000000_batch_sensitivity':[{k:t[k] for k in ['batches','observed_max_absolute_z','joint_reference_95','reference_tail_fraction']} for t in target['dependence_tiers']],
 'conditional_sine_centering_biases':b['production_fixed_span_benchmarks'],'counterexample':cb['counterexample'],
 'scope':'arithmetic and analytic checks passed; no asserted statistical acceptance/rejection of a zeta model'}
 (HERE/'validation_receipt.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['counterexample','conditional_sine_centering_biases']},indent=2))
if __name__=='__main__':main()
