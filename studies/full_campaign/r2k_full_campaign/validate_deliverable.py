from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,sys
from flint import arb
HERE=Path(__file__).resolve().parent

def main():
 s=json.loads((HERE/'summary_results.json').read_text());assert s['status']=='complete' and len(s['cases'])==30
 assert s['summaries']['all']['all_deterministic_targets_pass']
 total=0
 for n0 in [100000,1000000,10000000,20000000,50000000]:
  z=json.loads((HERE/'data'/f'zeros_{n0}.json').read_text());e=json.loads((HERE/'data'/f'empirical_stable_{n0}.json').read_text())
  assert z['n0']==n0 and z['count']==20000 and len(z['ordinate_balls'])==20000
  assert all(float.fromhex(h)==x for h,x in zip(z['ordinate_float_hex'],z['ordinates']))
  assert all(a<b for a,b in zip(z['ordinates'],z['ordinates'][1:]))
  assert e['stable_high_precision_before_subtraction'] and len(e['rows'])==6
  for r in e['rows']:
   m=r['continuous'];mu=Q(m['mean_count']['exact']);m2=Q(m['second_moment_count']['exact']);v=Q(m['variance']['exact']);L=Q(r['L'])
   assert v==m2-mu*mu and Q(m['mean_squared_deviation_from_L']['exact'])==v+(mu-L)**2 and v>=0
  total+=z['count']
 panels=0;points=0;targets=0
 for tag in ['', 'heldout_']:
  c=json.loads((HERE/(tag+'certified_results.json')).read_text());t=json.loads((HERE/(tag+'shifted_tail_results.json')).read_text())
  assert c['status']==t['status']=='complete';assert all(arb(r['infinite_product_radius'])<arb('1e-6') for r in t['results'])
  panels+=len(c['panels']);points+=c['gauss_euler_evaluations'];targets+=len(t['results'])
 for name in ['proof_checks.json','native_checks.json','native_production_checks.json','height_identity_checks.json']:
  assert json.loads((HERE/name).read_text())['status']=='PASS'
 o=json.loads((HERE/'extended_origin_results.json').read_text())['continuous_height_certificate'];assert o['certified_c_le_3'] and o['certified_positive']
 a=json.loads((HERE/'continuous_height_results.json').read_text());assert a['status']=='complete' and len(a['results'])==30 and a['all_targets_pass']
 primary={(r['n0'],r['L']):r for r in s['cases']};ar=[]
 for r in a['results']:
  pr=primary[(r['n0'],r['L'])];model=arb(r['infinite_product_height_average_ball']);emp=Q(pr['empirical_variance_exact']);emp=arb(emp.numerator)/emp.denominator
  er=arb(pr['empirical_ordinate_unfolding_error_upper']);combined=(model.rad()+er).upper();assert combined<arb('1e-6')
  diff=model-arb(pr['model_infinite_ball']);ar.append({'n0':r['n0'],'L':r['L'],'continuous_height_model_ball':model.str(45,radius=True,more=True),
   'difference_from_primary_ball':diff.str(45,radius=True,more=True),'residual_ball':(emp+arb(0,er.upper())-model).str(45,radius=True,more=True),
   'nominal_z':float((emp-model).mid())/pr['nominal_block_se'],'combined_deterministic_radius':combined.str(45,radius=True,more=True)})
 out={'status':'PASS','ordinate_count':total,'exact_empirical_cases':30,'certified_quantile_targets':targets,'quantile_gauss_panels':panels,'quantile_euler_point_evaluations':points,
 'continuous_height_cases':30,'continuous_height_euler_point_evaluations':a['euler_evaluations'],'continuous_height_gauss_panels':len(a['panels']),
 'continuous_height_max_combined_radius':max(float(arb(r['combined_deterministic_radius']).upper()) for r in ar),
 'continuous_height_max_absolute_primary_shift':max(float(arb(r['difference_from_primary_ball']).abs_upper()) for r in ar),
 'continuous_height_max_absolute_nominal_z':max(abs(r['nominal_z']) for r in ar),'continuous_height_comparisons':ar,
 'primary_summary':s['summaries'],'extended_origin_height_interval':o['T_interval'],'extended_origin_ratio_upper':o['uniform_R2_over_x2_upper']}
 (HERE/'validation_receipt.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='continuous_height_comparisons'},indent=2))
if __name__=='__main__':main()
