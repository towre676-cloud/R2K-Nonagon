from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,platform
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ctx.prec=128

def qa(s):
 q=Q(s);return arb(q.numerator)/q.denominator

def main():
 root=json.loads((HERE/'ordinate_error_results.json').read_text());obm=json.loads((HERE/'overlap_uncertainty_results.json').read_text())
 cases=[]
 for tag,kind in [('', 'replication'),('heldout_','held_out')]:
  np=HERE/(tag+'numerical_results.json');cp=HERE/(tag+'shifted_tail_results.json')
  if not(np.exists() and cp.exists()):continue
  num=json.loads(np.read_text());cert=json.loads(cp.read_text());assert num['status']==cert['status']=='complete'
  for h in num['heights']:
   for r in h['rows']:
    cr=[x for x in cert['results'] if x['n0']==h['n0'] and x['L']==r['L']];assert len(cr)==3
    model=sum((arb(x['infinite_product_variance_ball']) for x in cr),arb(0))/3
    finite=sum((arb(x['finite_product_variance_ball']) for x in cr),arb(0))/3
    er=next(x for hh in root['heights'] if hh['n0']==h['n0'] for x in hh['rows'] if x['L']==r['L'])
    empirical=qa(r['continuous']['variance']['exact']);eb=arb(er['continuous_variance_error_bound'])
    residual=empirical+arb(0,eb.upper())-model
    ob=next(x for hh in obm['heights'] if hh['n0']==h['n0'] for x in hh['rows'] if x['L']==r['L'])
    combined=(model.rad()+eb).upper();L=arb(r['L']);z=2*arb.pi()*L
    gue=L*(1-2*z.si()/arb.pi())+(z.log()+arb.const_euler()-z.ci()+1-z.cos())/arb.pi()**2
    cases.append({'kind':kind,'n0':h['n0'],'L':r['L'],
     'model_infinite_ball':model.str(40,radius=True,more=True),'model_finite_ball':finite.str(40,radius=True,more=True),
     'empirical_variance_exact':r['continuous']['variance']['exact'],'empirical_variance_float':r['continuous']['variance']['float'],
     'empirical_ordinate_unfolding_error_upper':eb.upper().str(40,radius=True,more=True),
     'combined_deterministic_radius':combined.str(40,radius=True,more=True),'deterministic_target_pass':bool(combined<arb('1e-6')),
     'residual_ball':residual.str(40,radius=True,more=True),'residual_midpoint':float(residual.mid()),
     'primary_nominal_z':float(residual.mid())/r['nominal_block_se'],'nominal_block_se':r['nominal_block_se'],
     'nominal_5block_z':float(residual.mid())/r['block_count_sensitivity']['5']['nominal_se'],
     'nominal_20block_z':float(residual.mid())/r['block_count_sensitivity']['20']['nominal_se'],
     'obm_profiles':[{**p,'nominal_z':float(residual.mid())/p['nominal_obm_se']} for p in ob['profiles']],
     'anchor_grid_shift':r['anchor_grid_shift'],'mean_bias_squared':r['continuous']['mean_count_bias_squared']['float'],
     'height_quantile_model_spread':r['height_model_spread'],
     'Simpson_vs_finite_certificate_error_bound':(arb(r['BK_height_avg'])-finite).abs_upper().str(40,radius=True,more=True),
     'GUE_exact_ball':gue.str(40,radius=True,more=True),'GUE_Simpson_error_bound':(arb(r['GUE_numeric'])-gue).abs_upper().str(40,radius=True,more=True)})
  
 def stats(rows):
  return {'cases':len(rows),'max_absolute_primary_nominal_z':max(abs(x['primary_nominal_z']) for x in rows),
  'max_model_radius':max(float(arb(x['model_infinite_ball']).rad().upper()) for x in rows),
  'max_empirical_ordinate_error':max(float(arb(x['empirical_ordinate_unfolding_error_upper']).upper()) for x in rows),
  'max_combined_deterministic_radius':max(float(arb(x['combined_deterministic_radius']).upper()) for x in rows),
  'all_deterministic_targets_pass':all(x['deterministic_target_pass'] for x in rows),
  'max_absolute_anchor_grid_shift':max(abs(x['anchor_grid_shift']) for x in rows),
  'max_Simpson_vs_finite_error':max(float(arb(x['Simpson_vs_finite_certificate_error_bound']).upper()) for x in rows),
  'max_GUE_Simpson_error':max(float(arb(x['GUE_Simpson_error_bound']).upper()) for x in rows)}
 out={'status':'complete' if len(cases)==30 else 'production_complete_heldout_pending','cases':cases,'summaries':{}}
 for kind in ['replication','held_out']:
  rr=[x for x in cases if x['kind']==kind]
  if rr:out['summaries'][kind]=stats(rr)
 out['summaries']['all']=stats(cases)
 out['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*') if p.suffix in ['.py','.c']}
 out['python']=platform.python_version();(HERE/'summary_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['summaries'],indent=2))
if __name__=='__main__':main()
