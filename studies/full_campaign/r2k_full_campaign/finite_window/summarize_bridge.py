from pathlib import Path
from fractions import Fraction as Q
import json,math,numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;ctx.prec=128

def qa(x):q=Q(x);return arb(q.numerator)/q.denominator

def main():
 dep=json.loads((HERE/'dependence_results.json').read_text());refs=np.load(HERE/'multiplier_reference.npz');old=json.loads((BASE/'continuous_height_results.json').read_text());fixed=json.loads((HERE/'fixed_continuous_height_results.json').read_text());pair=json.loads((HERE/'pair_closure_results.json').read_text());out=[]
 for w in dep['windows']:
  n=w['n0'];kind=w['geometry'];rows=[]
  if kind=='fixed_span':data=json.loads((HERE/f'fixed_data_{n}.json').read_text());rs=data['rows'];models=fixed['results'];N=data['point_count'];S=Q(18000)
  else:
   data=json.loads((BASE/'data'/f'empirical_stable_{n}.json').read_text());rs=next(h['rows'] for h in pair['heights'] if h['n0']==n);models=old['results'];N=len(data['unfolded_float_hex']);S=Q(float.fromhex(data['unfolded_float_hex'][-1]))
  residual=[]
  for r in rs:
   L=Q(r['L']);model=arb(next(m['infinite_product_height_average_ball'] for m in models if m['n0']==n and m['L']==float(L)));v=qa(r['variance_exact']);t=arb(r['translation_variance_ball']);rho=Q(N)/S
   comp=t+qa(L*(rho-1)-L*L*(rho*rho-1))
   er=arb(r['empirical_coordinate_error_bound']) if kind=='fixed_span' else arb(next(rr['empirical_ordinate_unfolding_error_upper'] for rr in json.loads((BASE/'summary_results.json').read_text())['cases'] if rr['n0']==n and rr['L']==float(L)))
   rb=v+arb(0,er.upper())-model;residual.append(float(rb.mid()))
   rows.append({'L':float(L),'empirical_variance_exact':r['variance_exact'],'continuous_height_model_ball':model.str(45,radius=True,more=True),'residual_ball':rb.str(45,radius=True,more=True),
    'combined_deterministic_radius':(model.rad()+er).upper().str(45,radius=True,more=True),'mean_count_minus_L':str(Q(r['mean_exact'])-L),
    'translation_variance_ball':r['translation_variance_ball'],'density_compensated_translation_ball':comp.str(45,radius=True,more=True),
    'anchor_minus_translation_ball':(v-t).str(45,radius=True,more=True),
    'anchor_minus_compensated_translation_ball':(v-comp).str(45,radius=True,more=True)})
  tiers=[];residual=np.array(residual)
  for t in w['tiers']:
   z=residual/np.array(t['nominal_se']);mx=float(np.max(np.abs(z)));reference=refs[t['reference_npz_key']]
   tiers.append({'batches':t['batches'],'nominal_z':z.tolist(),'nominal_se':t['nominal_se'],'observed_max_absolute_z':mx,
    'joint_reference_95':t['reference_95'],'reference_tail_fraction':float(np.mean(reference>=mx)),
    'reference_mode':t['reference_mode'],'scope':'conditional multiplier diagnostic; no asserted finite-sample significance level'})
  for i,r in enumerate(rows):r['anchor_influence_10batch_z']=next(t['nominal_z'][i] for t in tiers if t['batches']==10)
  out.append({'geometry':kind,'n0':n,'point_count':N,'span_exact':str(S),'sample_density_exact':str(Q(N)/S),'rows':rows,'dependence_tiers':tiers})
 result={'status':'PASS','windows':out,'all_fixed_target_bounds_below_1e_minus_6':all(arb(r['combined_deterministic_radius'])<arb('1e-6') for w in out if w['geometry']=='fixed_span' for r in w['rows'])}
 (HERE/'finite_window_summary.json').write_text(json.dumps(result,indent=2))
 for w in out:
  t=next(t for t in w['dependence_tiers'] if t['batches']==10);print(w['geometry'],w['n0'],'N',w['point_count'],'max |z|',t['observed_max_absolute_z'],'joint95',t['joint_reference_95'],'reference tail',t['reference_tail_fraction'])
if __name__=='__main__':main()
