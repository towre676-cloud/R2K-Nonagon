#!/usr/bin/env python3
"""Run all Part X campaigns, validate their result gates, write a single receipt."""
from pathlib import Path
import sys,os,subprocess,time,json,hashlib,platform,importlib.metadata
HERE=Path(__file__).resolve().parent
JOBS=['frame','frame_formalism','occupation','ground_submodularity','forests','dynamic','certification','nested_boxes','origin_cauchy','design','frequency_for_pairs','claude_uniform_rank','lossless_passivity','basis_policy','compensated_edge']
def require(condition,message):
 if not condition:raise RuntimeError(message)
def validate(name,d):
 if name in ('frame','occupation'):require(d['verdict']=='PASS',name+' verdict')
 elif name=='frame_formalism':require(d['total_assertion_groups']==175,'exact frame checks incomplete')
 elif name=='ground_submodularity':
  from fractions import Fraction
  require(d['grounding_submodularity_inequalities']==333,'ground inequalities incomplete')
  require(Fraction(d['minimum_exact_gap'])>0,'ground margin')
  require(Fraction(d['edge_counterexample']['submodularity_gap_for_reduction'])<0,'edge boundary control')
 elif name=='forests':
  require(d['exact_forests']['all_inverse_entries_and_six_clamp_identities_exact'],'forest identities')
  require(d['envelope_scan']['all_checks_passed'] and d['envelope_scan']['edge_directions']==672,'edge envelopes')
 elif name=='dynamic':
  for key in ['exact_clamp_all_pass','finite_J_all_exact_nonsingular','finite_leading_coefficient_all_exact','imaginary_finite_all_full_rank']:
   require(d['summary'][key],key)
  require(d['summary']['exact_clamp_networks']==21,'dynamic count')
 elif name=='certification':
  require(d['records_matched']==49 and not d['failures'],'archived contour campaign')
  f=d['full_first_zero_verifier']
  require(f['status']=='CERTIFIED' and f['segments_recomputed']==172,'first-zero analytic replay')
  require(d['zero_pivot_rejected'],'zero-pivot guard')
 elif name=='nested_boxes':
  require(d['parent_count']==5 and d['nested_parent_inertia']==[5,0],'spatial inertia')
  require(d['direct_parent_moments_overlap_child_sum'] and d['transformed_cross_blocks_contain_zero'],'moment gluing')
  for v in d['independent_full_child_verifiers']+[d['independent_full_union_verifier']]:
   require(v['status']=='CERTIFIED' and v['segments_recomputed']==v['segments_checked'],'full nested contour')
 elif name=='origin_cauchy':
  c=d['continuous_height_range_infinite']
  require(d['precision_bits']==192 and d['prime_cutoff']==2000000,'origin parameters')
  require(c['certified_positive'] and c['certified_c_le_3'],'uniform infinite-product origin bound')
 elif name=='design':
  require(d['required_checks']['passed']==d['required_checks']['count'],'design checks')
  require(sum(c['submodularity_elementary_checks'] for c in d['cases'])==34560,'design count')
  require(all(abs(c['greedy_D_efficiency_vs_optimal']-1)<1e-10 for c in d['cases']),'exhaustive optima')
 elif name=='frequency_for_pairs':
  require(d['summary']['exact_all_pass'] and d['summary']['exact_networks']==40,'frequency physical recoveries')
  require(d['summary']['observation_only_API_all_pass'] and d['summary']['negative_controls_all_rejected'],'observation-only acceptance gate')
  require(d['nongeneric_example']['rank_with_8_frequencies']==7,'column rank gate boundary')
 elif name=='claude_uniform_rank':
  require(d['all_full_rank'] and len(d['results'])==11,'uniform algebraic witnesses')
  require(d['results'][-1]['rank']==168,'n12 exact rank')
 elif name=='lossless_passivity':
  require(d['all_expected_decisions_agree'] and d['total_decision_cases']==27,'exact PR decisions')
  require(d['cases']['even_positive_root_counterexample']['positive_real'],'even tangency acceptance')
 elif name=='basis_policy':
  require(d['all50_inertias_decide'] and d['fallback_selected']==1 and d['structured_selected']==49,'basis fallback')
  fixture=HERE/'ten_zero_contour_fixture.json';meta=json.loads((HERE/'ten_zero_fixture_provenance.json').read_text())
  require(hashlib.sha256(fixture.read_bytes()).hexdigest()==meta['fixture_sha256'],'stress fixture hash')
  require(meta['full_verifier']['status']=='CERTIFIED' and meta['full_verifier']['count']==10,'stress fixture analytic evidence')
 elif name=='compensated_edge':
  require(d['random_checks']['all_physical'] and d['random_checks']['finite_controls']==600,'compensated physical controls')
  require(d['exact_fixture']['identity_exact'],'exact compensated optimum')
  require(d['estimation']['finite_risk_bounded_extension'] and d['estimation']['clipped_count']==0,'bounded estimator')
 else:raise RuntimeError('Unknown job')
def main():
 require(not sys.flags.optimize,'Run without -O: campaign assertions are part of the numerical checks.')
 logs=HERE/'replay_logs';logs.mkdir(exist_ok=True)
 receipt=dict(schema='nonagon.part_x.replay.v1',python=platform.python_version(),platform=platform.platform(),
  dependencies={p:importlib.metadata.version(p) for p in ['numpy','scipy','sympy','mpmath','python-flint']},campaigns=[],verdict='RUNNING')
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 for name in JOBS:
  start=time.perf_counter()
  result=subprocess.run([sys.executable,str(HERE/(name+'.py'))],cwd=HERE,env=env,capture_output=True,text=True)
  (logs/(name+'.log')).write_text(result.stdout+result.stderr)
  row=dict(program=name+'.py',returncode=result.returncode,elapsed_seconds=time.perf_counter()-start)
  receipt['campaigns'].append(row)
  require(result.returncode==0,name+' failed; see replay_logs/'+name+'.log')
  path=HERE/(name+'_results.json');d=json.loads(path.read_text());validate(name,d)
  row.update(verdict='PASS',output=path.name,output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256((HERE/(name+'.py')).read_bytes()).hexdigest())
  print(name+': PASS',flush=True)
 receipt['verdict']='PASS';(HERE/'replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('All fifteen Part X numerical campaigns passed.')
if __name__=='__main__':main()
