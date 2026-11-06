"""Select interpretable gate on repaired prior development only, then fresh test."""
import argparse,csv,datetime,hashlib,json,os,sys,time,warnings
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import gate
from gate import adaptive,core
ROOT=Path(__file__).resolve().parent
ADAPTIVE=ROOT.parent/'adaptive'

def h(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def score(row):return min(float(row['heldout_nrmse']),1e6) if row['heldout_nrmse'] is not None else 1e6

def select():
    rows=json.loads((ADAPTIVE/'development_raw.json').read_text());arrays=np.load(ADAPTIVE/'development_observations.npz')
    cases={}
    for row in rows:
        cid=row['case_id']
        if cid not in cases:
            cluster=row['cluster'];y=arrays[cluster+'_candidate_clean']*(1+row['noise']*arrays[cluster+'_standardized_noise'])
            features=gate.pilot_features(adaptive.POOL_W[adaptive.PILOT],y[adaptive.PILOT],row['order'],adaptive.BAND)
            cases[cid]=dict(case_id=cid,cluster=cluster,noise=row['noise'],order=row['order'],features=features,
                            pilot_observation_sha256=core.digest(adaptive.POOL_W[adaptive.PILOT],y[adaptive.PILOT]),outcomes={})
        if row['strategy'] in ['log_uniform','adaptive_disagreement']:cases[cid]['outcomes'][row['strategy']]=score(row)
    (ROOT/'development_features.json').write_text(json.dumps(list(cases.values()),indent=2))
    noisy=[r for r in cases.values() if r['noise']>0]
    clusters=sorted({r['cluster'] for r in noisy});rng=np.random.default_rng(391173)
    samples=rng.integers(0,len(clusters),(5000,len(clusters)))
    results=[]
    for rule in gate.candidate_rules():
        selected=[r for r in noisy if gate.matches(r['features'],rule)]
        selected_clusters=len({r['cluster'] for r in selected})
        base=[];policy=[]
        for cluster in clusters:
            rr=[r for r in noisy if r['cluster']==cluster]
            base.append(np.mean([min(r['outcomes']['log_uniform'],1.) for r in rr]))
            policy.append(np.mean([min(r['outcomes']['adaptive_disagreement' if gate.matches(r['features'],rule) else 'log_uniform'],1.) for r in rr]))
        base=np.array(base);policy=np.array(policy)
        ratio=float(np.mean(policy)/np.mean(base));boot=policy[samples].mean(axis=1)/base[samples].mean(axis=1)
        ci=np.quantile(boot,[.025,.975])
        eligible=rule['kind']=='always_log' or (selected_clusters>=4 and ratio<=.95 and ci[1]<=1.)
        results.append(dict(rule=rule,selected_positive_cases=len(selected),selected_system_clusters=selected_clusters,
                            capped_mean_error_ratio=ratio,dev_ratio_ci95=ci.tolist(),eligible=bool(eligible)))
    eligible=[r for r in results if r['eligible']]
    chosen=min(eligible,key=lambda r:r['capped_mean_error_ratio'])
    record=dict(rule=chosen['rule'],development_selection=chosen,
                selection='Small fixed rule set, >=4 selected systems, >=5% capped-mean risk reduction and bootstrap upper<=1 on prior DEVELOPMENT only; otherwise always_log. Development uncertainty is a selection heuristic, not confirmation.')
    (ROOT/'all_rules_development_scores.json').write_text(json.dumps(results,indent=2));(ROOT/'selected_rule.json').write_text(json.dumps(record,indent=2))
    protocol=dict(title='Observed-pilot adaptive acquisition gate',version=1,confirmation=dict(seed=69377123,systems_per_scenario=8),
                  primary='Mean clean heldout NRMSE ratio gated / always_log over positive-noise trials; paired whole-system bootstrap',
                  secondary='Geometric mean error, scenario cells, tail and runtime; all three strategies identical64uniquequery budget with33commonpilot',
                  selection=record['selection'],no_confirmation_for_selection=True,freeze='Local hash freeze before new policy confirmation, not independent preregistration',
                  features='Observed noisy pilot only; estimated pole gap/residue ratio/damping/training residual/model disagreement; no hidden scenario label, true poles or clean response',
                  source='Repaired adaptive development seed61739311, no original invalid auto-solver data',query_budget=64,common_pilot=33)
    (ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2))
    frozen=dict(frozen_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256={name:h(ROOT/name) for name in ['gate.py','experiment.py','protocol.json','selected_rule.json','all_rules_development_scores.json']},
                adaptive_sha256={name:h(ADAPTIVE/name) for name in ['acquisition.py','reference_core.py','protocol.json','development_raw.json']})
    (ROOT/'freeze.json').write_text(json.dumps(frozen,indent=2));print(json.dumps(record,indent=2))

def verify():
    frozen=json.loads((ROOT/'freeze.json').read_text())
    for name,wanted in frozen['sha256'].items():
        if h(ROOT/name)!=wanted:raise ValueError('policy file changed: '+name)
    for name,wanted in frozen['adaptive_sha256'].items():
        if h(ADAPTIVE/name)!=wanted:raise ValueError('adaptive dependency changed: '+name)

def confirmation():
    verify();spec=json.loads((ROOT/'protocol.json').read_text())['confirmation'];rule=json.loads((ROOT/'selected_rule.json').read_text())['rule'];rows=[];allarrays={};tstart=time.time()
    for scenario in adaptive.SCENARIOS:
      for i in range(spec['systems_per_scenario']):
        case=adaptive.system(scenario,i,spec['seed']);cluster=f'{scenario}-{i:03d}'
        allarrays[cluster+'_candidate_clean']=case['clean'];allarrays[cluster+'_heldout_clean']=case['heldout'];allarrays[cluster+'_standardized_noise']=case['noise']
        for eta in adaptive.NOISES:
          y=case['clean']*(1+eta*case['noise']);cid=f'{cluster}-eta{eta:g}'
          ft=time.perf_counter();selected_strategy,features=gate.choose_strategy(adaptive.POOL_W[adaptive.PILOT],y[adaptive.PILOT],case['order'],adaptive.BAND,rule);feature_seconds=time.perf_counter()-ft
          for strategy in ['always_log','always_adaptive','gated']:
            actual='log_uniform' if strategy=='always_log' else ('adaptive_disagreement' if strategy=='always_adaptive' else selected_strategy)
            requested=[]
            def query(indices):
                if len(set(requested+indices))!=len(requested)+len(indices):raise ValueError('duplicatequery')
                requested.extend(indices)
                if len(requested)>64:raise ValueError('budgetexceeded')
                return y[np.array(indices)].copy()
            row=dict(case_id=cid,cluster=cluster,scenario=scenario,noise=eta,strategy=strategy,actual_acquisition=actual,status='ok',features=features,
                     pool_observation_sha256=core.digest(adaptive.POOL_W,y),physical_solver_audit=case['solve_audit'])
            t0=time.perf_counter()
            with warnings.catch_warnings(record=True) as ws:
              warnings.simplefilter('always')
              try:
                model,indices,values,metadata=adaptive.acquire(actual,query,case['order'],spec['seed']+i*13+adaptive.SCENARIOS.index(scenario)*1001)
                row['seconds']=time.perf_counter()-t0+(feature_seconds if strategy=='gated' else 0)
                pred=model['predict'](1j*adaptive.TEST_W)
                if not np.all(np.isfinite(pred)):raise ValueError('nonfiniteprediction')
                row.update(heldout_nrmse=float(np.linalg.norm(pred-case['heldout'])/np.linalg.norm(case['heldout'])),observations=len(requested),selected_indices=indices,
                           metadata=metadata,selected_observation_sha256=core.digest(np.array(indices),values))
              except Exception as exc:row.update(status='failed',exception=f'{type(exc).__name__}: {exc}',heldout_nrmse=None,observations=len(requested),seconds=time.perf_counter()-t0)
              row['warnings']=[f'{w.category.__name__}: {w.message}' for w in ws]
            rows.append(row)
        print(cluster,len(rows),'fits',round(time.time()-tstart,1),'seconds',flush=True)
    (ROOT/'confirmation_raw.json').write_text(json.dumps(rows,indent=2))
    with (ROOT/'confirmation_raw.csv').open('w',newline='') as f:
        fields=sorted(set().union(*(r.keys() for r in rows)));w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()})
    allarrays['candidate_omega']=adaptive.POOL_W;allarrays['heldout_omega']=adaptive.TEST_W;np.savez_compressed(ROOT/'confirmation_observations.npz',**allarrays)
    summaries=[];comparisons=[];rng=np.random.default_rng(718931)
    for sc in ['ALL']+adaptive.SCENARIOS:
      for strategy in ['always_log','always_adaptive','gated']:
        rr=[r for r in rows if r['noise']>0 and r['strategy']==strategy and (sc=='ALL' or r['scenario']==sc)];e=np.array([score(r) for r in rr])
        summaries.append(dict(scenario=sc,strategy=strategy,trials=len(rr),failures=sum(r['status']!='ok' for r in rr),mean_nrmse=float(e.mean()),geometric_nrmse=float(10**np.log10(np.maximum(e,1e-14)).mean()),
                              median_nrmse=float(np.median(e)),p90_nrmse=float(np.quantile(e,.9)),median_seconds=float(np.median([r['seconds'] for r in rr])),adaptive_selected=sum(r['actual_acquisition']=='adaptive_disagreement' for r in rr)))
        if strategy=='always_log':continue
        aa={r['case_id']:r for r in rr};bb={r['case_id']:r for r in rows if r['strategy']=='always_log' and r['noise']>0 and (sc=='ALL' or r['scenario']==sc)}
        clusters=sorted({r['cluster'] for r in rr});a=[];b=[]
        for cluster in clusters:
            cids=[cid for cid in aa if aa[cid]['cluster']==cluster]
            a.append(np.mean([score(aa[cid]) for cid in cids]));b.append(np.mean([score(bb[cid]) for cid in cids]))
        a=np.array(a);b=np.array(b);resamples=rng.integers(0,len(a),(5000,len(a)));boot=a[resamples].mean(axis=1)/b[resamples].mean(axis=1)
        comparisons.append(dict(scenario=sc,strategy=strategy,comparator='always_log',mean_error_ratio=float(a.mean()/b.mean()),ratio_ci95=np.quantile(boot,[.025,.975]).tolist(),clusters=len(a),trials=len(rr)))
    result=dict(selected_rule=rule,summary=summaries,comparisons=comparisons,feature_selection_used_development_only=True,all_inputs_hash_frozen=True)
    (ROOT/'results.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['select','confirmation']);args=p.parse_args()
    select() if args.action=='select' else confirmation()
