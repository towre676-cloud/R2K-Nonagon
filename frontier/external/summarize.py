#!/usr/bin/env python3
"""All-device paired summaries; no model fitting or test-dependent selection."""
import collections,csv,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def run():
    protocol=json.loads((ROOT/'protocol.json').read_text())
    receipt=json.loads((ROOT/'confirmation_results.json').read_text())
    records=receipt['results'];devices=protocol['confirmation_devices'];rng=np.random.default_rng(protocol['bootstrap_seed'])
    cap=protocol['mandatory_failure_sensitivity']['selected_failure_error'];floor=protocol['mandatory_failure_sensitivity']['successful_error_floor']
    indexed={(r['device'],r['scenario'],r['budget'],r['method']):r for r in records}
    def err(r):return float(np.clip(r['test_error'],floor,cap)) if r['status']=='ok' else cap
    stats=[];pairs=[]
    for case in protocol['cases']:
        scenario,budget=case['scenario'],case['budget']
        for method in protocol['methods']:
            rows=[indexed[(d,scenario,budget,method)] for d in devices]
            vals=np.array([err(r) for r in rows])
            stats.append({'scenario':scenario,'budget':budget,'method':method,'device_count':len(vals),'failures':sum(r['status']!='ok' for r in rows),'median_capped_error':float(np.median(vals)),'mean_capped_error':float(np.mean(vals)),'p90_capped_error':float(np.quantile(vals,.9)),'max_capped_error':float(vals.max())})
        for baseline in ['loewner','vector_fitting','cubic_log_frequency']:
            a=np.array([err(indexed[(d,scenario,budget,'positive_real_relaxation')]) for d in devices]);b=np.array([err(indexed[(d,scenario,budget,baseline)]) for d in devices]);draw=rng.integers(0,len(devices),size=(10000,len(devices)))
            logs=np.log(a/b);difference=a-b
            ratios=np.exp(logs[draw].mean(axis=1));diffs=difference[draw].mean(axis=1)
            common=[d for d in devices if indexed[(d,scenario,budget,'positive_real_relaxation')]['status']=='ok' and indexed[(d,scenario,budget,baseline)]['status']=='ok']
            commonlog=[np.log(indexed[(d,scenario,budget,'positive_real_relaxation')]['test_error']/indexed[(d,scenario,budget,baseline)]['test_error']) for d in common]
            pairs.append({'scenario':scenario,'budget':budget,'numerator':'positive_real_relaxation','denominator':baseline,'n_devices':len(devices),'analysis':'mandatory all-device capped-failure sensitivity','geometric_error_ratio':float(np.exp(logs.mean())),'ratio_bootstrap_95_unadjusted':np.quantile(ratios,[.025,.975]).tolist(),'mean_error_difference':float(difference.mean()),'difference_bootstrap_95_unadjusted':np.quantile(diffs,[.025,.975]).tolist(),'nnls_win_count':int(np.sum(a<b)),'baseline_win_count':int(np.sum(a>b)),'common_success_devices':len(common),'common_success_geometric_ratio':float(np.exp(np.mean(commonlog))) if common else None})
    overall=[]
    for baseline in ['loewner','vector_fitting','cubic_log_frequency']:
        device_log=[]
        for d in devices:
            device_log.append(np.mean([np.log(err(indexed[(d,c['scenario'],c['budget'],'positive_real_relaxation')])/err(indexed[(d,c['scenario'],c['budget'],baseline)])) for c in protocol['cases']]))
        v=np.array(device_log);draw=rng.integers(0,len(devices),size=(10000,len(devices)));r=np.exp(v[draw].mean(axis=1))
        overall.append({'numerator':'positive_real_relaxation','denominator':baseline,'n_devices':len(devices),'tasks_per_device':len(protocol['cases']),'geometric_error_ratio':float(np.exp(v.mean())),'bootstrap_95_unadjusted':np.quantile(r,[.025,.975]).tolist(),'label':'Descriptive aggregate across all4 fixed tasks, clustered by device; no task cherry-picking.'})
    candidates=json.loads((ROOT/'confirmation_candidates.json').read_text())
    warning_counts=collections.Counter(w for r in candidates for w in r.get('warnings',[]))
    result={'scope':'Protected ID07–54confirmation; onebatch/oneSoC/oneinstrument; separatecells, not independentlaboratories. Intervals unadjusted for multiplecontrasts and conditional on shared measurement calibration.','protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),'confirmation_receipt_sha256':hashlib.sha256((ROOT/'confirmation_results.json').read_bytes()).hexdigest(),'device_count':len(devices),'case_count':receipt['case_count'],'selected_result_count':len(records),'candidate_count':receipt['candidate_count'],'candidate_failures':receipt['candidate_failure_count'],'selected_failures':receipt['selected_failure_count'],'per_task_method':stats,'paired_comparisons':pairs,'overall_device_clustered_descriptive':overall,'candidate_warning_counts':dict(warning_counts),'limitations':['Prediction targets stored measuredpoints, not noiseless physicaltruth.','NNLS is existing positive-real relaxationprior, not proof of uniqueNonagon advantage.','The model can be certified passive by construction; the physicalbattery outside measuredfrequencies cannot.','No tuning or hyperparametergrid edits fromconfirmationoutcomes.']}
    (ROOT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (ROOT/'scores.csv').open('w',newline='') as f:
        fields=['device','scenario','budget','method','status','test_error','validation_error'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r.get(k) for k in fields} for r in records)
    for r in stats:print(r['scenario'],r['budget'],r['method'],r['median_capped_error'])
    for r in overall:print('OVERALL',r['denominator'],r['geometric_error_ratio'],r['bootstrap_95_unadjusted'])
if __name__=='__main__':run()
