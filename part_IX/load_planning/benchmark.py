#!/usr/bin/env python3
"""Frozen prospective comparison of fixed and pilot-selected finite loads.

Both arms spend the identical pilot plus final canonical acquisition budget.
Only externally bounded intervals enter the planner and production inverse.
All circuits and instrument readings are synthetic; no apparatus is connected.
"""
from __future__ import annotations

import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from random import Random
import sys

HERE=Path(__file__).resolve().parent
CORE=HERE.parent
sys.path.insert(0,str(CORE))
from bounded_intervention import Interval,infer_bounds,canonical_plan,load_key
from planner import recommend,calibration_radius


def require(ok,message):
    if not ok:raise RuntimeError(message)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def encoded(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()


def solve(a,b):
    """Exact forward elimination; independent of reconstruction."""
    n=len(a);w=[[F(v) for v in row]+[F(b[i])] for i,row in enumerate(a)]
    for j in range(n):
        p=next(i for i in range(j,n) if w[i][j])
        w[j],w[p]=w[p],w[j]
        d=w[j][j];w[j]=[v/d for v in w[j]]
        for i in range(j+1,n):
            d=w[i][j];w[i]=[x-d*y for x,y in zip(w[i],w[j])]
    x=[F(0)]*n
    for i in range(n-1,-1,-1):x[i]=w[i][-1]-sum(w[i][j]*x[j] for j in range(i+1,n))
    return x


def response(k,loads):
    a=[row[:] for row in k]
    for j,t in loads.items():a[j][j]+=t
    return solve(a,[F(i==0) for i in range(len(k))])[0]


def fixture(n,family,scale,seed):
    rng=Random(seed);k=[[F(0) for j in range(n)] for i in range(n)]
    for i in range(n):k[i][i]=scale*F(rng.randint(2,4),4)
    edges=[]
    for i in range(n):
        for j in range(i+1,n):
            if family=='complete' or j==i+1:
                w=scale*F(rng.randint(2,8),4)
                k[i][i]+=w;k[j][j]+=w;k[i][j]-=w;k[j][i]-=w
                edges.append([i,j,str(w)])
    return {'id':f'n{n}_{family}_s{str(scale).replace("/","_")}',
            'n':n,'family':family,'scale':str(scale),'seed':seed,
            'K':[[str(v) for v in row] for row in k],'edges':edges}


def freeze():
    path=HERE/'protocol.json'
    if path.exists():raise FileExistsError('Protocol already frozen; refuse replacement')
    fixtures=[]
    for n in (3,4,6):
        for fi,family in enumerate(('path','complete')):
            for si,scale in enumerate(map(F,('1/16','1/4','1','4','16'))):
                fixtures.append(fixture(n,family,scale,202610312000+n*100+fi*10+si))
    hardware={'conductance_scale_S':'1/1000','current_A':'1/10000',
              'minimum_shunt_resistance_ohm':'20','maximum_shunt_resistance_ohm':'20000',
              'simultaneous_shunt_channels':2,'candidate_base_loads':['1/16','1/4','1','4','16']}
    budgets=[
        dict(id='ideal_reference',voltage_error_V='0',drift_error_V='0',
             current_relative_error='0',shunt_relative_error='0',relay_series_max_ohm='0'),
        dict(id='precision_calibration',voltage_error_V='1/100000000',drift_error_V='1/1000000000',
             current_relative_error='1/1000000',shunt_relative_error='1/1000000',relay_series_max_ohm='1/200'),
        dict(id='readout_limited',voltage_error_V='1/100000',drift_error_V='1/1000000',
             current_relative_error='1/1000000',shunt_relative_error='1/1000000',relay_series_max_ohm='1/200')]
    protocol={'schema':'nonagon.finite-load-benchmark-protocol.v2',
        'frozen_before_measurement_generation':True,'fixtures':fixtures,'budgets':budgets,
        'hardware':hardware,'pilot_base_load':'1','minimum_positive_edge_gap':'1/32',
        'prospective_policy':'planner uses common pilot scalar intervals and external calibration only; final readings generated after strength selection',
        'arms':['fixed_pilot_strength','pilot_selected_strength'],
        'fair_budget':'both arms acquire identical 2n−1 pilot readings and a complete n(n+1)/2 final plan; only repeated baseline intervals are intersected; for a repeated loaded setting use the narrower individually valid interval',
        'error_model':'independent bounded voltage, source-gain, shunt-calibration and relay errors; declared drift included in every scalar interval',
        'truth_seed':'202610312000 + 100*n + 10*family_index + scale_index',
        'reading_seed':'202610309900 + 100*fixture_index + budget_index',
        'hardware_status':'synthetic calibration benchmark; no physical apparatus connected',
        'primary_metrics':['full_matrix_enclosed','full_graph_disclosed','truth_containment','unresolved'],
        'fixture_sha256':hashlib.sha256(encoded(fixtures)).hexdigest(),
        'protocol_note':'An initial development protocol was superseded before this prospective campaign: independently calibrated repeated loaded readings need not enclose the same physical response and cannot be intersected. Corrected policy uses new fixture seeds. Initial protocol is retained for provenance.'}
    path.write_text(json.dumps(protocol,indent=2)+'\n')
    print(json.dumps({'protocol':str(path),'fixture_count':len(fixtures),'sha256':digest(path)}))


def make_reading(k,nominal,hw,rng,index,count):
    scale=F(hw['conductance_scale_S']);current=F(hw['current_A'])
    erel=F(hw['shunt_relative_error']);relaymax=F(hw['relay_series_max_ohm'])
    loads={};branches=[]
    for j,t in sorted(nominal.items()):
        delta=erel*F(rng.randint(-4,4),4)
        relay=relaymax*F(rng.randint(0,4),4)
        shunt=t*scale*(1+delta)
        actual=shunt/(1+relay*shunt)/scale
        loads[j]=actual
        branches.append({'site':j,'nominal':str(t),'actual':str(actual),
                         'relative_shunt_deviation':str(delta),'relay_ohm':str(relay)})
    gain=F(hw['current_relative_error'])*F(rng.randint(-4,4),4)
    actual_current=current*(1+gain)
    h=response(k,loads);v=actual_current*h/scale
    err=F(hw['voltage_error_V'])*F(rng.randint(-4,4),4)
    drift=F(hw['drift_error_V'])*F(2*index-count,max(1,count))
    measured=v+err+drift
    ve=F(hw['voltage_error_V'])+F(hw['drift_error_V'])
    crel=F(hw['current_relative_error'])
    lo=(measured-ve)*scale/(current*(1+crel))
    hi=(measured+ve)*scale/(current*(1-crel))
    require(lo>0,'Declared voltage envelope cannot certify positive readout')
    box=Interval(lo,hi)
    require(box.contains(h),'Voltage/current interval lost exact truth')
    record={'nominal_loads':{str(j):str(t) for j,t in nominal.items()},'branches':branches,
            'actual_current_A':str(actual_current),'true_response':str(h),
            'measured_voltage_V':str(measured),'readout_error_V':str(err),'drift_V':str(drift),
            'normalized_response_interval':box.json()}
    return record,box


def intersect(a,b):
    return Interval(max(a.lo,b.lo),min(a.hi,b.hi))


def score(k,bounds,gap):
    if bounds['status']!='enclosed':
        return {'matrix_enclosed':False,'full_graph':False,'disclosed_edges':0,
                'truth_checks':0,'max_K_interval_width':None}
    n=len(k);checks=0
    cols=[solve(k,[F(i==j) for i in range(n)]) for j in range(n)]
    g=[[cols[j][i] for j in range(n)] for i in range(n)]
    for key,truth in [('K_intervals',k),('G_intervals',g)]:
        for i,row in enumerate(truth):
            for j,v in enumerate(row):
                lo,hi=map(F,bounds[key][i][j]);require(lo<=v<=hi,'True matrix escaped interval')
                checks+=1
    disclosed=0
    for edge in bounds['edges']:
        i,j=edge['sites_zero_based'];w=-k[i][j];lo,hi=map(F,edge['weight_interval'])
        require(lo<=w<=hi,'True edge escaped interval');checks+=1
        status=edge['status']
        require(status!='inconsistent_with_gap','True fixture rejected gap class')
        require(status!='present' or w>0,'False positive edge')
        require(status!='absent_under_gap' or w==0,'False absence')
        disclosed+=status in ('present','absent_under_gap')
    return {'matrix_enclosed':True,'full_graph':disclosed==n*(n-1)//2,
            'disclosed_edges':disclosed,'truth_checks':checks,
            'max_K_interval_width':str(max(F(hi)-F(lo) for row in bounds['K_intervals'] for lo,hi in row))}


def run(output_dir):
    output_dir.mkdir(parents=True,exist_ok=True)
    protocol_path=HERE/'protocol.json';p=json.loads(protocol_path.read_text())
    require(hashlib.sha256(encoded(p['fixtures'])).hexdigest()==p['fixture_sha256'],'Frozen fixture hash differs')
    # Exclude exact matrix duplication with previously frozen engineering and
    # production fixtures, even when the generation seed differs.
    old=set()
    for oldpath in (CORE/'engineering_tests/instrument_fixtures').glob('*.json'):
        doc=json.loads(oldpath.read_text());old.add(hashlib.sha256(encoded(doc['K_normalized'])).hexdigest())
    for oldpath in (CORE/'data/cases').glob('*.certificate.json'):
        doc=json.loads(oldpath.read_text());old.add(hashlib.sha256(encoded(doc['recovered']['K'])).hexdigest())
    for f in p['fixtures']:
        require(hashlib.sha256(encoded(f['K'])).hexdigest() not in old,'Fresh fixture duplicates frozen matrix')
    rows=[];summary=defaultdict(Counter);selected=Counter();raw_count=0;unique_raw_count=0;truth_checks=0
    data_path=output_dir/'benchmark_cases.jsonl';example_written=False
    with data_path.open('w') as out:
        for fi,f in enumerate(p['fixtures']):
            n=f['n'];k=[[F(v) for v in row] for row in f['K']]
            for bi,budget in enumerate(p['budgets']):
                hw=dict(p['hardware'],**{key:value for key,value in budget.items() if key!='id'})
                t0=F(p['pilot_base_load']);pilotplan=canonical_plan(n,t0,2*t0)[:2*n-1]
                rng=Random(202610309900+100*fi+bi)
                pilot=[];pilotobs={}
                for ix,loads in enumerate(pilotplan):
                    record,box=make_reading(k,loads,hw,rng,ix,len(pilotplan))
                    pilot.append(record);pilotobs[load_key(loads)]=box
                pilotdoc={'schema':'nonagon.finite-load-pilot.v1','network_order':n,
                          'input_and_readout_site_zero_based':0,'pilot_base_load':str(t0),
                          'hardware':hw,'readings':[
                              {'nominal_loads':r['nominal_loads'],
                               'normalized_response_interval':r['normalized_response_interval']} for r in pilot]}
                choice=recommend(pilotdoc);chosen=F(choice['selected_base_load'])
                selected[budget['id'],str(chosen),choice['status']]+=1
                pilot_load_error=calibration_radius(hw,t0)
                for rec in pilot:
                    for branch in rec['branches']:
                        nominal=F(branch['nominal']);actual=F(branch['actual'])
                        require(nominal*(1-pilot_load_error)<=actual<=nominal*(1+pilot_load_error),
                                'Pilot shunt escaped declared calibration')
                        truth_checks+=1
                # Ground truth is consulted only after the recommendation has
                # been fixed, to test the pilot's conditional parameter boxes.
                g0=solve(k,[F(i==0) for i in range(n)])
                for single in choice['pilot_single_enclosures']:
                    j=single['site_zero_based'];gj=solve(k,[F(i==j) for i in range(n)])
                    qlo,qhi=map(F,single['q_interval']);blo,bhi=map(F,single['b_interval'])
                    require(qlo<=g0[j]*g0[j]<=qhi,'Pilot q enclosure lost truth')
                    require(blo<=gj[j]<=bhi,'Pilot b enclosure lost truth')
                    truth_checks+=2
                unique_raw_count+=len(pilot);truth_checks+=len(pilot)
                if not example_written and budget['id']=='precision_calibration':
                    (output_dir/'example_pilot.json').write_text(json.dumps(pilotdoc,indent=2)+'\n')
                    (output_dir/'example_recommendation.json').write_text(json.dumps(choice,indent=2)+'\n')
                    example_written=True
                arms=[]
                for arm,t in [('fixed_pilot_strength',t0),('pilot_selected_strength',chosen)]:
                    obs={};records=[];plan=canonical_plan(n,t,2*t)
                    # Equal prospective random error sequence for the same
                    # scalar ordinal; actual signals depend on selected loads.
                    final_rng=Random(202610310900+100*fi+bi)
                    for ix,loads in enumerate(plan):
                        rec,box=make_reading(k,loads,hw,final_rng,ix,len(plan))
                        records.append(rec);key=load_key(loads)
                        if key==():
                            obs[key]=intersect(pilotobs[key],box)
                        elif key in pilotobs and pilotobs[key].radius<box.radius:
                            # These individually valid intervals refer to
                            # different independently calibrated shunts. Select
                            # either full reading; do not intersect them.
                            obs[key]=pilotobs[key]
                        else:obs[key]=box
                    e=calibration_radius(hw,t)
                    bounds=infer_bounds(n,obs,t,2*t,gap=F(p['minimum_positive_edge_gap']),
                                        relative_load_error=e)
                    scored=score(k,bounds,F(p['minimum_positive_edge_gap']))
                    loadchecks=0
                    for rec in records:
                        for branch in rec['branches']:
                            nominal=F(branch['nominal']);actual=F(branch['actual'])
                            require(nominal*(1-e)<=actual<=nominal*(1+e),'True shunt escaped load calibration')
                            loadchecks+=1
                    raw_count+=len(pilot)+len(records)
                    unique_raw_count+=len(records)
                    truth_checks+=scored['truth_checks']+loadchecks+len(records)
                    row={'fixture':f['id'],'budget':budget['id'],'arm':arm,
                         'base_load':str(t),'status':bounds['status'],'reason':bounds.get('reason'),
                         'pilot_status':choice['status'],'pilot_scalar_count':len(pilot),
                         'final_scalar_count':len(records),'total_scalar_count':len(pilot)+len(records),
                         'relative_load_envelope':str(e),**scored}
                    rows.append(row);arms.append({'case':row,'readings':records,'bounds':bounds})
                    s=summary[budget['id'],arm];s['cases']+=1;s['matrix_enclosed']+=scored['matrix_enclosed']
                    s['full_graph']+=scored['full_graph'];s['disclosed_edges']+=scored['disclosed_edges']
                    s['possible_edges']+=n*(n-1)//2;s['scalar_readings']+=row['total_scalar_count']
                require(arms[0]['case']['total_scalar_count']==arms[1]['case']['total_scalar_count'],'Unequal budget')
                out.write(json.dumps({'fixture':f,'budget':budget,'hardware':hw,
                                      'pilot_readings':pilot,'pilot_document':pilotdoc,
                                      'recommendation':choice,'arms':arms},separators=(',',':'))+'\n')
    pairs=[]
    for budget in p['budgets']:
        outcomes=Counter()
        for f in p['fixtures']:
            arms=[r for r in rows if r['fixture']==f['id'] and r['budget']==budget['id']]
            fixed,planned=arms
            outcomes['planned_only_matrix_enclosed']+=planned['matrix_enclosed'] and not fixed['matrix_enclosed']
            outcomes['fixed_only_matrix_enclosed']+=fixed['matrix_enclosed'] and not planned['matrix_enclosed']
            outcomes['planned_only_full_graph']+=planned['full_graph'] and not fixed['full_graph']
            outcomes['fixed_only_full_graph']+=fixed['full_graph'] and not planned['full_graph']
            if fixed['matrix_enclosed'] and planned['matrix_enclosed']:
                wf=F(fixed['max_K_interval_width']);wp=F(planned['max_K_interval_width'])
                outcomes['paired_enclosures']+=1
                outcomes['planned_width_smaller']+=wp<wf;outcomes['fixed_width_smaller']+=wf<wp
                outcomes['equal_width']+=wf==wp
        pairs.append(dict(budget=budget['id'],**dict(outcomes)))
    receipt={'schema':'nonagon.finite-load-benchmark-results.v1','status':'PASSED',
        'hardware_status':p['hardware_status'],'protocol_sha256':digest(protocol_path),
        'fixture_count':len(p['fixtures']),'fresh_fixture_matrix_disjointness_checks':len(p['fixtures']),
        'prospective_cases':len(p['fixtures'])*len(p['budgets']),
        'arm_cases':len(rows),'scalar_readings_counting_shared_pilot_per_arm':raw_count,
        'unique_generated_scalar_readings':unique_raw_count,
        'truth_containment_checks':truth_checks,'false_edge_calls':0,'escaped_truth':0,
        'summary':[dict(budget=b,arm=a,**dict(s)) for (b,a),s in sorted(summary.items())],
        'selection_counts':[{'budget':b,'base_load':t,'status':status,'count':count}
                            for (b,t,status),count in sorted(selected.items())],
        'paired_comparisons':pairs,
        'source_sha256':{name:digest(HERE/name) for name in ('planner.py','benchmark.py')},
        'production_interval_source_sha256':digest(CORE/'bounded_intervention.py'),
        'records_sha256':digest(data_path),
        'limits':['First-order planner score is a diagnostic; all matrix/graph claims use the unchanged nonlinear interval gate.',
                  'A single globally shared strength preserves the production fixed-t,2t,pair-t schema.',
                  'No averaging or square-root-of-repeat-count error reduction is used for bounded errors.',
                  'Full graph absence requires the explicitly declared minimum positive edge gap1/32.',
                  'Synthetic data and external calibration envelopes do not constitute hardware qualification.']}
    (output_dir/'benchmark_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    return receipt


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--freeze',action='store_true')
    ap.add_argument('--output-dir',type=Path,default=HERE/'results')
    args=ap.parse_args()
    if args.freeze:freeze()
    else:run(args.output_dir)
