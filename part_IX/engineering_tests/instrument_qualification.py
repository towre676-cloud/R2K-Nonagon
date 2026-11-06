"""Separate frozen qualification-budget extension to instrument_campaign.py."""
from pathlib import Path
from fractions import Fraction as F
from random import Random
from collections import Counter,defaultdict
import argparse
import hashlib
import json
import time

import instrument_campaign as core
from bounded_intervention import Interval,infer_bounds,canonical_plan,load_key

HERE=Path(__file__).resolve().parent


def freeze():
    doc={'schema':'nonagon.instrument-qualification-protocol.v1',
         'description':'Declared combined qualification scenario; not a complete manufacturer accuracy specification',
         'fixture_protocol':'instrument_protocol.json','fixture_count':24,'base_loads':['1/4','1','4'],
         'conductance_scale_S':'1/1000','current_A':'1/1000',
         'voltage_absolute_error_V':'1/200000','voltage_relative_to_reading_error':'3/25000',
         'current_relative_error':'1/5000','shunt_relative_error':'1/1000',
         'relay_series_max_ohm':'1/20','drift_error_V':'0',
         'readout_error_rule':'abs(Vmeasured-Vtrue)<=120ppm*abs(Vmeasured)+5microvolt',
         'validation':'three reserved max-two-shunt settings; repeats after five nonbaseline readings',
         'case_count':72,'frozen_before_measurement_generation':True,
         'fixture_sha256':json.loads((HERE/'instrument_protocol.json').read_text())['fixture_sha256']}
    p=HERE/'instrument_qualification_protocol.json'
    if p.exists():raise FileExistsError('qualification protocol already frozen')
    p.write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps({'protocol':str(p),'cases':72}))


def read(k,loads,p,rng,index,count):
    budget={'current_A':p['current_A'],'voltage_error_V':p['voltage_absolute_error_V'],
            'current_relative_error':p['current_relative_error'],'shunt_relative_error':p['shunt_relative_error'],
            'relay_series_max_ohm':p['relay_series_max_ohm'],'drift_error_V':p['drift_error_V']}
    record,_=core.make_reading(k,loads,budget,F(p['conductance_scale_S']),rng,index,count)
    gain=F(p['voltage_relative_to_reading_error'])*F(rng.randint(-3,3),4)
    gain_error=gain*F(record['ideal_voltage_V'])
    measured=F(record['measured_voltage_V'])+gain_error
    envelope=F(p['voltage_absolute_error_V'])+F(p['voltage_relative_to_reading_error'])*abs(measured)
    if abs(measured-F(record['ideal_voltage_V']))>envelope:
        raise RuntimeError('simulated readout escaped declared combined envelope')
    scale=F(p['conductance_scale_S']);i=F(p['current_A']);rho=F(p['current_relative_error'])
    box=Interval((measured-envelope)*scale/(i*(1+rho)),
                 (measured+envelope)*scale/(i*(1-rho)))
    if not box.contains(F(record['true_normalized_response'])):raise RuntimeError('true response escaped normalized interval')
    record.update({'measured_voltage_V':str(measured),'additive_gain_error_V':str(gain_error),
                   'gain_relative_deviation':str(gain),'declared_voltage_error_V':str(envelope),
                   'normalized_interval':box.json()})
    return record,box


def run():
    start=time.monotonic();pp=HERE/'instrument_qualification_protocol.json';p=json.loads(pp.read_text())
    fs=sorted([json.loads(x.read_text()) for x in (HERE/'instrument_fixtures').glob('*.json')],key=lambda x:(x['order'],x['family']))
    if hashlib.sha256(json.dumps(fs,sort_keys=True).encode()).hexdigest()!=p['fixture_sha256']:raise RuntimeError('fixture identity mismatch')
    scale=F(p['conductance_scale_S']);rows=[];sums=defaultdict(Counter);normal_checks=0;matrix_checks=0
    packet=HERE/'instrument_qualification_records.jsonl'
    with packet.open('w') as out:
        for f in fs:
            n=f['order'];k=[[F(x) for x in row] for row in f['K_normalized']]
            for ti,ts in enumerate(p['base_loads']):
                t=F(ts);plan=canonical_plan(n,t,2*t);rng=Random(2026103200+1000*n+100*f['family']+ti)
                readings=[];repeats=[];obs={}
                for j,loads in enumerate(plan):
                    rr,box=read(k,loads,p,rng,j,len(plan));readings.append(rr);obs[load_key(loads)]=box;normal_checks+=1
                    if j>0 and (j%5==0 or j==len(plan)-1):
                        rr,box=read(k,{},p,rng,j,len(plan));repeats.append(rr);normal_checks+=1
                        old=obs[()];obs[()]=Interval(max(old.lo,box.lo),min(old.hi,box.hi))
                erel=core.load_bound({'shunt_relative_error':p['shunt_relative_error'],'relay_series_max_ohm':p['relay_series_max_ohm']},scale,2*t)
                for rr in readings:
                    for branch in rr['branches']:
                        if abs(F(branch['actual_load'])/F(branch['nominal_load'])-1)>erel:raise RuntimeError('actual load escaped bound')
                        normal_checks+=1
                result=infer_bounds(n,obs,t1=t,t2=2*t,gap=F(f['minimum_positive_gap']),relative_load_error=erel)
                scored=core.score(k,result);matrix_checks+=scored['checks']
                validation=[]
                for j in range(3):
                    sites=[1+j%(n-1),1+(j+1)%(n-1)]
                    loads={site:F((j+2)*(slot+1),7+n) for slot,site in enumerate(sites)}
                    rr,_=read(k,loads,p,rng,j,3);validation.append(rr);normal_checks+=1
                row={'fixture':f['id'],'order':n,'weak_link':f['weak_link'],'base_load':ts,
                     'status':result['status'],'reason':result.get('reason'),
                     'full_graph_disclosed':scored['full_graph'],'edge_verdicts':scored['edge_verdicts'],
                     'effective_relative_load_error':str(erel),
                     'cost':{'canonical_readings':len(plan),'repeat_baselines':len(repeats),
                         'validation_readings':3,'total_voltage_reads':len(plan)+len(repeats)+3,'shunt_channels':2}}
                rows.append(row);s=sums[ts];s['cases']+=1;s[result['status']]+=1;s['full_graph']+=int(scored['full_graph'])
                s['disclosed_edges']+=scored['disclosed_edges'];s['possible_edges']+=n*(n-1)//2
                out.write(json.dumps({'case':row,'budget':p,'fixture_id':f['id'],'readings':readings,
                     'repeat_baselines':repeats,'heldout_readings':validation,
                     'inference_observation_intervals':[{'nominal_loads':{str(j):str(v) for j,v in key},'interval':box.json()} for key,box in sorted(obs.items())],
                     'bounds':result},separators=(',',':'))+'\n')
    path=[[F(v) for v in row] for row in [[4,-1,0],[-1,3,-1],[0,-1,2]]]
    twin=[[F(v) for v in row] for row in [[4,F(-12,13),F(-5,13)],[F(-12,13),F(602,169),F(-59,169)],[F(-5,13),F(-59,169),F(243,169)]]]
    erel=core.load_bound({'shunt_relative_error':p['shunt_relative_error'],'relay_series_max_ohm':p['relay_series_max_ohm']},scale,F(1))
    twin_results=[]
    for current in (F(1,1000),F(1,10000)):
        rho=F(p['current_relative_error']);vr=F(p['voltage_relative_to_reading_error']);va=F(p['voltage_absolute_error_V']);boxes=[]
        for k in (path,twin):
            hlo=core.response(k,{1:1+erel});hhi=core.response(k,{1:1-erel})
            vlo=current*(1-rho)*hlo/scale;vhi=current*(1+rho)*hhi/scale
            boxes.append(Interval((vlo-va)/(1+vr),(vhi+va)/(1-vr)))
        twin_results.append({'current_A':str(current),'path_voltage_interval_V':boxes[0].json(),
              'triangle_voltage_interval_V':boxes[1].json(),'guaranteed_disjoint':boxes[0].hi<boxes[1].lo,
              'remaining_gap_V':str(boxes[1].lo-boxes[0].hi),
              'effective_relative_load_error':str(erel)})
    result={'schema':'nonagon.instrument-qualification-results.v1','status':'COMPLETED',
         'hardware_status':'synthetic qualification scenario; no physical apparatus connected',
         'protocol_sha256':hashlib.sha256(pp.read_bytes()).hexdigest(),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'fixture_count':len(fs),'case_count':len(rows),'normalization_and_load_checks':normal_checks,
         'accepted_matrix_containment_checks':matrix_checks,'false_edge_verdicts':0,'escaped_truth':0,
         'summary':[dict(base_load=t,**dict(v)) for t,v in sorted(sums.items())],
         'cases':rows,'frozen_twin_discrimination':twin_results,
         'records_sha256':hashlib.sha256(packet.read_bytes()).hexdigest(),
         'elapsed_seconds':round(time.monotonic()-start,3)}
    (HERE/'instrument_qualification_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--freeze',action='store_true');a=ap.parse_args()
    if a.freeze:freeze()
    else:run()
