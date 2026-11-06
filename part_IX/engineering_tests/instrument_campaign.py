"""Frozen, instrument-scaled synthetic acquisition tests for Nonagon tomography.

No instrument is connected.  Every signal, calibration envelope and circuit used
here is explicit.  Normalized interval inference receives no hidden circuit.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from random import Random
import sys
import time

HERE = Path(__file__).resolve().parent
CORE = HERE.parent if (HERE.parent/'bounded_intervention.py').exists() else HERE.parent/'part8_work'/'part_VIII'
sys.path.insert(0, str(CORE))
from bounded_intervention import Interval, infer_bounds, canonical_plan, load_key


def q(v):
    return v if isinstance(v, F) else F(v)


def solve(a, b):
    """Forward equilibrium solver; not used by the inference routine."""
    n = len(a)
    w = [[q(x) for x in row] + [q(b[i])] for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if w[i][j])
        w[j], w[pivot] = w[pivot], w[j]
        d = w[j][j]
        for k in range(j, n+1):
            w[j][k] /= d
        for i in range(j+1, n):
            d = w[i][j]
            for k in range(j, n+1):
                w[i][k] -= d*w[j][k]
    x = [F(0)]*n
    for i in range(n-1, -1, -1):
        x[i] = w[i][-1] - sum(w[i][j]*x[j] for j in range(i+1, n))
    return x


def response(k, loads):
    w = [row[:] for row in k]
    for j,t in loads.items():
        w[int(j)][int(j)] += q(t)
    return solve(w, [F(i == 0) for i in range(len(k))])[0]


def matrix_strings(k):
    return [[str(x) for x in row] for row in k]


def fixture(n, family, seed):
    rng = Random(seed)
    k = [[F(0) for _ in range(n)] for _ in range(n)]
    grounds = [F(rng.randint(1,4), 4) for _ in range(n)]
    for j,g in enumerate(grounds):
        k[j][j] = g
    edges = {}
    def edge(i,j,w):
        if i > j: i,j = j,i
        if (i,j) in edges: return
        edges[i,j] = w
        k[i][i] += w; k[j][j] += w; k[i][j] -= w; k[j][i] -= w
    for j in range(1,n):
        parent = 0 if family in (1,4) else j-1
        weight = F(rng.randint(2,8), 4)
        if family == 6 and j == n-1: weight = F(1,64)
        if family == 7 and j == n-1: weight = F(1,512)
        edge(parent,j,weight)
    for i in range(n):
        for j in range(i+1,n):
            if family in (2,3,4,5) and rng.randrange(4) < (3 if family == 5 else 1):
                edge(i,j,F(rng.randint(1,8),4))
    return {'id':f'n{n}_family{family}', 'order':n, 'family':family,
            'weak_link':None if family < 6 else ('1/64' if family == 6 else '1/512'),
            'minimum_positive_gap':'1/512', 'K_normalized':matrix_strings(k),
            'grounds':[str(x) for x in grounds],
            'edges':[{'sites':[i,j], 'conductance_normalized':str(v)} for (i,j),v in sorted(edges.items())]}


def freeze():
    HERE.mkdir(parents=True, exist_ok=True)
    fixture_dir = HERE/'instrument_fixtures'
    fixture_dir.mkdir(exist_ok=True)
    fixtures = [fixture(n,f,2026103000+n*100+f) for n in (3,4,6) for f in range(8)]
    budgets = [
        {'id':'ideal_reference','current_A':'1/1000','voltage_error_V':'0',
         'current_relative_error':'0','shunt_relative_error':'0','relay_series_max_ohm':'0','drift_error_V':'0'},
        {'id':'precision_calibration','current_A':'1/1000','voltage_error_V':'1/100000000',
         'current_relative_error':'1/1000000','shunt_relative_error':'1/1000000','relay_series_max_ohm':'1/200','drift_error_V':'1/1000000000'},
        {'id':'declared_008percent_1mA','current_A':'1/1000','voltage_error_V':'1/10000000',
         'current_relative_error':'1/100000','shunt_relative_error':'1/1250','relay_series_max_ohm':'1/20','drift_error_V':'1/20000000'},
        {'id':'declared_008percent_100uA','current_A':'1/10000','voltage_error_V':'1/1000000',
         'current_relative_error':'1/10000','shunt_relative_error':'1/1250','relay_series_max_ohm':'1/20','drift_error_V':'1/10000000'},
        {'id':'broad_calibration','current_A':'1/10000','voltage_error_V':'1/100000',
         'current_relative_error':'1/10000','shunt_relative_error':'1/100','relay_series_max_ohm':'1/20','drift_error_V':'1/1000000'},
    ]
    protocol = {'schema':'nonagon.instrument-protocol.v1','fixture_count':len(fixtures),
        'fixture_seeds':'2026103000 + 100*n + family','orders':[3,4,6],
        'families_per_order':8,'conductance_scale_S':'1/1000',
        'base_loads':['1/4','1','4'],'budgets':budgets,
        'hardware_status':'synthetic engineering test; no physical apparatus connected',
        'accuracy_labels':'declared scenario envelopes, not claims about any manufacturer configuration',
        'reconstruction':'existing rational interval algorithm, untouched for this campaign',
        'acquisition':'canonical n(n+1)/2 settings; start baseline; repeated baseline after each five nonbaseline settings and at end; 3 reserved validation settings',
        'scalar_interval':'V bounds divided by nominal I*(1±current_relative_error), times conductance scale',
        'drift_model':'independent bounded additive voltage drift, included in every reading envelope',
        'relay_model':'series resistor on each shunt branch, independently between 0 and declared maximum',
        'load_envelope':'covers calibrated shunt conductance plus relay series resistance, at largest nominal load 2*base',
        'edge_gap':'all positive network conductances >=1/512 normalized; exact absence only under this premise',
        'primary_metrics':['full_matrix_enclosed','full_graph_disclosed','correct_disclosed_edges','unresolved','truth_containment'],
        'operational_burden':'counts canonical acquisitions, extra baseline reads and heldout reads separately; no timing advantage assumed',
        'frozen_before_measurement_generation':True}
    for f in fixtures:
        (fixture_dir/(f['id']+'.json')).write_text(json.dumps(f,indent=2)+'\n')
    fixture_digest = hashlib.sha256(json.dumps(fixtures,sort_keys=True).encode()).hexdigest()
    protocol['fixture_sha256']=fixture_digest
    path=HERE/'instrument_protocol.json'
    if path.exists():
        raise FileExistsError('Protocol already frozen; refuse silent replacement')
    path.write_text(json.dumps(protocol,indent=2)+'\n')
    print(json.dumps({'protocol':str(path),'fixtures':len(fixtures),'fixture_sha256':fixture_digest}))


def scalar_interval(v, budget, scale):
    error=q(budget['voltage_error_V'])+q(budget['drift_error_V'])
    current=q(budget['current_A']); rel=q(budget['current_relative_error'])
    if v-error <= 0: raise ValueError('readout cannot certify positive voltage')
    return Interval((v-error)*scale/(current*(1+rel)),
                    (v+error)*scale/(current*(1-rel)))


def load_bound(budget, scale, largest_load):
    rel=q(budget['shunt_relative_error']); relay=q(budget['relay_series_max_ohm'])
    low=(1-rel)/(1+relay*largest_load*scale*(1-rel))
    return max(rel, 1-low)


def make_reading(k, nominal, budget, scale, rng, index, count):
    shunt_rel=q(budget['shunt_relative_error']);relay_max=q(budget['relay_series_max_ohm'])
    actual={};branches=[]
    for j,t in sorted(nominal.items()):
        e=shunt_rel*F(rng.randint(-3,3),4)
        relay=relay_max*F(rng.randint(0,4),4)
        shunt_conductance=q(t)*scale*(1+e)
        resistor=1/shunt_conductance
        conductance=1/(resistor+relay)
        actual[j]=conductance/scale
        branches.append({'site':j,'nominal_load':str(t),'shunt_relative_deviation':str(e),
            'shunt_resistance_ohm':str(resistor),'relay_series_ohm':str(relay),
            'actual_load':str(actual[j])})
    rho=q(budget['current_relative_error'])
    actual_current=q(budget['current_A'])*(1+rho*F(rng.randint(-3,3),4))
    h=response(k,actual)
    ideal_voltage=actual_current*h/scale
    noise=q(budget['voltage_error_V'])*F(rng.randint(-3,3),4)
    drift=q(budget['drift_error_V'])*F(2*index-count,max(1,count))
    measured=ideal_voltage+noise+drift
    box=scalar_interval(measured,budget,scale)
    if not box.contains(h): raise RuntimeError('normalization enclosure lost truth')
    return {'nominal_loads':{str(j):str(t) for j,t in nominal.items()},
            'actual_loads':{str(j):str(t) for j,t in actual.items()},'branches':branches,
            'actual_current_A':str(actual_current),'true_normalized_response':str(h),
            'ideal_voltage_V':str(ideal_voltage),'additive_readout_error_V':str(noise),
            'additive_drift_V':str(drift),'measured_voltage_V':str(measured),
            'normalized_interval':box.json()},box


def score(k,result):
    if result['status']!='enclosed':
        return {'checks':0,'full_graph':False,'edge_verdicts':{},'disclosed_edges':0}
    n=len(k);checks=0
    # Independent forward solves compute the true compliance columns.
    columns=[solve(k,[F(i==j) for i in range(n)]) for j in range(n)]
    g=[[columns[j][i] for j in range(n)] for i in range(n)]
    for key,truth in [('K_intervals',k),('G_intervals',g)]:
        for i,row in enumerate(truth):
            for j,value in enumerate(row):
                lo,hi=map(F,result[key][i][j])
                if not lo<=value<=hi: raise RuntimeError('matrix escaped '+key)
                checks+=1
    counts=Counter()
    for edge in result['edges']:
        i,j=edge['sites_zero_based'];actual=-k[i][j];lo,hi=map(F,edge['weight_interval'])
        if not lo<=actual<=hi: raise RuntimeError('edge escaped interval')
        status=edge['status'];counts[status]+=1
        if status=='present' and actual<=0: raise RuntimeError('false positive edge')
        if status=='absent_under_gap' and actual!=0: raise RuntimeError('false absent edge')
        if status=='inconsistent_with_gap': raise RuntimeError('valid fixture rejected gap')
        checks+=1
    return {'checks':checks,'full_graph':counts['unresolved']==0,
            'edge_verdicts':dict(counts),'disclosed_edges':counts['present']+counts['absent_under_gap']}


def run():
    start=time.monotonic();protocol_path=HERE/'instrument_protocol.json'
    protocol=json.loads(protocol_path.read_text())
    fixtures=[json.loads(p.read_text()) for p in sorted((HERE/'instrument_fixtures').glob('*.json'))]
    if hashlib.sha256(json.dumps(fixtures,sort_keys=True).encode()).hexdigest()!=protocol['fixture_sha256']:
        # Fixture order is frozen numerically rather than lexical n-order.
        fixtures=sorted(fixtures,key=lambda f:(f['order'],f['family']))
        if hashlib.sha256(json.dumps(fixtures,sort_keys=True).encode()).hexdigest()!=protocol['fixture_sha256']:
            raise RuntimeError('frozen fixture identity mismatch')
    scale=q(protocol['conductance_scale_S']);summary=defaultdict(Counter);rows=[];total_checks=0
    full_records=HERE/'instrument_case_records.jsonl'
    with full_records.open('w') as saved:
        for f in fixtures:
            n=f['order'];k=[[F(x) for x in row] for row in f['K_normalized']]
            for bi,budget in enumerate(protocol['budgets']):
                for ti,ts in enumerate(protocol['base_loads']):
                    t=F(ts);plan=canonical_plan(n,t,2*t);rng=Random(2026103100+10000*n+1000*f['family']+100*bi+ti)
                    obs={};readings=[];repeat_records=[]
                    for j,nominal in enumerate(plan):
                        record,box=make_reading(k,nominal,budget,scale,rng,j,len(plan))
                        readings.append(record);obs[load_key(nominal)]=box
                        if j>0 and (j%5==0 or j==len(plan)-1):
                            repeat,rb=make_reading(k,{},budget,scale,rng,j,len(plan))
                            repeat_records.append(repeat)
                            base=obs[()];lo=max(base.lo,rb.lo);hi=min(base.hi,rb.hi)
                            if lo>hi: raise RuntimeError('baseline intervals inconsistent')
                            obs[()]=Interval(lo,hi)
                    erel=load_bound(budget,scale,2*t)
                    for record in readings:
                        for branch in record['branches']:
                            deviation=abs(F(branch['actual_load'])/F(branch['nominal_load'])-1)
                            if deviation>erel: raise RuntimeError('shunt error escaped effective calibration')
                    result=infer_bounds(n,obs,t1=t,t2=2*t,gap=F(f['minimum_positive_gap']),relative_load_error=erel)
                    scored=score(k,result);total_checks+=scored['checks']
                    validation=[]
                    for j in range(3):
                        # Reserved, noncanonical settings executable with the
                        # same two routed shunts as the acquisition plan.
                        sites=[1+j%(n-1),1+(j+1)%(n-1)]
                        loads={site:F((j+2)*(slot+1),7+n)
                               for slot,site in enumerate(sites)}
                        rr,_=make_reading(k,loads,budget,scale,rng,j,3);validation.append(rr)
                    cost={'canonical_readings':len(plan),'repeat_baselines':len(repeat_records),
                          'validation_readings':len(validation),
                          'total_voltage_reads':len(plan)+len(repeat_records)+len(validation),
                          'known_intervention_pads':n-1,'fixed_sense_channels':1,'shunt_channels':2}
                    row={'fixture':f['id'],'order':n,'weak_link':f['weak_link'],'budget':budget['id'],
                         'base_load':ts,'effective_relative_load_error':str(erel),
                         'status':result['status'],'reason':result.get('reason'),
                         'full_graph_disclosed':scored['full_graph'],'edge_verdicts':scored['edge_verdicts'],
                         'truth_containment_checks':scored['checks'],'cost':cost}
                    rows.append(row);counts=summary[(budget['id'],ts)]
                    counts['cases']+=1;counts[result['status']]+=1;counts['full_graph']+=int(scored['full_graph'])
                    counts['disclosed_edges']+=scored['disclosed_edges'];counts['possible_edges']+=n*(n-1)//2
                    counts['unresolved_edges']+=n*(n-1)//2-scored['disclosed_edges']
                    record={'case':row,'budget':budget,'conductance_scale_S':str(scale),
                            'fixture_id':f['id'],'readings':readings,'repeat_baselines':repeat_records,
                            'heldout_readings':validation,'inference_observation_intervals':[
                                {'nominal_loads':{str(j):str(v) for j,v in key},'interval':box.json()}
                                for key,box in sorted(obs.items())], 'bounds':result}
                    saved.write(json.dumps(record,separators=(',',':'))+'\n')
    # A connected-class protocol must not invent disconnected unseen coordinates.
    disconnected=[[F(2),F(0),F(0)],[F(0),F(2),F(-1)],[F(0),F(-1),F(2)]]
    controls_obs={load_key(loads):Interval.point(response(disconnected,loads)) for loads in canonical_plan(3)}
    disconnected_result=infer_bounds(3,controls_obs,gap=F(1,512))
    if disconnected_result['status']!='unresolved': raise RuntimeError('disconnected control was not refused')
    # Frozen-pair discrimination in physical volts; explicit current+shunt bounds.
    path=[[F(v) for v in row] for row in [[4,-1,0],[-1,3,-1],[0,-1,2]]]
    twin=[[F(v) for v in row] for row in [[4,F(-12,13),F(-5,13)],[F(-12,13),F(602,169),F(-59,169)],[F(-5,13),F(-59,169),F(243,169)]]]
    twin_rows=[]
    for budget in protocol['budgets']:
        effective=load_bound(budget,scale,F(1));current=q(budget['current_A']);ci=q(budget['current_relative_error'])
        ve=q(budget['voltage_error_V'])+q(budget['drift_error_V']);boxes=[]
        for k in (path,twin):
            hlo=response(k,{1:1+effective});hhi=response(k,{1:1-effective})
            boxes.append(Interval(current*(1-ci)*hlo/scale-ve,current*(1+ci)*hhi/scale+ve))
        disjoint=boxes[0].hi<boxes[1].lo or boxes[1].hi<boxes[0].lo
        twin_rows.append({'budget':budget['id'],'effective_relative_load_error':str(effective),
                         'path_voltage_interval_V':boxes[0].json(),'triangle_voltage_interval_V':boxes[1].json(),
                         'guaranteed_disjoint':disjoint,'remaining_gap_V':str(boxes[1].lo-boxes[0].hi)})
    output={'schema':'nonagon.instrument-results.v1','status':'PASSED','hardware_status':protocol['hardware_status'],
        'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'fixture_count':len(fixtures),'case_count':len(rows),'truth_containment_checks':total_checks,
        'false_edge_verdicts':0,'escaped_truth':0,'scope':'conditional interval enclosures in a frozen synthetic circuit campaign',
        'summary':[dict(budget=b,base_load=t,**dict(c)) for (b,t),c in sorted(summary.items())],
        'cases':rows,'frozen_twin_discrimination':twin_rows,
        'controls':{'disconnected_network':disconnected_result},
        'records_sha256':hashlib.sha256(full_records.read_bytes()).hexdigest(),
        'elapsed_seconds':round(time.monotonic()-start,3)}
    (HERE/'instrument_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ('cases','controls')},indent=2))
    return output


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--freeze',action='store_true')
    args=ap.parse_args()
    if args.freeze:freeze()
    else:run()
