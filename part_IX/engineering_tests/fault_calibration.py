"""Prospective calibrated-shunt extension of the frozen catalog benchmark.

Uses fresh fixtures, explicit resistance/contact uncertainty and independent
per-acquisition actual shunts. The observer sees only finite prediction boxes,
readings and their externally declared absolute error bound.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import random
import time

import fault_benchmark as fb

HERE = Path(__file__).resolve().parent
PROTOCOL = {
    "schema": "nonagon.engineering.calibrated_fault_protocol.v1",
    "relationship_to_first_campaign": "Prospective extension with new confirmation seed, frozen before these additional cases. First readout-only campaign remains separate.",
    "confirmation_seed": 735826,
    "acquisition_seed": 182673,
    "confirmation_fixtures": 20,
    "specimens_per_fixture": 5,
    "dimension": 6,
    "catalog": "same preregistered five single-fault categories and healthy alternative as fault_protocol.json",
    "normalized_conductance_scale_S": "1/1000",
    "source_current_A": "1/1000",
    "normalized_compliance_equals_voltage_V": True,
    "resistor_relative_error": "1/1250",
    "resistor_absolute_error_ohm": "7/100",
    "relay_series_resistance_ohm": ["0", "1/20"],
    "calibration_scope": "declared simulated resistance envelope, not a claim to cover every operating specification of a product",
    "current_and_gain_error": "included in external normalized scalar error envelope, not separately randomized",
    "scalar_error_bounds": ["1/1000000", "1/100000", "1/10000"],
    "scalar_error_units": "V for the stated1mA source and1mS normalization",
    "actual_shunts": "Each active shunt in each acquisition independently varies inside its declared resistance and contact envelope; exact rational values retained for scoring only",
    "query_pool_and_policy": "frozen fault_benchmark definitions; strongest fixed plans maximize newly disjoint catalog prediction boxes before observing fault data",
    "budgets": [1, 6, 12, 21],
    "scope": "finite single-fault catalog with known healthy reference; no physical apparatus connected",
}


def load_bounds(query):
    low, high = {}, {}
    for j, nominal in query['loads'].items():
        t = F(nominal); resistance = F(1000)/t
        error = resistance*F(1,1250)+F(7,100)
        low[j] = str(F(1000)/(resistance+error+F(1,20)))
        high[j] = str(F(1000)/(resistance-error))
    return low, high


def signal_bounds(G, query):
    if not query['loads']:
        v=fb.prediction(G,query);return v,v
    low,high=load_bounds(query)
    # Static fixed-port response decreases in every positive grounding load.
    return fb.prediction(G,{**query,'loads':high}),fb.prediction(G,{**query,'loads':low})


def actual_query(query, rng):
    loads={}
    for j,nominal in query['loads'].items():
        R=F(1000)/F(nominal);error=R*F(1,1250)+F(7,100)
        resistance=R+error*F(rng.randint(-4,4),5)+F(rng.randint(0,5),100)
        loads[j]=str(F(1000)/resistance)
    return {**query,'loads':loads}


def survivors_after(boxes, measured, eta, survivors):
    return [i for i in survivors if boxes[i][0]-eta<=measured<=boxes[i][1]+eta]


def reading(boxes, actual, eta, survivors):
    lo,hi=actual-eta,actual+eta;choices={lo,hi,actual}
    for i in survivors:
        for v in [boxes[i][0]-eta,boxes[i][1]+eta]:
            if lo<=v<=hi:choices.add(v)
    return max(choices,key=lambda m:(len(survivors_after(boxes,m,eta,survivors)),abs(m-actual),-m))


def plan(query_ids,table,eta):
    count=len(next(iter(table.values())));pairs=list(fb.combinations(range(count),2));masks={}
    for q in query_ids:
        mask=0;boxes=table[q]
        for index,(i,j) in enumerate(pairs):
            gap=max(boxes[i][0]-boxes[j][1],boxes[j][0]-boxes[i][1])
            if gap>2*eta:mask|=1<<index
        masks[q]=mask
    remaining=(1<<len(pairs))-1;available=list(query_ids);selected=[]
    for _ in range(min(21,len(available))):
        q=max(available,key=lambda q:((masks[q]&remaining).bit_count(),masks[q].bit_count(),-query_ids.index(q)))
        selected.append(q);remaining&=~masks[q];available.remove(q)
    return selected


def robust_separation(table,queries,true,eta,ids):
    nearest=None
    for other in range(len(ids)):
        if other==true:continue
        gap=max(max(table[q][true][0]-table[q][other][1],table[q][other][0]-table[q][true][1]) for q in queries)-2*eta
        item=(gap,ids[other])
        if nearest is None or item<nearest:nearest=item
    return {'margin':str(nearest[0]),'nearest_alternative':nearest[1],
            'all_errors_identifiable':nearest[0]>0}


def main():
    started=time.monotonic();protocol_text=json.dumps(PROTOCOL,indent=2,sort_keys=True)+'\n'
    p=HERE/'fault_calibration_protocol.json'
    if p.exists() and p.read_text()!=protocol_text:raise RuntimeError('Frozen protocol mismatch')
    p.write_text(protocol_text);digest=hashlib.sha256(protocol_text.encode()).hexdigest()
    canonical,pool=fb.intervention_queries();scans=fb.conventional_queries();query_map={q['id']:q for q in pool+scans}
    canonical_ids=[q['id'] for q in canonical];pool_ids=[q['id'] for q in pool];scan_ids=[q['id'] for q in scans]
    generator=random.Random(PROTOCOL['confirmation_seed']);acquisition=random.Random(PROTOCOL['acquisition_seed'])
    fixtures=[];trials=[];summary={};containment_checks=0
    for fi in range(20):
        K=fb.fixture(generator);catalog=fb.catalog(K);ids=[c['id'] for c in catalog];Gs=[fb.inverse(c['K']) for c in catalog]
        table={qid:[signal_bounds(G,q) for G in Gs] for qid,q in query_map.items()}
        true_indices={kind:generator.choice([i for i,c in enumerate(catalog) if c['type']==kind]) for kind in fb.PROTOCOL['fault_types']}
        fdoc={'id':f'calibration_fixture_{fi:02d}','healthy_K':fb.encode_matrix(K),
              'catalog':[{**{k:v for k,v in c.items() if k!='K'},'K':fb.encode_matrix(c['K'])} for c in catalog],
              'prediction_intervals':{q:[[str(a),str(b)] for a,b in boxes] for q,boxes in table.items()},'plans':{}}
        # The physical acquisition for a specimen/query is shared across policies
        # and error levels, preventing per-policy generation advantages.
        actuals={}
        for kind,true in true_indices.items():
            actuals[kind]={}
            for qid,q in query_map.items():
                aq=actual_query(q,acquisition);signal=fb.prediction(Gs[true],aq)
                low,high=table[qid][true]
                if not low<=signal<=high:raise ValueError('actual load escaped physical prediction interval')
                actuals[kind][qid]=(aq['loads'],signal);containment_checks+=1
        for eta_text in PROTOCOL['scalar_error_bounds']:
            eta=F(eta_text)
            plans={'unloaded_port':['baseline'],'canonical_intervention':canonical_ids,
                   'greedy_fixed_intervention':plan(pool_ids,table,eta),'natural_reciprocal_scan':scan_ids,
                   'greedy_fixed_reciprocal_scan':plan(scan_ids,table,eta)}
            fdoc['plans'][eta_text]={name:{'queries':queries,'costs_at_budget':{str(b):fb.costs(queries[:b],query_map) for b in PROTOCOL['budgets'] if b<=len(queries)}} for name,queries in plans.items()}
            for kind,true in true_indices.items():
                for policy,queries in plans.items():
                    survivors=list(range(len(ids)));history=[];checkpoints={}
                    for count,qid in enumerate(queries,1):
                        loads,actual=actuals[kind][qid]
                        measured=reading(table[qid],actual,eta,survivors)
                        survivors=survivors_after(table[qid],measured,eta,survivors)
                        if true not in survivors:raise ValueError('true fault removed by calibrated uncertainty')
                        history.append({'query':qid,'actual_loads':loads,'actual':str(actual),'measured':str(measured),
                                        'error':str(measured-actual),'survivors':[ids[i] for i in survivors]})
                        if count in PROTOCOL['budgets']:
                            robust=robust_separation(table,queries[:count],true,eta,ids)
                            cp={'status':'identified' if len(survivors)==1 else 'unresolved','survivors':len(survivors),
                                'identified':ids[survivors[0]] if len(survivors)==1 else None,**robust}
                            checkpoints[str(count)]=cp
                            key=(eta_text,policy,count)
                            s=summary.setdefault(key,{'eta':eta_text,'policy':policy,'budget':count,'cases':0,'identified':0,'unresolved':0,'wrong':0,'all_error_box_identifiable':0,'by_fault_type':{}})
                            s['cases']+=1;s[cp['status']]+=1;s['all_error_box_identifiable']+=robust['all_errors_identifiable']
                            if cp['identified'] is not None and cp['identified']!=ids[true]:s['wrong']+=1
                            st=s['by_fault_type'].setdefault(kind,{'cases':0,'identified':0,'all_error_box_identifiable':0})
                            st['cases']+=1;st['identified']+=cp['status']=='identified';st['all_error_box_identifiable']+=robust['all_errors_identifiable']
                    trials.append({'fixture':fdoc['id'],'fault_type':kind,'true_candidate':ids[true],'eta':eta_text,
                                   'policy':policy,'readings':history,'checkpoints':checkpoints})
        fixtures.append(fdoc);print(f'Calibrated fixture{fi+1}/20 complete',flush=True)
    if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Protocol changed during run')
    result={'schema':'nonagon.engineering.calibrated_fault_results.v1','protocol_sha256':digest,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'confirmation_fixtures':20,
            'distinct_confirmation_specimens':100,'trial_count':len(trials),'physical_signal_containment_checks':containment_checks,
            'query_definitions':query_map,'summaries':list(summary.values()),'fixtures':fixtures,'trials':trials,
            'elapsed_seconds':round(time.monotonic()-started,3)}
    output=HERE/'fault_calibration_results.json';output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(output),'trials':len(trials),'physical_containment_checks':containment_checks,
                     'summary_at_10microvolt':[s for s in result['summaries'] if s['eta']=='1/100000' and s['budget'] in[6,21]],
                     'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
