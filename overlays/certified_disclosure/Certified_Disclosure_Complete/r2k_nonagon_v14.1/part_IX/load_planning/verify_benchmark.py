#!/usr/bin/env python3
"""Read-only independent audit of saved finite-load benchmark records.

Imports no planner, benchmark or interval inference. Exact rational forward
predictions use the inverse and one/two-site Woodbury formulas. This audit checks
recorded containment and decisions; it is not a new proof of interval propagation.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def key(loads):return tuple(sorted((int(j),F(t)) for j,t in loads.items()))
def asbox(xs):return tuple(map(F,xs))
def point(x):return(F(x),F(x))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    vals=[x*y for x in a for y in b];return(min(vals),max(vals))
def div(a,b):
    if b[0]<=0<=b[1]:raise ArithmeticError('Denominator intersects zero')
    return mul(a,(1/b[1],1/b[0]))
def positive(a):
    if a[0]<=0:raise ArithmeticError('Not strictly positive')
    return a

def invert(k):
    n=len(k);w=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(k)]
    for j in range(n):
        p=next(i for i in range(j,n) if w[i][j]);w[j],w[p]=w[p],w[j]
        pivot=w[j][j];w[j]=[x/pivot for x in w[j]]
        for i in range(n):
            if i==j:continue
            value=w[i][j];w[i]=[a-value*b for a,b in zip(w[i],w[j])]
    return [row[n:] for row in w]

def response(g,loads):
    active=sorted(loads)
    if not active:return g[0][0]
    if len(active)==1:
        j=active[0];t=loads[j]
        return g[0][0]-t*g[0][j]**2/(1+t*g[j][j])
    if len(active)!=2:raise RuntimeError('More than two shunts in saved record')
    i,j=active;a=1/loads[i]+g[i][i];c=1/loads[j]+g[j][j];b=g[i][j]
    x,y=g[0][i],g[0][j]
    return g[0][0]-(c*x*x-2*b*x*y+a*y*y)/(a*c-b*b)

def plan(n,t):
    ans=[()]
    for j in range(1,n):ans.extend([((j,t),),((j,2*t),)])
    ans.extend(((i,t),(j,t)) for i in range(1,n) for j in range(i+1,n))
    return ans

def radius(hw,t):
    e,r,s=map(F,(hw['shunt_relative_error'],hw['relay_series_max_ohm'],hw['conductance_scale_S']))
    return max(e,1-(1-e)/(1+r*2*t*s*(1-e)))

class Audit:
    def __init__(self):self.checks=0;self.readings=0;self.branches=0;self.decisions=0
    def check(self,condition,message):
        self.checks+=1
        if not condition:raise RuntimeError(message)
    def reading(self,record,g,hw,e):
        scale,current,crel=map(F,(hw['conductance_scale_S'],hw['current_A'],hw['current_relative_error']))
        erel,relaymax=map(F,(hw['shunt_relative_error'],hw['relay_series_max_ohm']))
        nominal=dict(key(record['nominal_loads']));actual={}
        self.check(len(nominal)<=2,'More than two nominal shunts')
        self.check(len(record['branches'])==len(nominal),'Branch count mismatch')
        for branch in record['branches']:
            j=branch['site'];t=F(branch['nominal']);delta=F(branch['relative_shunt_deviation']);r=F(branch['relay_ohm'])
            self.check(j in nominal and nominal[j]==t and j not in actual,'Branch address mismatch')
            self.check(abs(delta)<=erel and 0<=r<=relaymax,'Calibration realization escaped declared interval')
            v=t*(1+delta)/(1+r*t*scale*(1+delta))
            self.check(F(branch['actual'])==v,'Actual shunt equation failed')
            self.check(t*(1-e)<=v<=t*(1+e),'Actual shunt escaped symmetric inference envelope')
            actual[j]=v;self.branches+=1
        h=response(g,actual);self.check(h==F(record['true_response']),'Independent exact forward prediction differs')
        ai=F(record['actual_current_A']);self.check(current*(1-crel)<=ai<=current*(1+crel),'Current escaped bound')
        err,drift=map(F,(record['readout_error_V'],record['drift_V']))
        self.check(abs(err)<=F(hw['voltage_error_V']) and abs(drift)<=F(hw['drift_error_V']),'Voltage realization escaped envelope')
        measured=F(record['measured_voltage_V']);self.check(measured==ai*h/scale+err+drift,'Raw voltage equation failed')
        ve=F(hw['voltage_error_V'])+F(hw['drift_error_V'])
        expected=((measured-ve)*scale/(current*(1+crel)),(measured+ve)*scale/(current*(1-crel)))
        box=asbox(record['normalized_response_interval'])
        self.check(box==expected and box[0]<=h<=box[1],'Normalized response interval differs or loses truth')
        self.readings+=1
        return key(record['nominal_loads']),box
    def recommendation(self,doc,rec,g):
        n=doc['network_order'];t0=F(doc['pilot_base_load']);hw=doc['hardware']
        obs={key(r['nominal_loads']):asbox(r['normalized_response_interval']) for r in doc['readings']}
        self.check(set(obs)==set(plan(n,t0)[:2*n-1]),'Incomplete pilot schema')
        eta=(F(hw['voltage_error_V'])+F(hw['drift_error_V']))*F(hw['conductance_scale_S'])/F(hw['current_A'])+F(hw['current_relative_error'])*obs[()][1]
        self.check(eta==F(rec['raw_response_error_bound']),'Planner response uncertainty differs')
        e=radius(hw,t0);a=(t0*(1-e),t0*(1+e));b=(2*t0*(1-e),2*t0*(1+e));singles=[]
        try:
            sep=positive(sub(b,a))
            for j in range(1,n):
                d1=positive(sub(obs[()],obs[((j,t0),)]));d2=positive(sub(obs[()],obs[((j,2*t0),)]))
                y1,y2=div(a,d1),div(b,d2);slope=div(sub(y2,y1),sep)
                intercept=positive(sub(y1,mul(a,slope)));qb=positive(div(point(1),intercept));bb=div(slope,intercept)
                if bb[1]<=0:raise ArithmeticError('b incompatible with class')
                bb=(max(F(0),bb[0]),bb[1]);singles.append((j,qb,bb))
        except ArithmeticError:
            self.check(rec['status']=='fallback' and F(rec['selected_base_load'])==t0,'Unresolved pilot did not use fallback')
        else:
            self.check(rec['status']=='recommended','Separated pilot unexpectedly fell back')
            candidates=sorted(set(map(F,hw['candidate_base_loads'])));scores=[]
            indexed={F(r['base_load']):r for r in rec['candidates']}
            self.check(set(indexed)==set(candidates),'Candidate list differs')
            for t in candidates:
                score=max(eta*(4/t+8*bb[1]+4*t*bb[1]**2)/qb[0]+3*radius(hw,t) for j,qb,bb in singles)
                row=indexed[t];self.check(F(row['worst_site_local_diagnostic'])==score,'Independent objective differs')
                self.check(F(row['relative_load_envelope'])==radius(hw,t),'Candidate calibration radius differs')
                scores.append((score,t!=t0,t))
            chosen=min(scores)[2];self.check(F(rec['selected_base_load'])==chosen,'Finite minimax choice differs')
        actual={r['site_zero_based']:r for r in rec['pilot_single_enclosures']}
        self.check(set(actual)=={j for j,_,_ in singles},'Pilot enclosure sites differ')
        for j,qb,bb in singles:
            self.check(asbox(actual[j]['q_interval'])==qb and asbox(actual[j]['b_interval'])==bb,'Pilot local interval arithmetic differs')
            self.check(qb[0]<=g[0][j]**2<=qb[1] and bb[0]<=g[j][j]<=bb[1],'Pilot enclosure lost physical truth')
        chosen=F(rec['selected_base_load']);self.check(rec['total_scalar_count_with_pilot']==2*n-1+n*(n+1)//2,'Planner budget differs')
        if 'nominal_configuration_sequence' in rec:
            self.check([key(r['loads']) for r in rec['nominal_configuration_sequence']]==plan(n,chosen),'Recommended acquisition sequence differs')
        self.decisions+=1
    def bounds(self,bounds,k,g,expected,gap):
        saved={key(r['nominal_loads']):asbox(r['normalized_response_interval']) for r in bounds['acquisition_intervals']}
        self.check(saved==expected,'Production inference did not receive retained observation intervals')
        if bounds['status']!='enclosed':
            self.check(bounds['status']=='unresolved','Unknown inference verdict')
            return False,False,0,None
        n=len(k)
        for name,truth in [('K_intervals',k),('G_intervals',g)]:
            boxes=bounds[name];self.check(len(boxes)==n and all(len(row)==n for row in boxes),'Matrix dimensions differ')
            for i in range(n):
                for j in range(n):
                    lo,hi=asbox(boxes[i][j]);self.check(lo<=truth[i][j]<=hi,'Returned matrix enclosure loses truth')
        disclosed=0;seen=set()
        for edge in bounds['edges']:
            i,j=edge['sites_zero_based'];self.check(0<=i<j<n and (i,j) not in seen,'Edge address duplicate or invalid');seen.add((i,j))
            w=-k[i][j];lo,hi=asbox(edge['weight_interval']);status=edge['status']
            self.check(0<=lo<=w<=hi,'Edge interval loses truth')
            self.check(status in ('present','absent_under_gap','unresolved'),'Invalid graph decision')
            if status=='present':self.check(lo>0 and w>0,'False positive edge');disclosed+=1
            if status=='absent_under_gap':self.check(hi<gap and w==0,'False absent edge');disclosed+=1
        self.check(len(seen)==n*(n-1)//2,'Incomplete edge set')
        width=max(F(hi)-F(lo) for row in bounds['K_intervals'] for lo,hi in row)
        return True,disclosed==n*(n-1)//2,disclosed,width


def run(records,receipt_path,output):
    audit=Audit();receipt=json.loads(receipt_path.read_text());protocol=json.loads((HERE/'protocol.json').read_text())
    audit.check(sha(records)==receipt['records_sha256'],'Saved record hash differs')
    audit.check(sha(HERE/'protocol.json')==receipt['protocol_sha256'],'Protocol hash differs')
    for name,digest in receipt['source_sha256'].items():audit.check(sha(HERE/name)==digest,'Producer source hash differs: '+name)
    audit.check(sha(HERE.parent/'bounded_intervention.py')==receipt['production_interval_source_sha256'],'Production interval source changed')
    fixtures={f['id']:f for f in protocol['fixtures']};budgets={b['id']:b for b in protocol['budgets']}
    summary=defaultdict(Counter);selection=Counter();case_rows=[];seen=set();shared_raw=0;charged_raw=0;cache={}
    for line in records.open():
        packet=json.loads(line);f=packet['fixture'];budget=packet['budget'];hw=packet['hardware'];n=f['n'];k=[list(map(F,row)) for row in f['K']]
        audit.check(f==fixtures[f['id']] and budget==budgets[budget['id']],'Packet fixture or budget is not frozen')
        audit.check(hw==dict(protocol['hardware'],**{k:v for k,v in budget.items() if k!='id'}),'Packet hardware does not match frozen calibration budget')
        ident=(f['id'],budget['id']);audit.check(ident not in seen,'Duplicate prospective case');seen.add(ident)
        if f['id'] not in cache:cache[f['id']]=invert(k)
        g=cache[f['id']];gap=F(protocol['minimum_positive_edge_gap'])
        for i in range(n):
            audit.check(sum(k[i])>=0 and k[i][i]>0,'Fixture grounding class differs')
            for j in range(n):
                audit.check(k[i][j]==k[j][i] and g[i][j]>0,'Fixture symmetry/inverse positivity differs')
                if i!=j:audit.check(k[i][j]<=0 and (k[i][j]==0 or -k[i][j]>=gap),'Fixture edge gap premise fails')
        t0=F(protocol['pilot_base_load']);pilot={}
        for record in packet['pilot_readings']:
            config,box=audit.reading(record,g,hw,radius(hw,t0));pilot[config]=box;shared_raw+=1
        audit.check(list(pilot)==plan(n,t0)[:2*n-1],'Pilot plan differs')
        doc=packet['pilot_document'];audit.check(doc['hardware']==hw and doc['network_order']==n,'Pilot contract differs')
        audit.check({key(r['nominal_loads']):asbox(r['normalized_response_interval']) for r in doc['readings']}==pilot,'Planner pilot is not saved measured pilot')
        recommendation=packet['recommendation'];audit.recommendation(doc,recommendation,g)
        selection[budget['id'],recommendation['selected_base_load'],recommendation['status']]+=1
        audit.check(len(packet['arms'])==2,'Two comparison arms required');arm_rows=[]
        for arm in packet['arms']:
            case=arm['case'];t=F(case['base_load']);e=radius(hw,t);expected={};observed=[]
            if case['arm']=='fixed_pilot_strength':audit.check(t==t0,'Fixed arm changed strength')
            else:audit.check(case['arm']=='pilot_selected_strength' and t==F(recommendation['selected_base_load']),'Planned arm uses wrong strength')
            for record in arm['readings']:
                config,box=audit.reading(record,g,hw,e);shared_raw+=1;observed.append(config)
                if config==():box=(max(pilot[()][0],box[0]),min(pilot[()][1],box[1]))
                elif config in pilot and pilot[config][1]-pilot[config][0]<box[1]-box[0]:box=pilot[config]
                expected[config]=box
            audit.check(observed==plan(n,t),'Final canonical plan differs')
            audit.check(F(case['relative_load_envelope'])==e and F(arm['bounds']['contract']['relative_load_error'])==e,'Inference calibration envelope differs')
            enc,graph,disclosed,width=audit.bounds(arm['bounds'],k,g,expected,gap)
            audit.check((case['matrix_enclosed'],case['full_graph'],case['disclosed_edges'])==(enc,graph,disclosed),'Case summary disagrees with bounds')
            audit.check(case['max_K_interval_width']==(None if width is None else str(width)),'Case width differs')
            count=2*n-1+n*(n+1)//2
            audit.check(case['pilot_scalar_count']==2*n-1 and case['final_scalar_count']==n*(n+1)//2 and case['total_scalar_count']==count,'Arm budget differs')
            charged_raw+=count;s=summary[budget['id'],case['arm']];s['cases']+=1;s['matrix_enclosed']+=enc;s['full_graph']+=graph;s['disclosed_edges']+=disclosed;s['possible_edges']+=n*(n-1)//2;s['scalar_readings']+=count
            case_rows.append(case);arm_rows.append(case)
        audit.check(arm_rows[0]['total_scalar_count']==arm_rows[1]['total_scalar_count'],'Unequal paired budget')
    expected_summary=[dict(budget=b,arm=a,**dict(s)) for (b,a),s in sorted(summary.items())]
    audit.check(expected_summary==receipt['summary'],'Aggregate summaries differ')
    expected_selection=[{'budget':b,'base_load':t,'status':st,'count':v} for (b,t,st),v in sorted(selection.items())]
    audit.check(expected_selection==receipt['selection_counts'],'Strength-selection summaries differ')
    paired=[]
    for budget in protocol['budgets']:
        outcomes=Counter()
        for f in protocol['fixtures']:
            fixed,planned=[r for r in case_rows if r['fixture']==f['id'] and r['budget']==budget['id']]
            outcomes['planned_only_matrix_enclosed']+=planned['matrix_enclosed'] and not fixed['matrix_enclosed']
            outcomes['fixed_only_matrix_enclosed']+=fixed['matrix_enclosed'] and not planned['matrix_enclosed']
            outcomes['planned_only_full_graph']+=planned['full_graph'] and not fixed['full_graph']
            outcomes['fixed_only_full_graph']+=fixed['full_graph'] and not planned['full_graph']
            if fixed['matrix_enclosed'] and planned['matrix_enclosed']:
                wf,wp=map(F,(fixed['max_K_interval_width'],planned['max_K_interval_width']))
                outcomes['paired_enclosures']+=1;outcomes['planned_width_smaller']+=wp<wf
                outcomes['fixed_width_smaller']+=wf<wp;outcomes['equal_width']+=wf==wp
        paired.append(dict(budget=budget['id'],**dict(outcomes)))
    audit.check(paired==receipt['paired_comparisons'],'Paired outcome summaries differ')
    audit.check(len(seen)==receipt['prospective_cases']==len(fixtures)*len(budgets),'Prospective case count differs')
    audit.check(len(case_rows)==receipt['arm_cases'],'Arm case count differs')
    audit.check(charged_raw==receipt['scalar_readings_counting_shared_pilot_per_arm'],'Charged scalar count differs')
    audit.check(shared_raw==receipt['unique_generated_scalar_readings'],'Generated unique scalar count differs')
    result={'schema':'nonagon.finite-load-benchmark-independent-verification.v1','status':'passed','exact_predicates':audit.checks,
            'prospective_cases':len(seen),'arm_cases':len(case_rows),'pilot_only_strength_decisions':audit.decisions,
            'unique_raw_readings_replayed':audit.readings,'independent_actual_shunt_checks':audit.branches,
            'charged_scalar_readings':charged_raw,'false_edge_calls':0,'escaped_truth':0,
            'records_sha256':sha(records),'producer_receipt_sha256':sha(receipt_path),'checker_sha256':sha(__file__),
            'scope':'Independent exact raw physics, calibration, pilot objective, retained acquisition intervals, matrix/edge truth containment and aggregate counts; no reimplementation of the complete interval propagation proof.'}
    output.write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--records',type=Path,default=HERE/'results/benchmark_cases.jsonl')
    ap.add_argument('--receipt',type=Path,default=HERE/'results/benchmark_receipt.json')
    ap.add_argument('--output',type=Path,default=HERE/'results/independent_verification.json')
    a=ap.parse_args();print(json.dumps(run(a.records,a.receipt,a.output),indent=2))
