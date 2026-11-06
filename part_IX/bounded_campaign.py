"""Small exact-rational containment campaign for bounded intervention inference."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
from collections import Counter
import json
from bounded_intervention import (Interval, inverse, infer_bounds, observation_intervals, canonical_plan, load_key)

HERE=Path(__file__).resolve().parent

def oracle_for(k):
    def observe(loads):
        w=[row[:] for row in k]
        for j,v in loads.items():w[int(j)][int(j)]+=F(v)
        return inverse(w)[0][0]
    return observe


def test_containment(k,result):
    if result['status']!='enclosed':return 0
    g=inverse(k);count=0
    for key,truth in [('G_intervals',g),('K_intervals',k)]:
        for i,row in enumerate(truth):
            for j,value in enumerate(row):
                a,b=map(F,result[key][i][j])
                if not a<=value<=b:raise RuntimeError('True matrix escaped '+key)
                count+=1
    for edge in result['edges']:
        i,j=edge['sites_zero_based'];actual=-k[i][j]
        a,b=map(F,edge['weight_interval'])
        if not a<=actual<=b:raise RuntimeError('True edge escaped interval')
        if edge['status']=='present' and not actual>0:raise RuntimeError('False presence')
        if edge['status']=='absent_under_gap' and actual!=0:raise RuntimeError('False absence')
        if edge['status']=='inconsistent_with_gap':raise RuntimeError('Valid gap marked inconsistent')
        count+=1
    return count


def network(n,rng):
    k=[[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):k[i][i]=F(rng.randint(1,4),4)
    def edge(i,j):
        w=F(rng.randint(1,5),4)
        k[i][i]+=w;k[j][j]+=w;k[i][j]-=w;k[j][i]-=w
    seen=set()
    for i in range(1,n):
        j=rng.randrange(i);edge(i,j);seen.add((j,i))
    for i in range(n):
        for j in range(i+1,n):
            if (i,j) not in seen and rng.random()<F(1,4):edge(i,j)
    return k


def run():
    path=[[F(v) for v in row] for row in [[4,-1,0],[-1,3,-1],[0,-1,2]]]
    twin=[[F(v) for v in row] for row in [[4,F(-12,13),F(-5,13)],
              [F(-12,13),F(602,169),F(-59,169)],
              [F(-5,13),F(-59,169),F(243,169)]]]
    campaigns=[];checks=0
    etas=[F(0),F(1,10**12),F(1,10**9),F(1,10**6),F(1,10000)]
    representative=None
    for name,k in [('v7_path',path),('v7_triangle',twin)]:
        for eta in etas:
            result=infer_bounds(3,observation_intervals(3,oracle_for(k),eta),gap=F(1,4))
            # Triangle minimum edge is 59/169, exceeding the declared 1/4.
            checks+=test_containment(k,result)
            campaigns.append({'network':name,'eta':str(eta),'status':result['status'],
                'edge_verdicts':dict(Counter(x['status'] for x in result.get('edges',[]))),
                'reason':result.get('reason')})
            if name=='v7_path' and eta==F(1,10**9):representative=result
    rng=Random(20261030);networks=0;statuses=Counter();edge_counts=Counter();calibrated_counts=Counter();noisy_checks=0
    for n in range(2,7):
        for _ in range(5):
            k=network(n,rng);networks+=1
            for eta in [F(0),F(1,10**12),F(1,10**9),F(1,10**6)]:
                result=infer_bounds(n,observation_intervals(n,oracle_for(k),eta),gap=F(1,4))
                statuses[(str(eta),result['status'])]+=1
                checks+=test_containment(k,result)
                edge_counts.update(x['status'] for x in result.get('edges',[]))
            # Hidden physical load values and reported readings are independently
            # displaced within their externally supplied calibration envelopes.
            for rel in [F(1,10**9),F(1,10**6),F(1,100)]:
                eta=F(1,10**9);obs={}
                for nominal in canonical_plan(n):
                    actual={j:v*(1+rel*F(rng.randint(-3,3),4)) for j,v in nominal.items()}
                    center=oracle_for(k)(actual)+eta*F(rng.randint(-3,3),4)
                    obs[load_key(nominal)]=Interval(center-eta,center+eta)
                result=infer_bounds(n,obs,gap=F(1,4),relative_load_error=rel)
                calibrated_counts[(str(rel),result['status'])]+=1
                got=test_containment(k,result);checks+=got;noisy_checks+=got
    p=oracle_for(path)({1:F(1)});q=oracle_for(twin)({1:F(1)})
    separation=abs(p-q)
    sep_rows=[]
    for eta in [F(1,10000),F(1,1000),F(1,400)]:
        pbox=Interval(p-eta,p+eta);qbox=Interval(q-eta,q+eta)
        disjoint=pbox.hi<qbox.lo or qbox.hi<pbox.lo
        if disjoint!=(2*eta<separation):raise RuntimeError('Separation threshold failed')
        sep_rows.append({'eta':str(eta),'guaranteed_disjoint':disjoint})
    calibrated_separation=[]
    for rel in [F(1,100),F(1,1000),F(1,1000000)]:
        plo=oracle_for(path)({1:1+rel});phi=oracle_for(path)({1:1-rel})
        qlo=oracle_for(twin)({1:1+rel});qhi=oracle_for(twin)({1:1-rel})
        margin=(qlo-phi)/2
        if margin<=0:raise RuntimeError('Twin calibration intervals overlap')
        calibrated_separation.append({'relative_load_error':str(rel),
            'path_prediction_interval':[str(plo),str(phi)],
            'triangle_prediction_interval':[str(qlo),str(qhi)],
            'strict_additional_response_error_threshold':str(margin)})
    controls={}
    for label,func in [('reversed_interval',lambda:Interval(F(2),F(1))),
                       ('zero_denominator',lambda:F(1)/Interval(F(-1),F(1))),
                       ('negative_radicand',lambda:Interval(F(-1),F(1)).sqrt()),
                       ('binary_float_input',lambda:Interval.point(0.1))]:
        try:func()
        except (ValueError,TypeError,ZeroDivisionError):controls[label]='rejected'
        else:raise RuntimeError('Alteration not rejected '+label)
    out={'format':'nonagon-bounded-campaign/1','arithmetic':'fractions with exact outward dyadic square-root bounds',
         'measurement_contract':'simultaneous normalized-response intervals; exact known labels; optional independently varying relative load-calibration intervals',
         'network_count_random':networks,'random_orders':[2,3,4,5,6],
         'random_evaluations':sum(statuses.values()),
         'independently_noisy_load_and_response_evaluations':sum(calibrated_counts.values()),
         'noisy_calibration_status_counts':[{'relative_load_error':e,'status':s,'count':v} for (e,s),v in sorted(calibrated_counts.items())],
         'noisy_calibration_containment_checks':noisy_checks,
         'random_status_counts':[{'eta':eta,'status':s,'count':v} for (eta,s),v in sorted(statuses.items())],
         'random_edge_verdicts':dict(edge_counts),'matrix_and_edge_containment_checks':checks,
         'escaped_truth_or_false_edge_verdicts':0,'twin_cases':campaigns,
         'twin_one_load_separation':{'path':str(p),'triangle':str(q),'absolute_separation':str(separation),
                 'strict_per_reading_error_threshold':str(separation/2),'checks':sep_rows},
         'twin_discrimination_with_uncertain_shunts':calibrated_separation,
         'adversarial_controls':controls,'representative_path_eta_1e9':representative}
    (HERE/'bounded_campaign.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('representative_path_eta_1e9','twin_cases')},indent=2))
    return out

if __name__=='__main__':run()
