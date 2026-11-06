#!/usr/bin/env python3
"""Independent exact replay of the finite-library intervention campaign.

No planner, reconstruction or cofactor-model functions are imported. The checker
uses the independent Gaussian solver and verifies predicted data, complete
survivor histories, bounded adversarial observations, query choices and summaries.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

from verify_intervention import solve, rational


def check_campaign(document):
    count = 0
    def require(condition, reason):
        nonlocal count
        count += 1
        if not condition:
            raise ValueError(reason)
    require(document.get('schema') == 'nonagon.active_interventions.finite_library.v1', 'schema')
    candidates, queries, trials = document['candidates'], document['feasible_queries'], document['trials']
    require(len(candidates) == document['candidate_count'] == 21, 'candidate count')
    require(len(queries) == document['feasible_query_count'] == 15, 'query count')
    require(len(trials) == 189, 'trial count')
    ids = [candidate['id'] for candidate in candidates]
    qids = [query['id'] for query in queries]
    require(len(set(ids)) == len(ids) and len(set(qids)) == len(qids), 'unique identifiers')
    feasible = {q['id']:q for q in queries}
    predictions = {}
    baseK = [[rational(x) for x in row] for row in candidates[0]['K']]
    baseC = [[rational(x) for x in row] for row in candidates[0]['C']]
    def transform(A,Q):
        return [[sum(Q[k][i]*A[k][ell]*Q[ell][j] for k in range(3) for ell in range(3))
                 for j in range(3)] for i in range(3)]
    for i,candidate in enumerate(candidates):
        K = [[rational(x) for x in row] for row in candidate['K']]
        C = [[rational(x) for x in row] for row in candidate['C']]
        half = rational(candidate['half_angle_parameter'])
        require(half == F(i,100), 'candidate half-angle generation')
        cs,sn = (1-half*half)/(1+half*half),2*half/(1+half*half)
        Q = [[F(1),F(0),F(0)],[F(0),cs,sn],[F(0),-sn,cs]]
        require(transform(baseK,Q)==K and transform(baseC,Q)==C,'candidate exact matrix similarities')
        baseline=solve(K,[F(1),F(0),F(0)])[0]
        require(baseline==F(5,18)==rational(candidate['baseline_response']),'identical baseline')
        for j in range(3):
            require(sum(K[j])>=0 and sum(C[j])>=0,'candidate nonnegative grounds')
            for k in range(3):
                require(K[j][k]==K[k][j] and C[j][k]==C[k][j], 'candidate reciprocity')
                if j!=k:require(K[j][k]<=0 and C[j][k]<=0,'candidate nonnegative edge coefficients')
        predictions[candidate['id']]={}
        for query in queries:
            A=[row[:] for row in K]
            for node,load in query['loads'].items():
                require(str(int(node))==node and 0<=int(node)<3 and rational(load)>0,'known positive load')
                A[int(node)][int(node)]+=rational(load)
            predicted=solve(A,[F(1),F(0),F(0)])[0]
            require(predicted==rational(candidate['predictions'][query['id']]),'independent predicted response')
            predictions[candidate['id']][query['id']]=predicted

    def compatible(observations,eta):
        return [cid for cid in ids if all(abs(predictions[cid][qid]-value)<=eta for qid,value in observations)]
    def score(remaining,qid,eta):
        values=sorted(predictions[cid][qid] for cid in remaining)
        endpoints={v-eta for v in values}|{v+eta for v in values}
        worst=max(sum(abs(v-y)<=eta for v in values) for y in endpoints)
        separated=sum(values[j]-values[i]>2*eta for i in range(len(values)) for j in range(i+1,len(values)))
        minimum=min(values[j]-values[j-1] for j in range(1,len(values)))
        return (-worst,separated,minimum,values[-1]-values[0])

    observation_count=0
    seen_trials=set()
    by_group={}
    for trial in trials:
        policy,truth,eta=trial['policy'],trial['true_candidate'],rational(trial['eta'])
        require(policy in ('fixed','fixed_best_first','active') and truth in ids and eta>=0,'trial contract')
        key=(policy,truth,eta)
        require(key not in seen_trials,'duplicate trial')
        seen_trials.add(key)
        remaining=ids[:]; observations=[]; used=set()
        ranked=sorted(qids,key=lambda qid:score(ids,qid,eta),reverse=True) if policy=='fixed_best_first' else qids
        require(len(trial['trace'])<=4,'query budget')
        for index,step in enumerate(trial['trace']):
            observation_count+=1
            require(len(remaining)>1,'query after identification')
            query=step['query'];qid=query['id']
            require(qid in feasible and query==feasible[qid] and qid not in used,'feasible unused intervention')
            if policy=='active':
                best=max((q for q in qids if q not in used),key=lambda q:score(remaining,q,eta))
                selected_score=score(remaining,best,eta)
                require(qid==best and -selected_score[0]<len(remaining),'active minimax query selection')
                require(step['planner_worst_case_survivors']==-selected_score[0],'stored minimax bound')
            else:
                require(qid==next(q for q in ranked if q not in used),'fixed query order')
                require(step['planner_worst_case_survivors'] is None,'fixed trace metadata')
            value=rational(step['observed_response']);actual=predictions[truth][qid]
            error=value-actual
            require(abs(error)<=eta and error==rational(step['oracle_error']),'absolute oracle error bound')
            lo,hi=actual-eta,actual+eta
            endpoints={lo,hi,actual}
            for cid in remaining:
                prediction=predictions[cid][qid]
                endpoints.update(x for x in (prediction-eta,prediction+eta) if lo<=x<=hi)
            survivor_count=lambda y:sum(abs(predictions[cid][qid]-y)<=eta for cid in remaining)
            greatest=max(survivor_count(y) for y in endpoints)
            expected_observation=min(y for y in endpoints if survivor_count(y)==greatest)
            require(value==expected_observation,'bounded adversary maximizes survivors with declared tie break')
            observations.append((qid,value));used.add(qid)
            remaining=compatible(observations,eta)
            require(truth in remaining,'true candidate retained')
            require(step['step']==index+1 and step['survivors']==remaining,'full survivor-history replay')
        require(trial['intervention_queries_used']==len(trial['trace']),'stored query count')
        require(trial['survivors']==remaining,'final survivor list')
        identified=len(remaining)==1
        require(trial['status']==('identified' if identified else 'unresolved'),'final verdict')
        require(trial['identified_correctly'] is (identified and remaining[0]==truth),'identification truth')
        by_group.setdefault((eta,policy),[]).append(trial)
    for summary in document['summaries']:
        eta=rational(summary['absolute_componentwise_response_error'])
        for policy,statistics in summary['policies'].items():
            records=by_group[eta,policy]
            require(len(records)==21,'complete truth population per policy/error')
            require(statistics['identified_within_budget']==sum(r['identified_correctly'] for r in records),'summary identifications')
            require(statistics['unresolved']==sum(r['status']=='unresolved' for r in records),'summary unresolved')
            require(statistics['wrong_identifications']==0,'summary wrong decisions')
            require(rational(statistics['mean_queries_used'])==F(sum(len(r['trace']) for r in records),21),'summary mean query count')
            for row in statistics['by_equal_budget']:
                budget=row['intervention_budget']
                sizes=[len(r['trace'][min(budget,len(r['trace']))-1]['survivors']) if r['trace'] else len(ids) for r in records]
                require(row['identified']==sum(x==1 for x in sizes),'equal-budget identification count')
                require(rational(row['mean_survivors'])==F(sum(sizes),21),'equal-budget mean survivors')
                require(row['maximum_survivors']==max(sizes),'equal-budget largest survivor set')
    require(len(by_group)==9,'all policy/error combinations')
    return {'format':'nonagon-independent-active-receipts/1','status':'CERTIFIED','valid':True,
            'checks':count,'candidate_predictions':len(candidates)*len(queries),'trials':len(trials),
            'intervention_observations':observation_count,'wrong_identifications':0,
            'scope':'Independent rational Gaussian solves, exact candidate similarities, bounded oracle errors and worst-case choices, full survivor histories, active/fixed query choices, and all retained summary counts. Finite declared library only.',
            'optimized_python':not __debug__}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:report=check_campaign(json.loads(args.campaign.read_text()))
    except (ValueError,KeyError,TypeError,ZeroDivisionError) as exc:
        report={'status':'REJECTED','valid':False,'reason':str(exc)}
    text=json.dumps(report,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')
    raise SystemExit(0 if report['valid'] else 1)


if __name__=='__main__':main()
