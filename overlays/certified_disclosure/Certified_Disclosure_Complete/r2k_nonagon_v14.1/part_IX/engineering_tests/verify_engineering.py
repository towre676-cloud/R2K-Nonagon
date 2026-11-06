"""Independent checks of the Nonagon engineering test records.

The verifier does not import the tomography, interval inference, SPICE harness,
fault-planning or instrument-campaign implementations. It uses a separate exact
forward-elimination/back-substitution solver for circuit predictions and checks
the retained engineering truth and acquisition envelopes. Simulation records do
not establish calibration or performance of physical hardware.
"""
from __future__ import annotations

import argparse
import csv
import ast
import json
import hashlib
import math
import re
from itertools import combinations
from fractions import Fraction as F
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent


def core_path():
    return HERE.parent if (HERE.parent/'hardware_plan.py').exists() else HERE.parent/'part8_work'/'part_VIII'


class CheckFailure(ValueError):
    pass


class Checks:
    def __init__(self):
        self.count = 0

    def require(self, condition, message):
        self.count += 1
        if not condition:
            raise CheckFailure(message)


def q(value):
    if isinstance(value, bool):
        raise CheckFailure('Boolean is not a circuit quantity')
    return F(str(value)) if isinstance(value, float) else F(value)


def gaussian(matrix, right):
    """Exact forward elimination plus reverse substitution; no inverse routine."""
    n = len(matrix)
    if n == 0 or len(right) != n or any(len(row) != n for row in matrix):
        raise CheckFailure('Invalid independent linear system')
    a = [[q(x) for x in row] for row in matrix]
    b = [q(x) for x in right]
    for k in range(n):
        pivot = next((r for r in range(k, n) if a[r][k] != 0), None)
        if pivot is None:
            raise CheckFailure('Independent circuit system is singular')
        a[k], a[pivot] = a[pivot], a[k]
        b[k], b[pivot] = b[pivot], b[k]
        for i in range(k+1, n):
            if a[i][k]:
                factor = a[i][k] / a[k][k]
                for j in range(k+1, n):
                    a[i][j] -= factor * a[k][j]
                b[i] -= factor * b[k]
                a[i][k] = F(0)
    x = [F(0)] * n
    for i in reversed(range(n)):
        x[i] = (b[i]-sum(a[i][j]*x[j] for j in range(i+1, n))) / a[i][i]
    return x


def response(matrix, loads=None, source=0, readout=0):
    a = [[q(x) for x in row] for row in matrix]
    for j, amount in (loads or {}).items():
        a[int(j)][int(j)] += q(amount)
    return gaussian(a, [F(i == source) for i in range(len(a))])[readout]


def full_inverse(matrix):
    n = len(matrix)
    columns = [gaussian(matrix, [F(i == j) for i in range(n)]) for j in range(n)]
    return [[columns[j][i] for j in range(n)] for i in range(n)]


def within(interval, value):
    return q(interval[0]) <= q(value) <= q(interval[1])


def read_json(path):
    return json.loads(Path(path).read_text())


def check_interval_result(checks, matrix, result, prefix):
    """Check retained truth containment and the independent inversion witness."""
    n = len(matrix)
    checks.require(result.get('status') in ('enclosed', 'unresolved'), prefix+': known bounded status')
    if result['status'] != 'enclosed':
        return {'status': 'unresolved'}
    gtrue = full_inverse(matrix)
    gb = result['G_intervals']; kb = result['K_intervals']
    checks.require(len(gb) == n and len(kb) == n, prefix+': matrix order')
    for i in range(n):
        checks.require(len(gb[i]) == n and len(kb[i]) == n, prefix+': row order')
        for j in range(n):
            checks.require(within(gb[i][j], gtrue[i][j]), prefix+': compliance truth contained')
            checks.require(within(kb[i][j], matrix[i][j]), prefix+': conductance truth contained')
    witness = result['inverse_witness']
    mid = [[(q(v[0])+q(v[1]))/2 for v in row] for row in gb]
    rad = [[(q(v[1])-q(v[0]))/2 for v in row] for row in gb]
    independent = full_inverse(mid)
    supplied = [[q(x) for x in row] for row in witness['midpoint_inverse']]
    checks.require(supplied == independent, prefix+': independent midpoint inverse')
    norm = max(sum(abs(x) for x in row) for row in independent)
    enorm = max(sum(row) for row in rad)
    rho = norm*enorm
    checks.require(q(witness['midpoint_inverse_norm_inf']) == norm, prefix+': inverse norm')
    checks.require(q(witness['compliance_radius_norm_inf']) == enorm, prefix+': radius norm')
    checks.require(q(witness['rho']) == rho and rho < 1, prefix+': Neumann margin')
    error = norm*rho/(1-rho)
    checks.require(q(witness['stiffness_entry_error']) == error, prefix+': Neumann error')
    for edge in result['edges']:
        i, j = edge['sites_zero_based']; weight = -q(matrix[i][j])
        checks.require(within(edge['weight_interval'], weight), prefix+': edge truth contained')
        if edge['status'] == 'present':
            checks.require(weight > 0, prefix+': present edge is present')
        elif edge['status'] == 'absent_under_gap':
            checks.require(weight == 0, prefix+': absent edge is absent')
        elif edge['status'] == 'inconsistent_with_gap':
            gap = q(result['contract']['minimum_positive_edge_weight'])
            checks.require(0 < weight < gap, prefix+': gap contradiction is real')
        else:
            checks.require(edge['status'] == 'unresolved', prefix+': valid edge status')
    return {'status': 'enclosed'}


def plan_self_checks(checks):
    """Independently audit the retained v8 hardware-plan scaling and costs."""
    old = core_path()/'hardware_plan.py'
    for n in (3, 6):
        with (old.parent/'data'/f'hardware_{n}node_acquisition.csv').open() as f:
            rows=list(csv.DictReader(f))
        checks.require(len(rows) == n*(n+1)//2, 'hardware canonical configuration count')
        seen = set()
        for row in rows:
            checks.require(row['source_pad'] == row['sense_pad'] == 'node_1', 'hardware fixed port')
            checks.require(q(row['current_A']) == F(1, 10000), 'hardware default current')
            checks.require(q(row['conductance_scale_S']) == F(1, 1000), 'hardware conductance scale')
            loads = []
            for label in ('A', 'B'):
                node = row['shunt_'+label+'_pad']
                resistor = row['shunt_'+label+'_resistance_ohm']
                if node:
                    checks.require(node != 'node_1' and 1 < int(node[5:]) <= n, 'hardware known hidden pad')
                    checks.require(q(resistor) in (F(1000), F(500)), 'hardware normalized unit shunt')
                    loads.append((node, q(resistor)))
                else:
                    checks.require(not resistor, 'hardware unused channel disconnected')
            key = tuple(sorted(loads))
            checks.require(key not in seen, 'hardware configurations unique'); seen.add(key)
            checks.require(len(loads) <= 2, 'hardware at most two active shunts')
    # Exhaust all voltage/current calibration corners without trusting the helper's algebra.
    row = {'readout_V':'1/40','voltage_absolute_error_V':'1/1000000',
           'current_A':'1/10000','current_relative_error':'1/10000',
           'conductance_scale_S':'1/1000'}
    # Extract only the unit-normalization helper, avoiding any tomography import.
    tree=ast.parse(old.read_text())
    function=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='voltage_interval')
    isolated=ast.Module(body=[function],type_ignores=[])
    namespace={'F':F};exec(compile(isolated,str(old),'exec'),namespace)
    interval=namespace['voltage_interval'](row)
    for sv in (-1, 1):
        for si in (-1, 1):
            v = F(1,40)+sv*F(1,1000000)
            current = F(1,10000)*(1+si*F(1,10000))
            checks.require(within(interval, v*F(1,1000)/current), 'V/I calibration corner containment')
    # Separate the resistor setting resolution from its calibration tolerance.
    for resistance in (F(250), F(500), F(1000), F(4000)):
        rel = F(8,10000); relay = F(1,20)
        nominal_g = 1/resistance
        endpoints = [1/(resistance*(1-rel)), 1/(resistance*(1+rel)+relay)]
        envelope = max(abs(x/nominal_g-1) for x in endpoints)
        for rr in (resistance*(1-rel), resistance*(1+rel)):
            for rs in (F(0), relay):
                checks.require(abs((1/(rr+rs))/nominal_g-1) <= envelope,
                               'relay-plus-shunt effective conductance containment')
    return {'status':'passed','orders':[3,6],
            'scope':'Plan/scaling checks only; no instrument calibration established'}


def instrument_replay(checks,qualification=False):
    stem='instrument_qualification' if qualification else 'instrument'
    summary=read_json(HERE/(stem+'_results.json'))
    protocol_path=HERE/(stem+'_protocol.json')
    checks.require(hashlib.sha256(protocol_path.read_bytes()).hexdigest()==summary['protocol_sha256'],
                   'instrument frozen protocol identity')
    fixtures={p.stem:read_json(p) for p in (HERE/'instrument_fixtures').glob('*.json')}
    records=0;enclosed=0;packets=0;full_graph=0
    for line in (HERE/('instrument_qualification_records.jsonl' if qualification else 'instrument_case_records.jsonl')).open():
        record=json.loads(line);records+=1
        fixture=fixtures[record['fixture_id']]
        k=[[q(x) for x in row] for row in fixture['K_normalized']]
        n=len(k);budget=record['budget'];scale=q(record.get('conductance_scale_S',budget.get('conductance_scale_S')))
        effective=q(record['case']['effective_relative_load_error'])
        nominal_current=q(budget['current_A']);rho=q(budget['current_relative_error'])
        absolute_error=q(budget['voltage_absolute_error_V'] if qualification else budget['voltage_error_V'])
        drift=q(budget['drift_error_V'])
        shunt_error=q(budget['shunt_relative_error']);rmax=q(budget['relay_series_max_ohm'])
        for group in ('readings','repeat_baselines','heldout_readings'):
            for row in record[group]:
                packets+=1
                nominal={int(j):q(x) for j,x in row['nominal_loads'].items()}
                actual={int(j):q(x) for j,x in row['actual_loads'].items()}
                checks.require(len(actual)<=2,'instrument two routed shunts suffice including validation')
                checks.require(set(nominal)==set(actual),'instrument nominal/actual node set')
                for branch in row['branches']:
                    j=branch['site'];t=q(branch['nominal_load']);e=q(branch['shunt_relative_deviation'])
                    r=q(branch['shunt_resistance_ohm']);relay=q(branch['relay_series_ohm'])
                    checks.require(abs(e)<=shunt_error,'instrument shunt calibration bound')
                    checks.require(0<=relay<=rmax,'instrument relay bound')
                    checks.require(r==1/(t*scale*(1+e)),'instrument shunt resistance physical law')
                    checks.require(actual[j]==1/((r+relay)*scale),'instrument effective shunt conductance')
                    checks.require(q(branch['actual_load'])==actual[j],'instrument branch serialization')
                    # Effective canonical envelope is designed for strengths <=2t.
                    if group=='readings':
                        checks.require(abs(actual[j]/nominal[j]-1)<=effective,'instrument inferred load envelope')
                current=q(row['actual_current_A'])
                checks.require(nominal_current*(1-rho)<=current<=nominal_current*(1+rho),'instrument current envelope')
                h=response(k,actual)
                checks.require(h==q(row['true_normalized_response']),'independent actual circuit response')
                ideal=current*h/scale
                checks.require(ideal==q(row['ideal_voltage_V']),'instrument volts/current/conductance dimensions')
                noise=q(row['additive_readout_error_V']);d=q(row['additive_drift_V'])
                gain=q(row.get('additive_gain_error_V','0'))
                if qualification:
                    checks.require(gain==ideal*q(row['gain_relative_deviation']),'instrument gain error identity')
                checks.require(abs(noise)<=absolute_error and abs(d)<=drift,'instrument readout and drift envelopes')
                measured=ideal+noise+d+gain
                checks.require(measured==q(row['measured_voltage_V']),'instrument measured voltage serialization')
                verror=absolute_error+(q(budget['voltage_relative_to_reading_error'])*abs(measured) if qualification else 0)
                if qualification:
                    checks.require(q(row['declared_voltage_error_V'])==verror,'qualification measured-reading relative error rule')
                    checks.require(abs(noise+gain)<=verror,'qualification total readout error envelope')
                expected=[(measured-verror-drift)*scale/(nominal_current*(1+rho)),
                          (measured+verror+drift)*scale/(nominal_current*(1-rho))]
                checks.require([q(x) for x in row['normalized_interval']]==expected,'independent V/I normalization')
                checks.require(within(expected,h),'instrument actual response is enclosed')
        # Reconstruct only the observer's interval map, independently intersecting
        # baseline repeats; this does not rerun the interval inverse algorithm.
        expected_map={tuple(sorted((int(j),q(x)) for j,x in row['nominal_loads'].items())):
                      [q(x) for x in row['normalized_interval']] for row in record['readings']}
        for row in record['repeat_baselines']:
            old=expected_map[()];new=[q(x) for x in row['normalized_interval']]
            expected_map[()]=[max(old[0],new[0]),min(old[1],new[1])]
        for row in record['inference_observation_intervals']:
            key=tuple(sorted((int(j),q(x)) for j,x in row['nominal_loads'].items()))
            checks.require([q(x) for x in row['interval']]==expected_map[key],'observer receives retained interval only')
        result=check_interval_result(checks,k,record['bounds'],record['fixture_id'])
        enclosed+=result['status']=='enclosed'
        graph=record['bounds'].get('edges',[])
        disclosed=record['bounds']['status']=='enclosed' and all(e['status'] in ('present','absent_under_gap') for e in graph)
        full_graph+=disclosed
        checks.require(record['case']['full_graph_disclosed']==disclosed,'instrument full graph label')
        cost=record['case']['cost']
        checks.require(cost['canonical_readings']==n*(n+1)//2==len(record['readings']),'instrument acquisition cost')
        checks.require(cost['repeat_baselines']==len(record['repeat_baselines']),'instrument repeat cost')
        checks.require(cost['validation_readings']==len(record['heldout_readings'])==3,'instrument validation cost')
        checks.require(cost['total_voltage_reads']==sum(len(record[g]) for g in ('readings','repeat_baselines','heldout_readings')),
                       'instrument total scalar read cost')
    checks.require(records==summary['case_count']==(72 if qualification else 360),'instrument complete case count')
    return {'status':'passed','cases':records,'measurement_packets':packets,'enclosed_matrices':enclosed,
            'full_graph_cases':full_graph,'physical_hardware_measurements':0,
            'scope':'All retained actual circuit responses, input envelopes, inversion witnesses and truth containment replayed'}


def spice_replay(checks):
    doc=read_json(HERE/'spice_results.json');g=q(doc['conductance_scale_S']);current=q(doc['source_current_A'])
    count=0;voltages=0;max_error=0.0
    for case in doc['dc']['cases']:
        k=[[q(x) for x in row] for row in case['exact_K_normalized']];n=len(k)
        checks.require(case['acquisition_configurations']==n*(n+1)//2,'SPICE acquisition configuration count')
        for row in case['observations']:
            count+=1;loads=row['loads']
            checks.require(len(loads)<=2,'SPICE two-shunt hardware feasibility')
            a=[r[:] for r in k]
            for j,x in loads.items():a[int(j)][int(j)]+=q(x)
            expected=gaussian(a,[current/g]+[F(0)]*(n-1))
            checks.require([q(x) for x in row['exact_voltage_V']]==expected,'SPICE serialized exact voltages match independent solve')
            checks.require(q(row['exact_normalized_response'])==expected[0]*g/current,'SPICE response normalization')
            netlist=HERE/row['netlist'];raw=HERE/row['raw_log']
            checks.require(hashlib.sha256(netlist.read_bytes()).hexdigest()==row['netlist_sha256'],'SPICE netlist identity')
            matches=re.findall(r'v\(n(\d+)\)\s*=\s*([-+\deE.]+)',raw.read_text())
            parsed={int(j)-1:float(x) for j,x in matches}
            checks.require(len(parsed)==n,'SPICE raw log contains each node')
            for j,(observed,truth) in enumerate(zip(row['ngspice_voltage_V'],expected)):
                voltages+=1
                checks.require(float(observed)==parsed[j],'SPICE observation comes from retained external solver log')
                error=abs(float(observed)-float(truth));max_error=max(max_error,error)
                checks.require(error<=1e-11,'external SPICE voltage agrees with independent circuit solve')
    checks.require(count==doc['dc']['observation_count'],'SPICE total scalar configurations')
    # Independent transient check uses a matrix exponential at a sparse set of
    # retained times; no SPICE harness formula or eigenbasis is imported.
    import numpy as np
    from scipy.linalg import expm
    samples=0;max_transient_error=0.0
    for row in doc['transient']['cases']:
        k=np.array([[float(q(x)) for x in r] for r in row['exact_K_normalized']]);n=len(k)
        steady=np.array([float(x) for x in gaussian(row['exact_K_normalized'],[q(row['current_A'])/q(row['gscale_S'])]+[F(0)]*(n-1))])
        trace=np.loadtxt(HERE/row['raw_trace'],skiprows=1)
        checks.require(trace.shape[1]==n+1 and len(trace)==row['trace_rows'],'SPICE trace dimensions')
        indices=sorted(set(np.linspace(0,len(trace)-1,41,dtype=int)))
        for index in indices:
            t=trace[index,0]
            expected=steady-expm(-float(q(row['gscale_S']))*k*t/float(row['capacitance_each_F']))@steady
            error=float(np.max(np.abs(trace[index,1:]-expected)));max_transient_error=max(max_transient_error,error);samples+=1
            checks.require(error<=5e-7,'SPICE RC trace agrees with independent expm')
        tolerance=float(row['port_settling_tolerance_V'])
        later=trace[:,0]>=float(row['simulated_port_settling_time_s'])
        checks.require(float(np.max(np.abs(trace[later,1]-steady[0])))<=tolerance,'SPICE retained settling time stays within tolerance')
    return {'status':'passed','DC_configurations':count,'node_voltages':voltages,'max_voltage_error_V':max_error,
            'transient_samples':samples,'max_sampled_transient_error_V':max_transient_error,
            'physical_hardware_measurements':0}


def verify_fault_summaries(checks,doc):
    aggregated={}
    for trial in doc['trials']:
        for budget,checkpoint in trial['checkpoints'].items():
            key=(trial['eta'],trial['policy'],int(budget))
            row=aggregated.setdefault(key,{'cases':0,'identified':0,'unresolved':0,'wrong':0,
                                          'survivor_sum':0,'all_error_box_identifiable':0,'by_fault_type':{}})
            row['cases']+=1;row[checkpoint['status']]+=1;row['survivor_sum']+=checkpoint['survivors']
            row['wrong']+=checkpoint['identified'] is not None and checkpoint['identified']!=trial['true_candidate']
            row['all_error_box_identifiable']+=checkpoint.get('all_errors_identifiable',False)
            typ=row['by_fault_type'].setdefault(trial['fault_type'],{'cases':0,'identified':0,'all_error_box_identifiable':0})
            typ['cases']+=1;typ['identified']+=checkpoint['status']=='identified'
            typ['all_error_box_identifiable']+=checkpoint.get('all_errors_identifiable',False)
    seen=set()
    for row in doc['summaries']:
        key=(row['eta'],row['policy'],row['budget']);checks.require(key not in seen,'fault summary keys unique');seen.add(key)
        expected=aggregated[key]
        for field in ('cases','identified','unresolved','wrong','survivor_sum','all_error_box_identifiable'):
            if field in row:checks.require(row[field]==expected[field],'fault independently aggregated summary '+field)
        checks.require(set(row['by_fault_type'])==set(expected['by_fault_type']),'fault category summary keys')
        for typ,values in row['by_fault_type'].items():
            for field,value in values.items():checks.require(value==expected['by_fault_type'][typ][field],'fault category summary value')
    checks.require(seen==set(aggregated),'fault summary covers all checkpoints')
    return len(seen)


def fault_replay(checks):
    doc=read_json(HERE/'fault_results.json');queries=doc['query_definitions']
    checks.require(hashlib.sha256((HERE/'fault_protocol.json').read_bytes()).hexdigest()==doc['protocol_sha256'],
                   'fault frozen protocol identity')
    fixtures={f['id']:f for f in doc['fixtures']};tables={};predictions=0
    for fi,fixture in enumerate(doc['fixtures']):
        table={key:[] for key in queries}
        for candidate_index,candidate in enumerate(fixture['catalog']):
            k=candidate['K'];inverse_cache=full_inverse(k)
            for key,query in queries.items():
                if query['kind']=='scan':value=inverse_cache[query['source']][query['readout']]
                else:value=response(k,query['loads'],query['source'],query['readout'])
                table[key].append(value);predictions+=1
                checks.require(value==q(fixture['predictions'][key][candidate_index]),'independent fault candidate prediction')
        tables[fixture['id']]=table
        if (fi+1)%5==0:print(json.dumps({'verification_progress':'fault_predictions','fixtures':fi+1}),flush=True)
        for eta_text,plans in fixture['plans'].items():
            eta=q(eta_text)
            for policy,plan in plans.items():
                sequence=plan['queries']
                checks.require(len(sequence)==len(set(sequence)),'fault plan does not duplicate query budget')
                for budget_text,cost in plan['costs_at_budget'].items():
                    subset=sequence[:int(budget_text)];definitions=[queries[key] for key in subset]
                    expected={'scalar_readings':len(subset),'source_addresses':len({r['source'] for r in definitions}),
                              'readout_addresses':len({r['readout'] for r in definitions}),
                              'distinct_shunt_pads':len({j for r in definitions for j in r['loads']}),
                              'maximum_simultaneous_shunts':max([len(r['loads']) for r in definitions],default=0),
                              'configuration_count':len(subset),'source_block_minimum':len({r['source'] for r in definitions})}
                    for key,value in expected.items():checks.require(cost[key]==value,'fault fair scalar/address/configuration cost')
                if policy.startswith('greedy_fixed'):
                    pool=[key for key,r in queries.items() if (r['kind']=='scan')==('reciprocal' in policy)]
                    pairs=list(combinations(range(len(fixture['catalog'])),2));masks={}
                    for key in pool:
                        mask=0
                        for index,(i,j) in enumerate(pairs):
                            if abs(table[key][i]-table[key][j])>2*eta:mask|=1<<index
                        masks[key]=mask
                    remaining=(1<<len(pairs))-1;available=pool[:]
                    for actual in sequence:
                        expected=max(available,key=lambda key:((masks[key]&remaining).bit_count(),masks[key].bit_count(),-pool.index(key)))
                        checks.require(actual==expected,'fault frozen greedy plan uses catalog only')
                        available.remove(actual);remaining&=~masks[actual]
    measurements=0;checkpoints=0;wrong=0
    for trial in doc['trials']:
        fixture=fixtures[trial['fixture']];table=tables[trial['fixture']]
        ids=[c['id'] for c in fixture['catalog']];true=ids.index(trial['true_candidate']);eta=q(trial['eta'])
        checks.require(fixture['catalog'][true]['type']==trial['fault_type'],'fault specimen type')
        survivors=list(range(len(ids)))
        sequence=fixture['plans'][trial['eta']][trial['policy']]['queries']
        checks.require([r['query'] for r in trial['readings']]==sequence,'fault trial uses frozen query plan')
        for index,row in enumerate(trial['readings'],1):
            measurements+=1;values=table[row['query']];actual=values[true];measured=q(row['measured'])
            checks.require(actual==q(row['actual']) and measured-actual==q(row['error']),'fault actual signal/error identity')
            checks.require(abs(measured-actual)<=eta,'fault scalar error budget')
            survivors=[j for j in survivors if abs(values[j]-measured)<=eta]
            checks.require([ids[j] for j in survivors]==row['survivors'],'fault compatible candidate history')
            checks.require(true in survivors,'fault truth never excluded')
            checkpoint=trial['checkpoints'].get(str(index))
            if checkpoint is not None:
                checkpoints+=1;identified=ids[survivors[0]] if len(survivors)==1 else None
                checks.require(checkpoint['survivors']==len(survivors),'fault checkpoint survivor count')
                checks.require(checkpoint['identified']==identified,'fault checkpoint diagnosis')
                checks.require(checkpoint['status']==('identified' if identified else 'unresolved'),'fault checkpoint status')
                wrong+=identified is not None and identified!=trial['true_candidate']
    checks.require(len(doc['trials'])==doc['trial_count']==2000,'fault trial count')
    checks.require(wrong==0,'fault no incorrect accepted diagnosis')
    for record in doc['out_of_catalog_controls']:
        fixture=fixtures[record['fixture']];table=tables[record['fixture']]
        ids=[c['id'] for c in fixture['catalog']];survivors=list(range(len(ids)))
        for row in record['readings']:
            measured=response(record['K'],queries[row['query']]['loads'])
            checks.require(measured==q(row['measured']),'independent out-of-catalog circuit response')
            survivors=[j for j in survivors if table[row['query']][j]==measured]
            checks.require([ids[j] for j in survivors]==row['survivors'],'out-of-catalog survivors')
        checks.require(not survivors,'out-of-catalog control rejected')
    summary_tables=verify_fault_summaries(checks,doc)
    return {'status':'passed','candidate_predictions':predictions,'trials':len(doc['trials']),
            'measurement_decisions':measurements,'checkpoint_decisions':checkpoints,'wrong_identifications':wrong,
            'out_of_catalog_controls':len(doc['out_of_catalog_controls']),
            'summary_tables':summary_tables,
            'scope':'Complete exact candidate prediction, frozen-plan, cost and survivor replay'}


def calibration_fault_replay(checks):
    doc=read_json(HERE/'fault_calibration_results.json');queries=doc['query_definitions']
    protocol_path=HERE/'fault_calibration_protocol.json';protocol=read_json(protocol_path)
    checks.require(hashlib.sha256(protocol_path.read_bytes()).hexdigest()==doc['protocol_sha256'],
                   'calibrated-fault frozen protocol identity')
    scale=q(protocol['normalized_conductance_scale_S']);rrel=q(protocol['resistor_relative_error'])
    rabs=q(protocol['resistor_absolute_error_ohm']);relay_low,relay_high=map(q,protocol['relay_series_resistance_ohm'])
    fixtures={f['id']:f for f in doc['fixtures']};tables={};predictions=0
    bounds_by_query={}
    for key,query in queries.items():
        lower={};upper={}
        for j,t in query['loads'].items():
            nominal=1/(q(t)*scale);error=rrel*nominal+rabs
            lower[j]=1/((nominal+error+relay_high)*scale)
            upper[j]=1/((nominal-error+relay_low)*scale)
        bounds_by_query[key]=(lower,upper)
    for fi,fixture in enumerate(doc['fixtures']):
        table={key:[] for key in queries}
        for ci,candidate in enumerate(fixture['catalog']):
            k=candidate['K'];g=full_inverse(k)
            for key,query in queries.items():
                if query['kind']=='scan':lo=hi=g[query['source']][query['readout']]
                else:
                    lower,upper=bounds_by_query[key]
                    lo=response(k,upper);hi=response(k,lower)
                table[key].append((lo,hi));predictions+=1
                saved=doc['fixtures'][fi]['prediction_intervals'][key][ci]
                checks.require((lo,hi)==tuple(map(q,saved)),'independent calibrated candidate prediction interval')
        tables[fixture['id']]=table
        if (fi+1)%5==0:print(json.dumps({'verification_progress':'calibrated_fault_predictions','fixtures':fi+1}),flush=True)
        for eta_text,plans in fixture['plans'].items():
            eta=q(eta_text)
            for policy,plan in plans.items():
                sequence=plan['queries']
                for budget_text,cost in plan['costs_at_budget'].items():
                    subset=sequence[:int(budget_text)];definitions=[queries[key] for key in subset]
                    costs={'scalar_readings':len(subset),'source_addresses':len({r['source'] for r in definitions}),
                           'readout_addresses':len({r['readout'] for r in definitions}),
                           'distinct_shunt_pads':len({j for r in definitions for j in r['loads']}),
                           'maximum_simultaneous_shunts':max([len(r['loads']) for r in definitions],default=0),
                           'configuration_count':len(subset),'source_block_minimum':len({r['source'] for r in definitions})}
                    for field,value in costs.items():checks.require(cost[field]==value,'calibrated-fault resource costs')
                if not policy.startswith('greedy_fixed'):continue
                pool=[key for key,r in queries.items() if (r['kind']=='scan')==('reciprocal' in policy)]
                pairs=list(combinations(range(len(fixture['catalog'])),2));masks={}
                for key in pool:
                    mask=0
                    for index,(i,j) in enumerate(pairs):
                        a,b=table[key][i],table[key][j]
                        if max(a[0]-b[1],b[0]-a[1])>2*eta:mask|=1<<index
                    masks[key]=mask
                remaining=(1<<len(pairs))-1;available=pool[:]
                for actual in sequence:
                    expected=max(available,key=lambda key:((masks[key]&remaining).bit_count(),masks[key].bit_count(),-pool.index(key)))
                    checks.require(actual==expected,'calibrated-fault frozen catalog-only plan')
                    available.remove(actual);remaining&=~masks[actual]
    measurements=0;checkpoints=0;robust_positive=0;actual_cache={}
    for trial in doc['trials']:
        fixture=fixtures[trial['fixture']];table=tables[trial['fixture']]
        ids=[c['id'] for c in fixture['catalog']];true=ids.index(trial['true_candidate']);eta=q(trial['eta'])
        candidate=fixture['catalog'][true]
        checks.require(candidate['type']==trial['fault_type'],'calibrated fault specimen category')
        survivors=list(range(len(ids)));sequence=fixture['plans'][trial['eta']][trial['policy']]['queries']
        checks.require([r['query'] for r in trial['readings']]==sequence,'calibrated-fault trial uses frozen plan')
        for index,row in enumerate(trial['readings'],1):
            measurements+=1;key=row['query'];query=queries[key];lower,upper=bounds_by_query[key]
            loads=row['actual_loads'];checks.require(set(loads)==set(query['loads']),'calibrated actual load sites')
            for j,t in loads.items():checks.require(lower[j]<=q(t)<=upper[j],'calibrated actual shunt inside physical range')
            cache_key=(trial['fixture'],trial['true_candidate'],key)
            if cache_key not in actual_cache:
                actual_cache[cache_key]=(dict(loads),response(candidate['K'],loads,query['source'],query['readout']))
            stored_loads,actual=actual_cache[cache_key]
            checks.require(stored_loads==loads,'physical query shared fairly across policies/error levels')
            checks.require(actual==q(row['actual']),'independent calibrated physical signal')
            lo,hi=table[key][true]
            checks.require(lo<=actual<=hi,'calibrated signal inside candidate box')
            measured=q(row['measured']);checks.require(abs(measured-actual)<=eta and q(row['error'])==measured-actual,
                                                     'calibrated signal readout bound')
            survivors=[j for j in survivors if table[key][j][0]-eta<=measured<=table[key][j][1]+eta]
            checks.require([ids[j] for j in survivors]==row['survivors'],'calibrated survivor history')
            checks.require(true in survivors,'calibrated true candidate retained')
            checkpoint=trial['checkpoints'].get(str(index))
            if checkpoint is None:continue
            checkpoints+=1;identified=ids[survivors[0]] if len(survivors)==1 else None
            checks.require(checkpoint['identified']==identified and checkpoint['survivors']==len(survivors),'calibrated checkpoint diagnosis')
            checks.require(checkpoint['status']==('identified' if identified else 'unresolved'),'calibrated checkpoint status')
            alternatives=[]
            for other in range(len(ids)):
                if other==true:continue
                gap=max(max(table[qid][true][0]-table[qid][other][1],table[qid][other][0]-table[qid][true][1]) for qid in sequence[:index])-2*eta
                alternatives.append((gap,ids[other]))
            margin,nearest=min(alternatives)
            checks.require(q(checkpoint['margin'])==margin and checkpoint['nearest_alternative']==nearest,'independent all-error-box margin')
            checks.require(checkpoint['all_errors_identifiable']==(margin>0),'all-error-box diagnosis separated from selected noise outcome')
            if margin>0:
                robust_positive+=1
                checks.require(identified==trial['true_candidate'],'positive box margin implies actual correct identification')
    checks.require(len(doc['trials'])==doc['trial_count']==1500,'calibrated-fault complete trial count')
    summaries=verify_fault_summaries(checks,doc)
    return {'status':'passed','candidate_prediction_intervals':predictions,'trials':len(doc['trials']),
            'measurement_decisions':measurements,'checkpoint_decisions':checkpoints,
            'all_error_box_positive_checkpoints':robust_positive,'summary_tables':summaries,
            'scope':'Complete exact endpoint circuit bounds, physical load ranges, fair shared acquisitions, frozen plans, survivors and all-error-box margins'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'engineering_verification.json')
    args = parser.parse_args()
    checks = Checks()
    start=time.monotonic()
    required=['instrument_results.json','instrument_qualification_results.json','spice_results.json',
              'fault_results.json','fault_calibration_results.json','bench_import_controls.json']
    for name in required:checks.require((HERE/name).is_file(),'Required engineering stage is present: '+name)
    report = {'schema':'nonagon.engineering.independent-verification.v1',
              'status':'passed','physical_hardware_measurements':0,
              'scope':'Independent replay of retained simulation measurements and engineering envelopes',
              'hardware_plan':plan_self_checks(checks)}
    if (HERE/'instrument_results.json').exists():report['instrument_campaign']=instrument_replay(checks)
    if (HERE/'instrument_qualification_results.json').exists():report['instrument_qualification']=instrument_replay(checks,True)
    if (HERE/'spice_results.json').exists():report['spice']=spice_replay(checks)
    if (HERE/'fault_results.json').exists():report['fault_campaign']=fault_replay(checks)
    if (HERE/'fault_calibration_results.json').exists():report['calibrated_fault_campaign']=calibration_fault_replay(checks)
    bench=read_json(HERE/'bench_import_controls.json')
    checks.require(bench['status']=='passed' and len(bench['rejected'])==22,'bench importer control receipt')
    report['bench_import_controls']={'status':'passed','checks':bench['checks'],'malformed_inputs_rejected':len(bench['rejected']),
                                     'scope':'Synthetic CSV controls separately executed by bench_import_controls.py'}
    retained=[]
    for pattern in ('*.py','*_protocol.json','*_results.json','*_records.jsonl','bench_import_controls.json',
                    'bench_example_*.csv','instrument_fixtures/*.json','netlists/**/*.cir','netlists/**/*.log','netlists/**/*.dat'):
        retained.extend(HERE.glob(pattern))
    report['input_sha256']={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(retained)) if p.is_file()}
    report['checks'] = checks.count
    report['elapsed_seconds']=round(time.monotonic()-start,3)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
