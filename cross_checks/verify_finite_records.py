"""Independent exact replay of finite-load benchmark packets.

Imports no planner, acquisition generator, reconstruction or interval code.
Checks physical forward signals and truth containment of reported enclosures;
it is not a second implementation of the whole interval propagation proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def solve(matrix, rhs):
    n = len(matrix)
    a = [[F(x) for x in row] for row in matrix]
    b = list(map(F, rhs))
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        b[j], b[pivot] = b[pivot], b[j]
        for i in range(j+1, n):
            factor = a[i][j] / a[j][j]
            for k in range(j, n):
                a[i][k] -= factor*a[j][k]
            b[i] -= factor*b[j]
    x = [F(0)] * n
    for i in reversed(range(n)):
        x[i] = (b[i]-sum(a[i][j]*x[j] for j in range(i+1, n))) / a[i][i]
    return x


def run(records, benchmark_receipt, output):
    checks = 0
    def require(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(message)
    def within(box, value):
        return F(box[0]) <= value <= F(box[1])
    def key(loads):
        return tuple(sorted((int(j), F(t)) for j, t in loads.items()))
    def plan(n, t):
        return [()] + [((j, u),) for j in range(1,n) for u in (t,2*t)] + [
            ((i,t),(j,t)) for i in range(1,n) for j in range(i+1,n)]
    receipt = json.loads(benchmark_receipt.read_text())
    require(receipt['records_sha256'] == hashlib.sha256(records.read_bytes()).hexdigest(),
            'Benchmark record hash mismatch')
    cases = [json.loads(line) for line in records.read_text().splitlines()]
    measurements = arm_count = enclosures = 0
    for case in cases:
        fixture, hw = case['fixture'], case['hardware']
        k = [[F(v) for v in row] for row in fixture['K']]
        n = fixture['n']
        scale, current = F(hw['conductance_scale_S']), F(hw['current_A'])
        ve = F(hw['voltage_error_V'])+F(hw['drift_error_V'])
        gain, relative, relaymax = map(F, [hw['current_relative_error'],
            hw['shunt_relative_error'], hw['relay_series_max_ohm']])
        pilot = case['pilot_readings']
        doc = case['pilot_document']
        require(set(doc) == {'schema','network_order','input_and_readout_site_zero_based',
                            'pilot_base_load','hardware','readings'}, 'Pilot schema exposes extra data')
        require(doc['hardware'] == hw and doc['network_order'] == n,
                'Pilot hardware or labels changed')
        require(doc['readings'] == [{'nominal_loads': r['nominal_loads'],
                    'normalized_response_interval': r['normalized_response_interval']} for r in pilot],
                'Planner did not receive the retained pilot scalar intervals')
        t0 = F(doc['pilot_base_load'])
        require([key(r['nominal_loads']) for r in pilot] == plan(n,t0)[:2*n-1],
                'Pilot plan mismatch')
        recommendation = case['recommendation']
        chosen = F(recommendation['selected_base_load'])
        require(chosen in map(F, hw['candidate_base_loads']), 'Selected load outside finite hardware list')
        if recommendation['status'] == 'recommended':
            candidates = recommendation['candidates']
            best = min(candidates, key=lambda r:(F(r['worst_site_local_diagnostic']),
                F(r['base_load']) != t0,F(r['base_load'])))
            require(F(best['base_load']) == chosen, 'Selected load disagrees with disclosed objective')
            singles = {s['site_zero_based']:s for s in recommendation['pilot_single_enclosures']}
            for candidate in candidates:
                t = F(candidate['base_load'])
                e = max(relative, 1-(1-relative)/(1+relaymax*2*t*scale*(1-relative)))
                require(F(candidate['relative_load_envelope']) == e, 'Candidate calibration radius mismatch')
                scores = []
                for row in candidate['sites']:
                    s = singles[row['site_zero_based']]
                    qlo, bhi = F(s['q_interval'][0]), F(s['b_interval'][1])
                    eta = F(recommendation['raw_response_error_bound'])
                    objective = eta*(4/t+8*bhi+4*t*bhi*bhi)/qlo+3*e
                    require(F(row['total_local_relative_q_diagnostic']) == objective,
                            'Disclosed site objective mismatch')
                    scores.append(objective)
                require(F(candidate['worst_site_local_diagnostic']) == max(scores),
                        'Disclosed minimax objective mismatch')
        else:
            require(recommendation['status'] == 'fallback' and chosen == t0,
                    'Unresolved pilot changed the fallback')

        def replay_reading(row):
            nonlocal measurements
            measurements += 1
            nominal = {int(j):F(t) for j,t in row['nominal_loads'].items()}
            branches = row['branches']
            require({b['site'] for b in branches} == set(nominal), 'Shunt routing changed')
            actual = {}
            for branch in branches:
                j, t = branch['site'], nominal[branch['site']]
                delta, relay = F(branch['relative_shunt_deviation']), F(branch['relay_ohm'])
                require(abs(delta) <= relative and 0 <= relay <= relaymax, 'Shunt calibration escaped declared envelope')
                shunt = t*scale*(1+delta)
                actual[j] = shunt/(1+relay*shunt)/scale
                require(F(branch['nominal']) == t and F(branch['actual']) == actual[j],
                        'Physical series-relay shunt conversion mismatch')
            loaded = [r[:] for r in k]
            for j,t in actual.items():
                loaded[j][j] += t
            signal = solve(loaded, [F(i == 0) for i in range(n)])[0]
            require(signal == F(row['true_response']), 'Independent circuit response mismatch')
            actual_current = F(row['actual_current_A'])
            require(abs(actual_current/current-1) <= gain, 'Actual source current outside calibration')
            error, drift = F(row['readout_error_V']), F(row['drift_V'])
            require(abs(error) <= F(hw['voltage_error_V']) and abs(drift) <= F(hw['drift_error_V']),
                    'Actual voltage error outside calibration')
            voltage = actual_current*signal/scale+error+drift
            require(voltage == F(row['measured_voltage_V']), 'Voltage normalization mismatch')
            expected = [(voltage-ve)*scale/(current*(1+gain)),
                        (voltage+ve)*scale/(current*(1-gain))]
            require(list(map(F,row['normalized_response_interval'])) == expected,
                    'Submitted scalar interval differs from physical calibration conversion')
            require(within(row['normalized_response_interval'],signal), 'Scalar interval lost actual response')
        for row in pilot:
            replay_reading(row)
        arms = case['arms']
        require(len(arms) == 2, 'Expected fixed and planned arms')
        totals = []
        for index, arm in enumerate(arms):
            arm_count += 1
            row, bounds = arm['case'], arm['bounds']
            t = t0 if index == 0 else chosen
            require(F(row['base_load']) == t and [key(r['nominal_loads']) for r in arm['readings']] == plan(n,t),
                    'Prospective arm did not use its selected canonical plan')
            total = len(pilot)+len(arm['readings'])
            require(row['total_scalar_count'] == total and row['pilot_scalar_count'] == len(pilot)
                    and row['final_scalar_count'] == len(arm['readings']), 'Arm budget accounting mismatch')
            totals.append(total)
            for reading in arm['readings']:
                replay_reading(reading)
            if bounds['status'] == 'enclosed':
                enclosures += 1
                columns = [solve(k, [F(i == j) for i in range(n)]) for j in range(n)]
                g = [[columns[j][i] for j in range(n)] for i in range(n)]
                for name, truth in [('K_intervals',k),('G_intervals',g)]:
                    for i in range(n):
                        for j in range(n):
                            require(within(bounds[name][i][j],truth[i][j]), 'Matrix enclosure lost independent truth')
                for edge in bounds['edges']:
                    i,j = edge['sites_zero_based']; weight = -k[i][j]
                    require(within(edge['weight_interval'],weight), 'Edge enclosure lost truth')
                    require(edge['status'] != 'present' or weight > 0, 'False edge presence')
                    require(edge['status'] != 'absent_under_gap' or weight == 0, 'False edge absence')
            else:
                require(bounds['status'] == 'unresolved', 'Unknown inference status')
        require(totals[0] == totals[1], 'Unequal pilot-plus-final comparison budgets')
    require(len(cases) == receipt['prospective_cases'] and arm_count == receipt['arm_cases'],
            'Benchmark completeness mismatch')
    result = {'schema':'nonagon.finite-load-independent-record-replay.v1','status':'passed',
              'prospective_cases':len(cases),'arm_cases':arm_count,'physical_scalar_packets':measurements,
              'matrix_enclosures_checked':enclosures,'exact_predicates':checks,
              'records_sha256':hashlib.sha256(records.read_bytes()).hexdigest(),
              'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':__doc__.strip()}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    run(args.records,args.receipt,args.output)
