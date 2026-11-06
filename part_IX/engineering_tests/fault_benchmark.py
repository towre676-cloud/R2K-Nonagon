"""Frozen bounded-error fault diagnosis: interventions versus reciprocal scans.

This experiment identifies a member of an explicit single-fault catalog. It does
not infer an unrestricted physical graph and it does not claim hardware results.
All signals and compatibility decisions are rational; candidate predictions use
rank-one updates, allowing a separate Gaussian-solve replay to check them.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import random
import time

HERE = Path(__file__).resolve().parent
N = 6
PROTOCOL = {
    "schema": "nonagon.engineering.fault_protocol.v1",
    "question": "Identify a single specified fault in a known six-node passive fixture using calibrated scalar measurements.",
    "dimension": N,
    "development_seed": 438771,
    "confirmation_seed": 903671,
    "development_fixtures": 8,
    "confirmation_fixtures": 20,
    "confirmation_specimens_per_fixture": 5,
    "fault_types": ["edge_weak_5", "edge_weak_20", "edge_open", "unexpected_bridge", "ground_up_20"],
    "catalog_includes_healthy": True,
    "network_generator": "weighted six-cycle plus 0-3 extra chords, rational weights 1/4..2; all ground weights 1/4..1",
    "intervention_pool": "baseline; single hidden-site loads 1/4,1,2,4; equal pair loads 1/4,1,2,4; 61 distinct configurations",
    "canonical_intervention_order": "baseline, loads1/2 at each hidden site, pair load1;21 scalar configurations",
    "conventional_pool": "21 reciprocal compliance entries G_ij, 0<=i<=j<6; unit current at i and voltage at j",
    "policies": ["unloaded_port", "canonical_intervention", "greedy_fixed_intervention", "natural_reciprocal_scan", "greedy_fixed_reciprocal_scan"],
    "greedy_design": "Before any fault observation, repeatedly select query separating the largest number of still-unseparated catalog pairs by >2eta; ties follow fixed pool order. No test truth supplied.",
    "noise_bounds": ["0", "1/1000000", "1/100000", "1/10000"],
    "noise_units": "absolute normalized voltage/current compliance",
    "noise_oracle": "Choose an allowed reading maximizing current compatible-candidate count; tie choose largest absolute error, then smallest reading. Truth is used only by this acquisition oracle and scorer.",
    "budgets": [1, 6, 12, 21],
    "physical_scale": {"conductance_unit_S": "1/1000", "current_A": "1/10000", "voltage_error_V_per_normalized_eta": "1/10"},
    "load_calibration": "exact in this catalog comparison; a separate instrument campaign propagates shunt uncertainty",
    "comparison_contract": "Identical scalar error bounds and budgets. Reciprocal scans may move source and readout; fixed-port interventions keep both at site1 and move at most2 shunts. No timing/cost superiority claimed.",
    "confirmation_policy_changes": "none; code and protocol hashes retained",
}


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def solve(A, rhs):
    n = len(A)
    w = [[F(x) for x in row] + [F(rhs[i])] for i, row in enumerate(A)]
    for k in range(n):
        p = next(i for i in range(k, n) if w[i][k])
        w[k], w[p] = w[p], w[k]
        x = w[k][k]
        w[k] = [v/x for v in w[k]]
        for i in range(n):
            if i != k and w[i][k]:
                x = w[i][k]
                w[i] = [v-x*u for v, u in zip(w[i], w[k])]
    return [row[-1] for row in w]


def inverse(A):
    n = len(A)
    columns = [solve(A, [int(i == j) for i in range(n)]) for j in range(n)]
    return [list(row) for row in zip(*columns)]


def encode_matrix(A):
    return [[str(x) for x in row] for row in A]


def edge_update(K, i, j, change):
    out = [row[:] for row in K]
    out[i][i] += change; out[j][j] += change
    out[i][j] -= change; out[j][i] -= change
    return out


def physical(K):
    if any(K[i][j] != K[j][i] for i in range(N) for j in range(N)):
        return False
    if any(K[i][j] > 0 for i in range(N) for j in range(N) if i != j):
        return False
    if any(sum(row) < 0 for row in K) or not any(sum(row) > 0 for row in K):
        return False
    seen, todo = {0}, [0]
    while todo:
        i = todo.pop()
        for j in range(N):
            if K[i][j] < 0 and j not in seen:
                seen.add(j); todo.append(j)
    return len(seen) == N


def fixture(rng):
    edges = {}
    for i in range(N):
        j = (i+1) % N
        edges[tuple(sorted((i, j)))] = F(rng.randint(1, 8), 4)
    chords = [p for p in combinations(range(N), 2) if p not in edges]
    rng.shuffle(chords)
    for p in chords[:rng.randrange(4)]:
        edges[p] = F(rng.randint(1, 8), 4)
    K = [[F(0) for _ in range(N)] for _ in range(N)]
    for (i, j), w in sorted(edges.items()):
        K = edge_update(K, i, j, w)
    for i in range(N):
        K[i][i] += F(rng.randint(1, 4), 4)
    if not physical(K):
        raise ValueError("invalid fixture")
    return K


def catalog(K):
    out = [{"id": "healthy", "type": "healthy", "sites": [], "K": K}]
    for i, j in combinations(range(N), 2):
        weight = -K[i][j]
        if weight:
            for kind, change in [("edge_weak_5", -weight/20), ("edge_weak_20", -weight/5), ("edge_open", -weight)]:
                changed = edge_update(K, i, j, change)
                if physical(changed):
                    out.append({"id": f"{kind}_{i}_{j}", "type": kind, "sites": [i, j], "change": str(change), "K": changed})
        else:
            changed = edge_update(K, i, j, F(1, 4))
            out.append({"id": f"unexpected_bridge_{i}_{j}", "type": "unexpected_bridge", "sites": [i, j], "change": "1/4", "K": changed})
    for i in range(N):
        changed = [row[:] for row in K]
        change = sum(K[i])/5
        changed[i][i] += change
        out.append({"id": f"ground_up_20_{i}", "type": "ground_up_20", "sites": [i], "change": str(change), "K": changed})
    keys = [tuple(x for row in c["K"] for x in row) for c in out]
    if len(set(keys)) != len(keys) or not all(physical(c["K"]) for c in out):
        raise ValueError("catalog duplicates or nonphysical candidates")
    return out


def intervention_queries():
    baseline = {"id": "baseline", "kind": "shunt", "source": 0, "readout": 0, "loads": {}}
    canonical = [baseline]
    for j in range(1, N):
        canonical += [{"id": f"single_{j}_{t}", "kind": "shunt", "source": 0, "readout": 0, "loads": {str(j): str(t)}} for t in [1, 2]]
    canonical += [{"id": f"pair_{i}_{j}_1", "kind": "shunt", "source": 0, "readout": 0, "loads": {str(i): "1", str(j): "1"}} for i, j in combinations(range(1, N), 2)]
    pool = canonical[:]
    existing = {canonical_json(q["loads"]) for q in pool}
    for strength in [F(1, 4), F(1), F(2), F(4)]:
        for count in [1, 2]:
            for sites in combinations(range(1, N), count):
                loads = {str(i): str(strength) for i in sites}
                key = canonical_json(loads)
                if key not in existing:
                    existing.add(key)
                    pool.append({"id": "shunt_"+"_".join(map(str, sites))+"_"+str(strength), "kind": "shunt", "source": 0, "readout": 0, "loads": loads})
    return canonical, pool


def conventional_queries():
    return [{"id": f"G_{i}_{j}", "kind": "scan", "source": i, "readout": j, "loads": {}} for i in range(N) for j in range(i, N)]


def prediction(G, query):
    if query["kind"] == "scan":
        return G[query["source"]][query["readout"]]
    h = G[0][0]
    loads = [(int(j), F(t)) for j, t in query["loads"].items()]
    if not loads:
        return h
    j, t = loads[0]
    denom = 1+t*G[j][j]
    h -= t*G[0][j]**2/denom
    if len(loads) == 2:
        k, u = loads[1]
        cross = G[0][k]-t*G[0][j]*G[j][k]/denom
        diagonal = G[k][k]-t*G[j][k]**2/denom
        h -= u*cross**2/(1+u*diagonal)
    return h


def greedy_plan(query_ids, table, eta, budget=21):
    """Design consumes catalog predictions, not the selected fault or data."""
    candidates = len(next(iter(table.values())))
    pairs = list(combinations(range(candidates), 2))
    masks = {}
    for q in query_ids:
        values = table[q]
        mask = 0
        for k, (i, j) in enumerate(pairs):
            if abs(values[i]-values[j]) > 2*eta:
                mask |= 1 << k
        masks[q] = mask
    remaining = (1 << len(pairs))-1
    available = list(query_ids); out = []
    for _ in range(min(budget, len(available))):
        chosen = max(available, key=lambda q: ((masks[q] & remaining).bit_count(), masks[q].bit_count(), -query_ids.index(q)))
        out.append(chosen); remaining &= ~masks[chosen]; available.remove(chosen)
    return out


def compatible_ids(values, measurement, eta, survivors):
    """Observer receives no true label, circuit or perturbation information."""
    return [i for i in survivors if abs(values[i]-measurement) <= eta]


def adversarial_reading(values, true_index, eta, survivors):
    true = values[true_index]
    lo, hi = true-eta, true+eta
    choices = {lo, hi, true}
    for i in survivors:
        for value in [values[i]-eta, values[i]+eta]:
            if lo <= value <= hi:
                choices.add(value)
    return max(choices, key=lambda m: (len(compatible_ids(values, m, eta, survivors)), abs(m-true), -m))


def run_trial(table, plan, true_index, eta, candidate_ids):
    survivors = list(range(len(candidate_ids))); readings = []; checkpoints = {}
    for count, query in enumerate(plan, 1):
        measured = adversarial_reading(table[query], true_index, eta, survivors)
        exact = table[query][true_index]
        survivors = compatible_ids(table[query], measured, eta, survivors)
        if true_index not in survivors:
            raise ValueError("truth escaped its declared scalar error interval")
        readings.append({"query": query, "measured": str(measured), "actual": str(exact), "error": str(measured-exact), "survivors": [candidate_ids[i] for i in survivors]})
        if count in PROTOCOL["budgets"]:
            checkpoints[str(count)] = {"survivors": len(survivors), "status": "identified" if len(survivors) == 1 else "unresolved", "identified": candidate_ids[survivors[0]] if len(survivors) == 1 else None}
    return {"readings": readings, "checkpoints": checkpoints}


def costs(plan, query_map):
    queries = [query_map[q] for q in plan]
    return {"scalar_readings": len(queries), "source_addresses": len({q["source"] for q in queries}),
            "readout_addresses": len({q["readout"] for q in queries}),
            "distinct_shunt_pads": len({j for q in queries for j in q["loads"]}),
            "maximum_simultaneous_shunts": max([len(q["loads"]) for q in queries], default=0),
            "configuration_count": len(queries),
            "source_block_minimum": len({q["source"] for q in queries}),
            "timing_model": "not measured; repeated settings and real switching/settling costs are not assigned synthetic runtimes"}


def controls():
    values = [F(0), F(1), F(2)]
    checks = [compatible_ids(values, F(1), F(0), [0,1,2]) == [1],
              compatible_ids(values, F(1,2), F(1,2), [0,1,2]) == [0,1],
              compatible_ids(values, F(10), F(1), [0,1,2]) == []]
    if not all(checks):
        raise ValueError("compatibility controls failed")
    return len(checks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=HERE/'fault_results.json')
    args = ap.parse_args()
    started = time.monotonic(); HERE.mkdir(exist_ok=True)
    frozen_path = HERE/'fault_protocol.json'
    frozen_text = json.dumps(PROTOCOL, indent=2, sort_keys=True)+'\n'
    if frozen_path.exists() and frozen_path.read_text() != frozen_text:
        raise RuntimeError('Refusing to silently change the frozen protocol')
    frozen_path.write_text(frozen_text)
    protocol_hash = hashlib.sha256(frozen_text.encode()).hexdigest()
    canonical, pool = intervention_queries(); scans = conventional_queries()
    queries = {q['id']:q for q in pool+scans}
    if len(canonical) != 21 or len(pool) != 61 or len(scans) != 21:
        raise ValueError('Acquisition accounting failed')

    # Separate development fixtures exercise the frozen prediction machinery. No
    # confirmation score is used to select parameters, load grid or policy.
    dev_rng = random.Random(PROTOCOL['development_seed']); dev_checks = 0
    for _ in range(PROTOCOL['development_fixtures']):
        K=fixture(dev_rng); G=inverse(K)
        for q in pool[::7]:
            A=[row[:] for row in K]
            for j,t in q['loads'].items(): A[int(j)][int(j)]+=F(t)
            independent=solve(A,[1,0,0,0,0,0])[0]
            if independent != prediction(G,q):raise ValueError('Development update check failed')
            dev_checks += 1

    rng=random.Random(PROTOCOL['confirmation_seed']); fixtures=[]; trials=[]; summaries={}; out_of_catalog=[]
    for fixture_index in range(PROTOCOL['confirmation_fixtures']):
        K=fixture(rng); candidates=catalog(K); ids=[c['id'] for c in candidates]
        inverses=[inverse(c['K']) for c in candidates]
        table={qid:[prediction(G,q) for G in inverses] for qid,q in queries.items()}
        canonical_ids=[q['id'] for q in canonical]; pool_ids=[q['id'] for q in pool]; scan_ids=[q['id'] for q in scans]
        tuples=[tuple(table[q][i] for q in canonical_ids) for i in range(len(candidates))]
        if len(set(tuples)) != len(candidates):raise ValueError('Exact canonical catalog is not identifiable')
        selected={kind:rng.choice([i for i,c in enumerate(candidates) if c['type']==kind]) for kind in PROTOCOL['fault_types']}
        fdoc={'id':f'fixture_{fixture_index:02d}','healthy_K':encode_matrix(K),
              'catalog':[{**{k:v for k,v in c.items() if k!='K'},'K':encode_matrix(c['K'])} for c in candidates],
              'predictions':{q:[str(v) for v in vs] for q,vs in table.items()},'plans':{}}
        for eta_text in PROTOCOL['noise_bounds']:
            eta=F(eta_text)
            plans={'unloaded_port':['baseline'],'canonical_intervention':canonical_ids,
                   'greedy_fixed_intervention':greedy_plan(pool_ids,table,eta),
                   'natural_reciprocal_scan':scan_ids,
                   'greedy_fixed_reciprocal_scan':greedy_plan(scan_ids,table,eta)}
            fdoc['plans'][eta_text]={name:{'queries':plan,'costs_at_budget':{str(b):costs(plan[:b],queries) for b in PROTOCOL['budgets'] if b<=len(plan)}} for name,plan in plans.items()}
            for kind,true_index in selected.items():
                for policy,plan in plans.items():
                    result=run_trial(table,plan,true_index,eta,ids)
                    record={'fixture':fdoc['id'],'fault_type':kind,'true_candidate':ids[true_index],
                            'eta':eta_text,'policy':policy,**result}
                    trials.append(record)
                    for b,checkpoint in result['checkpoints'].items():
                        key=(eta_text,policy,b)
                        s=summaries.setdefault(key,{'eta':eta_text,'policy':policy,'budget':int(b),'cases':0,'identified':0,'unresolved':0,'wrong':0,'survivor_sum':0,'by_fault_type':{}})
                        s['cases']+=1;s[checkpoint['status']]+=1;s['survivor_sum']+=checkpoint['survivors']
                        if checkpoint['identified'] is not None and checkpoint['identified'] != ids[true_index]:s['wrong']+=1
                        ss=s['by_fault_type'].setdefault(kind,{'cases':0,'identified':0});ss['cases']+=1;ss['identified']+=checkpoint['status']=='identified'
        # A two-edge change is deliberately outside the single-fault catalog.
        edges=[(i,j) for i,j in combinations(range(N),2) if K[i][j]<0]
        hidden=[row[:] for row in K]
        for i,j in edges[:2]:hidden=edge_update(hidden,i,j,K[i][j]/5)
        Ghidden=inverse(hidden); survivors=list(range(len(candidates))); records=[]
        for qid in canonical_ids:
            m=prediction(Ghidden,queries[qid]);survivors=compatible_ids(table[qid],m,F(0),survivors)
            records.append({'query':qid,'measured':str(m),'survivors':[ids[i] for i in survivors]})
        if survivors:raise ValueError('Out-of-catalog negative control was silently assigned')
        out_of_catalog.append({'fixture':fdoc['id'],'K':encode_matrix(hidden),'declared_fault':'two edge weaknesses outside catalog','readings':records,'status':'inconsistent_with_catalog'})
        fixtures.append(fdoc)
        print(f'Fixture {fixture_index+1}/20 complete; {len(candidates)} candidates',flush=True)
    if hashlib.sha256(frozen_path.read_bytes()).hexdigest()!=protocol_hash:raise ValueError('Frozen protocol changed during confirmation')
    result={'schema':'nonagon.engineering.fault_results.v1','protocol_sha256':protocol_hash,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'exact rational synthetic circuit catalog; calibrated scalar-error enclosures, no apparatus measurements',
            'development_checks':dev_checks,'compatibility_controls':controls(),'confirmation_fixtures':len(fixtures),
            'distinct_confirmation_specimens':sum(5 for _ in fixtures),'trial_count':len(trials),
            'query_definitions':queries,'summaries':list(summaries.values()),'fixtures':fixtures,'trials':trials,
            'out_of_catalog_controls':out_of_catalog,'out_of_catalog_rejections':len(out_of_catalog),
            'elapsed_seconds':round(time.monotonic()-started,3)}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    concise=[{k:s[k] for k in ['eta','policy','budget','identified','unresolved','wrong']} for s in result['summaries'] if s['eta']=='1/100000' and s['budget'] in [6,21]]
    print(json.dumps({'output':str(args.output),'trials':len(trials),'confirmation_specimens':100,'selected_summaries':concise,'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
