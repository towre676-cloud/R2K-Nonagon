"""Audit recorded pilot inputs, gate decisions, budgets and frozen sources."""
import hashlib,json
from pathlib import Path
import numpy as np
import experiment,gate
from gate import adaptive,core

root=Path(__file__).resolve().parent
experiment.verify()
rows=json.loads((root/'confirmation_raw.json').read_text());a=np.load(root/'confirmation_observations.npz')
rule=json.loads((root/'selected_rule.json').read_text())['rule'];seen={}
for row in rows:
    indices=np.array(row['selected_indices']);cluster=row['cluster'];eta=row['noise']
    if len(indices)!=64 or len(set(indices))!=64:raise ValueError('unique64querybudgetfailed')
    if not set(adaptive.PILOT).issubset(set(indices)):raise ValueError('uncountedpilotinput')
    observed=a[cluster+'_candidate_clean']*(1+eta*a[cluster+'_standardized_noise'])
    if core.digest(adaptive.POOL_W,observed)!=row['pool_observation_sha256']:raise ValueError('sharednoisefailed')
    if core.digest(indices,observed[indices])!=row['selected_observation_sha256']:raise ValueError('selectedvaluesfailed')
    if row['case_id'] not in seen:
        seen[row['case_id']]=gate.pilot_features(adaptive.POOL_W[adaptive.PILOT],observed[adaptive.PILOT],
                                              2*len(adaptive.system(row['scenario'],int(cluster.rsplit('-',1)[1]),69377123)['truth']['m']),adaptive.BAND)
    expected=seen[row['case_id']]
    for key,value in expected.items():
        if key=='valid':
            if row['features'][key]!=value:raise ValueError('featurevaliditymismatch')
        elif isinstance(value,(int,float)):
            if not np.isclose(value,row['features'][key],rtol=1e-10,atol=1e-12):raise ValueError('pilotfeaturemismatch')
    if row['strategy']=='gated':
        actual='adaptive_disagreement' if gate.matches(expected,rule) else 'log_uniform'
        if row['actual_acquisition']!=actual:raise ValueError('decisionmismatch')
print(json.dumps(dict(rows_verified=len(rows),cases_verified=len(seen),pilot_in_every64_query_budget=True,
                     same_noisy_pilot_values_verified=True,features_recomputed=True,decisions_verified=True,
                     scientific_hashes_unchanged=True),indent=2))
