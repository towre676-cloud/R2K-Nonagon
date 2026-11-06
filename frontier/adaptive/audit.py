"""Read-only receipt verification after a completed confirmation run."""
import hashlib,json
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent
freeze=json.loads((root/'freeze.json').read_text())
for name,wanted in freeze['sha256'].items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=wanted:
        raise ValueError('scientific input changed: '+name)
rows=json.loads((root/'confirmation_raw.json').read_text())
arrays=np.load(root/'confirmation_observations.npz')
pool=arrays['candidate_omega'];heldout=arrays['heldout_omega']
if np.any(pool[:,None]==heldout[None,:]):raise ValueError('test/training frequency overlap')
def digest(*values):
    h=hashlib.sha256()
    for v in values:h.update(np.ascontiguousarray(v).tobytes())
    return h.hexdigest()
for row in rows:
    selected=np.array(row['selected_indices'])
    if len(set(selected))!=len(selected):raise ValueError('duplicate query')
    if len(selected)>64:raise ValueError('budget exceeded')
    if row['status']=='ok' and len(selected)!=64:raise ValueError('incomplete successful budget')
    cluster=row['cluster'];noise=row['noise']
    observed=arrays[cluster+'_candidate_clean']*(1+noise*arrays[cluster+'_standardized_noise'])
    if digest(pool,observed)!=row['candidate_observation_sha256']:raise ValueError('candidate noise table mismatch')
    if row['status']=='ok' and digest(selected,observed[selected])!=row['selected_observation_sha256']:
        raise ValueError('selected observations mismatch')
print(json.dumps({'rows_verified':len(rows),'successful':sum(r['status']=='ok' for r in rows),
                  'frozen_hashes_unchanged':True,'shared_noise_tables_and_queries_verified':True,
                  'protected_grid_disjoint':True,'budgets_and_uniqueness_verified':True},indent=2))
