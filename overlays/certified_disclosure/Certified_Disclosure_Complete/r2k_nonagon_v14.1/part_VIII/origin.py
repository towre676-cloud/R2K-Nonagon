#!/usr/bin/env python3
"""Apply the universal leading omitted-origin correction to archived production predictions.

This calculation does not estimate R2/x² in double precision and supplies no
uniform near-origin or Euler-product-tail certificate. No new zeros are computed.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / 'part_VI' / 'bk_multi_100000_1000000_10000000.json'
a = 1e-4
B = json.loads(ARCHIVE.read_text())
out = {'a': a, 'status': 'leading endpoint correction of archived values',
       'new_zero_runs': 0, 'uniform_remainder_certified': False,
       'conditional_remainder': 'If 0 <= R2(x) <= c*x^2 on (0,a], remainder <= 2*c*L*a^3/3',
       'source_sha256': hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(), 'rows': []}
for height in B['heights']:
    for row in height['rows']:
        L = row['L']; lead = 2*L*a-a*a
        corrected = row['BK_height_avg']-lead
        out['rows'].append({'n0': height['n0'], 'L': L,
            'leading_2La_minus_a2': lead, 'BK_old': row['BK_height_avg'],
            'BK_after_leading_correction': corrected, 'zeta': row['zeta'], 'se': row['se'],
            'residual_old': row['zeta']-row['BK_height_avg'],
            'residual_corrected': row['zeta']-corrected,
            'residual_in_archived_se': (row['zeta']-corrected)/row['se'],
            'GUE_old': row['GUE'], 'GUE_after_leading_correction': row['GUE']-lead})
out['max_abs_residual_in_archived_se'] = max(abs(r['residual_in_archived_se']) for r in out['rows'])
(HERE/'origin.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({'rows': len(out['rows']), 'status': out['status'],
    'max_abs_residual_in_archived_se': out['max_abs_residual_in_archived_se'],
    'uniform_remainder_certified': False}))
