#!/usr/bin/env python3
"""Post-confirmation descriptive check of a structural model restriction."""
import json,numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def run():
    rows=[]
    for i in range(1,55):
        d=np.genfromtxt(ROOT/'data'/f'ID{i:02}.csv',delimiter=',',skip_header=1);d=d[np.argsort(d[:,0])]
        q=int(np.argmin(d[:,1]));delta=float(d[-1,1]-d[q,1])
        rows.append({'device':f'ID{i:02}','phase':'development' if i<=6 else 'confirmation','minimum_real_frequency_hz':float(d[q,0]),'minimum_real_ohm':float(d[q,1]),'highest_frequency_hz':float(d[-1,0]),'highest_frequency_real_ohm':float(d[-1,1]),'high_frequency_increase_ohm':delta,'violates_nonincreasing_real_data':delta>0})
    out={'label':'Post-confirmation descriptive model-mismatch diagnosis; no model or scoring changes.','model_constraint':'For positive RC Debye terms with series ideal RLC, d Re Z(iomega)/domega = -sum_j 2 a_j tau_j^2 omega/(1+tau_j^2 omega^2)^2 <=0 for omega>0.','observed_nonmonotonic_confirmation_cells':sum(r['violates_nonincreasing_real_data'] for r in rows if r['phase']=='confirmation'),'confirmation_cells':48,'scope':'Applies to exact stored values. Unknown measurement uncertainty prevents inference of physical incompatibility or nonpassivity. More general positive-real RL relaxation terms can permit rising real part; not tested or selected here.','records':rows}
    (ROOT/'model_mismatch_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
