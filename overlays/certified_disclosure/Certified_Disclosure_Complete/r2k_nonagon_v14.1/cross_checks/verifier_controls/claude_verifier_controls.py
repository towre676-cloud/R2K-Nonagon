"""Adversarial interface checks of Claude's algebraic recovery helper.

This imports only the guarded verifier module, never the unguarded generator.
It does not modify the source release. Normal and optimized interpreters run in
isolated child processes; their recovered matrices receive an independent
Gaussian response replay and a physical-class check.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
TARGET=HERE.parents[1]/'research/claude_v8/part_VIII/verify_intervention.py'
T1,T2,TP=F(1,2),F(5),F(3,2)


def solve(a,b):
    n=len(a);m=[[F(v) for v in row]+[F(b[i])] for i,row in enumerate(a)]
    for k in range(n):
        p=next(i for i in range(k,n) if m[i][k]);m[k],m[p]=m[p],m[k]
        for i in range(k+1,n):
            c=m[i][k]/m[k][k]
            for j in range(k,n+1):m[i][j]-=c*m[k][j]
    x=[F(0)]*n
    for i in reversed(range(n)):x[i]=(m[i][-1]-sum(m[i][j]*x[j] for j in range(i+1,n)))/m[i][i]
    return x


def response(k,key):
    a=[row[:] for row in k]
    for j,t in key:a[j][j]+=t
    return solve(a,[F(i==0) for i in range(len(k))])[0]


def canonical(n):
    keys=[()]
    for j in range(1,n):keys.extend([((j,T1),),((j,T2),)])
    keys.extend(((i,TP),(j,TP)) for i in range(1,n) for j in range(i+1,n))
    return keys


def samples_from_g2(g00,g01_squared,g11):
    return {():F(g00),((1,T1),):g00-T1*g01_squared/(1+T1*g11),
            ((1,T2),):g00-T2*g01_squared/(1+T2*g11)}


def fixture(name):
    k=[[F(2),F(-1)],[F(-1),F(2)]];n=2
    measurements={key:response(k,key) for key in canonical(n)}
    if name=='valid_reference':pass
    elif name=='nonsquare_single_root':measurements=samples_from_g2(F(2,3),F(1,8),F(2,3))
    elif name=='indefinite_recovered_matrix':measurements=samples_from_g2(F(1),F(4),F(1))
    elif name=='negative_grounding_recovered_matrix':measurements=samples_from_g2(F(1,4),F(1,4),F(2))
    elif name=='missing_baseline':del measurements[()]
    elif name=='missing_single':del measurements[((1,T2),)]
    elif name=='unexpected_measurement':measurements[((77,F(1)),)]=F(123)
    elif name=='inconsistent_extra_load':measurements[((1,F(7)),)]=F(999)
    elif name=='nonsquare_pair_root':
        n=3;k=[[F(4),F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]]
        measurements={key:response(k,key) for key in canonical(n)}
        measurements[((1,TP),(2,TP))]+=F(1,1000000)
    else:raise ValueError(name)
    return n,measurements


def physical(k):
    n=len(k)
    if any(k[i][j]!=k[j][i] for i in range(n) for j in range(n)):return False
    if any(k[i][j]>0 for i in range(n) for j in range(n) if i!=j):return False
    if any(sum(row)<0 for row in k) or not any(sum(row)>0 for row in k):return False
    reached={0};changed=True
    while changed:
        changed=False
        for i in list(reached):
            for j in range(n):
                if k[i][j]<0 and j not in reached:reached.add(j);changed=True
    return len(reached)==n


def worker(name):
    spec=importlib.util.spec_from_file_location('claude_guarded_verifier',TARGET)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    n,measurements=fixture(name)
    result={'case':name,'optimized':not __debug__,'input':[
        {'loads':{str(j):str(t) for j,t in key},'response':str(value)} for key,value in measurements.items()]}
    try:
        recovered=module.recover(n,deepcopy(measurements),T1,T2,TP)
        canonical_keys=canonical(n)
        errors=[abs(response(recovered,key)-measurements[key]) for key in canonical_keys]
        result.update(outcome='returned',recovered_K=[[str(v) for v in row] for row in recovered],
                      canonical_forward_max_error=str(max(errors)),physical_grounded_network=physical(recovered),
                      unexpected_keys=[{str(j):str(t) for j,t in key} for key in measurements if key not in canonical_keys])
    except Exception as error:result.update(outcome='rejected',exception=type(error).__name__,message=str(error))
    print(json.dumps(result))


def main():
    global TARGET
    parser=argparse.ArgumentParser();parser.add_argument('--worker');parser.add_argument('--target',type=Path,default=TARGET);args=parser.parse_args()
    TARGET=args.target.resolve()
    if args.worker:worker(args.worker);return
    cases=['valid_reference','nonsquare_single_root','indefinite_recovered_matrix',
           'negative_grounding_recovered_matrix','missing_baseline','missing_single',
           'unexpected_measurement','inconsistent_extra_load','nonsquare_pair_root']
    rows=[]
    for optimized in (False,True):
        for name in cases:
            command=[sys.executable]+(['-O'] if optimized else [])+[str(Path(__file__).resolve()),'--worker',name,'--target',str(TARGET)]
            process=subprocess.run(command,cwd=HERE,text=True,capture_output=True,check=True)
            rows.append(json.loads(process.stdout))
    result={'schema':'nonagon.integration.claude-verifier-controls.v1',
            'source':str(TARGET),'source_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),
            'cases':len(cases),'runs':len(rows),'records':rows,
            'conclusion':'recover is an algebraic test helper, not a standalone certificate validator; preserve production checker and add explicit preflight/forward replay around cross-check calls'}
    (HERE/'claude_verifier_controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'cases':len(cases),'runs':len(rows),'optimized_inconsistent_returns':[
        row['case'] for row in rows if row['optimized'] and row['outcome']=='returned' and row['canonical_forward_max_error']!='0']}))


if __name__=='__main__':main()
