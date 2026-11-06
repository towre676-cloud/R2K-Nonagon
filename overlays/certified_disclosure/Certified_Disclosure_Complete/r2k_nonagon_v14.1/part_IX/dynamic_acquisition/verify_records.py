#!/usr/bin/env python3
"""Independent saved-record checks; imports no acquisition or inverse modules."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path


def verify(directory):
    checks=0
    def check(ok,message):
        nonlocal checks
        checks+=1
        if not ok: raise ValueError(message)

    def determinant(a):
        # Fraction elimination local to this verifier; Cramer's ratio below.
        a=[[F(x) for x in row] for row in a]; value=F(1); n=len(a)
        for j in range(n):
            pivot=next((i for i in range(j,n) if a[i][j]),None)
            if pivot is None: return F(0)
            if pivot!=j: a[pivot],a[j]=a[j],a[pivot]; value=-value
            value*=a[j][j]
            for i in range(j+1,n):
                ratio=a[i][j]/a[j][j]
                for k in range(j+1,n): a[i][k]-=ratio*a[j][k]
        return value

    def upperexp(q):
        return F(1,1<<(q.numerator//q.denominator))

    protocol=json.loads((directory/'protocol.json').read_text())
    records=json.loads((directory/'acquisition_records.json').read_text())
    receipt=json.loads((directory/'dynamic_receipt.json').read_text())
    check(receipt['records_sha256']==hashlib.sha256((directory/'acquisition_records.json').read_bytes()).hexdigest(),'Record digest')
    check(receipt['protocol_sha256']==hashlib.sha256((directory/'protocol.json').read_bytes()).hexdigest(),'Protocol digest')
    check(receipt['source_sha256']==hashlib.sha256((Path(__file__).parent/'replay.py').read_bytes()).hexdigest(),'Acquisition source digest')
    prior=protocol['certification_prior']
    check({k:prior[k] for k in ('beta','mass_min','mass_max','ground_min','stiffness_norm_inf_max')}==
          {'beta':'1/2','mass_min':'1','mass_max':'7/4','ground_min':'2','stiffness_norm_inf_max':'16'},'Frozen prior changed')
    beta=F(prior['beta']); mmin=F(prior['mass_min']); mmax=F(prior['mass_max'])
    ell=F(prior['ground_min'])/mmax; amplitude=1/(mmin*ell)
    horizon=F(protocol['horizon']); step=F(protocol['sample_interval'])
    sample=F(protocol['sample_absolute_error']); numeric=F(protocol['transformed_numerical_error'])
    scalar_count=0; matrix_count=0
    for record in records:
        f=record['fixture']; n=len(f['K'])
        k,c,m=[[[F(x) for x in row] for row in f[name]] for name in ('K','C','M')]
        check(all(c[i][j]==2*beta*m[i][j] for i in range(n) for j in range(n)),'Damping class')
        check(all(m[i][j]==0 for i in range(n) for j in range(n) if i!=j),'Diagonal mass')
        check(all(mmin<=m[i][i]<=mmax for i in range(n)),'Mass prior')
        check(all(sum(row)>=F(prior['ground_min']) for row in k),'Ground prior')
        check(max(sum(abs(x) for x in row) for row in k)<=F(prior['stiffness_norm_inf_max']),'Norm prior')
        check(all(k[i][j]==k[j][i] and (i==j or k[i][j]<0) for i in range(n) for j in range(n)),'Test stiffness class')
        expected=[()]
        for j in range(1,n): expected.extend([((j,F(1)),),((j,F(2)),)])
        expected.extend(((j,F(1)),(q,F(1))) for j in range(1,n) for q in range(j+1,n))
        seen={s:set() for s in (F(0),F(1),F(2))}
        for row in record['packets']:
            scalar_count+=1; s=F(row['s'])
            loads=tuple(sorted((int(j),F(t)) for j,t in row['loads'].items()))
            check(s in seen and loads not in seen[s] and loads in expected,'Canonical unique scalar setting')
            seen[s].add(loads)
            p=[[k[i][j]+s*c[i][j]+s*s*m[i][j] for j in range(n)] for i in range(n)]
            for j,t in loads:p[j][j]+=t
            exact=determinant([row[1:] for row in p[1:]])/determinant(p)
            check(exact==F(row['exact_reference']),'Independent Cramer response')
            lo,hi=map(F,row['interval']); center=F(row['center'])
            budget={q:F(v) for q,v in row['error_budget'].items()}
            w=budget['omega_lower']
            check(w>0 and w*w<=ell-beta*beta,'Frequency lower bound')
            upper=(F(prior['stiffness_norm_inf_max'])+max((t for _,t in loads),default=F(0)))/mmin
            check(budget['spectral_lower']==ell and budget['spectral_upper']==upper,'Prior-only spectrum')
            if s==0:
                tail=amplitude*(1+beta/w)*upperexp(beta*horizon); quadrature=F(0); observation=sample
            else:
                tail=amplitude*upperexp(s*horizon)*(1+(1+beta/w)*upperexp(beta*horizon))
                fourth=amplitude*(s**4+((s+beta)**2+upper-beta*beta)**2*(1+beta/w))
                quadrature=s*horizon*step**4*fourth/180; observation=s*horizon*sample
            check(budget['tail']==tail and budget['quadrature']==quadrature,'Analytic tail/Simpson budget')
            check(budget['sample']==observation and budget['numerical']==numeric,'Sample/numerical budget')
            radius=tail+quadrature+observation+numeric
            check(lo==center-radius and hi==center+radius and budget['total']==radius,'Submitted enclosure')
            check(lo<=exact<=hi,'Exact response containment')
        check(all(seen[s]==set(expected) for s in seen),'Complete acquisition plan')
        result=record['reconstruction']
        check(result['status']=='enclosed','Operator enclosure status')
        for name,truth in [('K',k),('C',c),('M',m)]:
            matrix_count+=1
            box=result[name+'_intervals']
            for i in range(n):
                for j in range(n):
                    lo,hi=map(F,box[i][j]);check(lo<=truth[i][j]<=hi,'Operator truth containment')
    check(scalar_count==96 and matrix_count==12 and len(records)==4,'Frozen campaign counts')
    check(receipt['real_laplace_scalar_settings']==scalar_count and receipt['operator_enclosures']==matrix_count,'Receipt counts')
    return dict(schema='nonagon.dynamic-step.independent.v1',status='PASSED',checks=checks,
                fixtures=len(records),scalar_settings=scalar_count,matrix_enclosures=matrix_count,
                imports_acquisition_or_inverse_code=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',nargs='?',type=Path,default=Path(__file__).parent)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify(args.directory);text=json.dumps(result,indent=2)
    if args.output:args.output.write_text(text+'\n')
    print(text)
