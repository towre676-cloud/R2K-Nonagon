#!/usr/bin/env python3
"""Part X: rigorous Arb comparison on archived contour moment balls.
PYTHONPATH=.cert_deps python experiments/certification.py
The archive inputs are contour evidence, not trusted point samples. Full analytic
recomputation is separately requested on the first-zero certificate.
"""
import sys, json, time, hashlib, math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.cert_deps'))
from flint import arb, acb, arb_mat, ctx
ROOT=Path(__file__).resolve().parents[1]
V9=ROOT if (ROOT/'part_VI').is_dir() else ROOT/'source/v9/r2k_nonagon_v9.0'
sys.path.insert(0,str(V9/'part_VI'))
import r2kverify

def sign(v):
    return 1 if v>0 else -1 if v<0 else 0

def mat(rows):return arb_mat(rows)
def sub(A,rr,cc):return mat([[A[i,j] for j in cc] for i in rr])
def nested(A):
    """Haynsworth, balanced split; single scalar pivots must exclude zero."""
    n=A.nrows()
    if n==1:
        q=sign(A[0,0]); assert q
        return [int(q>0),int(q<0)], [float(A[0,0].abs_lower())]
    k=n//2
    P=sub(A,range(k),range(k)); B=sub(A,range(k),range(k,n)); Q=sub(A,range(k,n),range(k,n))
    S=Q-B.transpose()*P.solve(B)
    i,p=nested(P); j,q=nested(S)
    return [i[0]+j[0],i[1]+j[1]],p+q

def moments(c):
    ks=c['moments_k']; totals=[acb(0) for _ in ks]
    for s in c['segments']:
        totals=[u+r2kverify.ball_from(v) for u,v in zip(totals,s['contrib'])]
    if c['symmetric_half_path']: return [x.imag/arb.pi() for x in totals]
    return [(x/(2*acb.pi()*acb(0,1))).real for x in totals]

def power_y(mu,K):
    # x^k = 2^-k sum_j binom(k,j) T_|k-2j|(x), y=(x+1)/2
    mx=[]
    for k in range(K+1):
        mx.append(sum((math.comb(k,j)*mu[abs(k-2*j)] for j in range(k+1)),arb(0))/arb(2)**k)
    return [sum((math.comb(k,j)*mx[j] for j in range(k+1)),arb(0))/arb(2)**k for k in range(K+1)]

def condensation(m,n):
    """All contiguous Hankel minors; n^2+O(n) arithmetic, zero denominator rejects.
    D_r^s=det(m[s+i+j],0<=i,j<r).
    """
    D={0:[arb(1)]*(2*n+1),1:list(m)}; ops=0
    out=[arb(1),m[0]]
    for r in range(1,n):
        vals=[]
        for s in range(2*(n-r)-1):
            den=D[r-1][s+2]; assert sign(den), 'zero/uncertain condensation pivot'
            vals.append((D[r][s]*D[r][s+2]-D[r][s+1]**2)/den);ops+=1
        D[r+1]=vals;out.append(vals[0])
    return out[1:],ops

def run_record(p):
    c=json.load(open(p));ctx.prec=max(128,c['method']['prec_bits'])
    r2kverify.integrity(c)
    mm=moments(c);n=int(round(float(mm[0].mid())))
    if n==0:return None
    if 2*n-2>=len(mm):return None
    H=mat([[(mm[i+j]+mm[abs(i-j)])/2 for j in range(n)] for i in range(n)])
    t=time.perf_counter(); flat=[sub(H,range(k),range(k)).det() for k in range(1,n+1)];tf=time.perf_counter()-t
    sg=[1]+[sign(v) for v in flat]
    if 0 in sg:return {'file':p.name,'status':'flat undecided'}
    neg=sum(a!=b for a,b in zip(sg,sg[1:]));fi=[n-neg,neg]
    t=time.perf_counter(); ni,piv=nested(H);tn=time.perf_counter()-t
    assert fi==ni
    yy=power_y(mm,2*n-2)
    t=time.perf_counter(); qd,ops=condensation(yy,n);tq=time.perf_counter()-t
    # T_k(2y-1): leading coefficient 2^(2k-1) for k>=1.
    # squared product = 2^(2(n-1)^2), including k=0 leading1.
    qdh=[v*arb(2)**(2*(k-1)**2) for k,v in enumerate(qd,1)]
    assert all(x.overlaps(y) for x,y in zip(flat,qdh))
    assert all(sign(x)==sign(y) for x,y in zip(flat,qdh))
    return {'file':p.name,'status':'matched','function':c['function'],'count':n,'flat_inertia':fi,'nested_inertia':ni,
            'qd_minor_signs':[sign(v) for v in qdh], 'flat_seconds':tf,'nested_seconds':tn,'condensation_seconds':tq,
            'condensation_updates':ops,'minimum_nested_pivot_abs_lower':min(piv),'flat_minor_balls':[v.str(20,radius=True) for v in flat],
            'condensation_minor_balls':[v.str(20,radius=True) for v in qdh], 'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    start=time.perf_counter();rows=[];fails=[]
    for p in sorted((V9/'part_VI/certs').glob('*.json')):
        c=json.load(open(p))
        if not isinstance(c,dict):continue
        if c.get('format') not in ('r2k-contour-certificate/1','r2k-contour-certificate/2'):continue
        try:
            rr=run_record(p)
            if rr:rows.append(rr)
        except Exception as e:fails.append({'file':p.name,'reason':str(e)})
    # Positive moment matrix at nodes -2,-1,1,2 has zero odd moments.
    # Raw condensation encounters mu_3=0 internally despite positive leading minors.
    symmetric=[arb(4) if k==0 else arb(2*(1+2**k) if k%2==0 else 0) for k in range(7)]
    zero_pivot=False
    try:condensation(symmetric,4)
    except AssertionError:zero_pivot=True
    fallback=[arb_mat([[symmetric[i+j] for j in range(k)] for i in range(k)]).det() for k in range(1,5)]
    assert all(sign(v)>0 for v in fallback)
    fresh=V9/'part_VI/certs/zero_1_b0p45.json'
    c=json.load(open(fresh));vf=r2kverify.verify(c,full=True,evaluator='series')
    result={'arithmetic':'python-flint 0.9.0 Arb, minimum128 bits','scope':'archived contour moment balls; one first-zero contour independently recomputed',
        'records':rows,'failures':fails,'records_matched':sum(x.get('status')=='matched' for x in rows),
        'off_axis_controls':[x['file'] for x in rows if x.get('flat_inertia',[0,0])[1]>0],
        'full_first_zero_verifier':vf,'zero_pivot_rejected':zero_pivot,'zero_pivot_exact_positive_fallback_minors':[str(v) for v in fallback],'elapsed_seconds':time.perf_counter()-start}
    (Path(__file__).resolve().parent/'certification_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()
