"""Symbolic cone proof of every trace span, plus exact convolution controls."""
from trace_derivation import PARTS,cumulant
from fractions import Fraction as Q
from collections import defaultdict
from pathlib import Path
import json,random

def symbolic():
    # j=k+d, k,d>=0. Every partial sum has two integer coefficients.
    freq=[(1,1),(-1,-1),(0,1),(0,-1)]
    coeff=defaultdict(Q)
    for part in PARTS:
        p=(0,0);points=[p]
        for block in part:
            p=tuple(p[q]+sum(freq[i][q] for i in block) for q in range(2));points.append(p)
        hi=[p for p in points if all(all(p[q]>=r[q] for q in range(2)) for r in points)]
        lo=[p for p in points if all(all(p[q]<=r[q] for q in range(2)) for r in points)]
        assert hi and lo,'needs another cone subdivision'
        span=tuple(hi[0][q]-lo[0][q] for q in range(2))
        coeff[span]+=Q((-1)**(len(part)-1),len(part))
    assert {p:c for p,c in coeff.items() if c}=={(1,1):Q(2),(1,2):Q(-1),(1,0):Q(-1)}
    return {str(p):str(c) for p,c in coeff.items() if c}

def convolution_controls():
    rng=random.Random(20261105);checks=0
    for N in [1,2,3,5,10]:
        K=12
        for _ in range(5):
            u=[Q(0)]+[Q(rng.randrange(1,10),rng.randrange(1,10)) for _ in range(K)]
            v=[Q(0)]+[Q(rng.randrange(1,10),rng.randrange(1,10)) for _ in range(K)]
            direct=sum((u[j]*v[k]*(cumulant(N,j,k)+(min(j,N)**2 if j==k else 0)) for j in range(1,K+1) for k in range(1,K+1)),Q(0))
            diag=sum((u[j]*v[j]*min(j,N)**2 for j in range(1,K+1)),Q(0))
            mx=2*sum((sum(u[1:m+1])*sum(v[1:m+1]) for m in range(1,N)),Q(0))
            hs=sum((u[j]*v[k]*max(N-j-k,0) for j in range(1,K+1) for k in range(1,K+1)),Q(0))
            ds=sum((u[j]*v[k]*max(N-abs(j-k),0) for j in range(1,K+1) for k in range(1,K+1)),Q(0))
            assert direct==diag+mx-hs-ds
            checks+=1
    return checks

if __name__=='__main__':
    out={'status':'PASS','all_75_spans_symbolically_proved_on_j_ge_k_cone':symbolic(),
         'exact_rational_convolution_controls':convolution_controls(),
         'scope':'mathematical identity plus finite arithmetic controls; no zeta distribution claim'}
    Path(__file__).with_name('trace_validation.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
