#!/usr/bin/env python3
"""Independent rational differential controls for VIII.2.4 and bounded-model controls.

The automatic differentiator derives the acquisition Jacobian from the rational
inverse itself. No production reconstruction or interval code is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import random

HERE = Path(__file__).resolve().parent

class Jet:
    def __init__(self, value, derivative=(F(0),)*3):
        self.v=F(value); self.d=tuple(derivative)
    def __add__(self, other):
        other=other if isinstance(other,Jet) else Jet(other)
        return Jet(self.v+other.v,tuple(a+b for a,b in zip(self.d,other.d)))
    __radd__=__add__
    def __neg__(self): return Jet(-self.v,tuple(-a for a in self.d))
    def __sub__(self, other): return self+-other if isinstance(other,Jet) else self+Jet(-other)
    def __rsub__(self, other): return Jet(other)+-self
    def __mul__(self, other):
        other=other if isinstance(other,Jet) else Jet(other)
        return Jet(self.v*other.v,tuple(a*other.v+self.v*b for a,b in zip(self.d,other.d)))
    __rmul__=__mul__
    def __truediv__(self, other):
        other=other if isinstance(other,Jet) else Jet(other)
        return Jet(self.v/other.v,tuple((a*other.v-self.v*b)/other.v**2 for a,b in zip(self.d,other.d)))

checks=0

def require(condition, message):
    global checks
    checks+=1
    if not condition: raise RuntimeError(message)


def solve(matrix,rhs):
    a=[list(map(F,row))+[F(v)] for row,v in zip(matrix,rhs)]
    for k in range(len(a)):
        pivot=next(i for i in range(k,len(a)) if a[i][k])
        a[k],a[pivot]=a[pivot],a[k]
        d=a[k][k]; a[k]=[x/d for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                d=a[i][k];a[i]=[x-d*y for x,y in zip(a[i],a[k])]
    return [row[-1] for row in a]


def response(K,loads):
    a=[row[:] for row in K]
    for i,t in loads.items(): a[i][i]+=t
    return solve(a,[1]+[0]*(len(a)-1))[0]


def main():
    u_values=[F(1,7),F(1,2),F(1),F(7,4),F(5)]
    for b in [F(1,5),F(1),F(13,4)]:
        for q in [F(1,17),F(2,7),F(1)]:
            for u in u_values:
                t=u/b; baseline=F(7)
                values=[baseline,baseline-t*q/(1+t*b),baseline-2*t*q/(1+2*t*b)]
                raw=[Jet(v,tuple(F(int(k==i)) for k in range(3))) for i,v in enumerate(values)]
                d1=raw[0]-raw[1];d2=raw[0]-raw[2]
                recovered=d1*d2/(2*t*(d2-d1))
                expected=((F(3,2)+2*u)/t,-2*(1+u)**2/t,(1+2*u)**2/(2*t))
                require(recovered.v==q,'exact rational inverse')
                require(recovered.d==expected,'automatic acquisition Jacobian')
                require(sum(map(abs,recovered.d))*t==4*(1+u)**2,'bounded rawreading l1 norm')
                expected_variance=F(13,2)/u**2+24/u+34+24*u+8*u**2
                require(sum(x*x for x in recovered.d)/b**2==expected_variance,'independent rawreading l2 variance')
                require(4/u+8+4*u>=16,'bounded optimum')
    # Positive derivative of the Gaussianstationarity equation gives unique root.
    left,right=0.5,1.5
    def stationarity(u): return -13/u**3-24/u**2+24+16*u
    for _ in range(100):
        mid=(left+right)/2
        if stationarity(mid)<0:left=mid
        else:right=mid
    optimum=(left+right)/2
    require(abs(optimum-0.9720604037153)<1e-13,'Gaussian optimum')
    rng=random.Random(938101)
    midpoint_checks=0
    for fixture in range(12):
        K=[[F(0)]*3 for _ in range(3)]
        for i,j in [(0,1),(1,2),(0,2)]:
            w=F(rng.randint(1,8),rng.randint(1,6))
            K[i][i]+=w;K[j][j]+=w;K[i][j]-=w;K[j][i]-=w
        for i in range(3):K[i][i]+=F(rng.randint(1,8),rng.randint(1,6))
        g=solve(K,[1,0,0]);q=g[2]**2;eta=F(1,10**6);delta=2*eta/q
        Kprime=[r[:] for r in K];Kprime[2][2]+=delta
        for _ in range(30):
            loads={j:F(rng.randint(0,40),rng.randint(1,8)) for j in [1,2]}
            h1,h2=response(K,loads),response(Kprime,loads)
            midpoint=(h1+h2)/2
            require(0<=h1-h2<=2*eta,'uniformresponse bound at control load')
            require(abs(midpoint-h1)<=eta and abs(midpoint-h2)<=eta,'common bounded transcript')
            midpoint_checks+=1
        estimate=(K[2][2]+Kprime[2][2])/2
        require(abs(estimate-K[2][2])==eta/q and abs(estimate-Kprime[2][2])==eta/q,'two-point error separation')
    receipt={'status':'PASSED','predicates':checks,'exact_differential_fixtures':45,
             'midpoint_control_loads':midpoint_checks,'gaussian_optimum_u':optimum,
             'bounded_optimum_u':'1','bounded_factor':'16',
             'scope':'First-order exact-load single-node extraction; bounded midpoint controls corroborate the separate mathematical proof.'}
    (HERE/'finite_design_theory.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()
