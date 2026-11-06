#!/usr/bin/env python3
"""Independent exact checks of the finite-load first-order objective.

This checker imports neither planner nor benchmark. Its dual-number derivative
is evaluated through the inverse as a rational function of three raw readings.
It does not establish the nonlinear interval inference algorithm.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent

class Dual:
    def __init__(self, value, derivative=(F(0),F(0),F(0))):
        self.value=F(value); self.derivative=tuple(map(F,derivative))
    @staticmethod
    def coerce(value):
        return value if isinstance(value,Dual) else Dual(value)
    def __add__(self, other):
        b=self.coerce(other)
        return Dual(self.value+b.value, [x+y for x,y in zip(self.derivative,b.derivative)])
    __radd__=__add__
    def __neg__(self): return Dual(-self.value,[-x for x in self.derivative])
    def __sub__(self, other): return self+-self.coerce(other)
    def __rsub__(self, other): return self.coerce(other)+-self
    def __mul__(self, other):
        b=self.coerce(other)
        return Dual(self.value*b.value,[x*b.value+self.value*y for x,y in zip(self.derivative,b.derivative)])
    __rmul__=__mul__
    def __truediv__(self, other):
        b=self.coerce(other)
        return Dual(self.value/b.value,[(x*b.value-self.value*y)/(b.value*b.value) for x,y in zip(self.derivative,b.derivative)])
    def __rtruediv__(self, other): return self.coerce(other)/self


def require(value, message):
    if not value: raise RuntimeError(message)


def run():
    checks=0; fixtures=0
    for t in [F(1,8),F(1,2),F(1),F(4),F(16)]:
        for u in [F(1,32),F(1,7),F(1,2),F(1),F(7,3),F(16)]:
            for q in [F(1,10000),F(1,7),F(2)]:
                b=u/t; h0=1+q/b
                raw=[h0,h0-t*q/(1+u),h0-2*t*q/(1+2*u)]
                hs=[Dual(x,[F(i==j) for j in range(3)]) for i,x in enumerate(raw)]
                d1=hs[0]-hs[1];d2=hs[0]-hs[2]
                recovered_q=d1*d2/(2*t*(d2-d1))
                recovered_b=(2*d1-d2)/(2*t*(d2-d1))
                require(recovered_q.value==q,'Inverse did not recover q');checks+=1
                require(recovered_b.value==b,'Inverse did not recover b');checks+=1
                claimed=((3+4*u)/(2*t),-2*(1+u)**2/t,(1+2*u)**2/(2*t))
                require(recovered_q.derivative==claimed,'Raw-reading gradient mismatch');checks+=3
                l1=sum(abs(x) for x in recovered_q.derivative)/q
                claimed_l1=(b/q)*(4/u+8+4*u)
                require(l1==claimed_l1,'Bounded-error coefficient mismatch');checks+=1
                l2=sum(x*x for x in recovered_q.derivative)/(q*q)
                claimed_l2=(b/q)**2*(F(13,2)/(u*u)+24/u+34+24*u+8*u*u)
                require(l2==claimed_l2,'Shared-baseline Gaussian coefficient mismatch');checks+=1
                # Equality gives the global minimum for every positive u.
                require(4/u+8+4*u-16==4*(u-1)**2/u,'Bounded optimum identity failed');checks+=1
                # Two drops share their baseline: Cov(d1,d2)=sigma^2.
                g1=d2.value*d2.value/(2*t*(d2.value-d1.value)**2)
                g2=-d1.value*d1.value/(2*t*(d2.value-d1.value)**2)
                covariance_variance=2*g1*g1+2*g2*g2+2*g1*g2
                require(covariance_variance==sum(x*x for x in recovered_q.derivative),'Drop covariance omitted');checks+=1
                # The two actual single loads have independent relative errors.
                # Differentiate through their physical forward drops, then the
                # unchanged nominal inverse; this measures calibration sensitivity.
                ta=Dual(t,[t,0,0]);tb=Dual(2*t,[0,2*t,0])
                dropa=ta*q/(1+ta*b);dropb=tb*q/(1+tb*b)
                q_cal=dropa*dropb/(2*t*(dropb-dropa))
                require(q_cal.value==q,'Calibration reference q moved');checks+=1
                require(q_cal.derivative==(2*q,-q,F(0)),'Relative shunt sensitivity mismatch');checks+=2
                require(sum(abs(x) for x in q_cal.derivative)==3*q,'Shunt worst-case coefficient mismatch');checks+=1
                fixtures+=1
    # Strict convexity f''=39/u^4+48/u^3+16 > 0 gives the unique optimum.
    lo,hi=F(9,10),F(11,10)
    def derivative_poly(u): return 16*u**4+24*u**3-24*u-13
    require(derivative_poly(lo)<0<derivative_poly(hi),'Gaussian root not bracketed');checks+=1
    for _ in range(128):
        mid=(lo+hi)/2
        if derivative_poly(mid)<0:lo=mid
        else:hi=mid
    require(F(9720603,10000000)<lo<hi<F(9720605,10000000),'Gaussian optimum moved');checks+=1
    # Treating the drops as independent double-counts common-baseline error.
    require(2*g1*g1+2*g2*g2!=covariance_variance,'Covariance negative control failed');checks+=1
    return {'schema':'nonagon.load-objective-independent-checks.v1','status':'passed',
            'fixtures':fixtures,'exact_predicates':checks,
            'bounded_optimum_u':'1','bounded_minimum_factor':'16','independent_relative_load_factor':'3',
            'gaussian_optimum_u_interval':[str(lo),str(hi)],
            'gaussian_optimum_u_decimal':float((lo+hi)/2),
            'derivative_source':'independent dual differentiation of rational inverse on raw h0,h(t),h(2t)',
            'scope':'exact finite-protocol first derivatives; nonlinear enclosures remain the separate production gate',
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=HERE/'independent_checks.json');a=ap.parse_args()
    result=run();a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
