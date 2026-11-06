#!/usr/bin/env python3
"""Uniform near-origin finite-production Euler product bound, using Arb.
A Cauchy envelope controls every uncomputed Taylor coefficient; no sampling.
"""
import sys,time,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'.optional-local-deps'))
import numpy as np
from flint import arb,acb,acb_series,ctx
ctx.prec=192;ctx.cap=40
N=20

def sieve(n):
 s=np.ones(n+1,dtype=bool);s[:2]=False
 for p in range(2,math.isqrt(n)+1):
  if s[p]:s[p*p::p]=False
 return np.flatnonzero(s)

def compute_products(Pmax):
 t=time.perf_counter();Z=acb_series([0,1],prec=N+1);A=acb_series([1],prec=N+1);B=acb_series([0],prec=N+1)
 for p in sieve(Pmax):
  p=int(p);lp=arb(p).log();E=(-Z*lp).exp()
  A=A*(1-(E-1)**2/arb(p-1)**2)
  B=B+lp**2/(p*(Z*lp).exp()-1)**2
 return A,B,time.perf_counter()-t

def envelope(T,rho):
 """Crude fully analytic Cauchy circle bound, valid also for infinite products.
 z-radius r = 2pi rho / log(T/2pi) <1/2. The containing box is
 used for regularized zeta derivative evaluation; radial elementary bounds for A,B.
 """
 L=(arb(T)/(2*arb.pi())).log();r=2*arb.pi()*rho/L
 box=acb(arb(0,r.upper()),arb(0,r.upper()))
 D=acb_series([acb(1)+box,acb(1)],prec=3).zeta(deflate=True)
 h=1+box*D[0];hp=D[0]+box*D[1];hpp=2*D[1]+2*box*D[2]
 H=h.abs_upper();Hlo=h.abs_lower();assert Hlo>0
 logpp=hpp.abs_upper()/Hlo+(hp.abs_upper()/Hlo)**2
 # n^(1-r)-1 >= n^(1-r) * (1-2^(r-1));
 # sum(log(n)^2 / n^(2-2r),n>=2) <= integral_1^infty + max term.
 # Safe monotone-integral bound using all integers intervals with n>=2:
 # f(n)<= integral_(n-1)^n [log(x+1)^2 / x^(2-2r)]dx
 # log(x+1)<=logx+log2 -> integral_1^infty exactly terms below.
 d=1-2*r; l2=arb(2).log()
 S=(2/d**3+2*l2/d**2+l2**2/d)/(1-arb(2)**(r-1))**2
 Aup=(r**2*S).exp();Bup=S
 Pup=(L*r).exp()*H**2*Aup
 M=1+((2+2*Pup)/r**2+2*logpp+2*Bup)/L**2
 return M,r,S

def row(T,A,B,rho=arb('0.01'),a=arb('0.0001'),c2_override=None):
 L=(arb(T)/(2*arb.pi())).log();alpha=2*arb.pi()/L
 Z=acb_series([0,1],prec=N+1)
 D=acb_series([1,1],prec=N+1).zeta(deflate=True);H=1+Z*D
 Hm=acb_series([H[k]*(-1)**k for k in range(N+1)],prec=N+1)
 P=(-L*Z).exp()*H*Hm*A
 logpp=H.log().derivative().derivative()
 # R is even and its constant coefficient is exactly zero by the displayed algebra.
 cs={}
 for k in range(2,N-2,2):
  q=2*(-P[k+2]+logpp[k]-B[k])/L**2*(acb(0,1)*alpha)**k
  assert q.imag.contains(0);cs[k]=q.real
 # Computed through x^16. Cauchy controls all even powers from x^18.
 M,r,S=envelope(T,rho)
 if c2_override is not None:cs[2]=c2_override
 maxk=max(cs);tail=(M*(a/rho)**(maxk+2)/(1-(a/rho)**2)).upper()
 # evenness permits geometric ratio(a/rho)^2. Terms x^0 and odd vanish exactly.
 dev=sum((cs[k].abs_upper()*a**(k-2) for k in cs if k>2),arb(0))+tail/a**2
 lower=cs[2].lower()-dev;upper=cs[2].upper()+dev
 return {'T':T if isinstance(T,(float,int)) else T.str(25,radius=True),'x2_coefficient_ball':cs[2].str(35,radius=True),'higher_coefficient_balls':{str(k):v.str(25,radius=True) for k,v in cs.items() if k>2},
  'uniform_R2_over_x2_lower':str(lower),'uniform_R2_over_x2_upper':str(upper),
  'cauchy_radius':str(rho),'cauchy_function_bound':str(M.upper()),'cauchy_remainder_abs_on_a':str(tail),'certified_positive':lower>0,'certified_c_le_3':upper<3,
  'variance_integral_remainder_bound_L16':str((arb(2)/3*upper*16*a**3).upper())}

def complete_euler_tail(A,B,Pmax,R=arb('0.1')):
    """Cauchy coefficients of the infinite tails, dominated by all integer n>P.
    Both analytic factors have real Taylor coefficients. No prime-density
    assumption or unverified prime-tail fit enters this envelope.
    """
    P=arb(Pmax);l=P.log();d=1-2*R
    J=P**(-d)*(l*l/d+2*l/d**2+2/d**3)
    SA=J/(1-1/P)**2
    umax=R**2*l*l*P**(-2+2*R)/(1-1/P)**2
    assert umax<1
    ML=R**2*SA/(1-umax)
    MB=J/(1-P**(-(1-R)))**2
    la=A.log(); aa=[];bb=[]
    for k in range(N+1):
        aa.append(la[k]+(acb(arb(0,(ML/R**k).upper())) if k>=2 else 0))
        bb.append(B[k]+acb(arb(0,(MB/R**k).upper())))
    return acb_series(aa,prec=N+1).exp(),acb_series(bb,prec=N+1),{'tail_circle_z_radius':str(R),'log_A_tail_sup_bound':str(ML.upper()),'B_tail_sup_bound':str(MB.upper()),'all_integer_integral_bound':str(J.upper())}

def continuous_range(A,B,Tlo,Thi):
    Tlo_exact=arb(str(Tlo));Thi_exact=arb(str(Thi))
    lo=(Tlo_exact/(2*arb.pi())).log();hi=(Thi_exact/(2*arb.pi())).log();L=lo.union(hi)
    Z=acb_series([0,1],prec=N+1)
    H=1+Z*acb_series([1,1],prec=N+1).zeta(deflate=True)
    Hm=acb_series([H[k]*(-1)**k for k in range(N+1)],prec=N+1)
    J=H*Hm*A;Q=J[4]-12*H.log()[4]+B[2]
    # d c2/dL = 8*pi²/L⁵ *(-J2 L²+3J3 L-4Q).
    numerator=(-J[2]*L**2+3*J[3]*L-4*Q).real
    assert numerator>0
    def c2(l):return (arb.pi()**2/3+4*arb.pi()**2*J[2]/l**2-8*arb.pi()**2*J[3]/l**3+8*arb.pi()**2*Q/l**4).real
    ci=c2(lo).union(c2(hi))
    rr=row(Tlo_exact.union(Thi_exact),A,B,c2_override=ci)
    rr['T_interval']=[Tlo,Thi];rr['endpoint_convention']='Exact decimal strings of the displayed archived endpoints';rr['c2_derivative_numerator_ball']=numerator.str(30,radius=True)
    rr['c2_derivative_numerator_lower']=str(numerator.lower());rr['c2_derivative_numerator_upper']=str(numerator.upper())
    rr['c2_coefficient_lower']=str(ci.lower());rr['c2_coefficient_upper']=str(ci.upper())
    return rr

def direct_check(T,x,Pmax):
    """Independent meromorphic formula, actual prime products at one point."""
    L=(arb(T)/(2*arb.pi())).log();r=2*arb.pi()*arb(x)/L;z=acb(0,r)
    Z=acb_series([1+z,1],prec=3).zeta()
    dlog=2*Z[2]/Z[0]-(Z[1]/Z[0])**2
    A=acb(1);B=acb(0)
    for p in sieve(Pmax):
        p=int(p);lp=arb(p).log();E=(-z*lp).exp()
        A=A*(1-(E-1)**2/arb(p-1)**2)
        B=B+lp**2/(p*(z*lp).exp()-1)**2
    term=dlog+(-L*z).exp()*(1+z).zeta()*(1-z).zeta()*A-B
    value=1+2*term.real/L**2
    return {'T':T,'x':x,'direct_finite_R2_ball':value.str(40,radius=True),'direct_finite_R2_over_x2_ball':(value/arb(x)**2).str(35,radius=True)}

def main():
 started=time.perf_counter()
 Pmax=int(sys.argv[1]) if len(sys.argv)>1 else 2000000
 A,B,dt=compute_products(Pmax)
 b=json.load(open(ROOT/'part_VI/bk_multi_100000_1000000_10000000.json'))
 Ts=[T for h in b['heights'] for T in h['T_quantiles']]
 AI,BI,tail=complete_euler_tail(A,B,Pmax)
 res={'precision_bits':ctx.prec,'prime_cutoff':Pmax,'prime_count':len(sieve(Pmax)),'product_seconds':dt,'scope':'finite production products plus independent all-integer envelopes for infinite Euler tails','rows':[row(T,A,B) for T in Ts],'infinite_euler_tail_envelope':tail,'infinite_product_rows':[row(T,AI,BI) for T in Ts],'continuous_height_range_infinite':continuous_range(AI,BI,min(Ts),max(Ts)),'independent_direct_checks':[direct_check(T,'0.0001',Pmax) for T in [Ts[1],Ts[4],Ts[7]]]}
 res['elapsed_seconds']=time.perf_counter()-started
 (Path(__file__).resolve().parent/'origin_cauchy_results.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
if __name__=='__main__':main()
