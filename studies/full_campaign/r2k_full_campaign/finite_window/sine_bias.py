"""Exact finite-span centering bias under the stationary sine-kernel DPP.
Closed form through the autocorrelation of the trapezoidal anchor weight.
This is a conditional benchmark, not an established finite-height zeta law.
"""
from fractions import Fraction as Q
from pathlib import Path
import math,json
from flint import arb,acb,ctx
HERE=Path(__file__).resolve().parent;ctx.prec=192

def qa(q):q=Q(q);return arb(q.numerator)/q.denominator
def poly_on(mid,D,L):
 c=[Q(0)]*4;S=D+L
 for i in range(3):
  for j in range(3):
   shift=S-i*D-j*L
   if mid+shift<=0:continue
   a=Q((-1)**(i+j)*math.comb(2,i)*math.comb(2,j),6*D*D)
   for k in range(4):c[k]+=a*math.comb(3,k)*shift**(3-k)
 return c

def I(k,x):
 if x==0:return arb(0)
 a=2*arb.pi();v=a*x
 if k==0:return -(1-v.cos())/x+a*v.si()
 if k==1:return v.log()+arb.const_euler()-v.ci()
 if k==2:return x-v.sin()/a
 return x*x/2-x*v.sin()/a+(1-v.cos())/(a*a)

def bias(S,L):
 S=Q(S);L=Q(L);D=S-L;assert D>=L>0
 cuts=sorted({Q(0),L,D-L,D,D+L});integral=arb(0)
 for lo,hi in zip(cuts,cuts[1:]):
  if hi==lo:continue
  c=poly_on((lo+hi)/2,D,L)
  integral+=sum((qa(v)*(I(k,qa(hi))-I(k,qa(lo))) for k,v in enumerate(c)),arb(0))
 ca0=qa(L*L/D-L**3/(3*D*D));B=ca0-integral/arb.pi()**2;assert B>0
 d=qa(D);l=qa(L)
 if L>=1:up=(l*l+2*l*l/arb.pi()**2*(d/l).log()+(l*l-arb(1)/3)/arb.pi()**4)/d**2
 else:up=(l*l+2*l*l/arb.pi()**2*d.log()+2*l*l*(1-l)/arb.pi()**2+2*l**3/(3*arb.pi()**4))/d**2
 assert B<up
 return B,up

def independent_fourier(S,L):
 l=qa(L);D=qa(S)-l;K=arb(40)
 def kernel(k):return l*l*(D*k).sinc_pi()**2*(l*k).sinc_pi()**2
 f=acb.integral(lambda k,_:2*k*kernel(k),0,1,abs_tol=arb('1e-20'))+acb.integral(lambda k,_:2*kernel(k),1,K,abs_tol=arb('1e-20'))
 tail=2/(3*arb.pi()**4*D**2*K**3)
 return f.real.union(f.real+tail)

def main():
 checks=[]
 for S,L in [(3,Q(1,2)),(5,1),(10,2)]:
  b,u=bias(S,L);f=independent_fourier(S,L);assert f.contains(b)
  checks.append({'S':S,'L':str(L),'closed_form':b.str(45,radius=True,more=True),'independent_Fourier_enclosure':f.str(45,radius=True,more=True),'status':'PASS'})
 rows=[]
 for L in [Q(1,2),Q(1),Q(2),Q(4),Q(8),Q(16)]:
  b,u=bias(18000,L);rows.append({'span':18000,'L':float(L),'expected_downward_bias_ball':b.str(45,radius=True,more=True),'elementary_upper_bound':u.upper().str(45,radius=True,more=True)})
 out={'status':'PASS','independent_controls':checks,'production_fixed_span_benchmarks':rows,
 'scope':'stationary unit-intensity sine DPP in a fixed external interval; E(anchor variance)=population number variance minus this nonnegative quantity'}
 (HERE/'sine_bias_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
