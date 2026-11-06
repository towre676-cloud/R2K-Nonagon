"""Certified continuous height average of the local BK variance model.
A Fubini reduction integrates the height variable analytically. This certifies
that averaging convention, not a finite-height local-stationarity theorem.
"""
from pathlib import Path
import sys,json,math,time,os
from flint import arb,acb,acb_series,ctx
import certified_campaign as cc
import origin_cauchy as oc
from native_euler import Euler
from shifted_tail import base_majorant
HERE=Path(__file__).resolve().parent;ctx.prec=128;ctx.cap=64;r0=arb(1)/64
FIXED=os.environ.get('R2K_HEIGHT_MODE')=='fixed'
RESULT=HERE/('finite_window/fixed_continuous_height_results.json' if FIXED else 'continuous_height_results.json')

def txt(x):return x.str(45,radius=True,more=True)
def save(p,d):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(d,indent=2));t.replace(p)
def K(q,r,L,lo,hi):
 def primitive(v):return (q*v).exp()*(L/q-r/(2*arb.pi())*(v/q-1/q**2))
 return primitive(hi)-primitive(lo)
def emom(j,lo,hi):
 def primitive(x):
  p=arb(1)
  for k in range(1,j+1):p=x**k-k*p
  return x.exp()*p
 return primitive(hi)-primitive(lo)

def main():
 start=time.time();models=[];zeros={}
 for n0 in [100000,1000000,10000000,20000000,50000000]:
  d=json.loads((HERE/'data'/f'zeros_{n0}.json').read_text());ta=arb(d['ordinate_balls'][0]);tb=arb(d['ordinate_balls'][-1]);la=(ta/(2*arb.pi())).log();last=(tb/(2*arb.pi())).log()
  for L0 in cc.LENGTHS:
   L=arb(L0)
   if FIXED:
    la=1+((arb(n0+1000)-arb(7)/8)/arb.const_e()).lambertw()
    lb=1+((arb(n0+19000)-L-arb(7)/8)/arb.const_e()).lambertw();D=18000-L
   else:
    lb=1+((last.exp()*(last-1)-L)/arb.const_e()).lambertw();D=lb.exp()*(lb-1)-la.exp()*(la-1)
   Rlo=2*arb.pi()*L/lb;Rhi=2*arb.pi()*L/la
   models.append({'n0':n0,'L':L0,'lo':la,'hi':lb,'D':D,'Rlo':Rlo,'Rhi':Rhi})
 E=Euler(2000000,128);td=json.loads((HERE/'tail_constants.json').read_text());S0=arb(td['S0_ball']);S2=arb(td['S2_ball'])
 oc.N=44;A,B,_=oc.compute_products(2000000);Z=acb_series([0,1],prec=45);H=1+Z*acb_series([1,1],prec=45).zeta(deflate=True)
 Hm=acb_series([H[k]*(-1)**k for k in range(45)],prec=45);J=H*Hm*A;dh=H.log().derivative().derivative()
 polynomials={0:[arb(0),arb(0),-arb(1)/2]}
 for k in range(2,42,2):
  sign=(-1)**(k//2);polynomials[k]=[-sign*J[k+2-j].real*((-1)**j)/math.factorial(j) for j in range(k+3)]
  polynomials[k][0]+=sign*(dh[k]-B[k]).real
 integrals=[]
 for m in models:
  lo,hi,D,L=m['lo'],m['hi'],m['D'],arb(m['L']);mom=[emom(j,lo,hi) for j in range(44)];s=arb(0)
  for k,poly in polynomials.items():
   a=sum((v*mom[j] for j,v in enumerate(poly)),arb(0));b=sum((v*mom[j+1] for j,v in enumerate(poly)),arb(0))
   s+=L*r0**(k+1)/(k+1)*a-r0**(k+2)/(2*arb.pi()*(k+2))*b
  s*=2/(arb.pi()*D)
  tlo=2*arb.pi()*lo.exp();thi=2*arb.pi()*hi.exp();M,_,_=oc.envelope(tlo.union(thi),arb('.1'));amin=2*arb.pi()/hi
  q=r0/(amin*arb('.1'));rem=M*q**42/(1-q*q)
  e=2/amin*(L*r0+r0*r0/(2*amin))*rem
  integrals.append(s+arb(0,e.upper()))
 cuts=[r0,arb(1)/32,arb(1)/16,arb(1)/8,arb(1)/4]
 for m in models:cuts.extend([m['Rlo'],m['Rhi']])
 mx=max(float(c.mid()) for c in cuts);cuts.extend(arb(k)/4 for k in range(2,math.ceil(mx*4)))
 cuts=sorted(cuts,key=lambda x:float(x.mid()));unique=[]
 for x in cuts:
  if not unique or float(x.mid())!=float(unique[-1].mid()):unique.append(x)
 cuts=unique;panels=[];nodecount=0
 checkpoint=RESULT
 if checkpoint.exists():
  old=json.loads(checkpoint.read_text())
  if old.get('prefix_integrals'):
   panels=old['panels'];integrals=[arb(x) for x in old['prefix_integrals']];nodecount=old['euler_evaluations']
 completed=len(panels)
 for idx,(a,b) in enumerate(zip(cuts,cuts[1:])):
  if idx<completed:continue
  c=(a+b)/2;h=(b-a)/2
  variable=any(float(b.mid())<=float(m['Rhi'].mid()) and float(c.mid())>=float(m['Rlo'].mid()) for m in models)
  R=(c*(arb('.25') if variable else arb('.8'))).min(arb('.25'));assert h<R
  q=(h/R).upper()
  _,meta=cc.disk_majorant(c,R);Db=arb(meta['D_upper']);Zb=arb(meta['zeta_product_upper']);Ab=arb(meta['A_upper']);Bb=arb(meta['B_upper'])
  rb=acb(arb(c.mid(),(R+c.rad()).upper()),arb(0,R.upper()));M=[];modes=[]
  for m in models:
   active=float(b.mid())<=float(m['Rhi'].mid());mode=float(c.mid())<float(m['Rlo'].mid());modes.append(mode)
   if not active:M.append(arb(0));continue
   cap=acb(m['hi']) if mode else 2*arb.pi()*m['L']/rb
   k1=K(acb(1),rb,m['L'],m['lo'],cap).abs_upper();kq=K(1-acb(0,1)*rb,rb,m['L'],m['lo'],cap).abs_upper()
   M.append((2/(arb.pi()*m['D'])*((Db+Bb)*k1+Zb*Ab*kq)).upper())
  n=8
  while not(4*h*max(M)*q**(2*n)/(1-q)<arb('1e-9')/(len(cuts)-1)):n+=1
  if n>100:raise RuntimeError('height quadrature budget exceeded')
  for k in range(n):
   rt,w=arb.legendre_p_root(n,k,weight=True);r=c+h*rt;AP,BP=E(r);zz=acb_series([acb(1,r),1],prec=3).zeta();Dlog=2*zz[2]/zz[0]-(zz[1]/zz[0])**2;C=zz[0]*zz[0].conjugate()*AP
   for j,m in enumerate(models):
    if not M[j]>0:continue
    cap=acb(m['hi']) if modes[j] else 2*arb.pi()*m['L']/r
    kernel=(Dlog-BP)*K(acb(1),r,m['L'],m['lo'],cap)+C*K(1-acb(0,r),r,m['L'],m['lo'],cap)
    integrals[j]+=h*w*2/(arb.pi()*m['D'])*kernel.real
  for j,m in enumerate(models):integrals[j]+=arb(0,(4*h*M[j]*q**(2*n)/(1-q)).upper())
  nodecount+=n;panels.append({'lo':txt(a),'hi':txt(b),'gauss_nodes':n,'cauchy_radius':txt(R),'majorant_max':txt(max(M))})
  save(checkpoint,{'status':'running','panels':panels,'prefix_integrals':[txt(x) for x in integrals],'euler_evaluations':nodecount})
  print('continuous height panel',idx+1,'/',len(cuts)-1,'nodes',n,'seconds',round(time.time()-start,1),flush=True)
 # Rigorous uniform-in-height infinite A-tail on the same finite windows.
 sigma=arb(1)/2;tcuts=[arb(0)]+[arb(1)/2**k for k in range(16,7,-1)]+[arb(k)/128 for k in range(1,65)]
 tcuts=sorted(tcuts,key=lambda x:float(x.mid()));leftbase=[(a,b,base_majorant(a,b,arb(0),arb(0),S0,S2)) for a,b in zip(tcuts,tcuts[1:]) if b>a]
 botbase=[(a,b,base_majorant(sigma,sigma,a,b,S0,S2)) for a,b in [(arb(0),r0)]+list(zip(cuts,cuts[1:]))]
 results=[]
 for j,m in enumerate(models):
  ell=m['lo'];amin=2*arb.pi()/m['hi'];L=arb(m['L'])
  left=sum((((b*b-a*a)/2)*(-ell*a).exp()*g for a,b,g in leftbase),arb(0))/amin
  bot=sum(((b-a)*(L+sigma/amin)*(-ell*sigma).exp()*g for a,b,g in botbase if float(b.mid())<=float(m['Rhi'].mid())),arb(0))
  right=arb(0)
  for a,b,g0 in leftbase:
   g=base_majorant(a,b,m['Rlo'],m['Rhi'],S0,S2);right+=(b*b-a*a)/2*(-ell*a).exp()*g/amin
  ae=(4/ell**2/amin*(left+bot+right)).upper();be=(S0/(2*arb.pi()**2)).upper();finite=L+integrals[j];full=finite+arb(0,(ae+be).upper())
  results.append({'n0':m['n0'],'L':m['L'],'ell_lower_ball':txt(m['lo']),'ell_upper_ball':txt(m['hi']),
    'finite_product_height_average_ball':txt(finite),'infinite_product_height_average_ball':txt(full),
    'A_tail_bound':txt(ae),'B_tail_bound':txt(be),'model_radius':txt(full.rad()),'radius_below_1e_minus_6':bool(full.rad()<arb('1e-6'))})
 save(checkpoint,{'status':'complete','panels':panels,'prefix_integrals':[txt(x) for x in integrals],'euler_evaluations':nodecount,'results':results,'all_targets_pass':all(r['radius_below_1e_minus_6'] for r in results),'elapsed_seconds':time.time()-start,
 'scope':('fixed coordinate spans; ' if FIXED else 'point-to-point spans; ')+'continuous average of local stationary BK variance; does not prove local stationarity of finite-height zeta correlations'})
if __name__=='__main__':main()
