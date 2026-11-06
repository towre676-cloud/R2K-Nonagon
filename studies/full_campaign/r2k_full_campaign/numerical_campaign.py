"""Fresh 3 x 20000-zero study, exact continuous anchors and stable BK curves.
Floating-point model quadrature is convergence evidence, not an interval proof.
"""
from pathlib import Path
import sys,json,time,hashlib,math,os
import numpy as np
from scipy.integrate import simpson
from flint import acb,acb_series,arb,ctx
BASE=Path(__file__).resolve().parent;ROOT=BASE.parent/'audit_work'
sys.path[:0]=[str(BASE.parent),str(ROOT/'part_IX'),str(ROOT/'part_X')]
from continuous_anchor_stats import statistics
from variance_quadrature import number_variance,sine_number_variance
import origin_cauchy as oc
ctx.prec=128;ctx.cap=40
LENGTHS=[.5,1.,2.,4.,8.,16.]
INDICES=[int(x) for x in os.environ.get('R2K_INDICES','100000,1000000,10000000').split(',')]
TAG=os.environ.get('R2K_CERT_TAG','')

def save(path,obj):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(path)
def stats(pts,L):
 s=statistics(pts,L);return {k:{'exact':str(v),'float':float(v)} for k,v in s.items()}
def get_zeros(n0):
 path=BASE/'data'/f'zeros_{n0}.json'
 if path.exists():return json.loads(path.read_text())
 start=time.time();zs=acb.zeta_zeros(n0,20000)
 # Retain enclosing balls and the exact binary64 ordinates used by the archive.
 vals=[float(z.imag.mid()) for z in zs]
 obj={'n0':n0,'count':20000,'precision_bits':ctx.prec,'elapsed_seconds':time.time()-start,
 'ordinate_balls':[z.imag.str(40,radius=True) for z in zs],'ordinates':vals,
 'ordinate_float_hex':[v.hex() for v in vals]}
 save(path,obj);print('zeros acquired',n0,round(obj['elapsed_seconds'],2),flush=True);return obj

def get_empirical(n0,zdata,stable=False):
 path=BASE/'data'/f"empirical_{'stable_' if stable else ''}{n0}.json"
 if path.exists():return json.loads(path.read_text())
 zs=np.array(zdata['ordinates']);u=zs/(2*np.pi)*np.log(zs/(2*np.pi*np.e))+7/8;u-=u[0]
 original=u.copy()
 if stable:
  twopi=2*arb.pi()
  def unfold(t):return t/twopi*(t/twopi).log()-t/twopi+arb(7)/8
  exact=[unfold(arb(t)) for t in zdata['ordinate_balls']];base=exact[0];u=np.array([float((x-base).mid()) for x in exact]);u[0]=0.
 rows=[]
 for L in LENGTHS:
  st=stats(u,L);anchors=np.arange(u[0],u[-1]-L,.01)
  counts=np.searchsorted(u,anchors+L)-np.searchsorted(u,anchors)
  blocks=[statistics(u[i*2000:(i+1)*2000],L) for i in range(10)]
  vs=[float(b['variance']) for b in blocks];se=float(np.std(vs,ddof=1)/np.sqrt(10))
  # Dependence sensitivity, with no claim of independent blocks.
  block_tiers={}
  for nb in [5,20]:
   k=20000//nb;vv=[float(statistics(u[i*k:(i+1)*k],L)['variance']) for i in range(nb)]
   block_tiers[str(nb)]={'block_variances':vv,'nominal_se':float(np.std(vv,ddof=1)/np.sqrt(nb))}
  rows.append({'L':L,'continuous':st,'grid_001_variance':float(np.var(counts)),
  'grid_001_mean':float(np.mean(counts)),'anchor_grid_shift':float(np.var(counts))-st['variance']['float'],
  'block_variances':vs,'nominal_block_se':se,'block_count_sensitivity':block_tiers})
  print('exact anchors',n0,L,flush=True)
 obj={'n0':n0,'unfolding':('Arb smooth Nbar(t)=t/(2pi)log(t/(2pi e))+7/8; subtract first in Arb, then round displacements to binary64' if stable else 'binary64 smooth Nbar(t)=t/(2pi)log(t/(2pi e))+7/8; subtract first value'),
 'stable_high_precision_before_subtraction':stable,'max_difference_from_original_unfolding':float(np.max(np.abs(u-original))),
 'unfolded_float_hex':[float(x).hex() for x in u],
 'T_quantiles':[float(np.quantile(zs,q)) for q in [.1,.5,.9]],'rows':rows}
 save(path,obj);return obj

class Evaluator:
 def __init__(self,Pmax):
  self.p=oc.sieve(Pmax).astype(float);self.lp=np.log(self.p);self.cut=int(np.searchsorted(self.p,200000,side='right'))
  self.A,self.B,_=oc.compute_products(Pmax);self.poly={}
 def coeff(self,T):
  ell=(arb(T)/(2*arb.pi())).log();alpha=2*arb.pi()/ell;Z=acb_series([0,1],prec=21)
  H=1+Z*acb_series([1,1],prec=21).zeta(deflate=True)
  Hm=acb_series([H[k]*(-1)**k for k in range(21)],prec=21)
  P=(-ell*Z).exp()*H*Hm*self.A;dh=H.log().derivative().derivative()
  return {k:float((2*(-P[k+2]+dh[k]-self.B[k])/ell**2*(acb(0,1)*alpha)**k).real.mid()) for k in range(2,18,2)}
 def value(self,x,T):
  ell=math.log(T/(2*math.pi));r=2*math.pi*x/ell
  if x<=.01:
   if T not in self.poly:self.poly[T]=self.coeff(T)
   main=sum(v*x**k for k,v in self.poly[T].items())
  else:main=None
  z=1j*r;phase=np.exp(-z*self.lp)
  logs=np.log1p(-(phase-1)**2/(self.p-1)**2)
  bterms=self.lp**2/(self.p/phase-1)**2
  aa=np.exp(np.sum(logs));bb=np.sum(bterms)
  al=np.exp(np.sum(logs[:self.cut]));bl=np.sum(bterms[:self.cut])
  S=acb_series([acb(1,r),1],prec=3).zeta();dlog=2*S[2]/S[0]-(S[1]/S[0])**2
  zz=S[0]*S[0].conjugate()*acb(0,-ell*r).exp()
  common=complex(dlog);prod=complex(zz)
  if main is None:main=1+2/ell**2*(common+prod*aa-bb).real
  low=1+2/ell**2*(common+prod*al-bl).real
  return main,low

def curve(n0,T,label,ev):
 path=BASE/'data'/f'curve_{n0}_{label}.npz'
 if path.exists():
  d=np.load(path);assert float(d['T'])==T,'cached height mismatch';return d['x'],d['main'],d['low']
 xs=np.linspace(0,16,1601);main=np.zeros(len(xs));low=main.copy();start=time.time()
 for i,x in enumerate(xs[1:],1):
  main[i],low[i]=ev.value(float(x),T)
  if i%400==0:print('curve',n0,label,i,'seconds',round(time.time()-start,1),flush=True)
 np.savez_compressed(path,x=xs,main=main,low=low,T=T)
 return xs,main,low

def main():
 start=time.time();emp=[]
 for n0 in INDICES:emp.append(get_empirical(n0,get_zeros(n0),stable=True))
 height_source=Path(os.environ['R2K_HEIGHTS_FILE']) if os.environ.get('R2K_HEIGHTS_FILE') else ROOT/'part_VI/bk_multi_100000_1000000_10000000.json'
 frozen=json.loads(height_source.read_text());frozen_T={h['n0']:h['T_quantiles'] for h in frozen['heights']}
 ev=Evaluator(2000000);result={'status':'running','scope':'fresh zeros, exact stored-dyadic anchor statistics, numerical full model quadrature','heights':[]}
 for e in emp:
  n0=e['n0'];curves=[curve(n0,T,str(i),ev) for i,T in enumerate(frozen_T[n0])];rows=[]
  for rr in e['rows']:
   L=rr['L'];bk=[number_variance(L,x,y) for x,y,_ in curves];bkavg=float(np.mean(bk))
   coarse=[number_variance(L,x[::2],y[::2]) for x,y,_ in curves]
   low=[number_variance(L,x,y) for x,_,y in curves]
   gue=number_variance(L,curves[0][0],1-np.sinc(curves[0][0])**2);ga=sine_number_variance(L)
   st=rr['continuous'];v=st['variance']['float'];se=rr['nominal_block_se']
   rows.append({**rr,'BK_quantiles':bk,'BK_height_avg':bkavg,'height_model_spread':max(bk)-min(bk),
   'BK_lowcut_height_avg':float(np.mean(low)),'prime_cutoff_shift':float(np.mean(low)-bkavg),
   'coarse_grid_shift':float(np.mean(coarse)-bkavg),'GUE_numeric':gue,'GUE_analytic':ga,'GUE_error':gue-ga,
   'residual':v-bkavg,'nominal_z':(v-bkavg)/se,'mean_square_residual':st['mean_squared_deviation_from_L']['float']-bkavg})
  result['heights'].append({'n0':n0,'T_quantiles':frozen_T[n0],'fresh_T_quantiles':e['T_quantiles'],'rows':rows})
  save(BASE/(TAG+'numerical_results.json'),result);print('height completed',n0,flush=True)
 result.update(status='complete',elapsed_seconds=time.time()-start)
 save(BASE/(TAG+'numerical_results.json'),result)
if __name__=='__main__':main()
