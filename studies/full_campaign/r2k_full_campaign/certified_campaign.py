"""Finite-product Gauss certificates plus rigorous infinite Euler-tail bounds.
All nodes, weights, point evaluations and analytic majorants use Arb.
The Gaussian remainder follows by exactness through degree 2n-1 and Cauchy.
"""
from pathlib import Path
import sys,json,time,math,ctypes,os
import numpy as np
from flint import arb,acb,acb_series,ctx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent/'audit_work'
sys.path.insert(0,str(ROOT/'part_X'))
import origin_cauchy as oc
from native_euler import Euler
ctx.prec=128;ctx.cap=64
P=2000000;Q=100000000;r0=arb(1)/64
LENGTHS=[.5,1.,2.,4.,8.,16.]
TAG=os.environ.get('R2K_CERT_TAG','')
RESULT=HERE/(TAG+'certified_results.json')
height_source=Path(os.environ['R2K_HEIGHTS_FILE']) if os.environ.get('R2K_HEIGHTS_FILE') else ROOT/'part_VI/bk_multi_100000_1000000_10000000.json'
archive=json.loads(height_source.read_text())
specs=[{'n0':h['n0'],'quantile':i,'T':T,'ell':(arb(T)/(2*arb.pi())).log()} for h in archive['heights'] for i,T in enumerate(h['T_quantiles'])]
for s in specs:s['alpha']=2*arb.pi()/s['ell']

def txt(x):return x.str(45,radius=True,more=True)
def val(x):return float(x.mid())
def save(path,data):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2));tmp.replace(path)

def tail_sums(E):
 path=HERE/'tail_constants.json'
 if path.exists():
  d=json.loads(path.read_text());return arb(d['S0_ball']),arb(d['S2_ball']),d
 ps=oc.sieve(Q).astype(np.uint64);ps=ps[ps>P];out=(ctypes.c_void_p*2)()
 f=E.lib.tail_sums;f.argtypes=[ctypes.POINTER(ctypes.c_uint64),ctypes.c_size_t,ctypes.POINTER(ctypes.c_void_p)];f.restype=ctypes.c_int
 assert f(ps.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),len(ps),out)==0
 v=[]
 for ptr in out:v.append(arb(ctypes.string_at(ptr).decode()));E.lib.release_string(ptr)
 # Every prime > Q is coprime to 2310. For decreasing f,
 # each allowed residue progression contributes <= f(Q)+W^-1 int_Q^infty f.
 W=2310;phi=480;q=arb(Q);density=arb(phi)/W;lp=q.log()
 j0=1/q;j2=(lp**2+2*lp+2)/q
 fac=(1-1/arb(P))**-2
 rem0=fac*(phi/q**2+density*j0)
 rem2=fac*(phi*lp**2/q**2+density*j2)
 S0=v[0]+arb(0,rem0.upper());S2=v[1]+arb(0,rem2.upper())
 # Positive tails: use nonnegative interval enclosure, not centered uncertainty.
 S0=v[0].union(v[0]+rem0);S2=v[1].union(v[1]+rem2)
 d={'cutoff':P,'exact_prime_sum_end':Q,'primes_summed':len(ps),'wheel':W,'wheel_totient':phi,
 'finite_S0':txt(v[0]),'finite_S2':txt(v[1]),'tail_S0_upper':txt(rem0.upper()),'tail_S2_upper':txt(rem2.upper()),
 'S0_ball':txt(S0),'S2_ball':txt(S2),'B_integrated_tail_bound':txt((S0/(2*arb.pi()**2)).upper())}
 save(path,d);return S0,S2,d

small=[int(p) for p in oc.sieve(1000)]
majorants={}
def euler_majorant(R):
 key=txt(R)
 if key in majorants:return majorants[key]
 d=1-2*R;assert d>0
 A=arb(1)
 for p in small:A*=1+(arb(p)**R+1)**2/arb(p-1)**2
 A*= (4*arb(1000)**(-d)/(d*(1-arb(1)/1000)**2)).exp()
 l2=arb(2).log();B=(2/d**3+2*l2/d**2+l2*l2/d)/(1-arb(2)**(R-1))**2
 majorants[key]=(A.upper(),B.upper());return majorants[key]

def disk_majorant(c,R):
 # Cover the entire square enclosing the Cauchy disk, proving zero exclusion
 # throughout it. Refinement is an enclosure operation, not point sampling.
 for parts in [1,2,4,8,16,32,64]:
  Dmax=arb(0);Zmax=arb(0);ok=True
  for j in range(parts):
   for k in range(parts):
    re=(c-R)+(2*j+1)*R/parts;im=-R+(2*k+1)*R/parts;rad=R/parts
    rb=acb(arb(re.mid(),(rad+re.rad()).upper()),arb(im.mid(),(rad+im.rad()).upper()));z=acb(0,1)*rb
    zp=acb_series([1+z,1],prec=3).zeta();zm=acb_series([1-z,1],prec=3).zeta()
    if not all(v.is_finite() for v in list(zp.coeffs())+list(zm.coeffs())) or not(zp[0].abs_lower()>0 and zm[0].abs_lower()>0):ok=False;break
    Dp=2*zp[2]/zp[0]-(zp[1]/zp[0])**2;Dm=2*zm[2]/zm[0]-(zm[1]/zm[0])**2
    Dmax=Dmax.max(Dp.abs_upper()).max(Dm.abs_upper());Zmax=Zmax.max(zp[0].abs_upper()*zm[0].abs_upper())
   if not ok:break
  if ok:
   A,B=euler_majorant(R)
   Ms=[(1+(2*Dmax+2*(s['ell']*R).exp()*Zmax*A+2*B)/s['ell']**2).upper() for s in specs]
   return Ms,{'square_subdivisions_per_axis':parts,'D_upper':txt(Dmax),'zeta_product_upper':txt(Zmax),'A_upper':txt(A),'B_upper':txt(B)}
 raise RuntimeError('Could not prove zero-free analytic Cauchy neighborhood')

def AP_real_bound(lo,hi):
 rb=lo.union(hi);z=acb(0,rb);A=arb(1)
 for p in small:
  phase=(-z*arb(p).log()).exp();A*= (1-(phase-1)**2/arb(p-1)**2).abs_upper()
 return (A*(arb(4)/999).exp()).upper()

def A_tail_interval(lo,hi,S0,S2):
 z=acb(0,lo.union(hi));H=1+z*acb_series([1+z,1],prec=1).zeta(deflate=True)[0]
 coeff=S2.upper() if not(lo>0) else S2.upper().min(4*S0.upper()/lo**2)
 return (H.abs_upper()**2*AP_real_bound(lo,hi)*coeff*(4*S0.upper()).exp()).upper()

def main():
 start=time.time();E=Euler(P,128);S0,S2,td=tail_sums(E)
 # One common r grid reuses the expensive Euler evaluation across all heights.
 endpoints={};cuts=[r0,arb(1)/32,arb(1)/16,arb(1)/8,arb(1)/4]
 for j,s in enumerate(specs):
  for L in LENGTHS:
   r=s['alpha']*L;cuts.append(r);endpoints.setdefault(float(r.mid()),[]).append((j,L))
 rmax=max(float(x.mid()) for x in cuts)
 cuts += [arb(k)/4 for k in range(2,math.ceil(rmax*4))]
 cuts=sorted(cuts,key=lambda a:float(a.mid()));unique=[]
 for x in cuts:
  if not unique or float(x.mid())!=float(unique[-1].mid()):unique.append(x)
 cuts=unique;alpha_min=min(s['alpha'].lower() for s in specs)
 print('panels',len(cuts)-1,'tail constants',td['S0_ball'],flush=True)
 # Origin polynomial, with a Cauchy bound for its omitted even coefficients.
 oc.N=44;AA,BB,dt=oc.compute_products(P);Z=acb_series([0,1],prec=45)
 H=1+Z*acb_series([1,1],prec=45).zeta(deflate=True);Hm=acb_series([H[k]*(-1)**k for k in range(45)],prec=45)
 dh=H.log().derivative().derivative();I0=[];I1=[];tailA=[]
 for s in specs:
  ell=s['ell'];alpha=s['alpha'];PP=(-ell*Z).exp()*H*Hm*AA
  cs={k:(2*(-PP[k+2]+dh[k]-BB[k])/ell**2*acb(0,1)**k).real for k in range(2,42,2)}
  u=sum((v*r0**(k+1)/(k+1) for k,v in cs.items()),arb(0));v=sum((v*r0**(k+2)/(k+2) for k,v in cs.items()),arb(0))
  M,_,_=oc.envelope(arb(s['T']),arb('.1'));q=r0/(alpha*arb('.1'));assert q<1
  rem=M*q**42/(1-q*q)
  I0.append(u+arb(0,(r0*rem).upper()));I1.append(v+arb(0,(r0*r0*rem/2).upper()))
  tailA.append(arb(0))
 # Bound A error on the origin interval separately.
 aorig=A_tail_interval(arb(0),r0,S0,S2)
 for j,s in enumerate(specs):tailA[j]=(r0*aorig).upper()
 panel_records=[];results=[];nodes_total=0
 checkpoint=RESULT
 if checkpoint.exists() and '--fresh' not in sys.argv:
  old=json.loads(checkpoint.read_text())
  if 'prefix_I0' in old:
   panel_records=old['panels'];results=old['results'];nodes_total=old['gauss_euler_evaluations'];I0=[arb(x) for x in old['prefix_I0']];I1=[arb(x) for x in old['prefix_I1']]
 completed=len(panel_records)
 for panel,(lo,hi) in enumerate(zip(cuts,cuts[1:])):
  if panel<completed:continue
  c=(lo+hi)/2;h=(hi-lo)/2;R=(c*arb('.8')).min(arb('.25'));assert c>R and h<R
  Ms,meta=disk_majorant(c,R);Mmax=max(Ms);q=(h/R).upper();n=8
  maxweight=2/alpha_min*(16+(c+R)/alpha_min)
  while not(4*h*Mmax*q**(2*n)/(1-q)*maxweight < arb('1e-9')/(len(cuts)-1)):n+=1
  if n>64:raise RuntimeError('quadrature budget exceeded')
  weights=[arb.legendre_p_root(n,k,weight=True) for k in range(n)]
  p0=[arb(0) for s in specs];p1=[arb(0) for s in specs]
  for root,w in weights:
   r=c+h*root;A,B=E(r);zz=acb_series([acb(1,r),1],prec=3).zeta()
   D=2*zz[2]/zz[0]-(zz[1]/zz[0])**2;C=zz[0]*zz[0].conjugate()*A
   for j,s in enumerate(specs):
    Rval=1+2/s['ell']**2*(D+acb(0,-s['ell']*r).exp()*C-B).real
    p0[j]+=h*w*Rval;p1[j]+=h*w*r*Rval
  nodes_total+=n
  at=A_tail_interval(lo,hi,S0,S2)
  for j,s in enumerate(specs):
   err=(4*h*Ms[j]*q**(2*n)/(1-q)).upper()
   I0[j]+=p0[j]+arb(0,err);I1[j]+=p1[j]+arb(0,((c+R)*err).upper())
   tailA[j]+=(hi-lo)*at
  # For each target, compute the sharper weighted A-tail prefix anew below.
  record={'lo':txt(lo),'hi':txt(hi),'cauchy_radius':txt(R),'gauss_nodes':n,'q_upper':txt(q),'R_majorant_max':txt(Mmax),'A_tail_unweighted_majorant':txt(at),**meta}
  panel_records.append(record)
  for j,L in endpoints.get(float(hi.mid()),[]):
   s=specs[j];alpha=s['alpha'];finite=L-L*L+2/alpha*(L*I0[j]-I1[j]/alpha)
   # Recompute weighted upper sum from certified interval majorants.
   weighted=L*r0*aorig
   for pr in panel_records:
    a=arb(pr['lo']);b=arb(pr['hi']);g=arb(pr['A_tail_unweighted_majorant']);weighted+=(b-a)*(L-a/alpha)*g
   aerr=(4/s['ell']**2/alpha*weighted).upper();berr=(S0/(2*arb.pi()**2)).upper()
   full=finite+arb(0,(aerr+berr).upper())
   results.append({'n0':s['n0'],'quantile':s['quantile'],'T':s['T'],'L':L,'finite_product_variance_ball':txt(finite),
   'finite_product_radius':txt(finite.rad()),'A_integrated_infinite_tail_bound':txt(aerr),'B_integrated_infinite_tail_bound':txt(berr),
   'infinite_product_variance_ball':txt(full),'infinite_product_radius':txt(full.rad())})
  data={'status':'running','prime_cutoff':P,'precision_bits':ctx.prec,'tail_constants':td,'panels':panel_records,'results':results,'gauss_euler_evaluations':nodes_total,'prefix_I0':[txt(v) for v in I0],'prefix_I1':[txt(v) for v in I1]}
  save(RESULT,data)
  print('panel',panel+1,'/',len(cuts)-1,'nodes',n,'targets',len(results),'seconds',round(time.time()-start,1),flush=True)
 data.update(status='complete',elapsed_seconds=time.time()-start)
 assert len(results)==len(specs)*len(LENGTHS)
 save(RESULT,data)
if __name__=='__main__':main()
