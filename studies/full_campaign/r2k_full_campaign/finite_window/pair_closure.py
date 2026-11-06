"""Exact finite-window pair/anchor closure and certified translation correction.
Geometric expansion of 1/(S-r) evaluates weighted pairs with exact moments.
"""
from pathlib import Path
from fractions import Fraction as Q
from bisect import bisect_left,bisect_right
import json,math,sys,time
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;ctx.prec=128
sys.path.insert(0,str(BASE.parent))
from continuous_anchor_stats import statistics

def qa(q):
 q=Q(q);return arb(q.numerator)/q.denominator

def grid(points,lengths):
 qs=[Q(float(x)) for x in points];ls=[Q(x) for x in lengths]
 e=max(q.denominator.bit_length()-1 for q in qs+ls);scale=1<<e
 return [int(q*scale) for q in qs],scale

def analyze(points,lengths,expected=None,span=None,eps=Q(0)):
 x,scale=grid(points,lengths);S=int(Q(span)*scale) if span is not None else x[-1];assert x[0]>=0 and x[-1]<=S
 maxlen=int(Q(max(lengths))*scale);margin=math.ceil(2*eps*scale)
 gaps=[]
 for i,a in enumerate(x):
  j=i+1
  while j<len(x) and x[j]-a<maxlen+margin+1:
   gaps.append(x[j]-a);j+=1
 gaps.sort();mom=[0]*7;at=0;rows=[]
 for L0 in lengths:
  L=Q(L0);li=int(L*scale);D=S-li;assert D>0
  while at<len(gaps) and gaps[at]<li:
   p=1;g=gaps[at]
   for k in range(7):mom[k]+=p;p*=g
   at+=1
  T=li*mom[0]-mom[1]
  left=sum(j*(li-a) for j,a in enumerate(x) if a<li)
  right=sum((len(x)-1-i)*(a-D) for i,a in enumerate(x) if a>D)
  W=T-left-right
  single_left=sum(li-a for a in x if a<li);single_right=sum(a-D for a in x if a>D)
  mu=Q(len(x)*li-single_left-single_right,D);var=mu+Q(2*W,D)-mu*mu
  if expected is not None:
   rr=next(r for r in expected if Q(r['L'])==L);assert var==Q(rr['continuous']['variance']['exact']);assert mu==Q(rr['continuous']['mean_count']['exact'])
  H=sum((Q(li*mom[k]-mom[k+1],S**(k+1)) for k in range(6)),Q(0))
  q=Q(li,S);rem=Q(T,S)*q**6/(1-q)
  Hb=qa(H).union(qa(H+rem));trans=qa(L-L*L)+2*Hb;delta=mu-L
  mean_correction=(1-2*L)*delta-delta*delta
  edge_correction=2*(qa(Q(W,D))-Hb)
  assert (trans+qa(mean_correction)+edge_correction).contains(qa(var))
  near=bisect_right(gaps,li+margin);dlo=Q(S,scale)-L-3*eps;assert dlo>0
  trans_error=2*near*(2*eps/dlo+L*eps/dlo**2)
  rows.append({'L':float(L),'mean_exact':str(mu),'variance_exact':str(var),'pair_count':mom[0],
   'unclipped_pair_weight_exact':str(Q(T,scale)),'clipped_pair_weight_exact':str(Q(W,scale)),
   'left_pair_loss_exact':str(Q(left,scale)),'right_pair_loss_exact':str(Q(right,scale)),
   'translation_variance_ball':trans.str(45,radius=True,more=True),'mean_correction_exact':str(mean_correction),
   'edge_correction_ball':edge_correction.str(45,radius=True,more=True),'geometric_tail_bound':str(2*rem),
   'translation_coordinate_error_bound':str(trans_error),'identity_exact':True})
 return rows

def selftest():
 import random
 rng=random.Random(3409)
 for trial in range(50):
  points=sorted(rng.sample(range(1,300),20));S=300;lengths=[1,4,20,100]
  rr=analyze(points,lengths,span=S)
  for r in rr:
   L=Q(r['L']);events=sorted({Q(0),Q(S)-L}|{Q(v) for x in points for v in [x,x-L] if 0<v<S-L})
   m1=Q(0);m2=Q(0)
   for a,b in zip(events,events[1:]):
    s=(a+b)/2;n=sum(s<=x<s+L for x in points);m1+=(b-a)*n;m2+=(b-a)*n*n
   assert Q(r['mean_exact'])==m1/(S-L);assert Q(r['variance_exact'])==m2/(S-L)-(m1/(S-L))**2
   direct=sum((Q(L-(b-a),S-(b-a)) for i,a in enumerate(points) for b in points[i+1:] if b-a<L),Q(0))
   truth=L-L*L+2*direct;assert arb(r['translation_variance_ball']).contains(qa(truth))
 return {'status':'PASS','independent_exact_cases':200}

def main():
 start=time.time();root=json.loads((BASE/'ordinate_error_results.json').read_text());model=json.loads((BASE/'continuous_height_results.json').read_text());hs=[]
 for n0 in [100000,1000000,10000000,20000000,50000000]:
  e=json.loads((BASE/'data'/f'empirical_stable_{n0}.json').read_text());u=[float.fromhex(x) for x in e['unfolded_float_hex']]
  epss=next(h['unfolded_coordinate_error_upper'] for h in root['heights'] if h['n0']==n0);eps=Q(math.nextafter(float(arb(epss).upper()),math.inf))
  rows=analyze(u,[.5,1,2,4,8,16],expected=e['rows'],eps=eps)
  for r in rows:
   m=arb(next(x['infinite_product_height_average_ball'] for x in model['results'] if x['n0']==n0 and x['L']==r['L']))
   t=arb(r['translation_variance_ball']);r['translation_minus_local_BK_ball']=(t-m).str(45,radius=True,more=True)
   r['empirical_minus_translation_ball']=(qa(Q(r['variance_exact']))-t).str(45,radius=True,more=True)
  hs.append({'n0':n0,'rows':rows});print('pair closure',n0,flush=True)
 out={'status':'PASS','selftest':selftest(),'heights':hs,'elapsed_seconds':time.time()-start,
 'scope':'exact algebra for existing point-to-point windows; stationary fixed-window unbiasedness does not automatically apply to these point-index-selected windows'}
 (HERE/'pair_closure_results.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
