"""Externally specified unfolded spans, not windows ending at selected points.
The common rule [n0+1000,n0+19000] was chosen after the original campaign;
these are exploratory geometry checks, not new preregistered held-out tests.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,math
from flint import arb,ctx
from pair_closure import analyze,grid
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;ctx.prec=128

def txt(a):return a.str(45,radius=True,more=True)
def event_moments(x,S,li):
 D=S-li;ev={0:0,D:0}
 for a in x:
  if 0<a<D:ev[a]=ev.get(a,0)-1
  if 0<a-li<D:ev[a-li]=ev.get(a-li,0)+1
 n=sum(0<a<=li for a in x);last=0;s1=0;s2=0
 for t,delta in sorted(ev.items()):
  s1+=(t-last)*n;s2+=(t-last)*n*n;n+=delta;last=t
 return Q(s1,D),Q(s2,D)-Q(s1,D)**2

def main():
 twopi=2*arb.pi();out=[]
 def unfold(t):return t/twopi*(t/twopi).log()-t/twopi+arb(7)/8
 for n0 in [100000,1000000,10000000,20000000,50000000]:
  z=json.loads((BASE/'data'/f'zeros_{n0}.json').read_text());a=arb(n0+1000);b=arb(n0+19000)
  first=unfold(arb(z['ordinate_balls'][0]));last=unfold(arb(z['ordinate_balls'][-1]));assert first<a and last>b
  pts=[];errors=[];source_indices=[]
  for j,t in enumerate(z['ordinate_balls']):
   y=unfold(arb(t))
   if y<a or y>b:continue
   assert y>a and y<b,'ambiguous fixed-window membership'
   x=y-a;f=float(x.mid());pts.append(f);source_indices.append(n0+j)
   errors.append((x-arb(f)).abs_upper())
  eps=max(errors);ep=Q(math.nextafter(float(eps),math.inf));Ls=[.5,1,2,4,8,16]
  rows=analyze(pts,Ls,span=18000,eps=ep);xi,scale=grid(pts,Ls);S=18000*scale
  for row in rows:
   L=Q(row['L']);li=int(L*scale);mu,v=event_moments(xi,S,li)
   assert mu==Q(row['mean_exact']) and v==Q(row['variance_exact'])
   near_margin=math.ceil(2*ep*scale);end=0;M=0
   for i,x in enumerate(xi):
    while end<len(xi) and xi[end]-x<=li+near_margin:end+=1
    M=max(M,end-i)
   threshold=2*sum(errors,arb(0))/(18000-arb(float(L)))
   row['maximum_count_bound']=M;row['empirical_coordinate_error_bound']=txt((4*M*threshold).upper())
  obj={'n0':n0,'absolute_unfolded_lower':n0+1000,'absolute_unfolded_upper':n0+19000,'span':18000,'point_count':len(pts),
   'points_hex':[p.hex() for p in pts],'source_zero_indices':source_indices,'max_coordinate_error':txt(eps),
   'lower_buffer':txt((a-first).lower()),'upper_buffer':txt((last-b).lower()),'rows':rows,'event_and_pair_moments_exactly_agree':True}
  (HERE/f'fixed_data_{n0}.json').write_text(json.dumps(obj,indent=2));out.append({'n0':n0,'point_count':len(pts)});print('fixed span',n0,len(pts),flush=True)
 (HERE/'fixed_preparation.json').write_text(json.dumps({'status':'PASS','windows':out,'all_exact_moment_crosschecks':30,'exploratory_not_preregistered':True},indent=2))
if __name__=='__main__':main()
