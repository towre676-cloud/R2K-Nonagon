"""Propagate certified ordinate balls to continuous-anchor variance.
For counts N,N' on normalized anchor y in [0,1], sum the movement of
all 2n threshold events. This bounds int|N-N'|; if both counts <= M,
|Var N - Var N'| <= 4M int|N-N'|. Exact stored-dyadic anchors are the reference.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,sys,math
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ctx.prec=128

def qa(x):
 q=Q(x);return arb(q.numerator)/q.denominator

def max_window(u,length):
 u=[Q(x) for x in u];length=Q(length);j=0;m=0
 for i,a in enumerate(u):
  while j<len(u) and u[j]-a<=length:j+=1
  m=max(m,j-i)
 return m

def one(n0):
 zp=HERE/'data'/f'zeros_{n0}.json';ep=HERE/'data'/f'empirical_{n0}.json'
 stable=HERE/'data'/f'empirical_stable_{n0}.json'
 if stable.exists():ep=stable
 if not(zp.exists() and ep.exists()):return None
 z=json.loads(zp.read_text());e=json.loads(ep.read_text());u=[float.fromhex(x) for x in e['unfolded_float_hex']]
 twopi=2*arb.pi()
 def unfold(t):return t/twopi*(t/twopi).log()-t/twopi+arb(7)/8
 vals=[unfold(arb(x)) for x in z['ordinate_balls']];base=vals[0];v=[x-base for x in vals];v[0]=arb(0)
 eps=max((x-qa(y)).abs_upper() for x,y in zip(v,u));rows=[]
 for r in e['rows']:
  L=Q(r['L']);D=Q(u[-1])-L;Db=v[-1]-qa(L);assert Db>0
  S=arb(0)
  for x,y in zip(v,u):
   for off in [Q(0),L]:S+=((x-qa(off))/Db-qa((Q(y)-off)/D)).abs_upper()
  Le=math.nextafter(float(L)+2*float(eps.upper()),math.inf)
  M=max_window(u,Le);bound=(4*M*S).upper()
  rows.append({'L':float(L),'maximum_count_bound':M,'threshold_movement_sum_upper':S.upper().str(40,radius=True,more=True),
    'continuous_variance_error_bound':bound.str(40,radius=True,more=True),'mean_count_error_bound':S.upper().str(40,radius=True,more=True)})
 return {'n0':n0,'unfolded_coordinate_error_upper':eps.str(40,radius=True,more=True),'rows':rows,
 'scope':'exact smooth unfolding of the certified ordinate balls versus the stored binary64 unfolding; no claim of an exactly stationary finite-height process'}

def main():
 result={'status':'complete','heights':[x for n in [100000,1000000,10000000,20000000,50000000] if (x:=one(n)) is not None]}
 (HERE/'ordinate_error_results.json').write_text(json.dumps(result,indent=2));print([(x['n0'],x['unfolded_coordinate_error_upper']) for x in result['heights']])
if __name__=='__main__':main()
