"""Independent controls for the two new numerical proof estimates."""
from fractions import Fraction as Q
from pathlib import Path
import sys,json,random
from flint import arb,ctx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from continuous_anchor_stats import statistics
ctx.prec=128

def check_gauss():
 c=arb(1)/2;h=arb(1)/2;R=arb(2);q=h/R;M=1/(5+c-R);n=8
 s=arb(0)
 for k in range(n):
  x,w=arb.legendre_p_root(n,k,weight=True);s+=h*w/(5+c+h*x)
 truth=(arb(6)/5).log();bound=4*h*M*q**(2*n)/(1-q)
 assert (s-truth).abs_upper()<bound
 return {'status':'PASS','function':'1/(5+x), [0,1], exact integral log(6/5)',
 'gauss_ball':s.str(40,radius=True),'exact_integral_ball':truth.str(40,radius=True),'cauchy_remainder_bound':bound.str(40,radius=True)}

def check_anchor_perturbation():
 rng=random.Random(691)
 worst=Q(0)
 for _ in range(100):
  u=[Q(0)]+[Q(x,7) for x in sorted(rng.sample(range(5,100),12))]
  v=[u[0]]+[x+Q(rng.randint(-2,2),10000) for x in u[1:]];L=Q(rng.randint(2,15),7)
  D=u[-1]-L;E=v[-1]-L
  S=sum((abs((b-off)/E-(a-off)/D) for a,b in zip(u,v) for off in [Q(0),L]),Q(0))
  actual=abs(statistics(u,L)['variance']-statistics(v,L)['variance']);bound=4*len(u)*S
  assert actual<=bound
  if bound:worst=max(worst,actual/bound)
 return {'status':'PASS','independent_exact_perturbed_point_sets':100,'largest_error_to_bound_ratio':float(worst)}

out={'gauss_remainder_control':check_gauss(),'empirical_perturbation_control':check_anchor_perturbation(),'status':'PASS'}
Path(__file__).with_name('proof_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
