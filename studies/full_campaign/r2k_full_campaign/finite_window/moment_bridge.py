"""Exact factorial-moment expansion and higher-order obstruction.
Classical moment partition calculus specialized to the anchor-variance functional.
No stochastic model for the zeta zeros is asserted.
"""
from fractions import Fraction as Q
from itertools import product,permutations,combinations
from pathlib import Path
import json,random,math
HERE=Path(__file__).resolve().parent

def expectation(law,fn):return sum((p*fn(bits) for bits,p in law.items()),Q(0))
def factorial_moment(law,indices):return expectation(law,lambda b:Q(all(b[i] for i in indices)))
def moment_formula(law,f,g,F=None,G=None):
 F=f if F is None else F;G=g if G is None else G;n=len(f)
 m1=sum((f[i]*factorial_moment(law,[i]) for i in range(n)),Q(0))+sum((g[i][j]*factorial_moment(law,[i,j]) for i,j in permutations(range(n),2)),Q(0))
 m2=sum((F[i]*factorial_moment(law,[i]) for i in range(n)),Q(0))+sum((G[i][j]*factorial_moment(law,[i,j]) for i,j in permutations(range(n),2)),Q(0))
 terms=[Q(0)]*4
 for i in range(n):terms[0]+=f[i]*F[i]*factorial_moment(law,[i])
 for i,j in permutations(range(n),2):terms[1]+=(f[i]*F[j]+2*f[i]*G[i][j]+2*F[i]*g[i][j]+2*g[i][j]*G[i][j])*factorial_moment(law,[i,j])
 for i,j,k in permutations(range(n),3):terms[2]+=(f[i]*G[j][k]+F[i]*g[j][k]+4*g[i][j]*G[i][k])*factorial_moment(law,[i,j,k])
 for i,j,k,l in permutations(range(n),4):terms[3]+=g[i][j]*G[k][l]*factorial_moment(law,[i,j,k,l])
 return m1,m2,terms,sum(terms,Q(0))-m1*m2

def functional(bits,f,g):return sum(f[i]*bits[i] for i in range(len(f)))+sum(g[i][j]*bits[i]*bits[j] for i,j in permutations(range(len(f)),2))
def law_from_weights(weights):
 z=sum(weights.values());return {k:Q(v,z) for k,v in weights.items()}
def run():
 bits=list(product([0,1],repeat=4));A={b:Q(3**sum(b),256) for b in bits};B={b:Q(3**sum(b)+(-1)**sum(b),256) for b in bits}
 assert sum(A.values())==sum(B.values())==1 and all(p>0 for p in B.values())
 for k in [1,2,3]:
  for subset in combinations(range(4),k):
   for pattern in product([0,1],repeat=k):
    assert expectation(A,lambda b:Q(tuple(b[i] for i in subset)==pattern))==expectation(B,lambda b:Q(tuple(b[i] for i in subset)==pattern))
 f=[Q(5,36)]*4;g=[[Q(0) if i==j else -Q(1,36) for j in range(4)] for i in range(4)]
 result=[]
 for label,law in [('independent_occupancy',A),('fourth_order_perturbation',B)]:
  mean=expectation(law,lambda b:functional(b,f,g));var=expectation(law,lambda b:(functional(b,f,g)-mean)**2)
  a,bb,terms,v=moment_formula(law,f,g);assert mean==a==bb and var==v
  result.append({'law':label,'mean_anchor_variance':str(mean),'sampling_variance_of_anchor_variance':str(var),'factorial_order_second_moment_terms':[str(t) for t in terms]})
 assert result[0]['mean_anchor_variance']=='11/48'
 assert result[0]['sampling_variance_of_anchor_variance']=='11/13824' and result[1]['sampling_variance_of_anchor_variance']=='1/1152'
 rng=random.Random(7249)
 for trial in range(20):
  n=5;bs=list(product([0,1],repeat=n));law=law_from_weights({b:rng.randint(1,20) for b in bs});f=[Q(rng.randint(-5,5),7) for _ in range(n)];F=[Q(rng.randint(-5,5),11) for _ in range(n)]
  g=[[Q(0) for _ in range(n)] for _ in range(n)];G=[[Q(0) for _ in range(n)] for _ in range(n)]
  for i in range(n):
   for j in range(i+1,n):g[i][j]=g[j][i]=Q(rng.randint(-3,3),13);G[i][j]=G[j][i]=Q(rng.randint(-3,3),17)
  m,M,t,cov=moment_formula(law,f,g,F,G)
  direct=expectation(law,lambda b:functional(b,f,g)*functional(b,F,G))-expectation(law,lambda b:functional(b,f,g))*expectation(law,lambda b:functional(b,F,G))
  assert cov==direct
 out={'status':'PASS','counterexample':{'circle_circumference':3,'nominal_intensity':1,'interval_length':'1/2','construction':'four sites spaced 3/4, independent jitter in [-1/16,1/16], common uniform rotation; occupancy laws differ only in fourth-order interaction',
 'identical_factorial_moment_measures_through_order':3,'results':result,'variance_difference':'1/13824'},'random_correlated_laws_covariance_checks':20,
 'scope':'general impossibility of deducing estimator sampling variance from only first/second/third-order correlations; not a claim that zeta zeros follow either example'}
 (HERE/'moment_bridge_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':run()
