"""Finite Feynman-Kac resolvent and occupation-time derivative checks."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad

def main():
 rng=np.random.default_rng(20161030);rows=[]
 for n in range(2,8):
  for k in range(3):
   W=np.zeros((n,n))
   for j in range(1,n):
    i=int(rng.integers(j));W[i,j]=W[j,i]=rng.uniform(.3,1.4)
   grounds=rng.uniform(.3,1.,n)
   K=np.diag(W.sum(1)+grounds)-W;D=np.diag(rng.uniform(.1,.7,n));A=K+D
   G=np.linalg.inv(A)
   # Generator -K kills at grounds. Feynman-Kac adds kill D.
   # Finite semigroup integral is tested against the resolvent; tail bounded
   # by exp(-lambda_min H)/lambda_min through the spectral norm.
   lam=np.linalg.eigvalsh(A)[0];T=40/lam
   integral,qe=quad(lambda t:expm(-A*t)[0,0],0,T,epsabs=2e-14,epsrel=2e-14,limit=120)
   tail=np.exp(-lam*T)/lam
   rel=abs(integral-G[0,0])/G[0,0]
   # Ordered products determine all mixed derivatives (occupation moments).
   # 3rd derivative checked independently by the 8 signed nearby resolvents.
   ids=list(range(min(n,3)));i,j,l=(ids*3)[:3];eps=2e-3
   from itertools import product, permutations
   derivative=0.
   for signs in product([-1,1],repeat=3):
    dd=np.zeros(n)
    for idx,sg in zip((i,j,l),signs):dd[idx]+=sg*eps
    derivative+=np.prod(signs)*np.linalg.inv(A+np.diag(dd))[0,0]/(8*eps**3)
   exact=-sum(G[0,p[0]]*G[p[0],p[1]]*G[p[1],p[2]]*G[p[2],0] for p in permutations((i,j,l)))
   # Mixed partial>=0 => grounding reduction submodular, checked in separate exact program.
   rows.append(dict(n=n,replicate=k,resolvent=G[0,0],semigroup_integral=integral,absolute_error=abs(integral-G[0,0]),relative_error=rel,analytic_tail_bound=tail,quadrature_error_estimate=qe,third_jet_exact=exact,third_jet_difference=derivative,third_jet_relative_error=abs(derivative-exact)/abs(exact)))
 assert max(x['relative_error'] for x in rows)<1e-12
 assert max(x['third_jet_relative_error'] for x in rows)<1e-3
 out=dict(schema='nonagon.part_x.occupation.v1',seed=20161030,cases=rows,summary=dict(cases=len(rows),max_resolvent_relative_error=max(x['relative_error'] for x in rows),max_third_jet_relative_error=max(x['third_jet_relative_error'] for x in rows)),verdict='PASS')
 Path(__file__).with_name('occupation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
