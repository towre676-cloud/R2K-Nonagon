"""Exact tests of global grounding-control submodularity and its edge boundary."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json, random
from frame_formalism import graph, inverse, add, diag, matmul, transpose

def outer(x): return [[a*b for b in x] for a in x]
def response(K,actions,S):
 A=[r[:] for r in K]
 for i in S:A=add(A,actions[i])
 return inverse(A)[0][0]
def check(K,actions):
 r=len(actions);V={}
 for mask in range(1<<r):V[mask]=response(K,actions,[j for j in range(r) if mask>>j&1])
 checks=0;worst=None
 for mask in V:
  for i in range(r):
   for j in range(i+1,r):
    if mask&(1<<i|1<<j):continue
    # Reduction f=h0-h is submodular iff h(A)+h(Aij)>=h(Ai)+h(Aj).
    gap=V[mask]+V[mask|1<<i|1<<j]-V[mask|1<<i]-V[mask|1<<j]
    checks+=1
    if worst is None or gap<worst[0]:worst=(gap,mask,i,j)
 return checks,worst,V

def main():
 rng=random.Random(10262026);total=0;min_gap=None;ground_networks=0
 for n in range(3,7):
  for trial in range(3):
   K=graph(n,rng);actions=[diag([F(1+j) if i==j else F(0) for i in range(n)]) for j in range(1,n)]
   count,worst,V=check(K,actions);assert worst[0]>=0
   total+=count;ground_networks+=1;min_gap=worst[0] if min_gap is None else min(min_gap,worst[0])
 edge_example=None
 for n in range(3,6):
  if edge_example:break
  for trial in range(100):
   K=graph(n,rng);pairs=list(combinations(range(n),2));actions=[]
   for u,v in pairs:
    b=[F(i==u)-F(i==v) for i in range(n)];actions.append(outer(b))
   count,worst,V=check(K,actions)
   if worst[0]<0:
    gap,mask,i,j=worst
    edge_example={'K':[[str(x) for x in row] for row in K],'edges':pairs,'A':[pairs[k] for k in range(len(pairs)) if mask>>k&1],'e':pairs[i],'f':pairs[j], 'submodularity_gap_for_reduction':str(gap),'gap_float':float(gap),'h_A':str(V[mask]),'h_Ae':str(V[mask|1<<i]),'h_Af':str(V[mask|1<<j]),'h_Aef':str(V[mask|1<<i|1<<j])}
    break
 assert edge_example
 out={'arithmetic':'exact rational','grounding_networks':ground_networks,'grounding_submodularity_inequalities':total,'minimum_exact_gap':str(min_gap),'edge_counterexample':edge_example}
 Path(__file__).with_name('ground_submodularity_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
