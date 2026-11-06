"""Exact intervention-frame compatibility tests; Python standard library only."""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import json, random

def transpose(A): return list(map(list,zip(*A)))
def matmul(A,B): return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def add(A,B): return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def sub(A,B): return [[a-b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def inverse(A):
 n=len(A); T=[list(map(F,row))+e for row,e in zip(A,eye(n))]
 for j in range(n):
  p=next(i for i in range(j,n) if T[i][j]);T[j],T[p]=T[p],T[j]
  d=T[j][j];T[j]=[v/d for v in T[j]]
  for i in range(n):
   if i!=j:
    d=T[i][j];T[i]=[a-d*b for a,b in zip(T[i],T[j])]
 return [r[n:] for r in T]
def det(A):
 if not A:return F(1)
 T=[list(map(F,row)) for row in A];p=F(1);n=len(T)
 for j in range(n):
  q=next((i for i in range(j,n) if T[i][j]),None)
  if q is None:return F(0)
  if q!=j:T[j],T[q]=T[q],T[j];p=-p
  d=T[j][j];p*=d
  for i in range(j+1,n):
   d0=T[i][j]/d
   for k in range(j+1,n):T[i][k]-=d0*T[j][k]
 return p
def minor(A,I,J=None):return [[A[i][j] for j in (I if J is None else J)] for i in I]
def diag(d):return [[x if i==j else F(0) for j in range(len(d))] for i,x in enumerate(d)]
def graph(n,rng):
 W=[[F(0) for _ in range(n)] for _ in range(n)]
 for i in range(n):
  for j in range(i):
   w=F(rng.randrange(1,8)) if i==j+1 or rng.random()<.55 else F(0)
   W[i][j]=W[j][i]=w
 return [[sum(W[i])+F(rng.randrange(1,5)) if i==j else -W[i][j] for j in range(n)] for i in range(n)]

def main():
 rng=random.Random(20261030); counts={k:0 for k in ['two_jet_exact','minor_polynomial_exact','frame_update_exact','elimination_exact','gauge_exact','third_jet_exact']}
 example=None
 for n in range(2,8):
  for sample in range(5):
   K=graph(n,rng);G=inverse(K);h=G[0][0];a=[G[i][0] for i in range(1,n)];S=minor(G,list(range(1,n)));r=n-1
   # The first and second jets recover the frame, choosing the positive M-matrix orientation.
   first=[-v*v for v in a];second=[[2*a[i]*S[i][j]*a[j] for j in range(r)] for i in range(r)]
   assert all(-first[i]==a[i]**2 for i in range(r))
   Grec=[[h]+a]+[[a[i]]+[second[i][j]/(2*a[i]*a[j]) for j in range(r)] for i in range(r)]
   assert Grec==G;counts['two_jet_exact']+=1
   th=[F(rng.randrange(1,6),rng.randrange(1,5)) for i in range(r)]
   H=inverse(add(K,diag([F(0)]+th)))[0][0]
   numerator=F(0);denominator=F(0)
   for k in range(r+1):
    for J in combinations(range(r),k):
     mon=F(1)
     for j in J:mon*=th[j]
     denominator+=det(minor(S,list(J)))*mon
     gram=[[h]+[a[j] for j in J]]+[[a[i]]+[S[i][j] for j in J] for i in J]
     numerator+=det(gram)*mon
   assert H==numerator/denominator;counts['minor_polynomial_exact']+=1
   R=inverse(add(S,diag([1/t for t in th])))
   a_col=[[v] for v in a]
   H2=h-matmul(matmul(transpose(a_col),R),a_col)[0][0]
   S2=sub(S,matmul(matmul(S,R),S));a2=sub(a_col,matmul(matmul(S,R),a_col))
   G2=inverse(add(K,diag([F(0)]+th)))
   assert H2==H and S2==minor(G2,list(range(1,n))) and a2==[[G2[i][0]] for i in range(1,n)]
   counts['frame_update_exact']+=1
   # Preserve the port and an intervention node while eliminating all other coordinates.
   if n>2:
    A=[0,1];J=list(range(2,n));B=minor(K,A,J)
    Kron=sub(minor(K,A),matmul(matmul(B,inverse(minor(K,J))),transpose(B)))
    assert inverse(Kron)==minor(G,A)
    assert inverse(add(Kron,diag([F(0),th[0]])))[0][0]==inverse(add(K,diag([F(0),th[0]]+[F(0)]*(n-2))))[0][0]
    counts['elimination_exact']+=1
   # Signature congruences fixing the port preserve all diagonal-intervention observations.
   signs=[F(1)]+[F(-1 if j%2 else 1) for j in range(1,n)]
   Kg=[[signs[i]*K[i][j]*signs[j] for j in range(n)] for i in range(n)]
   assert inverse(add(Kg,diag([F(0)]+th)))[0][0]==H
   counts['gauge_exact']+=1
   # The third Taylor coefficient for one combined direction matches the ordered-walk identity.
   D=diag([F(0)]+th);GD=matmul(G,D)
   coeff3=-matmul(matmul(matmul(GD,GD),GD),G)[0][0]
   walk=F(0)
   for i in range(r):
    for j in range(r):
     for k in range(r):walk-=a[i]*S[i][j]*S[j][k]*a[k]*th[i]*th[j]*th[k]
   assert coeff3==walk;counts['third_jet_exact']+=1
   if n==3 and sample==0:
    example={'K':[[str(x) for x in row] for row in K],'G':[[str(x) for x in row] for row in G], 'strengths':[str(x) for x in th], 'response':str(H),'numerator':str(numerator),'denominator':str(denominator)}
 out={'method':'exact rational arithmetic; no floating tolerances','seed':20261030,'network_sizes':list(range(2,8)),'networks_per_size':5,'checks':counts,'total_assertion_groups':sum(counts.values()),'example':example}
 path=Path(__file__).with_name('frame_formalism_results.json');path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
