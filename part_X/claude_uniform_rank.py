"""Exact uniform arrowhead witness for frequency-for-pairs, all n.

The symbolic proof is in claude_uniform_rank.md. The rational witness below
checks its square acquisition matrix at the actual half-integer sample grid.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from frame_formalism import inverse

def unknown_index(n):
 idx={}
 for kind in ['K','C']:
  for i in range(n):
   for j in range(i,n):idx[kind,i,j]=len(idx)
 for i in range(n):idx['M',i,i]=len(idx)
 return idx

def rational_rank(rows):
 A=[r[:] for r in rows];r=0
 for c in range(len(A[0])):
  p=next((i for i in range(r,len(A)) if A[i][c]),None)
  if p is None:continue
  A[r],A[p]=A[p],A[r];v=A[r][c]
  A[r]=[x/v for x in A[r]]
  for i in range(r+1,len(A)):
   if A[i][c]:
    v=A[i][c]
    A[i]=[a-v*b for a,b in zip(A[i],A[r])]
  r+=1
 return r

def witness(n):
 svals=[F(k,2) for k in range(n+2)];idx=unknown_index(n)
 K=[[F(0)]*n for _ in range(n)];C=[[F(0)]*n for _ in range(n)];M=[[F(0)]*n for _ in range(n)]
 for j in range(1,n):
  sj=svals[j];K[j][j]=-sj;C[j][j]=1-sj;M[j][j]=1;K[0][j]=K[j][0]=-1
 rows=[];special=[]
 for k,s in enumerate(svals):
  P=[[K[i][j]+s*C[i][j]+s*s*M[i][j] for j in range(n)] for i in range(n)]
  G=inverse(P);g=[G[i][0] for i in range(n)]
  if 1<=k<n:
   assert g==[F(-int(i==k)) for i in range(n)];special.append(k)
  else:
   assert g[0]
   for j in range(1,n):assert g[j]/g[0]==1/((s-svals[j])*(s+1))
  for i in range(n):
   row=[F(0)]*len(idx)
   for j in range(n):
    a,b=min(i,j),max(i,j)
    row[idx['K',a,b]]+=g[j];row[idx['C',a,b]]+=s*g[j]
   row[idx['M',i,i]]+=s*s*g[i]
   rows.append(row)
 r=rational_rank(rows);assert r==n*(n+2)
 return {'n':n,'rank':r,'unknowns':n*(n+2),'samples':[str(s) for s in svals],'coordinate_column_samples':special,'invertible_at_every_sample':True}

def main():
 t=time.perf_counter();rows=[witness(n) for n in range(2,13)]
 out={'method':'exact rational rank; uniform symbolic proof in companion memo','witness_scope':'nonphysical algebraic arrowhead witness; generic physical result follows by polynomial nonvanishing on a full-dimensional open family','results':rows,'all_full_rank':True,'elapsed_seconds':time.perf_counter()-t}
 p=Path(__file__).with_name('claude_uniform_rank_results.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
