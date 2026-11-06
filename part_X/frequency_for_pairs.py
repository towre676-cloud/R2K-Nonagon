#!/usr/bin/env python3
"""Independent frequency-for-pairs review and rank-gated exact recovery.
No imports from Claude's source tree. 40 exact tests; modular witnesses n=7..10;
raw-query Jacobian scout uses all information, including the otherwise unused Gjj.
Run: python review/frequency_for_pairs.py. Results written next to this program.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools,json,math,random
import numpy as np

def solve(A,b):
 n=len(A);B=[[F(v) for v in A[i]]+[F(b[i])] for i in range(n)]
 for k in range(n):
  p=next((i for i in range(k,n) if B[i][k]),None)
  if p is None:raise ValueError('rank gate: matrix is not full rank')
  B[k],B[p]=B[p],B[k];x=B[k][k];B[k]=[v/x for v in B[k]]
  for i in range(n):
   if i!=k and B[i][k]:
    x=B[i][k];B[i]=[a-x*b for a,b in zip(B[i],B[k])]
 return [r[-1] for r in B]
def inv(A):
 n=len(A);cols=[solve(A,[int(i==j) for i in range(n)]) for j in range(n)]
 return [[cols[j][i] for j in range(n)] for i in range(n)]
def rank(A):
 B=[[F(v) for v in r] for r in A];rr=0
 for k in range(len(B[0])):
  p=next((i for i in range(rr,len(B)) if B[i][k]),None)
  if p is None:continue
  B[rr],B[p]=B[p],B[rr];x=B[rr][k]
  for j in range(k,len(B[0])):B[rr][j]/=x
  for i in range(rr+1,len(B)):
   if B[i][k]:
    x=B[i][k]
    for j in range(k,len(B[0])):B[i][j]-=x*B[rr][j]
  rr+=1
 return rr

def network(n,family,rng):
 A=[[F(0)]*n for _ in range(n)]
 edges=[(j-1,j) for j in range(1,n)] if family=='path' else [(0,j) for j in range(1,n)] if family=='star' else list(itertools.combinations(range(n),2))
 for i,j in edges:
  v=F(rng.randint(1,9));A[i][i]+=v;A[j][j]+=v;A[i][j]-=v;A[j][i]-=v
 for i in range(n):A[i][i]+=rng.randint(1,7)
 return A

def pencil(K,C,M,s):return [[K[i][j]+s*C[i][j]+s*s*M[i][j] for j in range(len(K))] for i in range(len(K))]
def queries(n):return [{}]+[{j:F(t)} for j in range(1,n) for t in [1,3]]
def port(A,D):
 B=[r[:] for r in A]
 for j,t in D.items():B[j][j]+=t
 return solve(B,[1]+[0]*(len(A)-1))[0]
def column_from_data(n,y):
 g=[y[0]];diag=[y[0]]
 for j in range(1,n):
  d1=y[0]-y[2*j-1];d3=y[0]-y[2*j]
  q=(F(1)-F(1,3))/(1/d1-1/d3)
  a=math.isqrt(q.numerator);b=math.isqrt(q.denominator)
  if a*a!=q.numerator or b*b!=q.denominator:raise ValueError('invalid exact square or non-rational fixture')
  g.append(F(a,b));diag.append(q/d1-1)
 return g,diag

def linear_rows(n,gs,svals):
 pairs=list(itertools.combinations_with_replacement(range(n),2));idx={a:i for i,a in enumerate(pairs)};N=len(pairs);d=2*N+n;rows=[];rhs=[]
 for g,s in zip(gs,svals):
  for i in range(n):
   row=[F(0)]*d
   for j in range(n):
    a=idx[tuple(sorted((i,j)))];row[a]+=g[j];row[N+a]+=s*g[j]
   row[2*N+i]+=s*s*g[i];rows.append(row);rhs.append(F(i==0))
 return rows,rhs

def model_vector(K,C,M):
 n=len(K);pairs=list(itertools.combinations_with_replacement(range(n),2))
 return [A[i][j] for A in [K,C] for i,j in pairs]+[M[i][i] for i in range(n)]

def _exact(x):
 """External exact values may be Fraction, integer, or a rational string."""
 if isinstance(x,(float,np.floating)):raise ValueError('input gate: binary floats are not exact observations; supply Fraction or rational strings')
 try:return F(x)
 except (ValueError,TypeError,ZeroDivisionError) as exc:raise ValueError('input gate: invalid rational observation') from exc

def _full_rank_solve(rows,rhs,d):
 """Exact rectangular rank/consistency gate followed by back substitution."""
 B=[list(r)+[b] for r,b in zip(rows,rhs)];rr=0;piv=[]
 for k in range(d):
  p=next((i for i in range(rr,len(B)) if B[i][k]),None)
  if p is None:continue
  B[rr],B[p]=B[p],B[rr];z=B[rr][k]
  for j in range(k,d+1):B[rr][j]/=z
  for i in range(rr+1,len(B)):
   if B[i][k]:
    z=B[i][k]
    for j in range(k,d+1):B[i][j]-=z*B[rr][j]
  piv.append(k);rr+=1
 if any(all(v==0 for v in row[:d]) and row[d]!=0 for row in B):raise ValueError('regression consistency gate failed')
 if rr!=d:raise ValueError(f'rank gate failed: rank {rr} < {d} unknowns')
 x=[F(0)]*d
 for i in range(d-1,-1,-1):
  k=piv[i];x[k]=B[i][d]-sum(B[i][j]*x[j] for j in range(k+1,d))
 return x,rr

def _validate_physical(K,C,M):
 n=len(K)
 for name,A in [('K',K),('C',C)]:
  if any(A[i][j]!=A[j][i] for i in range(n) for j in range(n)):raise ValueError(f'physical gate: {name} is not symmetric')
  if any(A[i][j]>0 for i in range(n) for j in range(n) if i!=j):raise ValueError(f'physical gate: {name} has a positive off-diagonal')
  if any(sum(row)<0 for row in A):raise ValueError(f'physical gate: {name} has a negative ground')
 if not any(sum(row)>0 for row in K):raise ValueError('physical gate: K has no positive ground')
 seen={0}
 while True:
  nxt=seen|{j for i in seen for j in range(n) if K[i][j]<0}
  if nxt==seen:break
  seen=nxt
 if len(seen)!=n:raise ValueError('physical gate: spring graph is disconnected')
 if any(M[i][i]<=0 for i in range(n)) or any(M[i][j]!=0 for i in range(n) for j in range(n) if i!=j):raise ValueError('physical gate: masses must be positive and diagonal')
 return {'symmetric_K_C':True,'nonnegative_edge_weights':True,'nonnegative_grounds':True,'positive_K_ground':True,'connected_spring_graph':True,'positive_diagonal_mass':True}

def recover_from_observations(n,svals,records):
 """Recover a rational physical pencil using ONLY externally supplied raw readings.

 records[k] = {'baseline': h0, 'single_1': [h_j(1) for j=1..n-1],
               'single_3': [h_j(3) for j=1..n-1]}.
 svals must start at zero, strictly increase, and contain at least n+2 points.
 Inputs: Fraction/int/rational strings; extracted positive square roots must be
 rational (the exact rational-network backend). Binary floats are rejected. On failure
 this function raises ValueError naming the failed gate. No ground truth, seed,
 topology, hidden model, or saved expected answer enters this inverse.

 Return {'K': ..., 'C': ..., 'M': ..., 'certificate': ...}, all matrices exact.
 The final forward gate replays EVERY raw response, including the information
 in Gjj that is not used by the first-column linear regression.
 """
 if not isinstance(n,int) or n<1:raise ValueError('input gate: n must be a positive integer')
 sv=[_exact(s) for s in svals]
 if len(sv)<n+2 or len(records)!=len(sv):raise ValueError('input gate: require >= n+2 points and one raw record per point')
 if sv[0]!=0 or any(sv[i]>=sv[i+1] for i in range(len(sv)-1)):raise ValueError('input gate: Laplace points must start at zero and strictly increase')
 gs=[];flat_records=[]
 for rec in records:
  try:
   if len(rec['single_1'])!=n-1 or len(rec['single_3'])!=n-1:raise ValueError('input gate: wrong single-load record length')
   y=[_exact(rec['baseline'])]
   for j in range(n-1):y.extend([_exact(rec['single_1'][j]),_exact(rec['single_3'][j])])
  except (KeyError,TypeError) as exc:raise ValueError('input gate: malformed raw record') from exc
  if y[0]<=0:raise ValueError('positive-column gate: baseline must be positive')
  for j in range(1,n):
   if not (F(0)<y[2*j]<y[2*j-1]<y[0]):raise ValueError('positive-column gate: positive loads must strictly reduce compliance')
  try:g,diag=column_from_data(n,y)
  except (ValueError,ZeroDivisionError) as exc:raise ValueError('positive-column gate: invalid exact two-strength data') from exc
  if any(v<=0 for v in g+diag):raise ValueError('positive-column gate: recovered Green entries are not positive')
  gs.append(g);flat_records.append(y)
 rows,rhs=linear_rows(n,gs,sv);d=n*(n+2);x,rnk=_full_rank_solve(rows,rhs,d)
 pairs=list(itertools.combinations_with_replacement(range(n),2));N=len(pairs)
 K=[[F(0)]*n for _ in range(n)];C=[[F(0)]*n for _ in range(n)];M=[[F(0)]*n for _ in range(n)]
 for off,A in [(0,K),(N,C)]:
  for k,(i,j) in enumerate(pairs):A[i][j]=A[j][i]=x[off+k]
 for j in range(n):M[j][j]=x[2*N+j]
 physical=_validate_physical(K,C,M);replayed=0
 for k,s in enumerate(sv):
  A=pencil(K,C,M,s)
  for q,observed in zip(queries(n),flat_records[k]):
   if port(A,q)!=observed:raise ValueError(f'raw replay gate failed: frequency index {k}, load {q}')
   replayed+=1
 cert={'rank':rnk,'unknowns':d,'linear_equations':len(rows),'frequency_points':len(sv),'distinct_grounding_configurations':2*n-1,'real_scalar_observations':len(sv)*(2*n-1),'positive_first_column_gate':True,'physical_network_gate':physical,'all_raw_replay_gate':True,'forward_replayed_raw_readings':replayed,'replays_unused_Gjj_information':True,'maximum_exact_raw_residual':'0'}
 return {'K':K,'C':C,'M':M,'certificate':cert}

def _raw_records(K,C,M,svals):
 """Measurement simulator used by tests, deliberately outside the inverse API."""
 n=len(K);records=[]
 for s in svals:
  A=pencil(K,C,M,s);records.append({'baseline':port(A,{}),'single_1':[port(A,{j:F(1)}) for j in range(1,n)],'single_3':[port(A,{j:F(3)}) for j in range(1,n)]})
 return records

def raw_jacobian(K,C,M,svals):
 k,c,m=[np.array(x,dtype=float) for x in [K,C,M]];n=len(k);pairs=list(itertools.combinations_with_replacement(range(n),2));rows=[]
 for s in svals:
  s=float(s)
  P=k+s*c+s*s*m
  for D in queries(n):
   Pload=P+np.diag([float(D.get(j,0)) for j in range(n)])
   r=np.linalg.solve(Pload,np.eye(n)[:,0]);v=np.array([-(1 if i==j else 2)*r[i]*r[j] for i,j in pairs])
   rows.append(np.r_[v,s*v,-s*s*r*r])
 A=np.array(rows);sigma=np.linalg.svd(A,compute_uv=False)
 return {'rows':len(rows),'columns':A.shape[1],'numerical_rank':int(np.linalg.matrix_rank(A)),'minimum_singular_value':float(sigma[-1]),'condition':float(sigma[0]/sigma[-1])}

def qmod(x,p):
 x=F(x);return x.numerator*pow(x.denominator,-1,p)%p
def mod_solve(A,b,p):
 n=len(A);B=[[qmod(x,p) for x in A[i]]+[qmod(b[i],p)] for i in range(n)]
 for k in range(n):
  v=next((i for i in range(k,n) if B[i][k]),None)
  if v is None:raise ValueError('prime divides denominator')
  B[k],B[v]=B[v],B[k];u=pow(B[k][k],-1,p);B[k]=[(a*u)%p for a in B[k]]
  for i in range(n):
   if i!=k and B[i][k]:
    u=B[i][k];B[i]=[(a-u*b)%p for a,b in zip(B[i],B[k])]
 return [r[-1] for r in B]
def mod_rank(A,p):
 B=[[qmod(x,p) for x in r] for r in A];rr=0
 for k in range(len(B[0])):
  v=next((i for i in range(rr,len(B)) if B[i][k]),None)
  if v is None:continue
  B[rr],B[v]=B[v],B[rr];u=pow(B[rr][k],-1,p)
  for j in range(k,len(B[0])):B[rr][j]=B[rr][j]*u%p
  for i in range(rr+1,len(B)):
   if B[i][k]:
    u=B[i][k]
    for j in range(k,len(B[0])):B[i][j]=(B[i][j]-u*B[rr][j])%p
  rr+=1
 return rr

def raw_mod_rank(K,C,M,svals,p):
 n=len(K);pairs=list(itertools.combinations_with_replacement(range(n),2));rows=[]
 for s in svals:
  P=pencil(K,C,M,s)
  for D in queries(n):
   B=[r[:] for r in P]
   for j,t in D.items():B[j][j]+=t
   r=mod_solve(B,[1]+[0]*(n-1),p);sm=qmod(s,p)
   v=[(-(1 if i==j else 2)*r[i]*r[j])%p for i,j in pairs]
   rows.append(v+[(sm*a)%p for a in v]+[(-sm*sm*r[j]*r[j])%p for j in range(n)])
 return mod_rank(rows,p)

def main():
 rng=random.Random(10312026);out={'seed':10312026,'strengths':['1','3'],'exact':[],'modular':[],'raw_query_jacobian':[]}
 for trial in range(40):
  n=2+trial%5;family=['path','star','dense'][trial%3];K=network(n,family,rng);C=network(n,family,rng);M=[[F(rng.randint(1,7)) if i==j else F(0) for j in range(n)] for i in range(n)]
  svals=[F(k,2) for k in range(n+2)];gs=[]
  records=_raw_records(K,C,M,svals)
  api=recover_from_observations(n,svals,records)
  assert api['K']==K and api['C']==C and api['M']==M
  if trial==0:first_fixture=(K,C,M,svals,records)

  for s in svals:
   A=pencil(K,C,M,s);g,di=column_from_data(n,[port(A,q) for q in queries(n)])
   assert g==[r[0] for r in inv(A)] and di==[inv(A)[i][i] for i in range(n)]
   gs.append(g)
  rows,rhs=linear_rows(n,gs,svals)
  try:
   rec=solve(rows,rhs);ok=rec==model_vector(K,C,M);rnk=len(rows)
  except ValueError:rnk=rank(rows);ok=False
  out['exact'].append({'trial':trial,'n':n,'family':family,'rank':rnk,'unknowns':n*(n+2),'exact_recovery':ok,'frequencies':n+2,'configurations':2*n-1,'real_scalars':(2*n-1)*(n+2),'observation_only_certificate':api['certificate']})
  L=math.ceil(n*(n+2)/(2*n-1));sc=[F(k,2) for k in range(L)];J=raw_jacobian(K,C,M,sc)
  J['modular_rank']=raw_mod_rank(K,C,M,sc,2147483647);J['exact_full_rank_witness']=J['modular_rank']==n*(n+2)
  out['raw_query_jacobian'].append({'trial':trial,'n':n,'family':family,'frequencies':L,'distinct_configurations':2*n-1,**J})
  if trial%10==0:print('exact trial',trial,flush=True)
 # Modular rank is an exact nonvanishing witness over Q: all denominator pivots nonzero modulo p.
 prime=2147483647
 for n in range(7,11):
  for family in ['path','star','dense']:
   K=network(n,family,rng);C=network(n,family,rng);M=[[F(rng.randint(1,7)) if i==j else F(0) for j in range(n)] for i in range(n)];svals=[F(k,2) for k in range(n+2)];gs=[]
   for s in svals:gs.append(mod_solve(pencil(K,C,M,s),[1]+[0]*(n-1),prime))
   rows,rhs=linear_rows(n,gs,svals);r=mod_rank(rows,prime)
   out['modular'].append({'n':n,'family':family,'prime':prime,'rank':r,'unknowns':n*(n+2),'full_rank_witness':r==n*(n+2),'K':K,'C':C,'M':M,'svals':svals,'short_raw_frequencies':math.ceil(n*(n+2)/(2*n-1)),'short_raw_modular_rank':raw_mod_rank(K,C,M,[F(k,2) for k in range(math.ceil(n*(n+2)/(2*n-1)))],prime)})
 # Physical non-proportional damping example where FIRST-COLUMN linear recovery is deficient for all s.
 K=[[F(4),F(-3)],[F(-3),F(6)]];C=[[F(2),F(-1)],[F(-1),F(5)]];M=[[F(1),F(0)],[F(0),F(1)]]
 DK=[[F(1),F(-2)],[F(-2),F(4)]];DC=[[F(0),F(-1)],[F(-1),F(4)]];DM=[[F(0),F(0)],[F(0),F(1)]];tau=F(1,10)
 Kp=[[K[i][j]+tau*DK[i][j] for j in range(2)] for i in range(2)];Cp=[[C[i][j]+tau*DC[i][j] for j in range(2)] for i in range(2)];Mp=[[M[i][j]+tau*DM[i][j] for j in range(2)] for i in range(2)];sv=list(map(F,range(8)))
 gs=[[r[0] for r in inv(pencil(K,C,M,s))] for s in sv];gsp=[[r[0] for r in inv(pencil(Kp,Cp,Mp,s))] for s in sv];rows,rhs=linear_rows(2,gs,sv)
 out['nongeneric_example']={'K':K,'C':C,'M':M,'delta_K':DK,'delta_C':DC,'delta_M':DM,'tau':tau,'first_columns_equal':gs==gsp,'rank_with_8_frequencies':rank(rows),'unknowns':8,'unit_load_h_original_at_s1':port(pencil(K,C,M,F(1)),{1:F(1)}),'unit_load_h_perturbed_at_s1':port(pencil(Kp,Cp,Mp,F(1)),{1:F(1)}),'identity':'DeltaP(s)=u(s)u(s)^T, u(s)=[1,-(s+2)], and u(s)^T g(s)=0'}
 # Negative controls for the externally usable API.
 import copy
 controls={}
 def rejection(name,nn,ss,rr):
  try:recover_from_observations(nn,ss,rr)
  except ValueError as exc:controls[name]={'rejected':True,'reason':str(exc)}
  else:raise AssertionError(name+' corrupt fixture was accepted')
 kk,cc,mm,ss,rr=first_fixture;nn=len(kk)
 bad=copy.deepcopy(rr);bad[0]['single_1'][0]+=F(1,1000000)
 rejection('one_raw_value_changed',nn,ss,bad)
 # Two altered values retain EXACTLY the correct extracted first column.
 # This must reach and fail the final all-raw replay, not merely extraction.
 bad=copy.deepcopy(rr);at=1;A=pencil(kk,cc,mm,ss[at]);gg=[r[0] for r in inv(A)];v=[gg[1],-gg[0]]+[F(0)]*(nn-2)
 B=[[A[i][j]+F(1,100)*v[i]*v[j] for j in range(nn)] for i in range(nn)]
 assert [r[0] for r in inv(B)]==gg
 bad[at]['single_1'][0]=port(B,{1:F(1)});bad[at]['single_3'][0]=port(B,{1:F(3)})
 rejection('same_first_column_but_corrupt_unused_Gjj',nn,ss,bad)
 assert controls['same_first_column_but_corrupt_unused_Gjj']['reason'].startswith('raw replay gate failed')
 rejection('rank_deficient_physical_first_column_family',2,sv,_raw_records(K,C,M,sv))
 assert controls['rank_deficient_physical_first_column_family']['reason'].startswith('rank gate failed')
 out['api_negative_controls']=controls
 out['summary']={'exact_networks':40,'exact_all_pass':all(z['exact_recovery'] for z in out['exact']),'observation_only_API_all_pass':all(z['observation_only_certificate']['all_raw_replay_gate'] for z in out['exact']),'negative_controls_all_rejected':all(z['rejected'] for z in controls.values()),'modular_witnesses':len(out['modular']),'modular_all_fullrank':all(z['full_rank_witness'] for z in out['modular']),'raw_query_full_rank_count':sum(z['numerical_rank']==z['columns'] for z in out['raw_query_jacobian']),'raw_query_total':len(out['raw_query_jacobian']),'raw_query_exact_full_rank_witnesses':sum(z['exact_full_rank_witness'] for z in out['raw_query_jacobian']),'raw_query_minimum_point_counts':{n:math.ceil(n*(n+2)/(2*n-1)) for n in range(2,11)}}
 Path(__file__).with_name('frequency_for_pairs_results.json').write_text(json.dumps(out,indent=2,default=str));print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
