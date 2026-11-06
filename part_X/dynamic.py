#!/usr/bin/env python3
"""Part X: exact and binary64 damped-pencil acquisition experiments.
No saved data or ground-truth parameter is passed to either inverse.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools, json, math, platform, random
import numpy as np
import scipy
from scipy.optimize import least_squares
import mpmath as mp
OUT=Path(__file__).with_name('dynamic_results.json')

def exact_inv(A):
 n=len(A); B=[[F(x) for x in A[i]]+[F(i==j) for j in range(n)] for i in range(n)]
 for k in range(n):
  p=next(i for i in range(k,n) if B[i][k]); B[k],B[p]=B[p],B[k]
  z=B[k][k];B[k]=[x/z for x in B[k]]
  for i in range(n):
   if i!=k:
    z=B[i][k];B[i]=[x-z*y for x,y in zip(B[i],B[k])]
 return [r[n:] for r in B]
def sq(q):
 a=math.isqrt(q.numerator); b=math.isqrt(q.denominator)
 assert a*a==q.numerator and b*b==q.denominator
 return F(a,b)
def read(A,load):
 B=[r[:] for r in A]
 for j,t in load.items():B[j][j]+=t
 return exact_inv(B)[0][0]
def queries(n):
 return [{}]+[{j:t} for j in range(1,n) for t in (F(1),F(2))]+[{i:F(1),j:F(1)} for i,j in itertools.combinations(range(1,n),2)]
def exact_reconstruct(n,y):
 G=[[F(0)]*n for _ in range(n)];G[0][0]=y[0]; k=1
 for j in range(1,n):
  d1=y[0]-y[k];d2=y[0]-y[k+1];k+=2
  q=1/(2*(1/d1-1/d2));gjj=q/d1-1
  G[0][j]=G[j][0]=sq(q);G[j][j]=gjj
 for i,j in itertools.combinations(range(1,n),2):
  d=y[0]-y[k];k+=1;a=1+G[i][i];b=1+G[j][j]
  z=sq((d*a-G[0][i]**2)*(d*b-G[0][j]**2))
  G[i][j]=G[j][i]=(G[0][i]*G[0][j]-z)/d
 return exact_inv(G)
def net(n,family,rng):
 E={}
 for j in range(1,n):E[(j-1,j) if family=='path' else (0,j)]=F(rng.randint(3,8),3)
 if family=='dense':
  for i,j in itertools.combinations(range(n),2):E[(i,j)]=F(rng.randint(3,8),3)
 A=[[F(0)]*n for _ in range(n)]
 for (i,j),w in E.items():A[i][i]+=w;A[j][j]+=w;A[i][j]-=w;A[j][i]-=w
 for j in range(n):A[j][j]+=F(rng.randint(1,3),10*n)
 return A

def pencil(K,C,M,s):return [[K[i][j]+s*C[i][j]+s*s*M[i][j] for j in range(len(K))] for i in range(len(K))]
def connected_order(B):
 n=len(B);S=[0]
 while len(S)<n:
  S.append(next(j for j in range(n) if j not in S and any(B[j][i]!=0 for i in S)))
 return S

def clamp_data(A,order):
 return [exact_inv([[A[i][j] for j in order[:k]] for i in order[:k]])[0][0] for k in range(1,len(A)+1)]
def clamp_diagonals(offdiag,order,y):
 n=len(offdiag);R=[r[:] for r in offdiag];R[0][0]=1/y[0]
 for k in range(1,n):
  S=order[:k];j=order[k];Ai=exact_inv([[R[a][b] for b in S] for a in S]); b=[R[a][j] for a in S]
  v=[sum(Ai[a][q]*b[q] for q in range(k)) for a in range(k)]
  z=sum(b[a]*v[a] for a in range(k))
  assert v[0]!=0
  R[j][j]=z+v[0]**2/(y[k]-y[k-1])
 return R

def arr(A):return np.array(A,dtype=float)
def h(A,D=None):
 if D is None:D=np.zeros(len(A))
 return np.linalg.solve(A+np.diag(D),np.eye(len(A))[:,0])[0]
def float_protocol(A):
 return [h(A,np.array([float(e.get(j,0)) for j in range(len(A))])) for e in queries(len(A))]
def float_reconstruct(n,y,complex_=False):
 dtype=complex if complex_ else float
 G=np.zeros((n,n),dtype=dtype);G[0,0]=y[0];k=1
 for j in range(1,n):
  d1=y[0]-y[k];d2=y[0]-y[k+1];k+=2
  q=1/(2*(1/d1-1/d2));G[j,j]=q/d1-1
  G[0,j]=G[j,0]=np.sqrt(q)
 for i,j in itertools.combinations(range(1,n),2):
  d=y[0]-y[k];k+=1;a=1+G[i,i];b=1+G[j,j]
  z=np.sqrt((d*a-G[0,i]**2)*(d*b-G[0,j]**2))
  G[i,j]=G[j,i]=(G[0,i]*G[0,j]-z)/d
 return np.linalg.inv(G)
def diag_data(A,t):
 n=len(A);return np.array([h(A)]+[h(A,np.eye(n)[j]*t) for j in range(1,n)])
def diag_J(A,t):
 n=len(A);return np.array([-np.linalg.solve(A+(np.diag(np.eye(n)[l]*t) if l else 0),np.eye(n)[:,0])**2 for l in range(n)])
def fdet(A):
 B=[r[:] for r in A];v=F(1);n=len(B)
 for k in range(n):
  p=next((i for i in range(k,n) if B[i][k]),None)
  if p is None:return F(0)
  if p!=k:B[k],B[p]=B[p],B[k];v=-v
  v*=B[k][k]
  for i in range(k+1,n):
   a=B[i][k]/B[k][k]
   for j in range(k+1,n):B[i][j]-=a*B[k][j]
 return v

def main():
 rng=random.Random(80312026)
 R={'program':'dynamic.py','seed':80312026,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'arithmetic':['exact fractions','binary64','mpmath 80 decimal digits'],'counts':[],'exact_clamp':[],'float':[],'finite_local':[],'imaginary':[],'imaginary_finite_local':[]}
 for n in range(2,9):
  N=n*(n+1)//2
  R['counts'].append({'n':n,'old_real_readings':3*N,'dimension':n*(n+2),'two_real_plus_clamps':2*N+n,'one_complex_plus_clamps_real_scalars':2*N+n,'one_complex_plus_clamps_sampled_settings':N+n,'one_complex_plus_finite_static_distinct_grounding_configurations':N,'one_complex_plus_finite_static_sampled_settings':N+n,'growing_clamps_distinct_grounding_configurations':N+n-1})
  for family in ['path','star','dense']:
   K=net(n,family,rng);C=net(n,family,rng);M=[[F(i==j)*F(rng.randint(1,5),10) for j in range(n)] for i in range(n)]
   P1=pencil(K,C,M,F(1));P2=pencil(K,C,M,F(2));Q=queries(n)
   Kr=exact_reconstruct(n,[read(K,q) for q in Q]);P1r=exact_reconstruct(n,[read(P1,q) for q in Q])
   off2=[[F(0) if i==j else 2*P1r[i][j]-Kr[i][j] for j in range(n)] for i in range(n)]
   order=connected_order(off2);y=clamp_data(P2,order);P2r=clamp_diagonals(off2,order,y)
   Mr=[[(P2r[i][j]-2*P1r[i][j]+Kr[i][j])/2 for j in range(n)] for i in range(n)];Cr=[[P1r[i][j]-Kr[i][j]-Mr[i][j] for j in range(n)] for i in range(n)]
   assert Kr==K and Cr==C and Mr==M
   R['exact_clamp'].append({'n':n,'family':family,'readings':2*N+n,'exact_K':Kr==K,'exact_C':Cr==C,'exact_M':Mr==M,'min_clamp_gain':str(min(y[k]-y[k-1] for k in range(1,n)))})
   k,c,m=map(arr,(K,C,M));p1=k+c+m;p2=k+2*c+4*m
   kf=float_reconstruct(n,float_protocol(k));p1f=float_reconstruct(n,float_protocol(p1))
   off=2*p1f-kf;np.fill_diagonal(off,0)
   rec=off.copy();rec[0,0]=1/h(p2[np.ix_(order[:1],order[:1])]);yf=[]
   for z in range(1,n+1):yf.append(h(p2[np.ix_(order[:z],order[:z])]))
   for z in range(1,n):
    S=order[:z];j=order[z];a=rec[np.ix_(S,S)];b=rec[S,j];v=np.linalg.solve(a,b)
    rec[j,j]=b@v+v[0]**2/(yf[z]-yf[z-1])
   mf=(rec-2*p1f+kf)/2;cf=p1f-kf-mf
   er=lambda x,y:float(np.linalg.norm(x-y)/np.linalg.norm(y))
   R['float'].append({'n':n,'family':family,'K_relative_error':er(kf,k),'C_relative_error':er(cf,c),'M_relative_error':er(mf,m)})
   t=F(1);J=[]
   for l in range(n):
    B=[r[:] for r in P2]
    if l:B[l][l]+=t
    gi=exact_inv(B)
    J.append([-gi[j][0]**2 for j in range(n)])
   det=fdet(J);assert det!=0
   G=exact_inv(P2);g=[G[j][0] for j in range(n)]
   coeff=-F(2)**(n-1)*g[0]*math.prod(x*x for x in g[1:])*fdet(G)
   L=[[-x*x for x in g]]+[[2*g[j]*g[l]*G[j][l] for j in range(n)] for l in range(1,n)]
   assert fdet(L)==coeff
   jf=diag_J(p2,1.);sv=np.linalg.svd(jf,compute_uv=False)
   target=diag_data(p2,1.);offp=p2-np.diag(np.diag(p2));rowsum=-offp.sum(axis=1)
   sols=[];nruns=6
   # Positive ground parameterization: d=weighted degree+exp(x), all iterates physical.
   for zz in range(nruns):
    xx=np.log(np.diag(p2)-rowsum)+np.random.default_rng(100*n+zz).normal(0,.6,n)
    fit=least_squares(lambda x:(diag_data(offp+np.diag(rowsum+np.exp(x)),1.)-target)/target,x0=xx,xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=2000)
    sols.append({'relative_diagonal_error':float(np.linalg.norm(rowsum+np.exp(fit.x)-np.diag(p2))/np.linalg.norm(np.diag(p2))),'max_relative_data_residual':float(np.max(np.abs(fit.fun))),'success':bool(fit.success)})
   R['finite_local'].append({'n':n,'family':family,'t':1,'exact_J_determinant_nonzero':det!=0,'exact_J_determinant':str(det) if n<=3 else None,'leading_coefficient_identity':fdet(L)==coeff,'Jacobian_condition':float(sv[0]/sv[-1]),'multistart':sols})
   # Composite finite-load, imaginary-axis local inverse: all scalar quadratures counted.
   w=.1; Ap=k+1j*w*c-w*w*m; Ar=float_reconstruct(n,float_protocol(Ap),True)
   offk=Ar.real.copy();np.fill_diagonal(offk,0);deg=-offk.sum(axis=1)
   stat=diag_data(k,1.)
   fit=least_squares(lambda x:(diag_data(offk+np.diag(deg+np.exp(x)),1.)-stat)/stat,np.full(n,np.log(.1)),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=2500)
   kr=offk+np.diag(deg+np.exp(fit.x));cr=Ar.imag/w;mr=np.diag((np.diag(kr)-np.diag(Ar.real))/(w*w))
   pairs=list(itertools.combinations_with_replacement(range(n),2));jac=[]
   for load in queries(n):
    DD=np.diag([float(load.get(j,0)) for j in range(n)]);rr=np.linalg.solve(Ap+DD,np.eye(n)[:,0]);base=np.array([-rr[i]*rr[j]*(1 if i==j else 2) for i,j in pairs]);row=np.r_[base,1j*w*base,w*w*rr**2]
    jac.extend([row.real,row.imag])
   for l in range(n):
    DD=np.diag(np.eye(n)[l] if l else np.zeros(n));rr=np.linalg.solve(k+DD,np.eye(n)[:,0]);base=np.array([-rr[i]*rr[j]*(1 if i==j else 2) for i,j in pairs]);jac.append(np.r_[base,np.zeros(len(pairs)+n)])
   sj=np.linalg.svd(np.array(jac),compute_uv=False)
   R['imaginary_finite_local'].append({'n':n,'family':family,'omega':w,'t_static':1.,'real_scalar_readings':len(jac),'distinct_grounding_configurations':N,'configuration_frequency_samples':N+n,'parameter_dimension':n*(n+2),'Jacobian_min_singular_value':float(sj[-1]),'Jacobian_condition':float(sj[0]/sj[-1]),'Jacobian_full_numerical_rank':bool(np.linalg.matrix_rank(np.array(jac))==n*(n+2)),'K_relative_error':er(kr,k),'C_relative_error':er(cr,c),'M_relative_error':er(mr,m),'static_fit_max_relative_residual':float(np.max(np.abs(fit.fun)))})
   for w in [.001,.01,.03,.1,.3,1.,3.]:
    A=k+1j*w*c-w*w*m
    try:
     Ar=float_reconstruct(n,float_protocol(A),True)
     R['imaginary'].append({'n':n,'family':family,'omega':w,'relative_P_error':er(Ar,A),'physical_principal_branch_recovered':er(Ar,A)<1e-7})
    except Exception as ex:R['imaginary'].append({'n':n,'family':family,'omega':w,'error':str(ex)})
 # 80-digit imaginary phase branch + clamped static diagonal inverse, independent of Fraction inverse.
 mp.mp.dps=80
 hp=[]
 for n in [2,3,5,8]:
  K=net(n,'dense',rng);C=net(n,'dense',rng);M=[[F(i==j)*F(rng.randint(1,5),10) for j in range(n)] for i in range(n)]
  mm=lambda A:mp.matrix([[mp.mpf(x.numerator)/x.denominator for x in row] for row in A])
  k,c,m=map(mm,(K,C,M));w=mp.mpf('0.01');A=k+mp.j*w*c-w*w*m
  def mread(B,q):
   B=B.copy()
   for j,t in q.items():B[j,j]+=mp.mpf(t.numerator)/t.denominator
   return B**-1
  y=[mread(A,q)[0,0] for q in queries(n)];G=mp.matrix(n);G[0,0]=y[0];z=1
  for j in range(1,n):
   d1=y[0]-y[z];d2=y[0]-y[z+1];z+=2;q=1/(2*(1/d1-1/d2));G[j,j]=q/d1-1;G[0,j]=G[j,0]=mp.sqrt(q)
  for i,j in itertools.combinations(range(1,n),2):
   d=y[0]-y[z];z+=1;a=1+G[i,i];b=1+G[j,j]
   G[i,j]=G[j,i]=(G[0,i]*G[0,j]-mp.sqrt((d*a-G[0,i]**2)*(d*b-G[0,j]**2)))/d
  Ar=G**-1;kr=mp.matrix([[mp.re(Ar[i,j]) if i!=j else 0 for j in range(n)] for i in range(n)])
  # Dense graph order root,1,...,n−1.
  yy=[(k[:z,:z]**-1)[0,0] for z in range(1,n+1)];kr[0,0]=1/yy[0]
  for z in range(1,n):
   aa=kr[:z,:z];b=kr[:z,z:z+1];v=aa**-1*b;kr[z,z]=(b.T*v)[0]+v[0]**2/(yy[z]-yy[z-1])
  cr=mp.matrix([[mp.im(Ar[i,j])/w for j in range(n)] for i in range(n)]);mr=mp.matrix(n)
  for j in range(n):mr[j,j]=(kr[j,j]-mp.re(Ar[j,j]))/(w*w)
  err=lambda X,Y:mp.norm(X-Y)/mp.norm(Y)
  hp.append({'n':n,'omega':str(w),'P_relative_error':mp.nstr(err(Ar,A),12),'K_relative_error':mp.nstr(err(kr,k),12),'C_relative_error':mp.nstr(err(cr,c),12),'M_relative_error':mp.nstr(err(mr,m),12),'real_scalar_readings':n*(n+2)})
 R['imaginary_80dps']=hp
 # Fully printed three-node worked example, retained as exact fractions.
 K=[[F(3),F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]]
 C=[[F(2),F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]]
 M=[[F(i==j)*(i+1) for j in range(3)] for i in range(3)]
 P1=pencil(K,C,M,F(1));P2=pencil(K,C,M,F(2));order=[0,1,2]
 J=[]
 for l in range(3):
  B=[r[:] for r in P2]
  if l:B[l][l]+=1
  V=exact_inv(B);J.append([-V[j][0]**2 for j in range(3)])
 R['worked_example']={'K':K,'C':C,'M':M,'P1':P1,'P2':P2,'finite_queries':queries(3),'h_at_s0':[read(K,q) for q in queries(3)],'h_at_s1':[read(P1,q) for q in queries(3)],'growing_clamp_h_at_s2':clamp_data(P2,order),'finite_diagonal_h_at_s2':diag_data(arr(P2),1.).tolist(),'finite_diagonal_J_determinant':fdet(J),'count':15,'old_count':18}
 p2=arr(P2)
 R['finite_condition_scan']=[{'t':float(t),'condition_J':float(np.linalg.cond(diag_J(p2,t))),'det_J':float(np.linalg.det(diag_J(p2,t)))} for t in np.logspace(-8,8,33)]

 R['summary']={'imaginary_finite_all_full_rank':all(x['Jacobian_full_numerical_rank'] for x in R['imaginary_finite_local']),'imaginary_finite_max_K_error':max(x['K_relative_error'] for x in R['imaginary_finite_local']),'imaginary_finite_max_C_error':max(x['C_relative_error'] for x in R['imaginary_finite_local']),'imaginary_finite_max_M_error':max(x['M_relative_error'] for x in R['imaginary_finite_local']),'exact_clamp_networks':len(R['exact_clamp']),'exact_clamp_all_pass':all(x['exact_K'] and x['exact_C'] and x['exact_M'] for x in R['exact_clamp']),'finite_J_all_exact_nonsingular':all(x['exact_J_determinant_nonzero'] for x in R['finite_local']),'finite_leading_coefficient_all_exact':all(x['leading_coefficient_identity'] for x in R['finite_local']),'max_binary64_K_relative_error':max(x['K_relative_error'] for x in R['float']),'max_binary64_C_relative_error':max(x['C_relative_error'] for x in R['float']),'max_binary64_M_relative_error':max(x['M_relative_error'] for x in R['float']),'imaginary_principal_branch_success_by_omega':{str(w):sum(x.get('physical_principal_branch_recovered',False) for x in R['imaginary'] if x['omega']==w) for w in [.001,.01,.03,.1,.3,1.,3.]},'imaginary_total_each_omega':21,'finite_multistart_count':sum(len(x['multistart']) for x in R['finite_local']),'finite_multistart_worst_parameter_error':max(s['relative_diagonal_error'] for x in R['finite_local'] for s in x['multistart']),'finite_multistart_worst_residual':max(s['max_relative_data_residual'] for x in R['finite_local'] for s in x['multistart'])}
 OUT.write_text(json.dumps(R,indent=2,default=str));print(json.dumps(R['summary'],indent=2))
if __name__=='__main__':main()
