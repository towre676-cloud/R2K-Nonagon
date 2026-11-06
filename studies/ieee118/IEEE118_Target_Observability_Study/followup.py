#!/usr/bin/env python3
"""Explicit post-inspection follow-up: noise-oriented 180 rows plus known branch separator."""
import run as R
import numpy as np,json,time
from scipy.linalg import cholesky,solve_triangular
from scipy.optimize import least_squares
def main():
    start=time.perf_counter();buses,K,B,w,keys,nom,t,h0,pats=R.setup()
    saved=json.loads((R.OUT/'panels.json').read_text());panels=saved['panels'];pairs=saved['targets_reduced_indices']
    vt,_=R.target(w,B,pairs);yt,_=R.evaluate(w,B,pats,jac=False)
    sep=panels['historical']['180'][-1];out=[]
    for method in ('qr_norm','target_a','target_correlated'):
      sel=list(panels[method]['180'])
      if sep not in sel:sel.append(sep)
      z=np.zeros((len(sel),len(K)+1));z[:,0]=1
      for k,i in enumerate(sel):z[k,1+(pats[i][1][0] if pats[i][1] else 0)]=1
      L=cholesky(.5*np.eye(len(sel))+.25*z@z.T,lower=True)
      for sigma in (0,1e-7,1e-5):
       for trial in range(1 if sigma==0 else 5):
        rng=np.random.default_rng(700000+trial);noise=np.sqrt(.5)*rng.normal(size=len(sel))+.5*z@rng.normal(size=z.shape[1]);obs=yt[sel]+sigma*h0*noise
        cache={}
        def both(x):
         if 'x' not in cache or not np.array_equal(x,cache['x']):
          y,J=R.evaluate(np.exp(x),B,pats,sel);cache.update(x=x.copy(),f=solve_triangular(L,(y-obs)/h0,lower=True),j=solve_triangular(L,J/h0,lower=True))
         return cache['f'],cache['j']
        sol=least_squares(lambda x:both(x)[0],np.log(nom),jac=lambda x:both(x)[1],max_nfev=100,xtol=1e-11,ftol=1e-11,gtol=1e-11,bounds=(np.log(nom)-8,np.log(nom)+8))
        wh=np.exp(sol.x);v,_=R.target(wh,B,pairs)
        out.append(dict(method=method,budget=len(sel),added_index=sep,sigma=sigma,trial=trial,target_relative_rms=float(np.sqrt(np.mean(((v-vt)/vt)**2))),weight_relative_l2=float(np.linalg.norm(wh-w)/np.linalg.norm(w)),solver_success=bool(sol.success),nfev=int(sol.nfev),residual_norm=float(np.linalg.norm(sol.fun))))
       R.save('followup.json',{'scope':'Exploratory follow-up specified after primary results inspection; uses previously discovered branch-separating index, no global guarantee','rows':out,'seconds':time.perf_counter()-start})
       print(method,sigma,flush=True)
if __name__=='__main__':main()
