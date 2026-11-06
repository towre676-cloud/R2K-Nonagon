#!/usr/bin/env python3
"""Reproducible IEEE-118 known-topology design comparison. Run with BLAS threads=1."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json,time,hashlib,itertools,csv
from pathlib import Path
import numpy as np
from scipy.linalg import qr,solve_triangular,cholesky
from scipy.optimize import least_squares
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'))
import target_relative_panel as old
from fixtures import CASES
P=json.loads((ROOT/'protocol.json').read_text())
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def finite(x):return float(x) if np.isfinite(x) else None

def setup():
    # Independently parse the retained MATPOWER file, including tap ratios.
    import re
    text=(ROOT/'vendor/case118.m').read_text()
    block=re.search(r'mpc\.branch\s*=\s*\[(.*?)\];',text,re.S).group(1)
    edges=[]
    for line in block.splitlines():
        line=line.split('%')[0].strip().rstrip(';')
        if not line:continue
        a=line.split();tap=float(a[8]) or 1.
        if float(a[10]): edges.append((int(a[0]),int(a[1]),1/(float(a[3])*tap)))
    expected=CASES['case118']['edges']
    assert len(edges)==len(expected)
    assert all(a==d and b==e and abs(c-f)<1e-11 for (a,b,c),(d,e,f) in zip(edges,expected))
    assert hashlib.sha256((ROOT/'vendor/case118.m').read_bytes()).hexdigest()=='bc2e6f22b4b9e776572885ee4b50e4f4ab2ee0c5577e9126e86d906f14c4b5f7'
    buses,K,B,w,keys=old.edge_parameterization(CASES['case118'])
    nominal,scale,t,h0=old.nominal_design(B,K)
    pats=old.patterns(len(buses),t)
    return buses,K,B,w,keys,nominal,t,h0,pats

def evaluate(w,B,pats,indices=None,jac=True):
    # Vectorized rank-one / rank-two Woodbury evaluation, independent of old loop.
    if indices is None: indices=np.arange(len(pats))
    G=np.linalg.inv((B*w)@B.T);g=G[:,0]; ids=np.asarray(indices,int)
    X=np.repeat(g[:,None],len(ids),axis=1)
    for typ in ('single','pair'):
        cols=np.array([k for k,i in enumerate(ids) if pats[i][0]==typ],int)
        if not len(cols):continue
        i=np.array([pats[ids[k]][1][0] for k in cols]); a=np.array([pats[ids[k]][2][0] for k in cols])
        aa=G[i,i]+1/a
        if typ=='single': X[:,cols]-=G[:,i]*(g[i]/aa)[None,:]
        else:
            j=np.array([pats[ids[k]][1][1] for k in cols]);bb=G[j,j]+1/a;cc=G[i,j]
            det=aa*bb-cc*cc;u=(bb*g[i]-cc*g[j])/det;v=(aa*g[j]-cc*g[i])/det
            X[:,cols]-=G[:,i]*u[None,:]+G[:,j]*v[None,:]
    H=-(B.T@X).T**2*w if jac else None
    return X[0],H

def target(w,B,pairs):
    G=np.linalg.inv((B*w)@B.T); U=np.zeros((len(G),len(pairs)))
    for k,(i,j) in enumerate(pairs):U[i,k]=1;U[j,k]=-1
    X=G@U;v=np.sum(U*X,axis=0);J=-(B.T@X).T**2*w
    return v,J

def covariance(H,Z=None):
    # Stable square-root least-squares, never form H.T@H for inversion.
    p=H.shape[1]
    if Z is None:
        Q,R=qr(H,mode='economic');D=solve_triangular(R,Q.T)
        return D@D.T,D
    # R_noise=.5 I + .25 common + .25 group. QR includes nuisance prior.
    z=Z.shape[1]
    A=np.block([[H/np.sqrt(.5),Z/np.sqrt(.5)],[np.zeros((z,p)),2*np.eye(z)]])
    Q,R=qr(A,mode='economic');Ri=solve_triangular(R,np.eye(len(R)))
    C=Ri@Ri.T
    D=Ri[:p]@Q[:len(H),:].T/np.sqrt(.5)
    return C,D

def target_only_row(H,T,Z,method,budget,noise):
    """Target estimability does not require all parameter columns to be identifiable."""
    if noise=='correlated':
        L=cholesky(.5*np.eye(len(H))+.25*Z@Z.T,lower=True)
        A=solve_triangular(L,H,lower=True)
    else:A=H;L=None
    U,s,V=np.linalg.svd(A,full_matrices=True);rank=int(np.sum(s>s[0]*1e-12))
    gap=float(np.linalg.norm(T@V[rank:].T));ok=gap<1e-8
    out=dict(method=method,budget=budget,noise=noise,rank_rtol1e12=rank,condition=finite(s[0]/s[-1]),status='target_only_estimable_at_declared_tolerance' if ok else 'target_unresolved_at_declared_tolerance',target_kernel_gap=gap,target_rms_per_sigma=None,target_mc_rms_per_sigma=None,naive_iid_decoder_mc_rms_per_sigma=None,parameter_rms_per_sigma=None)
    if not ok:return out
    D=(T@V[:rank].T/s[:rank])@U[:,:rank].T
    if L is not None:D=solve_triangular(L.T,D.T,lower=False).T
    Ui,si,Vi=np.linalg.svd(H,full_matrices=False);ri=int(np.sum(si>si[0]*1e-12));Di=(T@Vi[:ri].T/si[:ri])@Ui[:,:ri].T
    rng=np.random.default_rng(430000+budget)
    if noise=='iid': E=rng.normal(size=(len(H),P['linear_trials']));var=np.sum(D*D)
    else:
        E=np.sqrt(.5)*rng.normal(size=(len(H),P['linear_trials']))+.5*Z@rng.normal(size=(Z.shape[1],P['linear_trials']))
        var=.5*np.sum(D*D)+.25*np.sum((D@Z)**2)
    out.update(target_rms_per_sigma=float(np.sqrt(var/len(T))),target_mc_rms_per_sigma=float(np.sqrt(np.mean((D@E)**2))),naive_iid_decoder_mc_rms_per_sigma=float(np.sqrt(np.mean((Di@E)**2))),decoder_target_identity_inf_residual=float(np.linalg.norm(D@H-T,ord=np.inf)))
    return out

def design(H,Z,T,seed,method,maxm):
    m,p=H.shape;selected=list(seed);used=np.zeros(m,bool);used[selected]=True
    if method=='qr_norm':
        order=np.argsort(-np.sum(H*H,axis=1));return selected+[int(i) for i in order if not used[i]][:maxm-p]
    aug=method=='target_correlated';A=np.column_stack([H,Z]) if aug else H
    s2=.5 if aug else 1.
    W=np.zeros((T.shape[0],A.shape[1]));W[:,:p]=T
    if method=='a_opt':W=np.eye(p)
    for k in range(p,maxm):
        if k==p or k%20==0:
            C,_=covariance(H[selected],Z[selected] if aug else None)
        V=A@C;den=s2+np.einsum('ij,ij->i',V,A)
        if method=='d_opt':score=den
        else:score=np.sum((V@W.T)**2,axis=1)/np.maximum(den,1e-30)
        score[used]=-np.inf;i=int(np.argmax(score));v=V[i]
        C-=np.outer(v,v)/den[i];C=(C+C.T)/2
        selected.append(i);used[i]=True
    return selected

def main():
    start=time.perf_counter();buses,K,B,w,keys,nom,t,h0,pats=setup();p=len(w);m=len(pats)
    rng=np.random.default_rng(P['seed']);pairs=[(i,(i+37)%len(K)) for i in range(0,len(K),10)]
    yn,Hn=evaluate(nom,B,pats);yt,Ht=evaluate(w,B,pats);Hn/=h0;Ht/=h0
    vn,Tn=target(nom,B,pairs);vt,Tt=target(w,B,pairs);Tn/=vn[:,None];Tt/=vt[:,None]
    # Every row has one group; common-reference plus first intervention node.
    Z=np.zeros((m,len(K)+1));Z[:,0]=1
    for i,pat in enumerate(pats):Z[i,1+(pat[1][0] if pat[1] else 0)]=1
    oldrec=json.loads((ROOT/'vendor/target_relative_results.json').read_text())['ieee118_adaptive']
    oldseed=oldrec['local_kernel_panel']['selected_indices'];adaptive=oldrec['adaptive_panel']['selected_indices']
    seed,rank,_,_=old.select_local_panel(Hn,p);assert rank==p
    # Preserve the historical panel as a separate comparator: QR may vary across BLAS versions.
    panels={};times={}
    for method in P['methods']:
        tick=time.perf_counter()
        if method=='random':seq=[0]+rng.choice(np.arange(1,m),378,replace=False).tolist()
        else:seq=design(Hn,Z,Tn,seed,method,379)
        panels[method]={str(k):seq[:k] for k in P['budgets'] if k<m}
        times[method]=time.perf_counter()-tick
        print('designed',method,round(times[method],2),flush=True)
    panels['historical']={'179':oldseed,'180':adaptive}
    panels['full']={str(m):list(range(m))}
    save('panels.json',{'panels':panels,'selection_seconds':times,'targets_reduced_indices':pairs,'targets_bus_pairs':[[buses[i],buses[j]] for i,j in pairs],'parameter_keys':keys,'baseline':h0,'t':t})
    # Numerical implementation checks against original loop and centered differences.
    checkids=rng.choice(m,20,replace=False);yy,JJ=evaluate(w,B,pats,checkids)
    yo,Jo=old.evaluate(w,B,pats,checkids,True)
    rel=np.max(np.abs(yy-yo))/h0;jrel=np.linalg.norm(JJ-Jo*w)/np.linalg.norm(JJ)
    cols=rng.choice(p,5,replace=False);fd=[]
    for j in cols:
        u=np.zeros(p);u[j]=1e-4
        yp,_=evaluate(w*np.exp(u),B,pats,checkids,False);ym,_=evaluate(w*np.exp(-u),B,pats,checkids,False)
        fd.append(np.linalg.norm((yp-ym)/2e-4-JJ[:,j])/max(np.linalg.norm(JJ[:,j]),1e-20))
    assert rel<1e-12 and jrel<1e-11 and max(fd)<1e-3
    checks={'vectorized_response_relative_max':rel,'vectorized_logjac_relative_error':jrel,'finite_difference_relative_errors':fd}
    save('checks.json',checks)
    rows=[]
    for method,pp in panels.items():
      for budget,sel in pp.items():
        H=Ht[sel];zs=Z[sel];sv=np.linalg.svd(H,compute_uv=False);rank=int(np.sum(sv>sv[0]*1e-12))
        for noise in ('iid','correlated'):
          if rank<p:
            rows.append(target_only_row(H,Tt,zs,method,int(budget),noise))
            continue
          C,D=covariance(H,zs if noise=='correlated' else None);Ct=C[:p,:p]
          iidC,iidD=covariance(H)
          Gnoise=np.eye(len(sel)) if noise=='iid' else None
          if noise=='iid': noise_draw=rng.normal(size=(len(sel),P['linear_trials']))
          else:noise_draw=np.sqrt(.5)*rng.normal(size=(len(sel),P['linear_trials']))+.5*zs@rng.normal(size=(zs.shape[1],P['linear_trials']))
          E=Tt@D@noise_draw;pred=float(np.sqrt(np.trace(Tt@Ct@Tt.T)/len(Tt)))
          assert np.isfinite(pred) and pred>0
          emp=float(np.sqrt(np.mean(E**2)))
          naive=Tt@iidD@noise_draw
          rows.append({'method':method,'budget':int(budget),'noise':noise,'rank_rtol1e12':rank,'condition':float(sv[0]/sv[-1]),'target_rms_per_sigma':pred,'target_mc_rms_per_sigma':emp,'naive_iid_decoder_mc_rms_per_sigma':float(np.sqrt(np.mean(naive**2))),'parameter_rms_per_sigma':float(np.sqrt(np.trace(Ct)/p))})
        print('evaluated',method,budget,flush=True)
    save('linear.json',rows)
    # Two-known-branch decision: one added reading, equal prior odds, Gaussian noise.
    wf=np.array(oldrec['saved_branch_weights_for_replay']['minimal_false_branch']);wb=np.array(oldrec['saved_branch_weights_for_replay']['minimal_true_branch'])
    yf,_=evaluate(wf,B,pats,jac=False);yb,_=evaluate(wb,B,pats,jac=False)
    sep=np.abs(yf-yb)/h0;sep[oldseed]=0;unused=np.setdiff1d(np.arange(m),oldseed)
    chosen=int(np.argmax(sep));qr_extra=next(i for i in panels['qr_norm']['180'] if i not in oldseed) if any(i not in oldseed for i in panels['qr_norm']['180']) else panels['qr_norm']['200'][-1]
    branch={'known_candidates_only':True,'historical_added_index':adaptive[-1],'max_absolute_gap_index':chosen,'max_relative_gap_index':oldrec['adaptive_branch_breaker']['index'],'normalized_max_gap':float(sep[chosen]),'branch_relative_weight_distance':float(np.linalg.norm(wf-wb)/np.linalg.norm(w)),'rows':[]}
    for sigma in (1e-7,1e-6,1e-5,1e-4):
      err=ndtr(-sep[unused]/(2*sigma))
      branch['rows'].append({'sigma':sigma,'max_gap_error':float(ndtr(-sep[chosen]/(2*sigma))),'historical_error':float(ndtr(-sep[adaptive[-1]]/(2*sigma))),'random_expected_error':float(err.mean()),'qr_extra_error':float(ndtr(-sep[qr_extra]/(2*sigma)))})
    save('branch.json',branch)
    # Every single channel deletion; rank and target-nullspace tests, no Gaussian prior.
    drops=[]
    for method in ('historical','target_a'):
      sel=panels[method]['180'];H=Ht[sel]
      for j in range(len(sel)):
        A=np.delete(H,j,axis=0);_,s,V=np.linalg.svd(A,full_matrices=True);rr=int(np.sum(s>s[0]*1e-12));N=V[rr:].T
        gap=float(np.linalg.norm(Tt@N))
        drops.append({'method':method,'removed_index':sel[j],'rank':rr,'target_nullspace_gap':gap,'target_recoverable_at_tolerance':bool(gap<1e-8)})
    save('dropout.json',drops)
    shifts=[]
    for sd in (.03,.1):
      for trial in range(5):
        ww=w*np.exp(rng.normal(0,sd,p));_,HH=evaluate(ww,B,pats);HH/=h0;vv,TT=target(ww,B,pairs);TT/=vv[:,None]
        for method in ('qr_norm','target_a','target_correlated'):
          for b in (180,250,379):
            sel=panels[method][str(b)];C,_=covariance(HH[sel]);metric=float(np.sqrt(np.trace(TT@C@TT.T)/len(TT)))
            shifts.append({'logweight_sd':sd,'trial':trial,'method':method,'budget':b,'target_rms_per_sigma':metric})
    save('shift.json',shifts)
    # Actual nonlinear estimation from nominal, not a truth-neighbourhood initialization.
    nonlinear=[]
    for method in P['nonlinear']['methods']:
      for budget in P['nonlinear']['budgets']:
        sel=panels[method][str(budget)];z=Z[sel]
        L=cholesky(.5*np.eye(len(sel))+.25*z@z.T,lower=True)
        for sigma in P['nonlinear']['sigma']:
          for trial in range(P['nonlinear']['trials']):
            rg=np.random.default_rng(900000+budget*100+trial)
            noise=np.sqrt(.5)*rg.normal(size=len(sel))+.5*z@rg.normal(size=z.shape[1]);obs=yt[sel]+sigma*h0*noise
            cache={}
            def both(x):
                if 'x' not in cache or not np.array_equal(x,cache['x']):
                    yy,JJ=evaluate(np.exp(x),B,pats,sel)
                    cache.update(x=x.copy(),f=solve_triangular(L,(yy-obs)/h0,lower=True),j=solve_triangular(L,JJ/h0,lower=True))
                return cache['f'],cache['j']
            tick=time.perf_counter()
            sol=least_squares(lambda x:both(x)[0],np.log(nom),jac=lambda x:both(x)[1],max_nfev=100,xtol=1e-11,ftol=1e-11,gtol=1e-11,bounds=(np.log(nom)-8,np.log(nom)+8))
            est=np.exp(sol.x);vhat,_=target(est,B,pairs)
            nonlinear.append({'method':method,'budget':budget,'sigma':sigma,'trial':trial,'target_relative_rms':float(np.sqrt(np.mean(((vhat-vt)/vt)**2))),'weight_relative_l2':float(np.linalg.norm(est-w)/np.linalg.norm(w)),'normalized_residual_norm':float(np.linalg.norm(sol.fun)),'nfev':int(sol.nfev),'solver_success':bool(sol.success),'seconds':time.perf_counter()-tick})
          save('nonlinear.json',nonlinear)
          print('nonlinear',method,budget,sigma,flush=True)
    save('receipt.json',{'status':'COMPLETED','seconds':time.perf_counter()-start,'protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),'matpower_sha256':hashlib.sha256((ROOT/'vendor/case118.m').read_bytes()).hexdigest(),'numpy':np.__version__,'linear_rows':len(rows),'nonlinear_fits':len(nonlinear),'dropout_checks':len(drops),'shift_rows':len(shifts),'limitations':['Exploratory simulations, no independent confirmation campaign','Local Jacobian uncertainty requires correct linearization','Only previously found nonlinear branches tested','Abstract correlated-noise model, not calibrated instrument data','Targets are 12 declared effective resistances; no AC voltage-estimation claim']})
if __name__=='__main__':main()
