#!/usr/bin/env python3
"""Target-relative measurement reduction for the MATPOWER DC benchmark.

This adapts PerfectJet's target-relative kernel criterion to the *differential*
of Nonagon's nonlinear finite-shunt response map.  For edge parameters w and
scalar responses F(w), local recovery of every supplied-topology edge weight is
certified at a nominal point when ker DF is {0}.  A rank-revealing row panel
therefore gives a locally sufficient sensor/intervention set.

Important scope: the target coordinates are the aggregate edges/grounds of the
supplied MATPOWER topology.  This is a known-topology calibration benchmark,
not blind topology discovery.  The inverse-function theorem gives local
uniqueness; a final adaptive reading is used to eliminate an observed remote
exact-fit branch, but that experiment is not a global injectivity theorem over
all positive weights.
"""
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.linalg import qr
from scipy.optimize import least_squares
import run_dc_reconstruction as R
from fixtures import CASES

HERE=Path(__file__).resolve().parent


def edge_parameterization(case):
    keep,K0=R.build_reduced(case); port=min(keep); buses,K=R.permute_port(keep,K0,port)
    idx={b:i for i,b in enumerate(buses)}; ref=case['ref']; acc={}
    for u,v,w in case['edges']:
        if u==ref or v==ref:
            key=('ground', v if u==ref else u)
        else:
            a,b=sorted((u,v)); key=('edge',a,b)
        acc[key]=acc.get(key,0.0)+float(w)
    keys=sorted(acc,key=str); weights=np.asarray([acc[k] for k in keys],float)
    B=np.zeros((len(buses),len(keys)),float)
    for c,k in enumerate(keys):
        if k[0]=='ground': B[idx[k[1]],c]=1.0
        else: B[idx[k[1]],c]=1.0; B[idx[k[2]],c]=-1.0
    Kcheck=(B*weights)@B.T
    if not np.allclose(Kcheck,K,rtol=1e-13,atol=5e-13):
        raise RuntimeError('aggregate edge parameterization does not reproduce reduced operator')
    return buses,K,B,weights,keys


def patterns(n,t):
    out=[('baseline',(),())]
    for j in range(1,n):
        out.append(('single',(j,),(t,)))
        out.append(('single',(j,),(2*t,)))
    for i in range(1,n):
        for j in range(i+1,n): out.append(('pair',(i,j),(t,t)))
    return out


def evaluate(weights,B,pats,indices=None,want_jacobian=False):
    if indices is None: indices=range(len(pats))
    K=(B*weights)@B.T; G=np.linalg.inv(K); g0=G[:,0]
    ys=[]; rows=[]
    for ii in indices:
        _,S,T=pats[ii]
        if not S:
            z=g0
        else:
            S=np.asarray(S,dtype=int); T=np.asarray(T,dtype=float)
            c=np.linalg.solve(G[np.ix_(S,S)]+np.diag(1/T),g0[S])
            z=g0-G[:,S]@c
        ys.append(float(z[0]))
        if want_jacobian:
            # d(e0^T A^-1 e0)/dw_e = -(b_e^T A^-1 e0)^2
            q=B.T@z; rows.append(-(q*q))
    return np.asarray(ys), (np.asarray(rows) if want_jacobian else None)


def select_local_panel(J,parameter_count):
    # Column scaling changes conditioning, not the kernel/rank condition.
    norms=np.linalg.norm(J,axis=0); norms[norms==0]=1.0; A=J/norms
    # Baseline is forced because it is also the pilot used to set the nominal scale.
    q=A[0]/np.linalg.norm(A[0]); residual=A-(A@q)[:,None]*q
    _,_,piv=qr(residual.T,mode='economic',pivoting=True,check_finite=False)
    selected=[0]+[int(i) for i in piv if i!=0][:parameter_count-1]
    block=A[selected]
    rank=int(np.linalg.matrix_rank(block,tol=1e-10)); sv=np.linalg.svd(block,compute_uv=False)
    return selected,rank,float(sv[-1]),float(sv[0]/sv[-1])


def solve_logweights(x0,B,pats,selected,yobs,max_nfev=240):
    yscale=max(abs(float(yobs[0])),float(np.ptp(yobs)),1e-12)
    def fun(x):
        y,_=evaluate(np.exp(x),B,pats,selected,False)
        return (y-yobs)/yscale
    def jac(x):
        w=np.exp(x); _,J=evaluate(w,B,pats,selected,True)
        return J*w[None,:]/yscale
    tic=time.perf_counter()
    sol=least_squares(fun,x0,jac=jac,method='trf',xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=max_nfev)
    return sol,float(np.linalg.norm(fun(sol.x))),time.perf_counter()-tic


def nominal_design(B,Ktrue):
    Gtrue=np.linalg.inv(Ktrue); Gunit=np.linalg.inv(B@B.T)
    # One baseline scalar supplies a global conductance scale for a unit-edge nominal graph.
    scale=float(Gunit[0,0]/Gtrue[0,0]); nominal=np.full(B.shape[1],scale)
    Gnom=np.linalg.inv((B*nominal)@B.T)
    t=float(1/np.median(np.diag(Gnom)[1:]))
    return nominal,scale,t,float(Gtrue[0,0])


def local_case_summary(name):
    buses,K,B,w,keys=edge_parameterization(CASES[name]); nominal,scale,t,h0=nominal_design(B,K)
    pats=patterns(len(buses),t); _,J=evaluate(nominal,B,pats,want_jacobian=True)
    selected,rank,smin,cond=select_local_panel(J,len(w))
    return {
      'case':name,'reduced_dimension':len(buses),'aggregate_parameter_count':len(w),
      'full_protocol_count':len(pats),'local_panel_count':len(selected),
      'full_to_local_reduction_factor':len(pats)/len(selected),'rank':rank,
      'normalized_sigma_min':smin,'normalized_condition_number':cond,
      'panel_types':{k:sum(pats[i][0]==k for i in selected) for k in ('baseline','single','pair')},
      'baseline_scale_nominal':scale,'load_strength_t':t,
    }


def run_ieee118():
    buses,K,B,wtrue,keys=edge_parameterization(CASES['case118']); p=len(wtrue)
    nominal,scale,t,h0=nominal_design(B,K); pats=patterns(len(buses),t)
    ytrue_all,_=evaluate(wtrue,B,pats); _,Jnom=evaluate(nominal,B,pats,want_jacobian=True)
    selected,rank,smin,cond=select_local_panel(Jnom,p)
    yobs=ytrue_all[selected]; xnom=np.log(nominal)

    # Deterministic second start chosen before branch inspection.
    rng=np.random.default_rng(20261030)
    rng.normal(0,0.25,p); rng.normal(0,0.50,p) # consume the frozen development sequence
    xfar=xnom+rng.normal(0,1.0,p)
    sol_a,res_a,sec_a=solve_logweights(xnom,B,pats,selected,yobs)
    sol_b,res_b,sec_b=solve_logweights(xfar,B,pats,selected,yobs)
    wa,wb=np.exp(sol_a.x),np.exp(sol_b.x)

    Ya,_=evaluate(wa,B,pats); Yb,_=evaluate(wb,B,pats)
    spread=np.abs(Ya-Yb)/np.maximum(np.maximum(np.abs(Ya),np.abs(Yb)),1e-12)
    spread[np.asarray(selected,dtype=int)]=-1.0
    added=int(np.argmax(spread)); selected_plus=selected+[added]; yplus=ytrue_all[selected_plus]
    fina,rfa,seca=solve_logweights(sol_a.x,B,pats,selected_plus,yplus)
    finb,rfb,secb=solve_logweights(sol_b.x,B,pats,selected_plus,yplus)
    wfa,wfb=np.exp(fina.x),np.exp(finb.x)

    # A third independent moderate start checks that the branch-breaking panel is not
    # only repairing the two warm starts used to choose the additional reading.
    rng2=np.random.default_rng(20261031); rng2.normal(0,0.15,p); xthird=xnom+rng2.normal(0,0.30,p)
    third,rft,sect=solve_logweights(xthird,B,pats,selected_plus,yplus)
    wft=np.exp(third.x)

    def err(w):
        return {'relative_l2':float(np.linalg.norm(w-wtrue)/np.linalg.norm(wtrue)),
                'max_relative':float(np.max(np.abs(w-wtrue)/wtrue))}
    return {
      'schema':'nonagon.target-relative-matpower.v1','case':'case118',
      'scope':'known-topology aggregate line/ground calibration; local target-identifiability theorem plus adaptive branch-disambiguation experiment',
      'not_claimed':'blind topology discovery or global injectivity over all positive edge weights',
      'reduced_dimension':len(buses),'aggregate_parameter_count':p,
      'aggregate_internal_edges':sum(k[0]=='edge' for k in keys),
      'aggregate_reference_ground_edges':sum(k[0]=='ground' for k in keys),
      'full_protocol_count':len(pats),'nominal_design':{'rule':'unit topology scaled from the one baseline response','scale':scale,'load_strength_t':t,'baseline_response':h0},
      'local_kernel_panel':{
        'count':len(selected),'rank':rank,'information_lower_bound_scalar_count':p,
        'hits_scalar_lower_bound':len(selected)==p,'normalized_sigma_min':smin,'normalized_condition_number':cond,
        'types':{k:sum(pats[i][0]==k for i in selected) for k in ('baseline','single','pair')},
        'full_to_panel_reduction_factor':len(pats)/len(selected),
        'selected_indices':selected,
      },
      'minimal_panel_multistart':[
        {'start':'baseline_scaled_uniform','residual_norm':res_a,'solve_seconds':sec_a,'nfev':int(sol_a.nfev),'error_to_truth':err(wa)},
        {'start':'frozen_lognormal_sigma_1','residual_norm':res_b,'solve_seconds':sec_b,'nfev':int(sol_b.nfev),'error_to_truth':err(wb)},
      ],
      'observed_remote_exact_fit_branch': bool(res_a < 1e-9 and res_b < 1e-9 and np.linalg.norm(wa-wb)/np.linalg.norm(wtrue)>1e-5),
      'adaptive_branch_breaker':{
        'selection_rule':'unmeasured candidate maximizing relative response disagreement between the two exact-fit branches',
        'index':added,'pattern':[pats[added][0],list(pats[added][1]),[float(x) for x in pats[added][2]]],
        'relative_branch_spread':float(spread[added]),'true_response':float(ytrue_all[added]),
        'false_branch_response_offset':float(Ya[added]-ytrue_all[added]),
        'true_branch_response_offset':float(Yb[added]-ytrue_all[added]),
      },
      'adaptive_panel':{
        'count':len(selected_plus),'full_to_panel_reduction_factor':len(pats)/len(selected_plus),
        'selected_indices':selected_plus,
        'warm_start_from_false_branch':{'residual_norm':rfa,'solve_seconds':seca,'nfev':int(fina.nfev),'error_to_truth':err(wfa)},
        'warm_start_from_true_branch':{'residual_norm':rfb,'solve_seconds':secb,'nfev':int(finb.nfev),'error_to_truth':err(wfb)},
        'independent_moderate_start':{'residual_norm':rft,'solve_seconds':sect,'nfev':int(third.nfev),'error_to_truth':err(wft)},
      },
      'saved_branch_weights_for_replay':{
        'minimal_false_branch':[float(x) for x in wa],
        'minimal_true_branch':[float(x) for x in wb],
        'adaptive_from_false':[float(x) for x in wfa],
        'adaptive_from_true':[float(x) for x in wfb],
        'adaptive_independent':[float(x) for x in wft],
      },
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,default=HERE/'target_relative_results.json'); args=ap.parse_args()
    out={
      'schema':'nonagon.target-relative-suite.v1',
      'theorem_transfer':{
        'source':'PerfectJet chapter 7 target-relative observability',
        'linear_statement':'ker(O) subset ker(G_target)',
        'nonagon_use':'apply to the Jacobian of the nonlinear response map; full column rank gives local edge-weight uniqueness by the inverse function theorem',
      },
      'local_panels':[local_case_summary(k) for k in ('case14','case30','case57','case118')],
      'ieee118_adaptive':run_ieee118(),
    }
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(args.output)
    x=out['ieee118_adaptive']; print('IEEE118',x['full_protocol_count'],'->',x['local_kernel_panel']['count'],'->',x['adaptive_panel']['count'])
    print('adaptive errors',x['adaptive_panel']['warm_start_from_false_branch']['error_to_truth'],x['adaptive_panel']['independent_moderate_start']['error_to_truth'])

if __name__=='__main__': main()
