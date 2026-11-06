#!/usr/bin/env python3
"""External validation of Nonagon finite diagonal-intervention reconstruction.

The frozen cases are external MATPOWER networks mapped to the incremental DC
susceptance operator.  The MATPOWER slack/reference bus is removed, turning
incident branches into positive grounds.  One remaining bus is a fixed scalar
observation/injection port; all other remaining buses receive calibrated
positive diagonal shunts.

This is a structural inverse benchmark, not an AC power-flow benchmark.
"""
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from fixtures import CASES

HERE = Path(__file__).resolve().parent


def build_reduced(case):
    buses=list(case['buses']); ref=case['ref']; keep=[b for b in buses if b!=ref]
    idx={b:i for i,b in enumerate(keep)}; n=len(keep)
    K=np.zeros((n,n),float)
    for u,v,w in case['edges']:
        if w <= 0: raise ValueError('strict benchmark requires positive DC branch weights')
        if u==ref and v==ref: continue
        if u==ref:
            K[idx[v],idx[v]] += w
        elif v==ref:
            K[idx[u],idx[u]] += w
        else:
            i,j=idx[u],idx[v]
            K[i,i]+=w; K[j,j]+=w; K[i,j]-=w; K[j,i]-=w
    return keep,K


def permute_port(keep,K,port_bus):
    p=keep.index(port_bus)
    order=[p]+[i for i in range(len(keep)) if i!=p]
    return [keep[i] for i in order],K[np.ix_(order,order)]


def response_pair(G,i,j,t):
    # i,j are hidden-coordinate indices in the permuted G (>=1)
    a=1.0/t+G[i,i]; b=1.0/t+G[j,j]; c=G[i,j]
    det=a*b-c*c
    vi,vj=G[0,i],G[0,j]
    delta=(b*vi*vi + a*vj*vj - 2*c*vi*vj)/det
    return G[0,0]-delta


def make_measurements(G,t):
    n=G.shape[0]; h0=float(G[0,0])
    singles=[]
    for j in range(1,n):
        g=G[0,j]**2; r=G[j,j]
        h1=h0-t*g/(1+t*r)
        h2=h0-2*t*g/(1+2*t*r)
        singles.append([h1,h2])
    pairs=[]
    for i in range(1,n):
        for j in range(i+1,n):
            pairs.append(response_pair(G,i,j,t))
    return h0,np.asarray(singles),np.asarray(pairs)


def reconstruct(h0,singles,pairs,t):
    n=singles.shape[0]+1
    G=np.empty((n,n),float); G.fill(np.nan); G[0,0]=h0
    g=np.empty(n,float); g[0]=h0*h0
    for j in range(1,n):
        d1=h0-singles[j-1,0]; d2=h0-singles[j-1,1]
        if not (d1>0 and d2>0) or abs(d2-d1) < np.finfo(float).eps*max(abs(h0),1.0):
            raise FloatingPointError(f'single-load degeneracy at coordinate {j}')
        inv_g=(1/d1-1/d2)/(1/t-1/(2*t))
        gj=1/inv_g
        r=(1/d1-inv_g/t)*gj
        if gj<=0 or r<=0 or not np.isfinite(gj+r):
            raise FloatingPointError(f'nonphysical single-load inverse at coordinate {j}')
        g[j]=gj; G[j,j]=r; G[0,j]=G[j,0]=math.sqrt(gj)
    k=0
    for i in range(1,n):
        for j in range(i+1,n):
            dl=h0-pairs[k]; k+=1
            if dl<=0: raise FloatingPointError(f'nonpositive pair response drop at {i},{j}')
            a=1/t+G[i,i]; b=1/t+G[j,j]
            disc=(dl*a-g[i])*(dl*b-g[j])
            scale=max(abs((dl*a-g[i])*(dl*b-g[j])),g[i]*g[j],1e-300)
            if disc < -5e-11*scale:
                raise FloatingPointError(f'negative pair discriminant at {i},{j}: {disc}')
            disc=max(0.0,disc)
            x=(G[0,i]*G[0,j]-math.sqrt(disc))/dl
            G[i,j]=G[j,i]=x
    if not np.all(np.isfinite(G)): raise FloatingPointError('nonfinite recovered Green matrix')
    return G


def response_from_G(G,inds,loads):
    if not inds: return float(G[0,0])
    S=np.asarray(inds,dtype=int); T=np.asarray(loads,dtype=float)
    A=G[np.ix_(S,S)] + np.diag(1/T)
    v=G[0,S]
    return float(G[0,0]-v@np.linalg.solve(A,v))


def topology_metrics(K,Khat):
    n=len(K); tri=np.triu_indices(n,1)
    true_w=-K[tri]; pred_w=-Khat[tri]
    mask=true_w>1e-10
    med=float(np.median(true_w[mask])); threshold=max(1e-10,1e-6*med)
    pmask=pred_w>threshold
    tp=int(np.sum(mask&pmask)); fp=int(np.sum(~mask&pmask)); fn=int(np.sum(mask&~pmask))
    precision=tp/(tp+fp) if tp+fp else 1.0; recall=tp/(tp+fn) if tp+fn else 1.0
    edge_rel=np.abs(pred_w[mask]-true_w[mask])/true_w[mask]
    nonedge=np.abs(pred_w[~mask])
    gt=K.sum(axis=1); gh=Khat.sum(axis=1); gmask=gt>1e-10; gpred=gh>threshold
    gtp=int(np.sum(gmask&gpred)); gfp=int(np.sum(~gmask&gpred)); gfn=int(np.sum(gmask&~gpred))
    return {
      'threshold':threshold,'true_internal_edges':int(mask.sum()),'predicted_internal_edges':int(pmask.sum()),
      'edge_precision':precision,'edge_recall':recall,
      'edge_relative_error_median':float(np.median(edge_rel)),'edge_relative_error_max':float(np.max(edge_rel)),
      'max_nonedge_leakage':float(np.max(nonedge)) if nonedge.size else 0.0,
      'positive_offdiagonal_violations':int(np.sum(Khat[tri]>threshold)),
      'true_reference_incident_edges':int(gmask.sum()),'predicted_reference_incident_edges':int(gpred.sum()),
      'reference_edge_precision':gtp/(gtp+gfp) if gtp+gfp else 1.0,
      'reference_edge_recall':gtp/(gtp+gfn) if gtp+gfn else 1.0,
      'reference_weight_relative_error_max':float(np.max(np.abs(gh[gmask]-gt[gmask])/gt[gmask])) if gmask.any() else 0.0,
    }


def heldout(G,Ghat,t,seed=0,count=100):
    rng=np.random.default_rng(seed); n=len(G); errors=[]; rel=[]
    for _ in range(count):
        m=int(rng.integers(1,min(7,n)))
        inds=sorted(rng.choice(np.arange(1,n),size=m,replace=False).tolist())
        loads=t*10**rng.uniform(-0.7,0.7,size=m)
        y=response_from_G(G,inds,loads); z=response_from_G(Ghat,inds,loads)
        errors.append(abs(y-z)); rel.append(abs(y-z)/max(abs(y),1e-300))
    return {'count':count,'max_abs_error':float(max(errors)),'rms_abs_error':float(np.sqrt(np.mean(np.square(errors)))),
            'max_relative_error':float(max(rel))}


def noisy_trial(base,singles,pairs,t,eta,rng,K,G,threshold_only=False):
    # independent componentwise bounded errors in [-eta,+eta]
    h0=base+rng.uniform(-eta,eta)
    s=singles+rng.uniform(-eta,eta,size=singles.shape)
    p=pairs+rng.uniform(-eta,eta,size=pairs.shape)
    try:
        Gh=reconstruct(h0,s,p,t); Kh=np.linalg.inv(Gh); Kh=(Kh+Kh.T)/2
        if not np.all(np.isfinite(Kh)): raise FloatingPointError('nonfinite K')
        tm=topology_metrics(K,Kh)
        return {'ok':True,'matrix_rel_fro':float(np.linalg.norm(Kh-K)/np.linalg.norm(K)),
                'edge_f1':float(2*tm['edge_precision']*tm['edge_recall']/max(tm['edge_precision']+tm['edge_recall'],1e-300)),
                'sign_violations':tm['positive_offdiagonal_violations']}
    except (FloatingPointError,np.linalg.LinAlgError,ValueError,ZeroDivisionError):
        return {'ok':False}


def run_case(name,noise_trials=12):
    keep,K0=build_reduced(CASES[name]); port=min(keep); buses,K=permute_port(keep,K0,port)
    eig=np.linalg.eigvalsh(K)
    if eig[0] <= 0: raise ValueError('reduced DC operator not positive definite')
    tic=time.perf_counter(); G=np.linalg.inv(K); inverse_s=time.perf_counter()-tic
    positivity=float(np.min(G));
    t=1.0/float(np.median(np.diag(G)[1:]))
    h0,singles,pairs=make_measurements(G,t)
    patterns=1+2*(len(K)-1)+(len(K)-1)*(len(K)-2)//2
    tic=time.perf_counter(); Gh=reconstruct(h0,singles,pairs,t); reconstruction_s=time.perf_counter()-tic
    Kh=np.linalg.inv(Gh); Kh=(Kh+Kh.T)/2
    matrix_rel=float(np.linalg.norm(Kh-K)/np.linalg.norm(K)); matrix_max=float(np.max(np.abs(Kh-K)))
    tm=topology_metrics(K,Kh); hv=heldout(G,Gh,t,seed=1000+len(K))
    # response dynamic range is the key finite-precision stress statistic
    all_drops=np.concatenate(([h0-singles[:,0].min(),], h0-singles[:,0], h0-singles[:,1], h0-pairs))
    positive_drops=all_drops[all_drops>0]
    noise=[]
    for rel_eta in [1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8]:
        eta=rel_eta*abs(h0); rng=np.random.default_rng(20261030+len(K)+int(-math.log10(rel_eta)))
        rows=[noisy_trial(h0,singles,pairs,t,eta,rng,K,G) for _ in range(noise_trials)]
        ok=[r for r in rows if r.get('ok')]
        noise.append({'relative_eta_to_baseline':rel_eta,'absolute_eta':eta,'trials':noise_trials,
                      'successful_reconstructions':len(ok),
                      'median_matrix_rel_fro':float(np.median([r['matrix_rel_fro'] for r in ok])) if ok else None,
                      'median_edge_f1':float(np.median([r['edge_f1'] for r in ok])) if ok else None,
                      'max_sign_violations':max((r['sign_violations'] for r in ok),default=None)})
    return {
      'case':name,'source_model':'MATPOWER active-branch incremental DC Bbus with slack removed',
      'original_buses':len(CASES[name]['buses']),'reference_bus':CASES[name]['ref'],'reduced_dimension':len(K),
      'fixed_port_bus':port,'active_branches_source':len(CASES[name]['edges']),
      'strict_positive_edge_model':True,'smallest_eigenvalue':float(eig[0]),'condition_number_K':float(eig[-1]/eig[0]),
      'green_min_entry':positivity,'load_strength_t':t,'load_rule':'t = 1 / median(G_jj), singles t and 2t, pairs t,t',
      'scalar_response_count':patterns,'distinct_load_patterns':patterns,
      'baseline_response':h0,'minimum_positive_protocol_drop':float(np.min(positive_drops)),
      'minimum_drop_over_baseline':float(np.min(positive_drops)/abs(h0)),
      'matrix_relative_frobenius_error':matrix_rel,'matrix_max_abs_error':matrix_max,
      'topology':tm,'heldout_multishunt_predictions':hv,
      'timing_seconds':{'dense_inverse_true':inverse_s,'explicit_G_reconstruction':reconstruction_s},
      'bounded_noise_study':noise,
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--noise-trials',type=int,default=12); ap.add_argument('--output',type=Path,default=HERE/'matpower_dc_results.json')
    args=ap.parse_args(); out={'schema':'nonagon.external.matpower_dc.v1','cases':[]}
    for name in ['case14','case30','case57','case118']:
        row=run_case(name,args.noise_trials); out['cases'].append(row)
        print(name, 'n=',row['reduced_dimension'],'relF=',f"{row['matrix_relative_frobenius_error']:.3e}",
              'edge P/R=',f"{row['topology']['edge_precision']:.3f}/{row['topology']['edge_recall']:.3f}",
              'min drop/base=',f"{row['minimum_drop_over_baseline']:.3e}")
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(args.output)

if __name__=='__main__': main()
