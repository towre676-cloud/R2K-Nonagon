"""Blinded, fixed-budget frequency-response identification challenge.

Fit functions receive only (s, noisy response, order, public band).
Hidden systems and clean held-out responses are confined to make_case/evaluate.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, time, warnings
from pathlib import Path
import numpy as np
import scipy
from scipy import linalg as la
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
BAND = (0.03, 8.0)
NSAMPLE = 64
NTEST = 257
SCENARIOS = ["ordinary", "colored_noise", "weak_link", "near_collision", "model_mismatch"]
NOISES = [0.0, 0.001, 0.01, 0.05]
METHODS = ["projected_loewner", "loewner_plus_vf", "balanced_plus_vf", "fixed_seed_vf"]

def digest(*arrays):
    h = hashlib.sha256()
    for x in arrays:
        h.update(np.ascontiguousarray(x).tobytes())
    return h.hexdigest()

def chain_response(s, m, k, c):
    n = len(m)
    K = np.zeros((n, n)); C = np.zeros((n, n))
    for j in range(n):
        K[j, j] += k[j]; C[j, j] += c[j]
        if j:
            K[j-1, j-1] += k[j]; K[j, j-1] -= k[j]; K[j-1, j] -= k[j]
            C[j-1, j-1] += c[j]; C[j, j-1] -= c[j]; C[j-1, j] -= c[j]
    M = np.diag(m); e = np.eye(n)[:, 0]
    H = np.array([la.solve(z*z*M + z*C + K, e, assume_a="gen")[0] for z in s])
    A = np.block([[np.zeros((n,n)), np.eye(n)], [-la.solve(M,K), -la.solve(M,C)]])
    return H, la.eigvals(A), K, C

def make_case(scenario, index, seed):
    """Truth-generating component. Fixed band/grid never use hidden poles."""
    rng = np.random.default_rng(np.random.SeedSequence([seed, SCENARIOS.index(scenario), index]))
    n = int(rng.integers(3, 6))
    m = rng.uniform(.7, 1.5, n); k = rng.uniform(1., 4., n); c = rng.uniform(.03, .16, n)
    if scenario == "weak_link":
        # A weak interior connection makes remote modes poorly visible at endpoint.
        k[n//2] *= .035; c[n//2] *= .035
    if scenario == "near_collision":
        # Two almost identical grounded blocks, weakly coupled. Still a genuine
        # chain, with extra grounding terms creating nearly coincident modes.
        n = 4
        m0 = rng.uniform(.8, 1.3, 2); k0 = rng.uniform(1., 3., 2)
        m = np.r_[m0, m0*(1+rng.uniform(-.002,.002,2))]
        k = np.r_[k0, .002, k0[1]]
        c = np.full(4, .015)
    train_w = np.geomspace(*BAND, NSAMPLE)
    # Fixed clean evaluation grid offset from the training lattice.
    test_w = np.exp(np.linspace(np.log(BAND[0]), np.log(BAND[1]), NTEST+2)[1:-1])
    s = 1j*train_w; st = 1j*test_w
    y, poles, K, C = chain_response(s, m,k,c)
    yt, _, _, _ = chain_response(st, m,k,c)
    if scenario == "near_collision":
        # Ground the far block like the near block: explicit extra spring/damper.
        K[2,2] += k0[0]; C[2,2] += c[0]
        e = np.eye(n)[:,0]; M=np.diag(m)
        y=np.array([la.solve(z*z*M+z*C+K,e)[0] for z in s])
        yt=np.array([la.solve(z*z*M+z*C+K,e)[0] for z in st])
        poles=la.eigvals(np.block([[np.zeros((n,n)),np.eye(n)],[-la.solve(M,K),-la.solve(M,C)]]))
    if scenario == "model_mismatch":
        # Small independent parasitic compliance, excluded from the supplied order.
        om = rng.uniform(1.,5.); damp=.025*om
        y += .025/(s*s+2*damp*s+om*om)
        yt += .025/(st*st+2*damp*st+om*om)
    z = (rng.standard_normal(NSAMPLE)+1j*rng.standard_normal(NSAMPLE))/np.sqrt(2)
    if scenario == "colored_noise":
        for j in range(1,NSAMPLE): z[j]=.9*z[j-1]+np.sqrt(1-.9**2)*z[j]
    return dict(s=s, st=st, clean=y, test_clean=yt, noise=z, poles=poles, order=2*n,
                truth=dict(m=m.tolist(),k=k.tolist(),c=c.tolist(),K=K.tolist(),C=C.tolist()))

def make_lr(s,y):
    # Conjugate augmentation carries no additional observation budget.
    ia=np.arange(0,len(s),2); ib=np.arange(1,len(s),2)
    sl=np.r_[s[ia],s[ia].conj()]; sr=np.r_[s[ib],s[ib].conj()]
    yl=np.r_[y[ia],y[ia].conj()]; yr=np.r_[y[ib],y[ib].conj()]
    L=(yl[:,None]-yr[None,:])/(sl[:,None]-sr[None,:])
    Ls=(sl[:,None]*yl[:,None]-sr[None,:]*yr[None,:])/(sl[:,None]-sr[None,:])
    return sl,sr,yl,yr,L,Ls

def projected(s,y,order,band,weights=None):
    sl,sr,yl,yr,L,Ls=make_lr(s,y)
    dl=np.ones(len(sl)); dr=np.ones(len(sr))
    if weights is not None:
        poles,residues=weights
        # Equalize sampled Cauchy-feature row energy, using the pilot estimated
        # residues only. This is a same-data heuristic inspired by residue balance,
        # not the archive's complex-plane adaptive sampling theorem.
        def wgt(nodes):
            energy=np.sum(np.abs(residues)[None,:]/np.maximum(np.abs(nodes[:,None]-poles[None,:])**2,1e-20),axis=1)
            w=1/np.sqrt(np.maximum(energy,1e-30)); w/=np.median(w)
            return np.clip(w,.1,10.)
        dl=wgt(sl); dr=wgt(sr)
    L=dl[:,None]*L*dr[None,:]; Ls=dl[:,None]*Ls*dr[None,:]
    U,sv,_=la.svd(np.c_[L,Ls],full_matrices=False)
    _,_,Vh=la.svd(np.r_[L,Ls],full_matrices=False)
    U=U[:,:order]; V=Vh.conj().T[:,:order]
    E=-U.conj().T@L@V; A=-U.conj().T@Ls@V
    B=U.conj().T@(dl*yl); C=(yr*dr)@V
    poles=la.eigvals(A,E)
    def predict(z):
        return np.array([C@la.solve(t*E-A,B) for t in z])
    return dict(poles=poles,predict=predict,metadata={"retained_singular_value":float(sv[order-1])})

def residue_fit(s,y,poles):
    sa=np.r_[s,s.conj()]; ya=np.r_[y,y.conj()]
    F=1/(sa[:,None]-poles[None,:])
    scale=np.linalg.norm(F,axis=0); scale=np.maximum(scale,1e-30)
    res=la.lstsq(F/scale[None,:],ya,cond=1e-12)[0]/scale
    return res

def stabilize_and_pair(poles,order,band):
    p=np.asarray(poles,dtype=complex)
    p=p[np.isfinite(p.real)&np.isfinite(p.imag)]
    if len(p)!=order: raise ValueError("nonfinite or wrong number of estimated poles")
    # Projection to stable poles is identical for both refinement methods.
    p=-np.maximum(np.abs(p.real),1e-5*band[0])+1j*p.imag
    # Conjugate pairing, including independent stable real poles.
    pos=list(p[p.imag>1e-8]); neg=list(p[p.imag < -1e-8]); real=list(p[np.abs(p.imag)<=1e-8].real)
    out=[]
    if pos and neg:
        a,b=linear_sum_assignment(abs(np.asarray(pos)[:,None]-np.asarray(neg)[None,:].conj()))
        usedp=set(a); usedn=set(b)
        for i,j in zip(a,b):
            q=(pos[i]+neg[j].conjugate())/2; out.extend([q,q.conjugate()])
        real += [q.real for j,q in enumerate(pos) if j not in usedp]
        real += [q.real for j,q in enumerate(neg) if j not in usedn]
    else: real += [q.real for q in pos+neg]
    out.extend(complex(q) for q in real)
    if len(out)!=order: raise ValueError("pole pairing failed")
    return np.asarray(out)

def vector_fit(s,y,order,band,initial,max_iter=8):
    """Standard barycentric vector-fitting iteration, same for both ablations.

    In each iteration solve linearized numerator/denominator residuals, relocate
    poles using a rank-one eigenproblem, project stable, then solve residues.
    Select by noisy training residual only, never by protected evaluation error.
    """
    poles=stabilize_and_pair(initial,order,band)
    sa=np.r_[s,s.conj()]; ya=np.r_[y,y.conj()]
    scale_y=max(np.linalg.norm(ya)/np.sqrt(len(ya)),1e-30); yn=ya/scale_y
    best=None
    for it in range(max_iter+1):
        res=residue_fit(s,y,poles)
        yp=np.sum(res[None,:]/(sa[:,None]-poles[None,:]),axis=1)
        loss=float(np.linalg.norm(yp-ya)/np.linalg.norm(ya))
        if not np.isfinite(loss):raise FloatingPointError("nonfinite vector-fit loss")
        if best is None or loss<best[0]:best=(loss,poles.copy(),res.copy(),it)
        if it==max_iter: break
        F=1/(sa[:,None]-poles[None,:]); Z=np.c_[F,-yn[:,None]*F]
        col=np.maximum(np.linalg.norm(Z,axis=0),1e-30)
        coef=la.lstsq(Z/col[None,:],yn,cond=1e-12)[0]/col
        new=la.eigvals(np.diag(poles)-np.ones((order,1))@coef[order:][None,:])
        poles=stabilize_and_pair(new,order,band)
    loss,poles,res,it=best
    return dict(poles=poles,predict=lambda z:np.sum(res[None,:]/(z[:,None]-poles[None,:]),axis=1),
                metadata={"training_nrmse":loss,"selected_iteration":it,"max_iterations":max_iter})

def fit(method,s,y,order,band):
    if method=="fixed_seed_vf":
        w=np.geomspace(band[0]*1.5,band[1]/1.5,order//2)
        p=np.r_[(-.05+1j)*w,(-.05-1j)*w]
        return vector_fit(s,y,order,band,p)
    pilot=projected(s,y,order,band)
    if method=="projected_loewner":return pilot
    if method=="loewner_plus_vf":return vector_fit(s,y,order,band,pilot["poles"])
    residues=residue_fit(s,y,pilot["poles"])
    balanced=projected(s,y,order,band,weights=(pilot["poles"],residues))
    return vector_fit(s,y,order,band,balanced["poles"])

def evaluate(model,case):
    pred=model["predict"](case["st"]); truth=case["test_clean"]
    if not np.all(np.isfinite(pred)):raise FloatingPointError("nonfinite heldout prediction")
    nrmse=float(np.linalg.norm(pred-truth)/np.linalg.norm(truth))
    floor=.01*np.sqrt(np.mean(abs(truth)**2))
    local=float(np.sqrt(np.mean(abs((pred-truth)/np.maximum(abs(truth),floor))**2)))
    p=model["poles"]; q=case["poles"]
    ii,jj=linear_sum_assignment(abs(p[:,None]-q[None,:])/np.maximum(abs(q[None,:]),1e-20))
    pole=float(np.mean(abs(p[ii]-q[jj])/np.maximum(abs(q[jj]),1e-20)))
    return dict(heldout_nrmse=nrmse,heldout_local_relative_rms=local,
                mean_relative_pole_error=pole,unstable_poles=int(np.sum(p.real>0)))

def run(split,per_cell,seed):
    rows=[]; observations={}; started=time.time()
    for sc in SCENARIOS:
      for j in range(per_cell):
        case=make_case(sc,j,seed); cluster=f"{sc}-{j:03d}"
        for eta in NOISES:
          y=case["clean"]*(1+eta*case["noise"]); dh=digest(case["s"],y)
          caseid=f"{cluster}-eta{eta:g}"
          observations[caseid+"_noisy"]=y
          observations[caseid+"_train_clean"]=case["clean"]
          observations[caseid+"_test_clean"]=case["test_clean"]
          observations[caseid+"_poles"]=case["poles"]
          for method in METHODS:
            row=dict(split=split,case_id=caseid,cluster=cluster,scenario=sc,noise=eta,
                     method=method,order=case["order"],observation_sha256=dh,
                     heldout_truth_sha256=digest(case["test_clean"]),status="ok")
            t0=time.perf_counter()
            with warnings.catch_warnings(record=True) as ws:
              warnings.simplefilter("always")
              try:
                model=fit(method,case["s"].copy(),y.copy(),case["order"],BAND)
                row["fit_seconds"]=time.perf_counter()-t0
                row.update(evaluate(model,case));row["metadata"]=model["metadata"]
              except Exception as exc:
                row.update(status="failed",exception=f"{type(exc).__name__}: {exc}",
                           fit_seconds=time.perf_counter()-t0,heldout_nrmse=None)
              row["warnings"]=[f"{w.category.__name__}: {w.message}" for w in ws]
            rows.append(row)
        print(f"{split} {cluster}: {len(rows)} fits elapsed {time.time()-started:.1f}s",flush=True)
    # Same-data assertions inspect every single case, including failed fits.
    for cid in {r["case_id"] for r in rows}:
        rr=[r for r in rows if r["case_id"]==cid]
        assert len(rr)==len(METHODS) and len({r["observation_sha256"] for r in rr})==1
    (ROOT/f"{split}_raw.json").write_text(json.dumps(rows,indent=2))
    keys=["split","case_id","cluster","scenario","noise","method","order","observation_sha256",
          "heldout_truth_sha256","status","heldout_nrmse","heldout_local_relative_rms","mean_relative_pole_error",
          "unstable_poles","fit_seconds","exception","warnings","metadata"]
    with (ROOT/f"{split}_raw.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for row in rows:w.writerow({k:(json.dumps(v) if isinstance(v,(list,dict)) else v) for k,v in row.items()})
    observations["train_omega"]=case["s"].imag;observations["test_omega"]=case["st"].imag
    np.savez_compressed(ROOT/f"{split}_observations.npz",**observations)
    return rows

def summarize(rows):
    def score(r):return min(max(float(r["heldout_nrmse"]),1e-14),1e6) if r["heldout_nrmse"] is not None else 1e6
    cells=[]
    for scenario in ["ALL"]+SCENARIOS:
      for eta in ["ALL", "POSITIVE"]+NOISES:
        subset=[r for r in rows if (scenario=="ALL" or r["scenario"]==scenario) and
                (eta=="ALL" or (eta=="POSITIVE" and r["noise"]>0) or r["noise"]==eta)]
        for method in METHODS:
          rr=[r for r in subset if r["method"]==method];v=np.array([score(r) for r in rr])
          cells.append(dict(scenario=scenario,noise=eta,method=method,n=len(rr),failures=sum(r["status"]!="ok" for r in rr),
                     median_nrmse=float(np.median(v)),p90_nrmse=float(np.quantile(v,.9)),
                     geometric_mean_nrmse=float(10**np.mean(np.log10(v))),
                     mean_capped_nrmse=float(np.mean(v)),median_fit_seconds=float(np.median([r["fit_seconds"] for r in rr])),
                     warning_trials=sum(bool(r["warnings"]) for r in rr)))
    pairs=[];rng=np.random.default_rng(778198)
    for competitor in ["projected_loewner","loewner_plus_vf","fixed_seed_vf"]:
      for scenario in ["ALL"]+SCENARIOS:
        subset=[r for r in rows if r["noise"]>0 and (scenario=="ALL" or r["scenario"]==scenario)]
        balanced={r["case_id"]:r for r in subset if r["method"]=="balanced_plus_vf"}
        other={r["case_id"]:r for r in subset if r["method"]==competitor}
        cluster={}
        wins=0
        for cid,r in balanced.items():
            d=np.log10(score(other[cid]))-np.log10(score(r))
            cluster.setdefault(r["cluster"],[]).append(d);wins+=score(r)<score(other[cid])
        # Cluster bootstrap respects the same system being reused at noise levels.
        a=np.array([np.mean(v) for v in cluster.values()]);n=len(a)
        boot=np.mean(a[rng.integers(0,n,(5000,n))],axis=1)
        ci=np.quantile(boot,[.025,.975])
        pairs.append(dict(scenario=scenario,comparison=f"balanced_plus_vf vs {competitor}",
                     positive_noise_trials=len(balanced),system_clusters=n,
                     mean_log10_error_improvement=float(np.mean(a)),
                     geometric_error_ratio_balanced_over_comparator=float(10**(-np.mean(a))),
                     ratio_ci95=[float(10**(-ci[1])),float(10**(-ci[0]))],
                     paired_win_fraction=float(wins/len(balanced))))
    return dict(primary="Positive-noise clean heldout transfer NRMSE; ratio balanced_plus_vf / loewner_plus_vf; cluster bootstrap",
                failure_penalty=1e6,error_floor=1e-14,summary=cells,paired_comparisons=pairs,
                same_observation_hash_per_case_verified=True,versions=dict(numpy=np.__version__,scipy=scipy.__version__))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("split",choices=["development","confirmation"]);args=p.parse_args()
    protocol=json.loads((ROOT/"protocol.json").read_text())
    spec=protocol[args.split]
    if args.split=="confirmation":
        frozen=json.loads((ROOT/"freeze.json").read_text())
        assert frozen["code_sha256"]==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"code changed after freeze"
        assert frozen["protocol_sha256"]==hashlib.sha256((ROOT/"protocol.json").read_bytes()).hexdigest(),"protocol changed after freeze"
    rows=run(args.split,spec["systems_per_scenario"],spec["seed"])
    (ROOT/("results.json" if args.split=="confirmation" else "development_results.json")).write_text(json.dumps(summarize(rows),indent=2))
