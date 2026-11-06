"""Closed-loop physical-frequency acquisition under a fixed observation budget.

The designer receives a noisy query callback, frequency pool, order and public
band. The callback has no route to truth, heldout frequencies or noise fraction.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, sys, time, warnings
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy import linalg as la
import reference_core as core

ROOT=Path(__file__).resolve().parent
BAND=(.03,8.); NPOOL=1025; NTEST=1024; BUDGET=64
SCENARIOS=["ordinary","colored_noise","weak_link","near_collision","model_mismatch"]
NOISES=[0.,.001,.01,.05]
POOL_W=np.geomspace(*BAND,NPOOL)
TEST_W=np.exp(np.log(BAND[0])+(np.arange(NTEST)+.37)/NTEST*np.log(BAND[1]/BAND[0]))
LOG_GRID=np.round(np.linspace(0,NPOOL-1,BUDGET)).astype(int)
PILOT=np.unique(np.r_[LOG_GRID[::2],LOG_GRID[-1]])
BASELINES=["log_uniform","linear_coverage","random","resonance_directed"]
ADAPTIVE=["adaptive_iopt","adaptive_disagreement"]

def system(scenario,index,seed):
    rng=np.random.default_rng(np.random.SeedSequence([seed,SCENARIOS.index(scenario),index]))
    n=int(rng.integers(3,6));m=rng.uniform(.7,1.5,n);k=rng.uniform(1.,4.,n);c=rng.uniform(.03,.16,n)
    if scenario=="weak_link":k[n//2]*=.035;c[n//2]*=.035
    if scenario=="near_collision":
        n=4;m0=rng.uniform(.8,1.3,2);k0=rng.uniform(1.,3.,2)
        m=np.r_[m0,m0*(1+rng.uniform(-.002,.002,2))]
        k=np.r_[k0,.002,k0[1]];c=np.array([.015,.015,.00003,.015])
    _,p,K,C=core.chain_response(1j*POOL_W[:1],m,k,c);M=np.diag(m)
    if scenario=="near_collision":K[2,2]+=k0[0];C[2,2]+=.015
    e=np.eye(n)[:,0]
    solve_audit=dict(evaluations=0,max_backward_residual=0.,max_numpy_relative_difference=0.)
    def response(w):
        values=[]
        for t in w:
            Q=(1j*t)**2*M+1j*t*C+K
            answer=la.solve(Q,e,assume_a="gen")
            relative_residual=la.norm(Q@answer-e)/(la.norm(Q)*la.norm(answer)+la.norm(e))
            if not np.isfinite(relative_residual) or relative_residual>1e-11:
                raise ValueError("physical solve residual check failed")
            independent=np.linalg.solve(Q,e)
            difference=la.norm(answer-independent)/max(la.norm(independent),1e-30)
            if not np.isfinite(difference) or difference>1e-9:
                raise ValueError("independent NumPy physical solve disagreement")
            solve_audit["evaluations"]+=1
            solve_audit["max_backward_residual"]=max(solve_audit["max_backward_residual"],float(relative_residual))
            solve_audit["max_numpy_relative_difference"]=max(solve_audit["max_numpy_relative_difference"],float(difference))
            values.append(answer[0])
        return np.array(values)
    clean=response(POOL_W);heldout=response(TEST_W)
    extra=None
    if scenario=="model_mismatch":
        om=rng.uniform(1.,5.);d=.025*om;extra=dict(omega=om,damping=d,numerator=.025)
        clean+=.025/((1j*POOL_W)**2+2*d*1j*POOL_W+om**2)
        heldout+=.025/((1j*TEST_W)**2+2*d*1j*TEST_W+om**2)
    z=(rng.standard_normal(NPOOL)+1j*rng.standard_normal(NPOOL))/np.sqrt(2)
    if scenario=="colored_noise":
        for j in range(1,NPOOL):z[j]=.95*z[j-1]+np.sqrt(1-.95**2)*z[j]
    return dict(clean=clean,heldout=heldout,noise=z,order=2*n,solve_audit=solve_audit,
                truth=dict(m=m.tolist(),k=k.tolist(),c=c.tolist(),K=K.tolist(),C=C.tolist(),extra=extra))

def sorted_fit(selected,values,order):
    idx=np.argsort(selected);s=1j*POOL_W[np.array(selected)[idx]];y=np.array(values)[idx]
    return core.fit("loewner_plus_vf",s,y,order,BAND),s,y

def coverage_score(selected):
    # Distance to observed frequencies in log-frequency coordinates.
    return np.min(abs(np.arange(NPOOL)[:,None]-np.array(selected)[None,:]),axis=1).astype(float)

def select_unique_targets(targets,selected,count):
    out=[];taken=set(selected)
    for w in targets:
        rank=np.argsort(abs(POOL_W-w))
        j=next(int(i) for i in rank if int(i) not in taken)
        out.append(j);taken.add(j)
        if len(out)==count:break
    return out

def information_batch(model,s,y,selected,count):
    """I-optimal linearized pole/residue design, with model-estimated noise.

    The covariance is a design surrogate, not a calibrated statistical interval.
    Every 4th acquisition is exploration outside this function.
    """
    p=model["poles"];r=core.residue_fit(s,y,p);sp=1j*POOL_W
    F=np.c_[1/(sp[:,None]-p[None,:]),r[None,:]/(sp[:,None]-p[None,:])**2]
    norms=np.maximum(np.sqrt(np.mean(abs(F)**2,axis=0)),1e-20);F/=norms
    pred=model["predict"](sp)
    noise_var=np.maximum(abs(pred),.02*np.sqrt(np.mean(abs(pred)**2)))**2
    noise_var/=np.median(noise_var)
    j=F[np.array(selected)]/np.sqrt(noise_var[np.array(selected)])[:,None]
    information=j.conj().T@j
    ridge=1e-7*max(float(np.trace(information).real/len(p)/2),1.)
    covariance=la.inv(information+ridge*np.eye(F.shape[1]))
    G=F.conj().T@F/NPOOL
    chosen=[];taken=set(selected)
    for _ in range(count):
        FC=F@covariance
        numerator=np.einsum('ij,jk,ik->i',FC,G,FC.conj()).real
        variance=np.einsum('ij,ij->i',FC,F.conj()).real
        score=np.maximum(numerator,0)/np.maximum(noise_var+variance,1e-30)
        score[list(taken)]=-np.inf
        ix=int(np.argmax(score));chosen.append(ix);taken.add(ix)
        v=covariance@F[ix].conj()
        den=noise_var[ix]+float((F[ix]@v).real)
        covariance-=np.outer(v,v.conj())/max(den,1e-30)
    return chosen

def acquire(strategy,query,order,random_seed):
    selected=PILOT.tolist();values=query(selected).tolist();design_warnings=[];rounds=0
    if strategy=="log_uniform":
        new=[int(i) for i in LOG_GRID if int(i) not in set(selected)]
        values+=query(new).tolist();selected+=new
    elif strategy=="linear_coverage":
        new=select_unique_targets(np.linspace(*BAND,BUDGET-len(selected)),selected,BUDGET-len(selected))
        values+=query(new).tolist();selected+=new
    elif strategy=="random":
        rng=np.random.default_rng(random_seed)
        available=np.array(sorted(set(range(NPOOL))-set(selected)))
        new=rng.choice(available,BUDGET-len(selected),replace=False).tolist()
        values+=query(new).tolist();selected+=new
    else:
      while len(selected)<BUDGET:
        count=min(4,BUDGET-len(selected));rounds+=1
        with warnings.catch_warnings(record=True) as ws:
          warnings.simplefilter("always")
          try:
            model,s,y=sorted_fit(selected,values,order)
            if strategy=="adaptive_iopt":new=information_batch(model,s,y,selected,max(count-1,0))
            else:
                yp=model["predict"](1j*POOL_W)
                if strategy=="adaptive_disagreement":
                    alt=core.fit("fixed_seed_vf",s,y,order,BAND)
                    other=alt["predict"](1j*POOL_W)
                    score=abs(yp-other)
                else:
                    # Resonance-directed exploitation weighted by distance to data.
                    score=abs(yp)
                new=[];taken=set(selected)
                for _ in range(max(count-1,0)):
                    dist=np.min(abs(np.arange(NPOOL)[:,None]-np.array(list(taken))[None,:]),axis=1)
                    # No duplicate queries; encourage resolving a region rather than
                    # concentrating an entire batch on adjacent pool points.
                    q=score*np.minimum(dist/4,1);q[list(taken)]=-np.inf
                    ix=int(np.argmax(q));new.append(ix);taken.add(ix)
            # One exploration frequency per batch, including all strategies.
            if len(new)<count:
                q=coverage_score(selected+new);q[selected+new]=-np.inf
                new.append(int(np.argmax(q)))
          except Exception as exc:
            design_warnings.append(f"round{rounds}: {type(exc).__name__}: {exc}; coverage fallback")
            new=[]
            for _ in range(count):
                q=coverage_score(selected+new);q[selected+new]=-np.inf
                new.append(int(np.argmax(q)))
          design_warnings.extend(f"round{rounds}: {w.category.__name__}: {w.message}" for w in ws)
        values+=query(new).tolist();selected+=new
    if len(selected)!=BUDGET or len(set(selected))!=BUDGET:raise ValueError("observation budget/uniqueness violation")
    model,s,y=sorted_fit(selected,values,order)
    return model,selected,np.array(values),dict(design_rounds=rounds,design_warnings=design_warnings)

def run(split,protocol):
    spec=protocol[split];methods=BASELINES+(ADAPTIVE if split=="development" else [protocol["selected_adaptive"]])
    rows=[];arrays={};started=time.time()
    for sc in SCENARIOS:
      for i in range(spec["systems_per_scenario"]):
        case=system(sc,i,spec["seed"]);cluster=f"{sc}-{i:03d}"
        arrays[cluster+"_candidate_clean"]=case["clean"];arrays[cluster+"_heldout_clean"]=case["heldout"]
        arrays[cluster+"_standardized_noise"]=case["noise"]
        for eta in NOISES:
          observed=case["clean"]*(1+eta*case["noise"])
          caseid=f"{cluster}-eta{eta:g}";poolhash=core.digest(POOL_W,observed)
          for method in methods:
            row=dict(split=split,case_id=caseid,cluster=cluster,scenario=sc,noise=eta,strategy=method,
                     order=case["order"],candidate_observation_sha256=poolhash,status="ok",physical_solver_audit=case["solve_audit"])
            requested=[]
            def query(indices):
                if any(j in requested for j in indices):raise ValueError("repeated noisy query")
                requested.extend(indices)
                if len(requested)>BUDGET:raise ValueError("observation budget exceeded")
                return observed[np.array(indices)].copy()
            t0=time.perf_counter()
            with warnings.catch_warnings(record=True) as ws:
              warnings.simplefilter("always")
              try:
                # The acquisition interface has no eta/scenario/truth access.
                model,selected,values,meta=acquire(method,query,case["order"],spec["seed"]+i*13+SCENARIOS.index(sc)*1001)
                row["fit_and_acquisition_seconds"]=time.perf_counter()-t0
                prediction=model["predict"](1j*TEST_W)
                if not np.all(np.isfinite(prediction)):raise FloatingPointError("nonfinite test prediction")
                row.update(heldout_nrmse=float(la.norm(prediction-case["heldout"])/la.norm(case["heldout"])),
                           selected_indices=selected,selected_observation_sha256=core.digest(np.array(selected),values),
                           observations=len(requested),metadata=meta)
              except Exception as exc:
                row.update(status="failed",exception=f"{type(exc).__name__}: {exc}",heldout_nrmse=None,
                           observations=len(requested),selected_indices=requested,fit_and_acquisition_seconds=time.perf_counter()-t0)
              row["warnings"]=[f"{w.category.__name__}: {w.message}" for w in ws]
            rows.append(row)
        print(f"{split} {cluster}: {len(rows)} fits, elapsed {time.time()-started:.1f}s",flush=True)
    for cid in {r["case_id"] for r in rows}:
        rr=[r for r in rows if r["case_id"]==cid]
        if len({r["candidate_observation_sha256"] for r in rr})!=1:raise ValueError("noise pool mismatch")
    (ROOT/f"{split}_raw.json").write_text(json.dumps(rows,indent=2))
    with (ROOT/f"{split}_raw.csv").open("w",newline="") as f:
        keys=sorted(set().union(*(r.keys() for r in rows)));w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})
    arrays["candidate_omega"]=POOL_W;arrays["heldout_omega"]=TEST_W
    np.savez_compressed(ROOT/f"{split}_observations.npz",**arrays)
    return rows

def summarize(rows,primary):
    methods=sorted({r["strategy"] for r in rows})
    def error(r):return min(max(float(r["heldout_nrmse"]),1e-14),1e6) if r["heldout_nrmse"] is not None else 1e6
    cells=[]
    for sc in ["ALL"]+SCENARIOS:
      for eta in ["ALL","POSITIVE"]+NOISES:
        for method in methods:
          rr=[r for r in rows if r["strategy"]==method and (sc=="ALL" or r["scenario"]==sc) and (eta=="ALL" or (eta=="POSITIVE" and r["noise"]>0) or eta==r["noise"])]
          e=np.array([error(r) for r in rr])
          cells.append(dict(scenario=sc,noise=eta,strategy=method,n=len(rr),failures=sum(r["status"]!="ok" for r in rr),
                            geometric_nrmse=float(10**np.mean(np.log10(e))),median_nrmse=float(np.median(e)),p90_nrmse=float(np.quantile(e,.9)),
                            mean_capped_nrmse=float(np.mean(e)),max_nrmse=float(np.max(e)),
                            median_seconds=float(np.median([r["fit_and_acquisition_seconds"] for r in rr])),
                            warning_trials=sum(bool(r["warnings"]) or bool(r.get("metadata",{}).get("design_warnings",[])) for r in rr)))
    pairs=[];rng=np.random.default_rng(771931)
    for strategy in [m for m in methods if m!="log_uniform"]:
      for sc in ["ALL"]+SCENARIOS:
        rr=[r for r in rows if r["noise"]>0 and (sc=="ALL" or r["scenario"]==sc)]
        a={r["case_id"]:r for r in rr if r["strategy"]==strategy};b={r["case_id"]:r for r in rr if r["strategy"]=="log_uniform"}
        clusters={};wins=0
        for cid,r in a.items():
            clusters.setdefault(r["cluster"],[]).append(np.log10(error(r)/error(b[cid])));wins+=error(r)<error(b[cid])
        v=np.array([np.mean(x) for x in clusters.values()]);n=len(v)
        boot=np.mean(v[rng.integers(0,n,(5000,n))],axis=1);ci=np.quantile(boot,[.025,.975])
        pairs.append(dict(scenario=sc,strategy=strategy,comparator="log_uniform",geometric_error_ratio=float(10**np.mean(v)),
                          ratio_ci95=[float(10**ci[0]),float(10**ci[1])],win_fraction=wins/len(a),trials=len(a),clusters=n))
    return dict(primary=f"{primary} vs log_uniform; positive-noise clean heldout NRMSE; clustered paired error ratio",summary=cells,
                comparisons=pairs,shared_candidate_noise_verified=True,query_budget=BUDGET,common_pilot=len(PILOT),
                versions=dict(numpy=np.__version__,scipy=core.scipy.__version__))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("split",choices=["development","confirmation"]);args=p.parse_args()
    protocol=json.loads((ROOT/"protocol.json").read_text())
    if args.split=="confirmation":
        freeze=json.loads((ROOT/"freeze.json").read_text())
        for filename,wanted in freeze["sha256"].items():
            actual=hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()
            if actual!=wanted:raise ValueError("frozen file changed: "+filename)
    rows=run(args.split,protocol)
    result=summarize(rows,protocol.get("selected_adaptive","not_yet_selected"))
    (ROOT/("results.json" if args.split=="confirmation" else "development_results.json")).write_text(json.dumps(result,indent=2))
