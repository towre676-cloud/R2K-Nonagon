#!/usr/bin/env python3
"""Replay fixed-catalog local D-optimal design and symmetry-opening tests.

All floating-point calculations use NumPy binary64; finite-difference checks
independently verify analytic physical-parameter response derivatives.
"""
from __future__ import annotations
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc


def network(n, edges, grounds):
    K = np.diag(np.asarray(grounds, dtype=float))
    columns = []
    names = []
    weights = []
    for u, v, w in edges:
        b = np.zeros(n)
        b[u], b[v] = 1., -1.
        K += w * np.outer(b,b)
        columns.append(b)
        weights.append(w)
        names.append(f"edge {u}-{v}")
    for j,w in enumerate(grounds):
        b=np.zeros(n); b[j]=1.
        columns.append(b);weights.append(w); names.append(f"ground {j}")
    return K,np.array(columns),np.array(weights),names


def catalog(n,t):
    configs=[("baseline",np.zeros(n))]
    for j in range(1,n):
        for fac in [1,2]:
            d=np.zeros(n);d[j]=fac*t
            configs.append((f"node {j}: {fac}t",d))
    for j,k in itertools.combinations(range(1,n),2):
        d=np.zeros(n);d[j]=t;d[k]=t
        configs.append((f"pair {j},{k}: t,t",d))
    return configs


def response(K,d):
    e=np.zeros(len(K));e[0]=1.
    return float(np.linalg.solve(K+np.diag(d),e)[0])


def logdet(J):
    sg,ld=np.linalg.slogdet(J)
    if sg != 1: raise ArithmeticError("nonpositive information")
    return float(ld)


def case(name,n,edges,grounds,budget=4):
    K,B,w,names=network(n,edges,grounds)
    F0=response(K,np.zeros(n))
    t=1./float(np.median(np.diag(np.linalg.inv(K))))
    configs=catalog(n,t)
    sigma=2e-4
    prior_sd=.2
    p=len(w)
    J0=np.eye(p)/prior_sd**2
    blocks=[];grads=[]; fd_errors=[]
    for label,d in configs:
        v=np.linalg.solve(K+np.diag(d),np.eye(n)[:,0])
        # eta is log weight: w_a(eta)=w_a*exp(eta_a).
        g=-w*(B@v)**2/F0
        grads.append(g)
        blocks.append(np.outer(g,g)/sigma**2)
        h=2e-5
        for a in range(p):
            Q=w[a]*np.outer(B[a],B[a])
            fd=(response(K+math.expm1(h)*Q,d)-response(K+math.expm1(-h)*Q,d))/(2*h*F0)
            fd_errors.append(abs(fd-g[a]))
    blocks=np.asarray(blocks);grads=np.asarray(grads)
    m=len(configs)
    J=J0.copy(); selected=[]; step_gains=[]
    for _ in range(budget):
        choices=[(logdet(J+blocks[i])-logdet(J),i) for i in range(m) if i not in selected]
        gain,i=max(choices,key=lambda z:(z[0],-z[1]))
        selected.append(i);step_gains.append(gain);J+=blocks[i]
    fg=logdet(J)-logdet(J0)
    all_sets=[]
    for comb in itertools.combinations(range(m),budget):
        val=logdet(J0+sum((blocks[i] for i in comb),start=np.zeros_like(J0)))-logdet(J0)
        all_sets.append((val,comb))
    best,opt=max(all_sets,key=lambda z:z[0])
    vals=np.array([a[0] for a in all_sets])
    # Check every elementary diminishing-return inequality Delta_i(S)>=Delta_i(S+j).
    infos=[]; lds=[]
    for mask in range(1<<m):
        A=J0.copy()
        for i in range(m):
            if mask>>i&1:A+=blocks[i]
        infos.append(A);lds.append(logdet(A))
    min_dr=math.inf; worst=None; dr_checks=0
    for mask in range(1<<m):
        absent=[i for i in range(m) if not mask>>i&1]
        for i,j in itertools.combinations(absent,2):
            # Symmetric mixed difference; both ordered choices give same gap.
            gap=lds[mask|1<<i]+lds[mask|1<<j]-lds[mask]-lds[mask|1<<i|1<<j]
            if gap<min_dr:min_dr=gap;worst=[mask,i,j]
            dr_checks+=1
    # Matrix determinant lemma gives all rank-one marginal gains.
    mdl_errors=[]
    trace_errors=[]
    for idx in range(m):
        g=grads[idx]/sigma
        direct=logdet(J0+blocks[idx])-logdet(J0)
        lemma=math.log1p(float(g@np.linalg.solve(J0,g)))
        mdl_errors.append(abs(direct-lemma))
        # d/d alpha logdet(J+alpha gg^T) = g^T (J+alpha gg^T)^-1 g.
        alpha=.37;hh=1e-4
        numeric=(logdet(J0+(alpha+hh)*blocks[idx])-logdet(J0+(alpha-hh)*blocks[idx]))/(2*hh)
        exact=float(g@np.linalg.solve(J0+alpha*blocks[idx],g))
        trace_errors.append(abs(numeric-exact))
    eig=np.linalg.eigvalsh(J)
    return {
        "name":name,"n":n,"parameters":names,"parameter_count":p,
        "K":K.tolist(),"response_baseline":F0,"common_strength_t":t,
        "coordinate_contract":"log physical weights w_a(eta)=w_a exp(eta_a), normalized response F/F0",
        "sigma_normalized_response":sigma,"prior_precision":25.0,
        "catalog_size":m,"budget":budget,
        "catalog":[{"label":c[0],"loads":c[1].tolist()} for c in configs],
        "selected_greedy":[configs[i][0] for i in selected],
        "selected_optimal":[configs[i][0] for i in opt],
        "greedy_step_logdet_gains":step_gains,
        "greedy_logdet_gain":fg,"optimal_logdet_gain":best,
        "greedy_D_efficiency_vs_optimal":math.exp((fg-best)/p),
        "uniform_subset_mean_logdet_gain":float(vals.mean()),
        "uniform_subset_p05_logdet_gain":float(np.quantile(vals,.05)),
        "uniform_subset_p95_logdet_gain":float(np.quantile(vals,.95)),
        "greedy_D_factor_vs_geometric_mean_uniform":math.exp((fg-float(vals.mean()))/p),
        "exhaustive_budget_subsets":len(all_sets),
        "submodularity_elementary_checks":dr_checks,
        "minimum_diminishing_return_gap":min_dr,"minimum_gap_case":worst,
        "max_gradient_abs_finite_difference_error":max(fd_errors),
        "max_matrix_determinant_lemma_abs_error":max(mdl_errors),
        "max_logdet_derivative_abs_finite_difference_error":max(trace_errors),
        "selected_information_eigenvalues":eig.tolist(),
        "selected_information_condition_number":float(eig[-1]/eig[0]),
        "selected_local_linear_covariance_trace":float(np.trace(np.linalg.inv(J))),
        "prior_covariance_trace":p*prior_sd**2,
        "gradient_matrix_rank":int(np.linalg.matrix_rank(grads)),
        "gradient_matrix_singular_values":np.linalg.svd(grads,compute_uv=False).tolist(),
    }


def symmetry_opening():
    # Diamond with equal paths 0-1-3 and 0-2-3. 1--2 is equipotential
    # at port 0 before a load; all shunt routing and calibration are explicit.
    n=4
    edges=[(0,1,1.),(0,2,1.),(1,3,1.),(2,3,1.),(1,2,.7)]
    grounds=[.2,.3,.3,.4]
    K,B,w,names=network(n,edges,grounds)
    edge_idx=names.index("edge 1-2")
    b=B[edge_idx]
    F0=response(K,np.zeros(n))
    records=[]
    for load in [0.,.1,.3,1.,3.,10.]:
        d=np.zeros(n);d[1]=load
        v=np.linalg.solve(K+np.diag(d),np.eye(n)[:,0])
        deriv=-float(b@v)**2
        h=1e-4
        fd=(response(K+h*np.outer(b,b),d)-response(K-h*np.outer(b,b),d))/(2*h)
        records.append({"load_node_1":load,"port_response":float(v[0]),"potential_difference_1_minus_2":float(b@v),"edge_derivative_physical":deriv,"edge_derivative_finite_difference":fd,"squared_derivative_Fisher_per_unit_variance":deriv**2})
    def rational_solve(A,b):
        Z=[list(map(Fraction,row))+[Fraction(rhs)] for row,rhs in zip(A,b)]
        for j in range(len(Z)):
            scale=Z[j][j];Z[j]=[x/scale for x in Z[j]]
            for k in range(len(Z)):
                if k==j:continue
                scale=Z[k][j];Z[k]=[x-scale*y for x,y in zip(Z[k],Z[j])]
        return [row[-1] for row in Z]
    A=[[Fraction('2.2'),-1,-1,0],[-1,3,Fraction('-.7'),-1],[-1,Fraction('-.7'),3,-1],[0,-1,-1,Fraction('2.4')]]
    exact=[]
    for ts in ['0','.1','.3','1','3','10']:
        t=Fraction(ts); AA=[row.copy() for row in A];AA[1][1]+=t
        v=rational_solve(AA,[1,0,0,0]);delta=v[1]-v[2]
        exact.append({'load':str(t),'response':str(v[0]),'potential_difference':str(delta),'physical_edge_derivative':str(-delta*delta)})
    return {"K":K.tolist(),"edge":"1-2","baseline_response":F0,"records":records,"exact_rational_records":exact,
            "principle":"a one-node ground breaks path-exchange symmetry and opens the exactly blind transverse edge"}


def bayesian_scalar_edge():
    """van Trees for log edge weight in an exactly soluble 3-node diamond.

    K(w)=[[3,-1,-1],[-1,2+w,-w],[-1,-w,2+w]].
    Only w=.7 exp(eta) is unknown, eta~N(0,.2^2). The unloaded
    response is exactly 1/2 for all eta, while load theta at node 2
    gives a rational response and opens scalar Fisher information.
    """
    sd=.2; sigma=2e-4; w0=.7; max_z=10.
    rows=[];forward_errors=[];gradient_errors=[]
    for w in np.geomspace(.01,100,21):
        K=np.array([[3.,-1.,-1.],[-1.,2+w,-w],[-1.,-w,2+w]])
        for theta in [0.,.1,.3,1.,3.,10.]:
            d=np.array([0.,0.,theta])
            direct=response(K,d)
            predicted=(4+4*w+theta*(2+w))/(8+8*w+theta*(5+3*w))
            forward_errors.append(abs(direct-predicted))
            v=np.linalg.solve(K+np.diag(d),np.array([1.,0.,0.]))
            direct_derivative=-2*w*(v[1]-v[2])**2
            predicted_derivative=-2*w*theta**2/(8+8*w+theta*(5+3*w))**2
            gradient_errors.append(abs(direct_derivative-predicted_derivative))
    if max(forward_errors)>1e-13 or max(gradient_errors)>1e-13:
        raise ArithmeticError('Bayesian family forward/gradient identity failed')
    for theta in [0.,.1,.3,1.,3.,10.,None]:
        def grad(z):
            w=w0*math.exp(sd*z)
            if theta is None:return -2*w/(5+3*w)**2
            return -2*w*theta**2/(8+8*w+theta*(5+3*w))**2
        def integrand(z):
            return grad(z)**2*math.exp(-z*z/2)/math.sqrt(2*math.pi)
        mean_g2,quad_error=quad(integrand,-max_z,max_z,epsabs=1e-16,epsrel=1e-12)
        bound=(1/30 if theta is None else theta**2/(2*(8+5*theta)*(8+3*theta)))
        tail_upper=bound**2*float(erfc(max_z/math.sqrt(2)))
        information=mean_g2/sigma**2
        rows.append({"load":"clamp" if theta is None else theta,
          "expected_squared_normalized_log_weight_gradient":mean_g2,
          "quadrature_reported_abs_error":quad_error,
          "omitted_tail_analytic_upper_bound":tail_upper,
          "uniform_gradient_absolute_bound":bound,
          "expected_data_Fisher_information_one_reading":information,
          "van_Trees_MSE_lower_bound_one_reading":1/(1/sd**2+information),
          "van_Trees_RMSE_lower_bound_one_reading":math.sqrt(1/(1/sd**2+information)),
          "van_Trees_MSE_lower_bound_four_repeats":1/(1/sd**2+4*information),
          "baseline_response_exact":"1/2"})
    return {"unknown":"eta=log(w/0.7)","prior":"Gaussian(0,0.2^2)","sigma_normalized_response":sigma,
      "model":"F_theta(w)=[4+4w+theta(2+w)]/[8+8w+theta(5+3w)]",
      "analytic_derivative":"partial_eta(F_theta/F0)=-2w theta^2/[8+8w+theta(5+3w)]^2",
      "quadrature":"SciPy adaptive Gauss-Kronrod, standard normal z on [-10,10], analytic uniform-gradient Gaussian-tail bound",
      "forward_identity_checks":len(forward_errors),"max_forward_identity_abs_error":max(forward_errors),
      "gradient_identity_checks":len(gradient_errors),"max_gradient_identity_abs_error":max(gradient_errors),
      "records":rows}


def main():
    cases=[
      case("four-node grounded path",4,[(0,1,1.),(1,2,.7),(2,3,1.4)],[.4,.2,.35,.6]),
      case("2x2 grounded grid",4,[(0,1,1.),(0,2,.8),(1,3,1.2),(2,3,.9)],[.3,.5,.2,.4]),
      case("asymmetric grounded diamond",4,[(0,1,1.),(0,2,.8),(1,3,.65),(2,3,1.3),(1,2,.45)],[.25,.4,.3,.55]),
    ]
    out={"schema":"nonagon.part-x.design.v1","arithmetic":"IEEE binary64, deterministic inputs; exact Fraction symmetry checks","cases":cases,"symmetry_opening":symmetry_opening(),"bayesian_scalar_edge":bayesian_scalar_edge()}
    required=[]
    for c in cases:
        required.extend([c["greedy_D_efficiency_vs_optimal"]<=1.+1e-10,
                         c["minimum_diminishing_return_gap"]>=-1e-8,
                         c["max_gradient_abs_finite_difference_error"]<1e-8,
                         c["max_matrix_determinant_lemma_abs_error"]<1e-9])
    if not all(required):raise ArithmeticError("design verification failed")
    out["required_checks"]={"count":len(required),"passed":int(sum(required))}
    path=Path(__file__).with_name("design_results.json")
    path.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"path":str(path),"checks":out["required_checks"],"cases":[{k:c[k] for k in ["name","greedy_logdet_gain","optimal_logdet_gain","greedy_D_efficiency_vs_optimal","minimum_diminishing_return_gap"]} for c in cases]},indent=2))

if __name__=="__main__":main()
