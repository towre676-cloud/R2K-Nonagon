"""Exact diagonal-congruence experiment and compensated far-endpoint floor.

K_t=T K T, T_vv=t, T_ii=1 otherwise. Baseline port data are invariant;
calibrated diagonal grounds break the congruence. This program verifies the
exact maximum possible response separation and physical network constraints.
"""
from pathlib import Path
from itertools import combinations
import json, math
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
RNG=np.random.default_rng(20311030)

def network(n,edges,grounds):
    K=np.diag(np.array(grounds,dtype=float))
    for i,j,w in edges:
        K[i,i]+=w; K[j,j]+=w; K[i,j]-=w; K[j,i]-=w
    return K

def response(K,D): return np.linalg.inv(K+np.diag(D))[0,0]

def random_checks():
    max_identity_error=0.; max_attained_error=0.; max_excess=-math.inf
    all_physical=True; finite_controls=0; rows=[]
    for case in range(20):
        n=3+case%6
        edges=[(i,i+1,float(RNG.uniform(.2,2))) for i in range(n-1)]
        for i,j in combinations(range(n),2):
            if j>i+1 and RNG.uniform()<.3: edges.append((i,j,float(RNG.uniform(.1,1.5))))
        grounds=RNG.uniform(.4,2,n); K=network(n,edges,grounds); G=np.linalg.inv(K)
        v=int(RNG.integers(1,n)); epsilon_max=min(grounds[i]/w for i,j,w in edges if j==v) if any(j==v for i,j,w in edges) else math.inf
        epsilon_max=min(epsilon_max,min((grounds[j]/w for i,j,w in edges if i==v),default=math.inf))
        t=1+min(.2,epsilon_max/2); T=np.eye(n); T[v,v]=t; Kt=T@K@T
        gp=Kt.sum(axis=1); all_physical=bool(all_physical and np.min(gp)>=-1e-14)
        S=G[0,v]**2/G[v,v]; envelope=S*(t-1)/(t+1)
        optimum=np.zeros(n); optimum[v]=t/G[v,v]
        attained=response(Kt,optimum)-response(K,optimum)
        max_attained_error=max(max_attained_error,abs(attained-envelope))
        for _ in range(30):
            D=np.zeros(n); D[1:]=10**RNG.uniform(-3,3,n-1)
            transformed=D.copy(); transformed[v]/=t*t
            shift=response(Kt,D)-response(K,D)
            identity=response(K,transformed)-response(K,D)
            max_identity_error=max(max_identity_error,abs(shift-identity))
            max_excess=max(max_excess,shift-envelope)
            finite_controls+=1
        rows.append(dict(n=n,vertex=v,t=t,epsilon_max=float(epsilon_max),
                         new_ground_min=float(np.min(gp)),S=S,envelope=envelope,attained=attained))
    assert max_identity_error<2e-14 and max_attained_error<2e-14 and max_excess<2e-14 and all_physical
    return dict(networks=20,finite_controls=finite_controls,all_physical=all_physical,
                max_congruence_identity_error=max_identity_error,max_attained_envelope_error=max_attained_error,
                max_random_envelope_excess=max_excess,rows=rows)

def exact_fixture():
    K=sp.Matrix([[4,-1,-2],[-1,4,-1],[-2,-1,5]])
    G=K.inv(); v=2; t=sp.Rational(6,5); T=sp.diag(1,1,t); Kt=T*K*T
    theta=t/G[v,v]; D=sp.diag(0,0,theta)
    attained=(Kt+D).inv()[0,0]-(K+D).inv()[0,0]
    S=G[0,v]**2/G[v,v]; envelope=sp.factor(S*(t-1)/(t+1))
    assert attained==envelope
    transformed=sp.diag(0,0,theta/t**2)
    assert (Kt+D).inv()[0,0]==(K+transformed).inv()[0,0]
    # Exact physical grounds and target edge transformation.
    old_grounds=K*sp.ones(3,1); new_grounds=Kt*sp.ones(3,1)
    assert all(z>=0 for z in new_grounds)
    assert -Kt[1,2]==t*(-K[1,2])
    return dict(K=str(K.tolist()),t=str(t),vertex=v,old_grounds=[str(z) for z in old_grounds],
                new_grounds=[str(z) for z in new_grounds],theta_optimal=str(theta),S=str(S),
                envelope=str(envelope),attained=str(attained),identity_exact=True)

def path_floor():
    rows=[]
    for n in (3,5,7,9,11):
        K=network(n,[(i,i+1,1.) for i in range(n-1)],[.5]*n); G=np.linalg.inv(K)
        u,v=n-2,n-1; w=1.; m=n*(n+1)//2; sigma=1e-12
        a=G[0,v]; r=G[v,v]; S=a*a/r; delta=2*sigma/(math.sqrt(m)*S)
        t=1+delta; T=np.eye(n); T[v,v]=t; Kt=T@K@T
        theta=np.zeros(n);theta[v]=t/r
        exact_sep=S*delta/(2+delta)
        tail=sigma*w*r/(math.sqrt(m)*a*a)
        assert np.min(Kt.sum(axis=1))>=0 and delta<=.5
        rows.append(dict(n=n,edge=[u,v],m=m,sigma=sigma,G0far=float(a),Gfarfar=float(r),wGfarfar=float(w*r),
                         congruence_delta=float(delta),admissible_delta_max=.5,new_ground_min=float(np.min(Kt.sum(axis=1))),
                         maximum_mean_separation=float(exact_sep),sigma_over_sqrt_m=float(sigma/math.sqrt(m)),
                         tail_error_threshold=float(tail),local_efficient_sd_along_gauge=float(2*tail),
                         gaussian_KL_upper_bound=float(m*exact_sep**2/(2*sigma**2))))
    return rows

def estimation_example():
    K=network(7,[(i,i+1,1.) for i in range(6)],[.5]*7); G=np.linalg.inv(K)
    v=6; w0=1.; F0=G[0,0]; S=G[0,v]**2/G[v,v]
    # The single intervention theta*Gvv=1 is locally Fisher-optimal for t=1.
    h0=F0-S/2; m=100;sigma=S*1e-4
    hbar=h0+sigma/math.sqrt(m)*RNG.standard_normal(80000)
    y=(F0-hbar)/S
    # Bounded extension makes the estimator defined with finite risk on every
    # Gaussian outcome, while preserving the local inverse around y=1/2.
    y_clipped=np.clip(y,.25,.75)
    clipped_count=int(np.count_nonzero(y_clipped!=y))
    t_est=np.sqrt(1/y_clipped-1); w_est=w0*t_est
    expected_sd=2*sigma*w0/(math.sqrt(m)*S)
    return dict(n=7,vertex=v,m=m,sigma=float(sigma),replicates=len(hbar),
                first_order_sd=float(expected_sd),empirical_sd=float(np.std(w_est,ddof=1)),
                ratio=float(np.std(w_est,ddof=1)/expected_sd),
                estimator='y=(F0-hbar)/S; t=sqrt(1/clip(y,1/4,3/4)-1); w=t*w0',
                clipped_count=clipped_count,clip_interval=[.25,.75],
                globally_defined=True,finite_risk_bounded_extension=True,
                estimator_range=[float(w0/math.sqrt(3)),float(w0*math.sqrt(3))],
                model_scope='one-parameter diagonal-congruence nuisance family')

def main():
    out=dict(random_checks=random_checks(),exact_fixture=exact_fixture(),path_floor=path_floor(),estimation=estimation_example())
    (HERE/'compensated_edge_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='path_floor'},indent=2)[:5000])
if __name__=='__main__':main()
