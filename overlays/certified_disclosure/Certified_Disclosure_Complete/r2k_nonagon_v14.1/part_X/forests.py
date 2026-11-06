"""Part X: positive forest observables and the exact edge-activation envelope.

Run with CODEX_PRIMARY_RUNTIME_PYTHON. Seeded binary64 checks are accompanied by
integer forest enumeration, exact rational identities, and 80-digit checks.
"""
from pathlib import Path
from itertools import combinations
import json, math
import numpy as np
import sympy as sp
import mpmath as mp

HERE=Path(__file__).resolve().parent
RNG=np.random.default_rng(20261030)

def network(n, edges, grounds, exact=False):
    K=sp.diag(*grounds) if exact else np.diag(np.array(grounds,dtype=float))
    for u,v,w in edges:
        K[u,u]+=w; K[v,v]+=w; K[u,v]-=w; K[v,u]-=w
    return K

def clamp_green(K, S):
    n=len(K); U=[i for i in range(n) if i not in S]
    H=np.zeros_like(K)
    if U: H[np.ix_(U,U)]=np.linalg.inv(K[np.ix_(U,U)])
    return H

def gamma(G,u,v):
    av=G[0,u]-G[0,v]*G[u,v]/G[v,v]
    au=G[0,v]-G[0,u]*G[u,v]/G[u,u]
    return max(av,au),av,au

def forest_components(n, edges, choose):
    parent=list(range(n)); w=1
    def find(i):
        while parent[i]!=i: i=parent[i]
        return i
    for j in choose:
        u,v,z=edges[j]; a,b=find(u),find(v)
        if a==b: return None
        parent[a]=b; w*=z
    return tuple(find(i) for i in range(n)),w

def forest_sums(n, edges, grounds):
    """Enumerate all spanning forests of the graph with added ground node n."""
    extended=edges+[(i,n,g) for i,g in enumerate(grounds) if g]
    records={c:[] for c in (1,2,3)}
    for c in records:
        for subset in combinations(range(len(extended)),n+1-c):
            f=forest_components(n+1,extended,subset)
            if f is not None: records[c].append(f)
    Z=sum(w for comp,w in records[1])
    Zij=[[sum(w for c,w in records[2] if c[i]==c[j] and c[i]!=c[n]) for j in range(n)] for i in range(n)]
    return Z,Zij,records

def relative(a,b):
    return float(np.max(np.abs(np.asarray(a,dtype=float)-np.asarray(b,dtype=float)))/max(1.,float(np.max(np.abs(np.asarray(b,dtype=float))))))

def exact_forest_test():
    n=5; edges=[(0,1,2),(0,2,3),(1,2,5),(1,3,7),(2,4,11),(3,4,13),(1,4,17)]; grounds=[1,2,3,5,7]
    K=network(n,edges,grounds,True); G=K.inv(); Z,Zij,rec=forest_sums(n,edges,grounds)
    assert Z==K.det()
    assert sp.Matrix(Zij)/Z==G
    rows=[]
    for u,v in combinations(range(1,n),2):
        denominator=sum(w for c,w in rec[2] if c[n]!=c[v])
        numerator=sum(w for c,w in rec[3] if c[n]!=c[v] and c[0]==c[u] and c[0] not in (c[n],c[v]))
        av=G[0,u]-G[0,v]*G[u,v]/G[v,v]
        assert av==sp.Rational(numerator,denominator)
        rows.append(dict(u=u,v=v,clamp_v_forest_numerator=numerator,clamp_v_forest_denominator=denominator,value=str(av)))
    negative=[]
    for I in combinations(range(n),2):
        for J in combinations(range(n),2):
            minor=G.extract(I,J).det()
            if minor<0: negative.append(dict(rows=list(I),cols=list(J),minor=str(minor)))
    assert negative
    return dict(nodes=n,edges=edges,grounds=grounds,tree_sum=Z,determinant=int(K.det()),two_forest_sums=Zij,
                forest_counts={str(c):len(r) for c,r in rec.items()},clamp_observables=rows,
                arbitrary_negative_minor=negative[0],all_inverse_entries_and_six_clamp_identities_exact=True)

def envelope_scan():
    """Compare the formula with all ideal-clamp subsets and random finite loads."""
    max_corner_error=0.; max_finite_violation=-math.inf; max_derivative_error=0.; cases=0; finite_count=0
    for n in range(3,10):
        for rep in range(8):
            edges=[(i,i+1,float(RNG.uniform(.2,2))) for i in range(n-1)]
            for i,j in combinations(range(n),2):
                if j>i+1 and RNG.uniform()<.35: edges.append((i,j,float(RNG.uniform(.05,2))))
            K=network(n,edges,RNG.uniform(.1,1.5,n)); G=np.linalg.inv(K)
            corner_H=[clamp_green(K,set(S)) for k in range(n) for S in combinations(range(1,n),k)]
            for u,v in combinations(range(1,n),2):
                gam,av,au=gamma(G,u,v)
                observed=max(abs(H[0,u]-H[0,v]) for H in corner_H)
                max_corner_error=max(max_corner_error,abs(observed-gam))
                assert gam>0
                for q in range(10):
                    loads=np.zeros(n); loads[1:]=10**RNG.uniform(-4,5,n-1)
                    H=np.linalg.inv(K+np.diag(loads)); contrast=H[0,u]-H[0,v]
                    max_finite_violation=max(max_finite_violation,abs(contrast)-gam)
                    delta=float(10**RNG.uniform(-5,1)); b=np.eye(n)[u]-np.eye(n)[v]
                    Hp=np.linalg.inv(K+np.diag(loads)+delta*np.outer(b,b))
                    drop=H[0,0]-Hp[0,0]; predicted=delta*contrast**2/(1+delta*(b@H@b))
                    max_derivative_error=max(max_derivative_error,abs(drop-predicted))
                    finite_count+=1
                cases+=1
    return dict(networks=56,edge_directions=cases,random_finite_controls=finite_count,max_all_corner_envelope_error=max_corner_error,
                max_finite_envelope_excess=max_finite_violation,max_rank_one_drop_absolute_error=max_derivative_error,
                all_checks_passed=bool(max_corner_error<2e-14 and max_finite_violation<2e-14 and max_derivative_error<2e-14))

def symmetric_example():
    w=sp.symbols('w',positive=True); K=sp.Matrix([[3,-1,-1],[-1,2+w,-w],[-1,-w,2+w]])
    G=sp.simplify(K.inv()); h=sp.simplify(G[0,0]); av=sp.simplify(G[0,1]-G[0,2]*G[1,2]/G[2,2])
    assert h==sp.Rational(1,2)
    assert sp.simplify(av-1/(5+3*w))==0
    rows=[]
    for ww in [.1,.3,1.,3.,10.]:
        A=np.array(K.subs(w,ww)).astype(float); H=np.linalg.inv(A); gam,*_=gamma(H,1,2)
        rows.append(dict(edge_weight=ww,baseline_response=float(H[0,0]),baseline_edge_contrast=float(H[0,1]-H[0,2]),
                         gamma=gam,clamped_response=float(clamp_green(A,{2})[0,0]),
                         ideal_fisher_per_reading_times_sigma2=gam**4))
    return dict(K='[[3,-1,-1],[-1,2+w,-w],[-1,-w,2+w]]',baseline='1/2',gamma=str(av),
                clamped_response=str(sp.simplify((K.extract([0,1],[0,1]).inv())[0,0])),rows=rows)

def weak_link_and_distance():
    weak=[]
    x=sp.symbols('x',positive=True)
    rational_K=network(4,[(0,1,x),(1,2,1),(1,3,1),(2,3,sp.Rational(7,10))],[1,1,sp.Rational(1,2),sp.Rational(1,2)],True)
    rational_G=rational_K.inv()
    rational_gamma=sp.factor(rational_G[0,2]-rational_G[0,3]*rational_G[2,3]/rational_G[3,3])
    assert sp.simplify(rational_gamma-5*x/(28+39*x))==0
    for w in np.logspace(-7,-1,7):
        # Port 0, hub 1, symmetric endpoints 2 and 3: baseline edge blindness.
        K=network(4,[(0,1,w),(1,2,1),(1,3,1),(2,3,.7)],[1,1,.5,.5]); G=np.linalg.inv(K)
        gam,av,au=gamma(G,2,3)
        weak.append(dict(w=float(w),G0u=float(G[0,2]),G0v=float(G[0,3]),baseline_contrast=float(G[0,2]-G[0,3]),
                         gamma=gam,gamma_over_w=gam/w,resolution_scale=gam**-2))
    distances=[]
    for depth in range(1,13):
        # A grounded stem ends in two equal deep endpoints; their mutual edge
        # is invisible at baseline at every depth.
        hub=depth; u=depth+1; v=depth+2; n=depth+3
        edges=[(i,i+1,1) for i in range(depth)]+[(hub,u,1),(hub,v,1),(u,v,.7)]
        K=network(n,edges,[.5]*n); G=np.linalg.inv(K); gam,av,au=gamma(G,u,v)
        exact_gamma=sp.Rational(30*2**depth,112*4**depth+23)
        assert abs(gam-float(exact_gamma))<2e-15
        distances.append(dict(depth=depth,n=n,G0u=float(G[0,u]),baseline_contrast=float(G[0,u]-G[0,v]),
                              gamma=gam,gamma_exact=str(exact_gamma),resolution_scale=gam**-2))
    slope=float(np.polyfit(np.log([r['w'] for r in weak[:5]]),np.log([r['resolution_scale'] for r in weak[:5]]),1)[0])
    return dict(weak_link=weak,distance=distances,weak_link_loglog_resolution_slope=slope,
                weak_link_gamma_exact=str(rational_gamma),distance_gamma_exact='30*2^d/(112*4^d+23)',
                distance_last_step_resolution_ratio=distances[-1]['resolution_scale']/distances[-2]['resolution_scale'])

def finite_strength_and_estimation():
    K=np.array([[3.,-1,-1],[-1,3,-1],[-1,-1,3]])
    G=np.linalg.inv(K); gam=gamma(G,1,2)[0]; rows=[]
    for beta in [.1,.3,1.,3.,9.,100.]:
        theta=beta/G[2,2]; H=np.linalg.inv(K+np.diag([0,0,theta])); q=H[0,1]-H[0,2]
        assert abs(q/gam-beta/(1+beta))<2e-14
        rows.append(dict(beta=beta,theta=theta,contrast=q,contrast_over_gamma=q/gam,
                         sd_inflation_vs_clamp=(gam/q)**2))
    # Nuisance-known local inversion is explicit for the clamped symmetric example.
    # h=(2+w)/(5+3w); w=(2-5h)/(3h-1).
    sigma=1e-6; m=100; w0=1.; h0=3/8
    readings=h0+(sigma/math.sqrt(m))*RNG.standard_normal(200000)
    clipped_readings=np.clip(readings,.35,.40)
    estimates=(2-5*clipped_readings)/(3*clipped_readings-1)
    empirical=float(np.std(estimates,ddof=1)); crlb=sigma/(math.sqrt(m)*gam**2)
    return dict(finite_grounding=rows,local_edge_estimation=dict(edge_weight=w0,sigma=sigma,readings_per_estimate=m,
           replicates=len(estimates),empirical_sd=empirical,first_order_optimal_sd=crlb,ratio=empirical/crlb,
           estimator='w=(2-5*clip(h,.35,.40))/(3*clip(h,.35,.40)-1)',clipped_count=int(np.count_nonzero(readings!=clipped_readings)),clip_response_interval=[.35,.40],mean_estimate=float(np.mean(estimates))))

def high_precision():
    mp.mp.dps=80; w=mp.mpf('0.0000001')
    K=mp.matrix([[1+w,-w,0,0],[-w,3+w,-1,-1],[0,-1,mp.mpf('2.2'),mp.mpf('-.7')],[0,-1,mp.mpf('-.7'),mp.mpf('2.2')]])
    G=K**-1; av=G[0,2]-G[0,3]*G[2,3]/G[3,3]
    U=[0,1,2]; R=mp.matrix([[K[i,j] for j in U] for i in U])**-1
    error=abs(av-R[0,2]); assert error<mp.mpf('1e-75')
    # Binary64 calculation compared with the same exact decimal network.
    Gn=np.linalg.inv(np.array(K.tolist(),dtype=float)); gn=gamma(Gn,2,3)[0]
    return dict(decimal_digits=80,w=str(w),gamma=mp.nstr(av,72),clamp_identity_absolute_error=mp.nstr(error,5),
                binary64_gamma=float(gn),binary64_relative_error=mp.nstr(abs(mp.mpf(gn)-av)/av,6))

def archived_family_prediction():
    """Exact forest scale plus analytic Jacobian predicts the archived raw-noise SD.

    Complex-step differentiates the complete inversion, preserving the one shared
    raw baseline. The archive supplies only its already generated Monte Carlo SD.
    """
    archive=HERE/'viii_stability_reference.json'
    old=json.loads(archive.read_text())
    def recover(data,n,ths):
        G=np.zeros((n,n),dtype=complex); G[0,0]=data[0]; pos=1
        for j in range(1,n):
            t1,t2=ths[j]; d1=data[0]-data[pos]; d2=data[0]-data[pos+1]; pos+=2
            squared=(1/t1-1/t2)/(1/d1-1/d2)
            G[j,j]=squared/d1-1/t1; G[0,j]=G[j,0]=np.sqrt(squared)
        for i,j in combinations(range(1,n),2):
            drop=data[0]-data[pos]; pos+=1
            a=1/ths[i][0]+G[i,i]; b=1/ths[j][0]+G[j,j]
            G[i,j]=G[j,i]=(G[0,i]*G[0,j]-np.sqrt((drop*a-G[0,i]**2)*(drop*b-G[0,j]**2)))/drop
        return np.linalg.inv(G)
    def prediction(K,sigma):
        G=np.linalg.inv(K); n=len(K); ths={j:(2**-.5/G[j,j],4/G[j,j]) for j in range(1,n)}
        interventions=[{}]
        for j in range(1,n): interventions.extend([{j:ths[j][0]},{j:ths[j][1]}])
        interventions.extend({i:ths[i][0],j:ths[j][0]} for i,j in combinations(range(1,n),2))
        h=[]
        for I in interventions:
            A=K.copy()
            for j,t in I.items(): A[j,j]+=t
            h.append(np.linalg.inv(A)[0,0])
        h=np.asarray(h,dtype=complex); derivative=[]
        for j in range(len(h)):
            pert=h.copy(); pert[j]+=1e-30j
            derivative.append(recover(pert,n,ths)[-1,-1].imag/1e-30)
        return sigma*np.linalg.norm(derivative)
    weak=[]; path=[]
    for row in old['weak_link']:
        w=row['w']; K=network(2,[(0,1,w)],[1,1]); pred=prediction(K,1e-9)
        weak.append(dict(w=w,forest_numerator=w,tree_sum=1+2*w,forest_G01=w/(1+2*w),
                         predicted_sd=pred,archived_sd=row['sd_K11_protocol'],archived_over_prediction=row['sd_K11_protocol']/pred))
    for row in old['distance']:
        n=row['n']; edges=[(i,i+1,1) for i in range(n-1)]; K=network(n,edges,[.5]*n)
        exact_K=network(n,edges,[sp.Rational(1,2)]*n,True); tree=exact_K.det(); pred=prediction(K,1e-12)
        path.append(dict(n=n,forest_numerator=1,tree_sum=str(tree),forest_G0far=str(1/tree),
                         predicted_sd=pred,archived_sd=row['sd_K_far'],archived_over_prediction=row['sd_K_far']/pred))
    return dict(weak_link=weak,path=path,prediction_method='forest ratio + full complex-step Jacobian with shared raw baseline')

def main():
    out=dict(exact_forests=exact_forest_test(),envelope_scan=envelope_scan(),symmetric_edge=symmetric_example(),
             hard_families=weak_link_and_distance(),finite_and_estimation=finite_strength_and_estimation(),high_precision=high_precision(),
             table_VIII_3_prediction=archived_family_prediction())
    (HERE/'forests_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k in ['envelope_scan','high_precision']},indent=2))
if __name__=='__main__': main()
