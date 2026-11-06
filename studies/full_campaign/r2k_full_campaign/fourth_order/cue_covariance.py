"""Exact CUE covariance identity, Arb finite Fourier sums, analytic tails.

This certifies an explicitly circular fixed-rank model, not a zeta null law.
Run from any directory; requires python-flint and numpy.
"""
from pathlib import Path
import json,time
from fractions import Fraction as Q
from flint import arb, arb_poly, arb_mat, ctx
import numpy as np

HERE=Path(__file__).resolve().parent
LS=[Q(1,2),Q(1),Q(2),Q(4),Q(8),Q(16)]
S=18000
K=180000
NS=[17999,18000,18001]

def enc(x):
    return {'ball':x.str(40,radius=True,more=True),'mid':float(x.mid()),'radius':float(x.rad())}

def periodic_variance(d,L):
    # Exact dyadic arithmetic, including artificial seam pairs.
    xx=[Q(float.fromhex(x)) for x in d['points_hex']]
    n=len(xx); total=Q(0);pairs=0
    for i,x in enumerate(xx):
        j=i+1
        while j<n and xx[j]-x<L:
            total+=L-(xx[j]-x);pairs+=1;j+=1
    for i,x in enumerate(xx):
        j=0
        while j<i and S-x+xx[j]<L:
            total+=L-(S-x+xx[j]);pairs+=1;j+=1
    mu=Q(n)*L/S
    v=mu+2*total/S-mu*mu
    return {'mean_exact':str(mu),'variance_exact':str(v),'variance':float(v),
            'unordered_close_pairs':pairs,
            'scope':'exact for stored dyadic points; circular seam introduced'}

def main():
    ctx.prec=128;t=time.time();pi=arb.pi()
    w=[];prefix=[]
    for L in LS:
        l=arb(L.numerator)/L.denominator
        a=[arb(0)]+[((arb(j)*l/S).sin_pi()/(pi*j))**2 for j in range(1,K+1)]
        w.append(a);p=[arb(0)]
        for j in range(1,max(NS)):p.append(p[-1]+a[j])
        prefix.append(p)
        print('weights',str(L),round(time.time()-t,2),flush=True)
    results={str(N):{'mean':[],'covariance':[[None]*6 for _ in LS]} for N in NS}
    for N in NS:
        for a,L in zip(w,LS):
            l=arb(L.numerator)/L.denominator/S
            # Exact infinite Fourier total, with min(j,N) split at N.
            mean=2*sum((a[j]*(j-N) for j in range(1,N)),arb(0))+N*l*(1-l)
            results[str(N)]['mean'].append(enc(mean))
    for i in range(6):
        for j in range(i,6):
            # product u(z) v(z^-1): lag d at coefficient K+d.
            corr=arb_poly(w[i])*arb_poly(list(reversed(w[j])))
            hank=arb_poly(w[i][:max(NS)])*arb_poly(w[j][:max(NS)])
            for N in NS:
                diag=sum((w[i][m]*w[j][m]*min(m,N)**2 for m in range(1,K+1)),arb(0))
                mx=2*sum((prefix[i][m]*prefix[j][m] for m in range(1,N)),arb(0))
                hs=sum((hank[m]*(N-m) for m in range(2,N)),arb(0))
                ds=N*corr[K]+sum(((corr[K+m]+corr[K-m])*(N-m) for m in range(1,N)),arb(0))
                finite=4*(diag+mx-hs-ds)
                # |kappa|<=N, supported on |j-k|<N. Bound omitted
                # rows/columns with weights <=1/(pi*j)^2; integral test.
                tail=4*arb(N)**2/(3*pi**4*K**3)*(1+4/(1-arb(N)/K)**2)
                whole=finite+arb(0,tail.upper())
                cell={'finite':enc(finite),'tail_bound_ball':enc(tail),
                      'tail_bound':float(tail.upper()),'enclosure':enc(whole),
                      'Gaussian_pairing_finite':enc(4*diag),
                      'fourth_cumulant_finite':enc(4*(mx-hs-ds))}
                results[str(N)]['covariance'][i][j]=cell
                results[str(N)]['covariance'][j][i]=cell
            print('covariance',i,j,round(time.time()-t,2),flush=True)
    empirical=[]
    for file in sorted((HERE.parent/'finite_window').glob('fixed_data_*.json')):
        d=json.loads(file.read_text());N=d['point_count'];model=results[str(N)]
        rows=[]
        for idx,L in enumerate(LS):
            p=periodic_variance(d,L);mean=model['mean'][idx]['mid'];var=model['covariance'][idx][idx]['enclosure']['mid']
            p.update({'L':float(L),'CUE_mean':mean,'residual':p['variance']-mean,'CUE_sd':var**.5,
                      'standardized_residual':(p['variance']-mean)/var**.5})
            rows.append(p)
        empirical.append({'n0':d['n0'],'rank':N,'rows':rows})
    out={'status':'complete','scope':'conditional CUE_N circle model; not finite-line sine DPP or a zeta theorem',
         'span':S,'cutoff':K,'precision_bits':ctx.prec,'lengths':[float(x) for x in LS],
         'models':results,'empirical_periodized':empirical,'elapsed_seconds':time.time()-t}
    for N in NS:
        mat=np.array([[c['enclosure']['mid'] for c in row] for row in results[str(N)]['covariance']])
        err=max(c['enclosure']['radius'] for row in results[str(N)]['covariance'] for c in row)
        eig=np.linalg.eigvalsh(mat)
        assert eig[0]>6*err
        results[str(N)]['positive_definite_margin_lower']=float(eig[0]-6*err)
        # This last eigenvalue margin uses floating eigensolver: a numerical
        # diagnostic only, not an interval spectral certificate.
        results[str(N)]['spectral_margin_scope']='floating eigensolver diagnostic'
        # Rigorous positivity by exact-dyadic preconditioning and Gershgorin.
        # Numerical Cholesky proposes a basis; Arb proves its invertibility
        # and positivity of the congruence for every matrix in the enclosure.
        A=[[arb(c['enclosure']['ball']) for c in row] for row in results[str(N)]['covariance']]
        prop=np.linalg.inv(np.linalg.cholesky(mat))
        P=[[arb(float(x)) for x in row] for row in prop]
        det=arb_mat(P).det();assert abs(det)>0
        B=[[sum((P[i][k]*A[k][l]*P[j][l] for k in range(6) for l in range(6)),arb(0)) for j in range(6)] for i in range(6)]
        margins=[B[i][i]-sum((abs(B[i][j]) for j in range(6) if i!=j),arb(0)) for i in range(6)]
        assert all(x>0 for x in margins)
        results[str(N)]['positivity_certificate']={'method':'dyadic congruence and interval Gershgorin',
                'preconditioner_hex':[[float(x).hex() for x in row] for row in prop],
                'determinant':enc(det),'Gershgorin_margins':[enc(x) for x in margins]}
    (HERE/'cue_covariance_results.json').write_text(json.dumps(out,indent=2))
    print('COMPLETE',round(time.time()-t,2),flush=True)

if __name__=='__main__':main()
