from pathlib import Path
import json,sys,math
import numpy as np
from flint import arb,acb,arb_mat,ctx
import full_matrix as f
HERE=Path(__file__).resolve().parent

def positive_certificate(A):
    mid=np.array([[float(x.mid()) for x in row] for row in A]);P=np.linalg.inv(np.linalg.cholesky(mid))
    P=[[arb(float(x)) for x in row] for row in P];det=arb_mat(P).det();assert abs(det)>0
    B=[[sum((P[i][k]*A[k][l]*P[j][l] for k in range(6) for l in range(6)),arb(0)) for j in range(6)] for i in range(6)]
    margins=[B[i][i]-sum((abs(B[i][j]) for j in range(6) if j!=i),arb(0)) for i in range(6)]
    assert all(x>0 for x in margins)
    return {'method':'exact-dyadic congruence and interval Gershgorin','preconditioner_hex':[[float(x.mid()).hex() for x in row] for row in P],
            'determinant_ball':f.txt(det),'positive_margins':[f.txt(x) for x in margins]}

def main():
    ctx.prec=192;d=json.loads((HERE/'full_matrix_results.json').read_text());assert d['status']=='complete'
    W=[[arb(c['ball']) for c in row] for row in d['influence_covariance']];V=[[arb(c['ball']) for c in row] for row in d['anchor_covariance']]
    checks=0;maxradius=arb(0)
    for i,l in enumerate(f.LS):
        for j,m in enumerate(f.LS[i:],i):
            p=json.loads((HERE/'matrix_cases'/('L%s_M%s.json'%(l,m))).read_text());assert p['status']=='complete'
            assert p['contract'][:3]==[l,m,f.S]
            prefix=sum((arb(r['weighted_integral']) for r in p['panels']),arb(0));near=2*prefix/((f.S-l)*(f.S-m))
            assert near.overlaps(arb(p['near']['ball']))
            err=sum((arb(p[k]['ball']) for k in ['spatial_error','separation_error','tail_error']),arb(0))
            replay=near+arb(p['far']['ball'])+arb(0,err.upper())
            assert replay.overlaps(W[i][j]) and W[i][j].overlaps(W[j][i])
            maxradius=maxradius.max(V[i][j].rad());checks+=1
    native=[]
    for l,m in [(.5,.5),(.5,16.),(2.,8.),(4.,4.),(8.,16.),(16.,16.)]:
        C=f.Counts(l,m)
        for h in [arb(0),arb('.3'),abs(arb(l)-m)/2+arb('.1'),(arb(l)+m)/2+arb('.1'),arb(64)]:
            a=C(h);b=C.python(h);assert a.overlaps(b),(l,m,h,a,b)
            native.append({'L':l,'M':m,'h':str(h),'native_ball':str(a),'independent_python_ball':str(b)})
    old=json.loads((HERE/'sine_fixed_line_results.json').read_text())
    assert W[3][3].overlaps(arb(old['influence_variance_enclosure']['ball']))
    # Gershgorin padding supplies an explicit Loewner upper covariance.
    M=[[arb(float(V[i][j].mid())) for j in range(6)] for i in range(6)]
    radius=[]
    for i in range(6):
        need=sum(((V[i][j].rad()+abs(V[i][j].mid()-M[i][j])).upper() for j in range(6)),arb(0)).upper()
        pad=arb(math.nextafter(float(need),math.inf));assert pad>=need;radius.append(pad)
    U=[[M[i][j]+(radius[i] if i==j else arb(0)) for j in range(6)] for i in range(6)]
    assert all(x.is_exact() for row in U for x in row)
    posW=positive_certificate(W);posV=positive_certificate(V);posU=positive_certificate(U)
    out={'status':'PASS','completed_distinct_entries':checks,'native_python_controls':native,
         'scalar_L4_independent_production_overlap':True,'maximum_anchor_entry_radius':f.txt(maxradius),
         'influence_positive_certificate':posW,'anchor_positive_certificate':posV,
         'Loewner_upper_covariance':[[f.enc(x) for x in row] for row in U],
         'Loewner_upper_covariance_exact_rational':[[str(x.fmpq()) for x in row] for row in U],
         'Loewner_upper_positive_certificate':posU,
         'simultaneous_model_bound':'P((V-EV)^T U^-1(V-EV)>=6/alpha)<=alpha, by Markov; no Gaussian assumption',
         'confidence_95_percent_quadratic_threshold':120,'scope':'conditional sine-DPP theorem for a fixed declared observation; no zeta or post-selection coverage claim'}
    (HERE/'full_matrix_validation.json').write_text(json.dumps(out,indent=2))
    print('PASS: 21 entries, 30 native/Python controls, all covariance positivity certificates; max radius',maxradius)

if __name__=='__main__':main()
