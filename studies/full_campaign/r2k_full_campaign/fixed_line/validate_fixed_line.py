"""Independent identities, quadrature controls, tail replay and acceptance."""
from pathlib import Path
import sys,json
from flint import arb,acb,arb_mat,ctx
HERE=Path(__file__).resolve().parent
import sine_covariance as m

def init():
    ctx.prec=256
    m.QX=[arb.legendre_p_root(m.NX,k,weight=True) for k in range(m.NX)]
    nw=m.nodes(arb(0),arb(m.L));m.XX=[x for x,w in nw];m.SW=[w.sqrt() for x,w in nw]
    K=arb_mat([[(m.XX[i]-m.XX[j]).sinc_pi()*m.SW[i]*m.SW[j] for j in range(m.NX)] for i in range(m.NX)])
    eye=arb_mat(m.NX,m.NX)
    for i in range(m.NX):eye[i,i]=1
    m.JJ=eye-2*K
    return K

def main():
    K=init();d=json.loads((HERE/'sine_fixed_line_results.json').read_text())
    exact=json.loads((HERE/'exact_controls.json').read_text());assert exact['status']=='PASS'
    assert d['status']=='complete' and (d['L'],d['span'],d['near_cutoff'])==(4,18000,64)
    # Recheck exact-contour arithmetic checkpoints rather than trusting a prefix.
    total=arb(0)
    for i,p in enumerate(d['panels']):
        assert arb(p['lo'])==arb(i)/2 and arb(p['hi'])==arb(i+1)/2
        total+=arb(p['weighted_integral_ball'])
    assert len(d['panels'])==128
    assert (2*total/m.D**2).overlaps(arb(d['near_influence_variance']['ball']))
    # At zero separation: ordinary fourth cumulant plus twice variance squared.
    K2=K*K;K3=K2*K;K4=K2*K2
    var=m.L-K2.trace();kap=m.L-7*K2.trace()+12*K3.trace()-6*K4.trace()
    assert m.overlap_F(arb(0)).overlaps(kap+2*var**2)
    assert m.overlap_F(arb(m.L)).overlaps(m.disjoint_F(arb(m.L)))
    # Analytic Fourier integral recurrence, independently adaptive-integrated.
    oldD=m.D;m.D=65
    integral_controls=[]
    for freq,power in [(0,2),(2,2),(2,8),(4,16),(4,32)]:
        z=m.harmonic_integrals(freq,power)[power];o=arb(freq)*arb.pi()
        numeric=acb.integral(lambda x,_:(acb(0,o)*x).exp()/x**power,64,65,abs_tol=arb('1e-65'),rel_tol=arb('1e-55'))
        assert numeric.overlaps(z)
        integral_controls.append({'harmonic_pi':freq,'power':power,'recurrence':str(z),'adaptive_integral':str(numeric)})
    m.D=oldD
    # Replay stored tail coefficients and compare pointwise with full kernels.
    coeff=d['tail_harmonic_coefficients'];checks=[]
    for h in [arb(64),arb(257)/4,arb(65),arb(96)]:
        poly=sum((acb(c['coefficient_real'],c['coefficient_imag'])*acb(0,arb(c['harmonic_pi'])*arb.pi()*h).exp()/h**c['power'] for c in coeff),acb(0))
        full=m.disjoint_F(h);c=arb(m.L)/(arb.pi()*(h-m.L));q=(arb(m.L)/h)**m.JORDER
        bound=(1+2*m.L)**2*c*c*(2*q+q*q)+4*c**4*((1+q)**4-1)
        assert abs(poly.real-full)<bound and poly.imag.contains(0)
        checks.append({'h':str(h),'full_kernel':str(full),'polynomial':str(poly),'absolute_error_bound':str(bound)})
    harmonic={o:m.harmonic_integrals(o,4*m.JORDER) for o in [-4,-2,0,2,4]}
    replay=sum((acb(c['coefficient_real'],c['coefficient_imag'])*(m.D*harmonic[c['harmonic_pi']][c['power']]-harmonic[c['harmonic_pi']][c['power']-1]) for c in coeff),acb(0))*2/m.D**2
    assert replay.real.overlaps(arb(d['far_influence_variance']['ball'])) and replay.imag.contains(0)
    W=arb(d['influence_variance_enclosure']['ball']);V=arb(d['anchor_variance_sampling_variance_enclosure']['ball'])
    assert W>0 and V>0 and V.rad()<arb('8e-10')
    Q=arb(d['centering_square_variance_upper']['ball']);B=arb(d['centering_bias']['ball'])
    assert Q.overlaps(2*B*B+(arb(m.L)/m.D)**2*B)
    circular=json.loads((HERE.parent/'fourth_order/cue_covariance_results.json').read_text())
    circle=arb(circular['models']['18000']['covariance'][3][3]['enclosure']['ball'])
    diff=V-circle
    out={'status':'PASS','exact_projection_count_checks':exact['exact_joint_count_covariance_checks'],
         'exact_projection_centering_checks':exact['exact_projection_block_identity_and_centering_bounds'],
         'near_panel_prefix_checks':len(d['panels']),'zero_and_touching_separation_controls':2,
         'adaptive_harmonic_integral_controls':integral_controls,'tail_full_kernel_controls':checks,
         'fixed_line_sampling_variance_enclosure':m.enc(V),'line_minus_CUE_variance_enclosure':m.enc(diff),
         'absolute_line_minus_CUE_variance_bound':str(abs(diff).upper()),
         'zeta_sampling_law_certified':False}
    (HERE/'validation_receipt.json').write_text(json.dumps(out,indent=2));print('PASS; fixed-line variance',V,'line minus CUE',diff)

if __name__=='__main__':main()
