"""Certified scalar covariance for the actual fixed-line anchor geometry.

Conditional unit-intensity sine projection DPP only. No zeta sampling-law claim.
All finite evaluation, spatial/h quadrature and long-distance tails are enclosed.
"""
from pathlib import Path
import sys,json,time,math
from collections import defaultdict
from flint import arb,acb,arb_mat,acb_mat,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'finite_window'))
from sine_bias import bias
L=4;S=18000;D=S-L;H=64;NX=48;NH=20;JORDER=8

def txt(x):return x.str(55,radius=True,more=True)
def enc(x):return {'ball':txt(x),'mid':float(x.mid()),'radius':float(x.rad())}
def integer_power(x,n):
    out=arb(1)
    for _ in range(n):out*=x
    return out
def trace_product(A,B):
    n=A.nrows()
    return sum((A[i,j]*B[j,i] for i in range(n) for j in range(n)),arb(0))

def nodes(lo,hi):
    c=(lo+hi)/2;half=(hi-lo)/2
    return [(c+half*r,half*w) for r,w in QX]

def overlap_F(h):
    # Partition the union into left-only, intersection, right-only pieces.
    # h strictly inside (0,L) in production Gaussian panels.
    pieces=[(arb(0),h,1,0),(h,arb(L),1,1),(arb(L),h+L,0,1)]
    data=[]
    for lo,hi,a,b in pieces:
        for x,w in nodes(lo,hi):data.append((x,w,a,b,a*b))
    n=len(data);sw=[q[1].sqrt() for q in data]
    K=arb_mat([[(data[i][0]-data[j][0]).sinc_pi()*sw[i]*sw[j] for j in range(n)] for i in range(n)])
    X=arb_mat([[K[i,j]*data[i][2] for j in range(n)] for i in range(n)])
    Y=arb_mat([[K[i,j]*data[i][3] for j in range(n)] for i in range(n)])
    Z=arb_mat([[K[i,j]*data[i][4] for j in range(n)] for i in range(n)])
    XX=X*X;YY=Y*Y;XY=X*Y
    cov=Z.trace()-XY.trace()
    kap=(Z.trace()-XY.trace()-2*trace_product(Z,Z)-2*trace_product(Z,X)-2*trace_product(Z,Y)
         +2*trace_product(X,YY)+2*trace_product(Y,XX)+8*trace_product(Z,XY)
         -4*trace_product(XX,YY)-2*trace_product(XY,XY))
    return kap+2*cov**2

def disjoint_F(h):
    C=arb_mat([[(XX[j]-XX[i]+h).sinc_pi()*SW[i]*SW[j] for j in range(NX)] for i in range(NX)])
    Z=JJ*C*JJ;T=C*C.transpose()
    f2=-sum((Z[i,j]*C[i,j] for i in range(NX) for j in range(NX)),arb(0))
    return f2+2*(T.trace()**2-trace_product(T,T))

def harmonic_integrals(omega,maxp):
    """J_p = integral_H^D exp(i omega h) h^-p dh, p>=1."""
    if omega==0:
        return {p:acb((arb(D)/H).log() if p==1 else (arb(H)**(1-p)-arb(D)**(1-p))/(p-1)) for p in range(1,maxp+1)}
    if omega<0:return {p:v.conjugate() for p,v in harmonic_integrals(-omega,maxp).items()}
    o=arb(omega)*arb.pi();u=o*H;v=o*D
    out={1:acb(v.ci()-u.ci(),v.si()-u.si())}
    eD=acb(0,v).exp();eH=acb(0,u).exp()
    for p in range(2,maxp+1):out[p]=(acb(0,o)*out[p-1]-(eD/arb(D)**(p-1)-eH/arb(H)**(p-1)))/(p-1)
    return out

def tail_polynomial():
    # C(h)=sum h^-r [E_(r,+) exp(i*pi*h)+E_(r,-) exp(-i*pi*h)].
    E={};pi=arb.pi();I=acb(0,1);n=NX
    for r in range(1,JORDER+1):
        mat=acb_mat([[integer_power(-(XX[j]-XX[i]),r-1)*SW[i]*SW[j]/(2*pi)*
                      acb((XX[j]-XX[i]).sin_pi(),-(XX[j]-XX[i]).cos_pi()) for j in range(n)] for i in range(n)])
        assert all(z.is_finite() for z in mat.entries()),('E',r)
        E[r,1]=mat;E[r,-1]=mat.conjugate()
    keys=list(E);Z={k:acb_mat(JJ)*E[k]*acb_mat(JJ) for k in keys}
    # Native matrix products compute all elementwise bilinear traces at once.
    ez=acb_mat([[Z[k][i,j] for i in range(n) for j in range(n)] for k in keys])
    ee=acb_mat([[E[k][i,j] for k in keys] for i in range(n) for j in range(n)])
    gram=ez*ee;assert all(z.is_finite() for z in gram.entries()),'gram2'
    coeff=defaultdict(acb);Ds={}
    for ii,(r,a) in enumerate(keys):
        for jj,(s,b) in enumerate(keys):
            coeff[r+s,a+b]-=gram[ii,jj]
            key=(r+s,a+b);prod=E[r,a]*E[s,b].transpose()
            Ds[key]=Ds[key]+prod if key in Ds else prod
    dk=list(Ds);tr=[Ds[k].trace() for k in dk]
    left=acb_mat([[Ds[k][i,j] for i in range(n) for j in range(n)] for k in dk])
    right=acb_mat([[Ds[k][j,i] for k in dk] for i in range(n) for j in range(n)])
    gram=left*right;assert all(z.is_finite() for z in gram.entries()),'gram4'
    for i,(p,a) in enumerate(dk):
        for j,(q,b) in enumerate(dk):coeff[p+q,a+b]+=2*(tr[i]*tr[j]-gram[i,j])
    assert all(z.is_finite() for z in coeff.values()),'coefficients'
    integrals={o:harmonic_integrals(o,4*JORDER) for o in [-4,-2,0,2,4]}
    assert all(z.is_finite() for vals in integrals.values() for z in vals.values()),'integrals'
    total=acb(0)
    for (p,o),c in coeff.items():total+=c*(D*integrals[o][p]-integrals[o][p-1])
    assert total.imag.contains(0)
    # Pointwise geometric remainder and Schatten norm inequalities.
    HH=arb(H);ll=arb(L);den=1-ll/HH;j=JORDER
    e2=(1+2*ll)**2*ll**2/(pi**2*den**2)*(2*ll**j/((j+1)*HH**(j+1))+ll**(2*j)/((2*j+1)*HH**(2*j+1)))
    e4=4*ll**4/(pi**4*den**4)*sum((math.comb(4,r)*ll**(j*r)/((3+j*r)*HH**(3+j*r)) for r in range(1,5)),arb(0))
    bound=2*(e2+e4)/D
    rows=[{'power':p,'harmonic_pi':o,'coefficient_real':txt(c.real),'coefficient_imag':txt(c.imag)} for (p,o),c in sorted(coeff.items())]
    return 2*total.real/D**2,bound,rows

def main():
    global QX,XX,SW,JJ
    ctx.prec=256;t=time.time()
    QX=[arb.legendre_p_root(NX,k,weight=True) for k in range(NX)]
    qh=[arb.legendre_p_root(NH,k,weight=True) for k in range(NH)]
    nw=nodes(arb(0),arb(L));XX=[x for x,w in nw];SW=[w.sqrt() for x,w in nw]
    K=arb_mat([[(XX[i]-XX[j]).sinc_pi()*SW[i]*SW[j] for j in range(NX)] for i in range(NX)])
    eye=arb_mat(NX,NX)
    for i in range(NX):eye[i,i]=1
    JJ=eye-2*K
    prefix=arb(0);records=[];checkpoint=HERE/'near_checkpoint.json'
    if '--resume' in sys.argv and checkpoint.exists():
        old=json.loads(checkpoint.read_text())
        assert old['precision_bits']==ctx.prec and old['NX']==NX and old['NH']==NH and old['H']==H
        assert (old['L'],old['S'])==(L,S)
        records=old['panels'];prefix=arb(old['prefix'])
    for a in range(len(records),2*H):
        lo=arb(a)/2;hi=lo+arb(1)/2;c=(lo+hi)/2;half=(hi-lo)/2;panel=arb(0)
        for x,w in qh:
            h=c+half*x;f=overlap_F(h) if a<2*L else disjoint_F(h)
            panel+=half*w*(D-h)*f
        prefix+=panel
        records.append({'lo':txt(lo),'hi':txt(hi),'weighted_integral_ball':txt(panel)})
        if a%16==15:print('h through',float(hi),round(time.time()-t,1),flush=True)
        checkpoint.write_text(json.dumps({'precision_bits':ctx.prec,'NX':NX,'NH':NH,'H':H,'L':L,'S':S,'prefix':txt(prefix),'panels':records}))
    near=2*prefix/D**2
    far,tail_error,coeff=tail_polynomial();print('tail complete',round(time.time()-t,1),flush=True)
    # Tensor Gauss error: one complex spatial variable occurs in <=2 kernel
    # edges; total domains have length <=L. All other weights remain positive.
    spatial=4*(2*arb.pi()*L).exp()*arb(2)**(-2*NX)*sum(arb(c*r)*arb(L)**r for r,c in [(1,1),(2,9),(3,16),(4,8)])
    # h-analyticity of the *quadrature polynomial*, using absolute complex
    # weights <=L+2rho on each moving A/B interval and entire sinc kernels.
    rho=arb(1);U=L+2*rho;e=(arb.pi()*rho).exp()
    MF=U*e+7*U**2*e**2+12*U**3*e**3+6*U**4*e**4+2*(U*e+U**2*e**2)**2
    q=arb(1)/4
    h_error=2/arb(D)**2*(2*H)*4*(arb(1)/4)*(D+H+rho)*MF*q**(2*NH)/(1-q)
    W=near+far+arb(0,(spatial+h_error+tail_error).upper())
    print('BUDGET',txt(near),txt(far),txt(spatial),txt(h_error),txt(tail_error),flush=True)
    assert W>0
    B,Bup=bias(S,L);m=arb(L)/D
    Qvar=2*B**2+m**2*B
    sdQ=Qvar.sqrt();sdW=W.sqrt()
    sdV=(sdW-sdQ).union(sdW+sdQ);assert sdV>0
    varV=sdV**2
    x=2*arb.pi()*L
    pop=L-2*L/arb.pi()*x.si()+(1-x.cos()+arb.const_euler()+x.log()-x.ci())/arb.pi()**2
    out={'status':'complete','scope':'stationary unit-intensity sine projection DPP on the line; not a zeta law',
         'span':S,'L':L,'anchor_length':D,'precision_bits':ctx.prec,'spatial_gauss_order':NX,
         'separation_gauss_order':NH,'near_cutoff':H,'tail_geometric_order':JORDER,
         'near_influence_variance':enc(near),'far_influence_variance':enc(far),
         'spatial_tensor_error_bound':enc(spatial),'separation_quadrature_error_bound':enc(h_error),
         'long_distance_polynomial_error_bound':enc(tail_error),'influence_variance_enclosure':enc(W),
         'centering_bias':enc(B),'centering_square_variance_upper':enc(Qvar),
         'centering_sd_perturbation_upper':enc(sdQ),'anchor_variance_sampling_variance_enclosure':enc(varV),
         'anchor_variance_sampling_sd_enclosure':enc(sdV),'population_number_variance':enc(pop),
         'expected_empirical_anchor_variance':enc(pop-B),'panels':records,
         'tail_harmonic_coefficients':coeff,'elapsed_seconds':time.time()-t}
    (HERE/'sine_fixed_line_results.json').write_text(json.dumps(out,indent=2))
    print('COMPLETE',enc(varV),enc(sdV),round(time.time()-t,1),flush=True)

if __name__=='__main__':main()
