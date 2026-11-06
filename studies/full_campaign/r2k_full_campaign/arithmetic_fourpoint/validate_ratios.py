from pathlib import Path
import sys,json,math,os
from flint import arb,acb,acb_series,ctx
from ratios_jet import Jet,local_arithmetic,correlation,rmt_controls,prod
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
sys.path.insert(0,str(BASE))
from native_euler import Euler

def series_ratio(num,den,K,n):
    def polynomial(roots):
        p=[Jet(n,1)]
        for r in roots:
            q=[Jet(n) for _ in range(len(p)+1)]
            for i,c in enumerate(p):q[i]+=c;q[i+1]-=c*r
            p=q
        return p
    top=polynomial(num);bottom=polynomial(den);out=[]
    for k in range(K+1):
        v=top[k] if k<len(top) else Jet(n)
        for j in range(1,min(k,len(bottom)-1)+1):v-=bottom[j]*out[k-j]
        out.append(v)
    return out

def local_controls():
    n=4;K=256;eta=arb('.01');out=[]
    for p in [2,3,5]:
        A=[];B=[]
        for i,v in enumerate(['.13','.47','-.83','-1.29']):
            a=Jet(n,acb(0,arb(v)));a.a[1<<i]=acb(1)
            (A if i<2 else B).append(a)
        C=[Jet(n,a.a[0]) for a in A];D=[Jet(n,a.a[0]) for a in B]
        # Include a nontrivial swap so the arithmetic value is not identically one.
        A=[A[0],-B[1]];B=[B[0],-C[1]]
        # C is fixed; the second B variable is deliberately an independent jet
        # in this control to exercise every mixed component.
        B[1].a[1<<1]=acb(-1)
        lp=arb(p).log();root=lambda a:(-(a+arb(1)/2)*lp).exp()
        u=[root(a) for a in A];v=[root(a) for a in B];c=[root(a) for a in C];d=[root(a) for a in D]
        f=series_ratio(c,u,K,n);g=series_ratio(d,v,K,n);ct=sum((x*y for x,y in zip(f,g)),Jet(n))
        one=Jet(n,1);norm=prod((1-a*b for a in u for b in v),one)*prod((1-a*b for a in c for b in d),one)
        norm/=prod((1-a*b for a in u for b in d),one)*prod((1-a*b for a in c for b in v),one)
        approx=ct*norm;res=local_arithmetic(A,B,C,D,p)
        q=arb(p)**(-arb(1)/2+eta);k=K+1
        c0=math.comb(k+3,k)**2
        ratio=arb(k+4)**2/(k+1)**2
        tail=c0*q**(2*k)/(1-q*q*ratio)*((1+q*q)/(1-q*q))**8/eta**4
        pad=acb(arb(0,tail.upper()),arb(0,tail.upper()))
        assert all(a.overlaps(b+pad) for a,b in zip(res.a,approx.a))
        out.append({'prime':p,'series_degree':K,'all_16_mixed_coefficients_match':True,'Cauchy_series_tail_upper':str(tail)})
    return out

def main():
    ctx.prec=192;ctx.cap=12;local=local_controls();rmt=rmt_controls();checks=[]
    # Separate elementary prime loops avoid the optional accelerator's
    # deliberately single-initialization global state.
    os.environ['R2K_PURE_PYTHON']='1'
    ell=(arb('9503636.19016305477494')/(2*arb.pi())).log();x=arb('.37');r=2*arb.pi()*x/ell
    for P in [31,97,379,997,10007]:
        value,_,_=correlation([arb(0),x],ell,'arithmetic',P)
        E=Euler(P,192);A,B=E(r);z=acb_series([acb(1,r),1],prec=3).zeta();dl=2*z[2]/z[0]-(z[1]/z[0])**2
        direct=1+2/ell**2*(dl+acb(0,-ell*r).exp()*z[0]*z[0].conjugate()*A-B).real
        assert value.overlaps(direct)
        checks.append({'prime_cutoff':P,'compiler':str(value),'independent_BK_pair_evaluator':str(direct)})
    out={'status':'PASS','RMT_determinant_controls':rmt,'local_Euler_residue_vs_series_controls':local,
         'finite_prime_BK_pair_controls':checks,'scope':'validates assembly and finite Euler differentiation; does not certify an infinite-prime four-point tail or ratios-conjecture remainder'}
    (HERE/'ratios_validation.json').write_text(json.dumps(out,indent=2));print('PASS',len(rmt),'RMT determinants, 48 local mixed coefficients, five independent BK pair controls')

if __name__=='__main__':main()
