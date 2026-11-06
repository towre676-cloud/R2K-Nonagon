"""Infinite-product BK witnesses excluding a scalar Hermitian DPP completion.
The obstruction concerns the declared analytic model, not actual zeta correlations.
"""
from pathlib import Path
import sys,json
from flint import arb,acb,acb_series,ctx
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
sys.path.insert(0,str(BASE))
from native_euler import Euler

def txt(x):return x.str(45,radius=True,more=True)
def main():
    ctx.prec=192;ctx.cap=8;P=2000000;E=Euler(P,192)
    tails=json.loads((BASE/'tail_constants.json').read_text());s0=arb(tails['S0_ball']);s2=arb(tails['S2_ball'])
    rows=[]
    for n0,xs in [(100000,'0.95'),(1000000,'0.97'),(10000000,'0.98'),(20000000,'0.98'),(50000000,'1.92')]:
        data=json.loads((BASE/'data'/('zeros_%s.json'%n0)).read_text());T=arb(data['ordinate_balls'][10000]);x=arb(xs)
        ell=(T/(2*arb.pi())).log();r=2*arb.pi()*x/ell;A,B=E(r)
        z=acb_series([acb(1,r),1],prec=3).zeta();dl=2*z[2]/z[0]-(z[1]/z[0])**2
        finite=1+2/ell**2*(dl+acb(0,-ell*r).exp()*z[0]*z[0].conjugate()*A-B).real
        # For every omitted prime, |(p^-ir-1)^2/(p-1)^2| <=4/(p-1)^2.
        logA=4*s0/(1-arb(4)/P**2)
        err=2/ell**2*(z[0].abs_upper()**2*A.abs_upper()*(logA.exp()-1)+s2)
        full=finite+arb(0,err.upper());assert full>1
        rows.append({'n0':n0,'T_ball':txt(T),'x':xs,'finite_product_R2':txt(finite),
                     'infinite_product_tail_bound':txt(err),'infinite_product_R2':txt(full),'R2_minus_one':txt(full-1)})
    out={'status':'PASS','prime_cutoff':P,'precision_bits':192,'witnesses':rows,
         'theorem':'A scalar Hermitian determinantal process with intensity one has rho2(x,y)=1-|K(x,y)|^2<=1; none can realize this BK R2 at the certified witnesses.',
         'scope':'excludes exact scalar Hermitian DPP completion of the declared local stationary BK model; does not exclude general processes, mixtures, non-Hermitian representations, or approximate asymptotic kernels; not a zeta theorem'}
    (HERE/'dpp_obstruction.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
