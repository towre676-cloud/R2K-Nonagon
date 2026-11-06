"""Exploratory coincidence/symmetry checks of the finite-prime local model.
These are point controls, not a uniform collision-removal or prime-tail proof.
"""
from pathlib import Path
import json,time
from flint import arb,ctx
from ratios_jet import correlation
HERE=Path(__file__).resolve().parent

def main():
    ctx.prec=192;ctx.cap=12;t=time.time();P=97
    ell=(arb('9503636.19016305477494')/(2*arb.pi())).log();xs=[arb(x) for x in ['0','.37','1.13','2.41']]
    base,_,_=correlation(xs,ell,'arithmetic',P);rows=[]
    for label,ys in [('translation',[x+arb('.23') for x in xs]),('permutation',[xs[i] for i in [2,0,3,1]])]:
        v,_,_=correlation(ys,ell,'arithmetic',P);assert v.overlaps(base)
        rows.append({'control':label,'value':v.str(45,radius=True,more=True),'matches_base':True})
    collision=[]
    for s in ['.1','.01','.001']:
        eps=arb(s);ys=[arb(0),eps,arb('1.13'),arb('2.41')]
        v,_,_=correlation(ys,ell,'arithmetic',P)
        collision.append({'separation':s,'finite_prime_rho4':v.str(45,radius=True,more=True),
                          'rho4_over_separation_squared':(v/(eps*eps)).str(45,radius=True,more=True)})
        print('collision',s,v,flush=True)
    out={'status':'complete','scope':'exploratory finite-prime local-ratios point evidence; no uniform coincidence theorem',
         'prime_cutoff':P,'base_value':base.str(45,radius=True,more=True),'symmetry_controls':rows,'collision_checks':collision,'seconds':time.time()-t}
    (HERE/'collision_diagnostics.json').write_text(json.dumps(out,indent=2));print('COMPLETE',flush=True)

if __name__=='__main__':main()
