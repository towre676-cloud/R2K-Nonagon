"""Rigorous integrated A-tail bound by downward contour displacement in r.
The A-tail difference cancels the zeta product pole at zero. On the left
vertical edge F(-it) is real, so only t/alpha contributes to the real
weighted integral. On the right edge the weight is exactly it/alpha.
"""
from pathlib import Path
import json,sys
from flint import arb,acb,acb_series,ctx
import certified_campaign as cc
HERE=Path(__file__).resolve().parent;ctx.prec=128;ctx.cap=64
P=cc.P;small=cc.small

def txt(x):return x.str(45,radius=True,more=True)

def base_majorant(tlo,thi,rlo,rhi,S0,S2):
 z=acb(tlo.union(thi),rlo.union(rhi))
 Hp=1+z*acb_series([1+z,1],prec=1).zeta(deflate=True)[0]
 Hm=1-z*acb_series([1-z,1],prec=1).zeta(deflate=True)[0]
 AP=arb(1)
 for p in small:
  e=(-z*arb(p).log()).exp();AP*=(1-(e-1)**2/arb(p-1)**2).abs_upper()
 AP*=(arb(4)/999).exp()
 K=(1+arb(P)**(-tlo))**2*S0.upper()
 lower=tlo*tlo+rlo*rlo
 coeff=S2.upper() if not(lower>0) else S2.upper().min(K/lower)
 return (Hp.abs_upper()*Hm.abs_upper()*AP*coeff*(4*S0.upper()).exp()).upper()

def main():
 td=json.loads((HERE/'tail_constants.json').read_text());S0=arb(td['S0_ball']);S2=arb(td['S2_ball']);sigma=arb(1)/2
 d=json.loads(cc.RESULT.read_text());assert d['status']=='complete'
 # Exact dyadic partition of both vertical edges; no integration sampling.
 tcuts=sorted({arb(0).str(),sigma.str(),*[ (arb(k)/128).str() for k in range(1,65)],*[ (arb(1)/2**k).str() for k in range(8,17)]},key=lambda s:float(arb(s).mid()))
 tcuts=[arb(t) for t in tcuts]
 leftbase=[]
 for a,b in zip(tcuts,tcuts[1:]):leftbase.append((a,b,base_majorant(a,b,arb(0),arb(0),S0,S2)))
 # Bottom-edge partition inherited from the Gauss certificate, with origin.
 rpanels=[(arb(0),cc.r0)]+[(arb(p['lo']),arb(p['hi'])) for p in d['panels']]
 bottom=[(a,b,base_majorant(sigma,sigma,a,b,S0,S2)) for a,b in rpanels]
 results=[]
 for record in d['results']:
  s=next(s for s in cc.specs if s['n0']==record['n0'] and s['quantile']==record['quantile']);ell=s['ell'];alpha=s['alpha'];L=arb(record['L']);R=alpha*L
  left=sum((((b*b-a*a)/2)*(-ell*a).exp()*g for a,b,g in leftbase),arb(0))/alpha
  bot=arb(0)
  for a,b,g in bottom:
   if float(b.mid())>float(R.mid()):break
   bot+=(b-a)*(L-a/alpha+sigma/alpha)*(-ell*sigma).exp()*g
  right=arb(0)
  for a,b in zip(tcuts,tcuts[1:]):
   g=base_majorant(a,b,R,R,S0,S2)
   right+=(b*b-a*a)/2*(-ell*a).exp()*g/alpha
  Aerr=(4/ell**2/alpha*(left+bot+right)).upper();Berr=(S0/(2*arb.pi()**2)).upper()
  finite=arb(record['finite_product_variance_ball']);full=finite+arb(0,(Aerr+Berr).upper())
  results.append({**record,'A_tail_previous_real_path_bound':record['A_integrated_infinite_tail_bound'],
   'A_integrated_infinite_tail_bound':txt(Aerr),'A_contour_left_bound':txt((4/ell**2/alpha*left).upper()),
   'A_contour_bottom_bound':txt((4/ell**2/alpha*bot).upper()),'A_contour_right_bound':txt((4/ell**2/alpha*right).upper()),
   'infinite_product_variance_ball':txt(full),'infinite_product_radius':txt(full.rad()),'radius_below_1e_minus_6':bool(full.rad()<arb('1e-6'))})
  print('shifted tail',s['n0'],s['quantile'],float(L),txt(Aerr),flush=True)
 out={'status':'complete','sigma':txt(sigma),'prime_cutoff':P,'results':results,
 'all_target_radii_below_1e_minus_6':all(r['radius_below_1e_minus_6'] for r in results),
 'scope':'rigorous full infinite-product BK variance at declared quantile heights; excludes finite-height model discrepancy'}
 (HERE/(cc.TAG+'shifted_tail_results.json')).write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
