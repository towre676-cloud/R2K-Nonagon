from pathlib import Path
import json
from flint import arb,acb,ctx
from continuous_height import K,emom
ctx.prec=128
lo=arb('9.3');hi=arb('9.6');L=arb(4);D=hi.exp()*(hi-1)-lo.exp()*(lo-1);a=2*arb.pi()*L/hi;b=2*arb.pi()*L/lo
# Independent adaptive integration of the elementary averaged kernel.
f1=lambda r,_:K(acb(1),r,L,lo,hi)
f2=lambda r,_:K(acb(1),r,L,lo,2*arb.pi()*L/r)
I=acb.integral(f1,0,a,abs_tol=arb('1e-24'),rel_tol=arb('1e-24'))+acb.integral(f2,a,b,abs_tol=arb('1e-24'),rel_tol=arb('1e-24'))
computed=L+2/(arb.pi()*D)*I.real
truth=L+2*L*L/D*(hi.ei()-lo.ei())
assert computed.overlaps(truth)
rows=[]
for j in [0,1,2,5,20,40]:
 v=acb.integral(lambda x,_:x.exp()*x**j,lo,hi,abs_tol=arb('1e-24'),rel_tol=arb('1e-24')).real
 w=emom(j,lo,hi);assert v.overlaps(w)
 rows.append({'degree':j,'adaptive_integral_overlaps_elementary_primitive':True})
# Limit of the normalized height kernel as the height interval shrinks.
e=arb(10);eps=arb('1e-8');r=arb(1);Dd=(e+eps).exp()*(e+eps-1)-(e-eps).exp()*(e-eps-1)
k=K(acb(1),r,L,e-eps,e+eps).real/Dd
expected=(L-e*r/(2*arb.pi()))/e
assert (k-expected).abs_upper()<arb('1e-12')
out={'status':'PASS','F_equal_one_height_average':{'computed':computed.str(40,radius=True,more=True),'independent_Ei_identity':truth.str(40,radius=True,more=True)},'exponential_moments':rows,'narrow_height_limit':'PASS'}
Path(__file__).with_name('height_identity_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
