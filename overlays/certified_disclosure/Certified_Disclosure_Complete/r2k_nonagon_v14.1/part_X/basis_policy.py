#!/usr/bin/env python3
"""Certified basis policy: structured route only when signs decide.
Default: replay retained contour moment enclosures. --regenerate-ten regenerates
(and fully verifies) the common 10-zero fixture before comparison.
"""
import sys,json,hashlib,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.cert_deps'))
EXPERIMENTS=ROOT/'part_X' if (ROOT/'part_X/certification.py').exists() else ROOT/'experiments'
sys.path.insert(0,str(EXPERIMENTS))
from flint import arb,arb_mat,ctx
from certification import moments,power_y,condensation,sign,sub
V9=ROOT if (ROOT/'part_VI').exists() else ROOT/'source/v9/r2k_nonagon_v9.0'
CLAUDE=ROOT/'source/claude_v10/r2k_nonagon_v10.0'
if not CLAUDE.exists():CLAUDE=V9
sys.path.insert(0,str(CLAUDE/'part_VI'))
import certgen,r2kverify
HERE=Path(__file__).resolve().parent
FIXTURE=HERE/'ten_zero_contour_fixture.json'

def prepare_fixture(regenerate=False):
 if regenerate or not FIXTURE.exists():
  c=certgen.generate('xi',10.,51.,.45,list(range(22)),True,mode='series',prec=80)
  FIXTURE.write_text(json.dumps(c,indent=1))
  v=r2kverify.verify(c,full=True,evaluator='series')
  assert v['status']=='CERTIFIED' and v['count']==10
  meta={'full_verifier':v,'fixture_sha256':hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
   'producer_sha256':hashlib.sha256((CLAUDE/'part_VI/certgen.py').read_bytes()).hexdigest(),
   'tmint_sha256':hashlib.sha256((CLAUDE/'part_VI/tmint.py').read_bytes()).hexdigest(),
   'independent_verifier_sha256':hashlib.sha256((CLAUDE/'part_VI/r2kverify.py').read_bytes()).hexdigest(),
   'common_input':'generator_moments, unchanged Claude10[10,51] test parameters'}
  (FIXTURE.parent/'ten_zero_fixture_provenance.json').write_text(json.dumps(meta,indent=2))

def decide(p,ten=False):
 c=json.load(open(p));ctx.prec=c['method']['prec_bits'] if ten else max(128,c['method']['prec_bits'])
 r2kverify.integrity(c)
 mu=[arb(x[0],x[1]) for x in c['generator_moments']] if ten else moments(c)
 n=int(round(float(mu[0].mid())))
 if n<1 or len(mu)<2*n-1:return None
 H=arb_mat([[(mu[i+j]+mu[abs(i-j)])/2 for j in range(n)] for i in range(n)])
 flat=[sub(H,range(k),range(k)).det() for k in range(1,n+1)]
 fs=[sign(v) for v in flat];assert all(fs)
 hd=[];failure=None;ops=0
 try:
  qd,ops=condensation(power_y(mu,2*n-2),n)
  hd=[v*arb(2)**(2*(k-1)**2) for k,v in enumerate(qd,1)]
 except AssertionError as e:failure=str(e)
 hs=[sign(v) for v in hd]
 eligible=len(hd)==n and all(hs) and hs==fs and all(a.overlaps(b) for a,b in zip(hd,flat))
 chosen=hd if eligible else flat
 ss=[1]+[sign(v) for v in chosen];neg=sum(a!=b for a,b in zip(ss,ss[1:]))
 rr={'file':p.name,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'zeros':n,
  'selected_method':'shifted_hankel' if eligible else 'direct_chebyshev_fallback','selected_inertia':[n-neg,neg],
  'flat_minor_signs':fs,'hankel_minor_signs':hs,'hankel_updates':ops,'hankel_failure':failure,
  'all_selected_signs_decide':all(sign(v) for v in chosen),'last_direct_minor':flat[-1].str(25,radius=True),
  'last_hankel_minor':hd[-1].str(25,radius=True) if len(hd)==n else None,
  'last_hankel_over_direct_radius':float(hd[-1].rad()/flat[-1].rad()) if len(hd)==n else None}
 return rr

def main():
 start=time.perf_counter();prepare_fixture('--regenerate-ten' in sys.argv)
 rows=[]
 for p in sorted((V9/'part_VI/certs').glob('*.json')):
  c=json.load(open(p))
  if not isinstance(c,dict) or c.get('format') not in ('r2k-contour-certificate/1','r2k-contour-certificate/2'):continue
  r=decide(p)
  if r:rows.append(r)
 assert len(rows)==49
 rows.append(decide(FIXTURE,ten=True))
 assert rows[-1]['selected_method']=='direct_chebyshev_fallback'
 assert rows[-1]['selected_inertia']==[10,0]
 result={'arithmetic':'Arb ball arithmetic; archived comparisons>=128bits, identical10zeroinput80bits',
  'policy':'use shifted Hankel only if every leading-minor sign decides, matches the direct reference, and balls overlap; otherwise choose direct Chebyshev',
  'scope':'comparison of retained moment enclosures;10zerofixture independently full verified on creation',
  'records':rows,'structured_selected':sum(x['selected_method']=='shifted_hankel' for x in rows),
  'fallback_selected':sum(x['selected_method']=='direct_chebyshev_fallback' for x in rows),
  'all50_inertias_decide':all(x['all_selected_signs_decide'] for x in rows),'elapsed_seconds':time.perf_counter()-start}
 assert result['structured_selected']==49 and result['fallback_selected']==1
 (HERE/'basis_policy_results.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()
