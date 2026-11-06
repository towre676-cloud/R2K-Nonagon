from pathlib import Path
import sys,json,copy,hashlib,importlib.util,difflib
from flint import arb,acb,ctx
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[2]/'original_sources/unpacked/r2k_nonagon_v15.0/part_VI'
if (HERE/'certs').exists() and (HERE/'r2kverify.py').exists():SOURCE=HERE
# The original source can also be supplied explicitly for an extracted package.
if len(sys.argv)>1:SOURCE=Path(sys.argv[1]).resolve()
import r2kverify_hardened as h

def seal(c):
    c=copy.deepcopy(c);c.pop('sha256',None);c['sha256']=hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest();return c

def main():
    rows=[]
    for p in sorted((SOURCE/'certs').glob('*.json')):
        try:c=json.loads(p.read_text())
        except Exception:continue
        if not isinstance(c,dict) or not str(c.get('format','')).startswith('r2k-contour-certificate/'):continue
        v=h.verify(c,False);assert v['status']=='CERTIFIED',(p,v)
        rows.append({'file':p.name,'count':v['count'],'segments_checked':v['segments_checked'],'mode':'light stored evidence only'})
    c=json.loads((SOURCE/'certs/V_e10.json').read_text());bad=[]
    def test(label,mutate,full=False,evaluator='auto'):
        b=copy.deepcopy(c);mutate(b);v=h.verify(seal(b),full,evaluator)
        assert v['status']=='REJECTED',(label,v)
        bad.append({'case':label,'reason':v['reason']})
    test('half flag integer',lambda b:b.update(symmetric_half_path=1))
    test('empty moments',lambda b:b.update(moments_k=[]))
    test('boolean moment index',lambda b:b.update(moments_k=[False,1,2,3,4]))
    test('negative disk radius',lambda b:b['segments'][0].update(R2='-0x1.0p-4'))
    test('infinite coordinate',lambda b:b['rectangle'].update(lo='inf'))
    test('negative recorded radius',lambda b:b['segments'][0]['contrib'][0].__setitem__(1,-1))
    test('nonfinite recorded ball',lambda b:b['segments'][0]['contrib'][0].__setitem__(0,'nan'))
    test('bad precision',lambda b:b['method'].update(prec_bits=0))
    test('boolean precision',lambda b:b['method'].update(prec_bits=True))
    test('invalid Taylor order',lambda b:b['method'].update(K=0))
    test('DFT below requested coefficients',lambda b:b['method'].update(NDFT=b['method']['K']+1),True)
    test('unsupported evaluator',lambda b:None,False,'invalid')
    test('unsupported function',lambda b:b.update(function='unknown'))
    for label,cert,totals in [
        ('nonreal full-contour moments',{'symmetric_half_path':False},[2*acb.pi()*acb(0,1)*acb(1,1)]),
        ('negative count',{'symmetric_half_path':False},[-2*acb.pi()*acb(0,1)])]:
        try:h.decide(cert,totals,[0]);raise AssertionError(label+' accepted')
        except h.Reject as e:bad.append({'case':label,'reason':str(e)})
    v=h.verify([]);assert v['status']=='REJECTED';bad.append({'case':'nonobject certificate','reason':v['reason']})
    assert h.decide({'symmetric_half_path':False},[2*acb.pi()*acb(0,1)],[0])['count']==1
    assert h.decide({'symmetric_half_path':False},[acb(0)],[0])['count']==0
    original=(SOURCE/'r2kverify.py').read_text();patched=(HERE/'r2kverify_hardened.py').read_text()
    (HERE/'r2kverify_hardening.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='a/part_VI/r2kverify.py',tofile='b/part_VI/r2kverify.py')))
    out={'status':'PASS','original_verifier_sha256':hashlib.sha256(original.encode()).hexdigest(),
         'hardened_verifier_sha256':hashlib.sha256(patched.encode()).hexdigest(),
         'retained_light_acceptances':rows,'adversarial_rejections':bad,'positive_and_empty_count_controls':2,
         'scope':'acceptance hardening and retained-evidence replay; light mode is not a fresh analytic proof'}
    (HERE/'hardening_checks.json').write_text(json.dumps(out,indent=2));print('PASS',len(rows),'retained light acceptances;',len(bad),'adversarial rejections')

if __name__=='__main__':main()
