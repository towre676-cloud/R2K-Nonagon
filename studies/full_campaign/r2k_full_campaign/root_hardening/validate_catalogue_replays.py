from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
    vh=hashlib.sha256((HERE/'r2kverify_hardened.py').read_bytes()).hexdigest()
    low=json.loads((HERE/'retained_full_replay.json').read_text());assert low['status']=='complete' and low['verifier_sha256']==vh
    results=dict(low['results']);assert len(results)==46
    for name,fn in [('V_e10','V_e10_replay.json'),('V_e12a','V_e12a_parallel_replay.json'),('V_e12b','V_e12b_parallel_replay.json')]:
        r=json.loads((HERE/fn).read_text());results[name+'.json']=r
        if 'contract' in r:
            assert r['contract']['verifier_sha256']==vh
            assert r['contract']['certificate_sha256']==hashlib.sha256((HERE/'certs'/(name+'.json')).read_bytes()).hexdigest()
    count=segments=0;rows=[]
    for p in sorted((HERE/'certs').glob('*.json')):
        c=json.loads(p.read_text())
        if not isinstance(c,dict) or not str(c.get('format','')).startswith('r2k-contour-certificate/'):continue
        r=results[p.name];assert r['status']=='CERTIFIED' and r['mode']=='full'
        assert r['segments_recomputed']==r['segments_checked']==len(c['segments']) and r['count']>=0
        assert r['imag_parts_contain_zero']
        if 'source_sha256' in r:assert r['source_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
        count+=1;segments+=r['segments_recomputed'];rows.append({'file':p.name,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'count':r['count'],'inertia':r.get('inertia'),'segments':r['segments_recomputed']})
    assert count==49 and segments==5446
    extra=json.loads((HERE/'extra_fixture_replays.json').read_text());assert extra['status']=='complete' and len(extra['results'])==2
    for name,r in extra['results'].items():
        p=HERE/'extra_fixtures'/name;c=json.loads(p.read_text())
        assert r['status']=='CERTIFIED' and r['mode']=='full' and r['segments_recomputed']==r['segments_checked']==len(c['segments'])
        assert r['source_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
    out={'status':'PASS','full_certificates':count,'fresh_segments':segments,'separate_fixture_replays':extra['results'],
         'total_full_replays_including_fixtures':51,'total_fresh_segments_including_fixtures':6434,'verifier_sha256':vh,'catalogue':rows,
         'scope':'validates completed analytic receipts and source identity; run the full replay scripts for independent recomputation'}
    (HERE/'catalogue_validation.json').write_text(json.dumps(out,indent=2));print('PASS 49 full certificates / 5446 fresh analytic segments')

if __name__=='__main__':main()
