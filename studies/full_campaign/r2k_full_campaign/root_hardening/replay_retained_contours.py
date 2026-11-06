from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing,json,time,hashlib,sys
HERE=Path(__file__).resolve().parent

def work(path):
    import r2kverify_hardened as h
    p=Path(path);t=time.time();c=json.loads(p.read_text());r=h.verify(c,full=True)
    r['seconds']=time.time()-t;r['source_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    return p.name,r

def main():
    files=[]
    for p in sorted((HERE/'certs').glob('*.json')):
        data=json.loads(p.read_text())
        if isinstance(data,dict) and str(data.get('format','')).startswith('r2k-contour-certificate/') and p.stem not in ['V_e10','V_e12a','V_e12b']:files.append(p)
    assert len(files)==46
    out={};dest=HERE/'retained_full_replay.json'
    if '--resume' in sys.argv and dest.exists():
        old=json.loads(dest.read_text());assert old['verifier_sha256']==hashlib.sha256((HERE/'r2kverify_hardened.py').read_bytes()).hexdigest()
        out=old['results']
        for p in files:
            if p.name in out:assert out[p.name]['source_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        jobs=[pool.submit(work,str(p)) for p in files if p.name not in out]
        for future in as_completed(jobs):
            name,r=future.result();assert r['status']=='CERTIFIED',(name,r);out[name]=r
            tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps({'status':'complete' if len(out)==len(files) else 'running',
                'verifier_sha256':hashlib.sha256((HERE/'r2kverify_hardened.py').read_bytes()).hexdigest(),'results':out},indent=2));tmp.replace(dest)
            print(name,r['status'],r['count'],r.get('inertia'),'seconds',round(r['seconds'],1),flush=True)
    assert len(out)==46
    print('ALL 46 RETAINED CONTOURS FRESHLY REPLAYED',sum(r['segments_recomputed'] for r in out.values()),'segments',flush=True)

if __name__=='__main__':main()
