from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing,json,hashlib,time
HERE=Path(__file__).resolve().parent

def work(path):
    import r2kverify_hardened as h
    p=Path(path);t=time.time();r=h.verify(json.loads(p.read_text()),full=True)
    r.update(seconds=time.time()-t,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest());return p.name,r

def main():
    out={}
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        jobs=[pool.submit(work,str(p)) for p in sorted((HERE/'extra_fixtures').glob('*.json'))]
        for f in as_completed(jobs):
            name,r=f.result();assert r['status']=='CERTIFIED',(name,r);out[name]=r
            (HERE/'extra_fixture_replays.json').write_text(json.dumps({'status':'complete' if len(out)==2 else 'running','results':out},indent=2))
            print(name,'CERTIFIED',r['count'],r.get('inertia'),r['segments_recomputed'],'segments',round(r['seconds'],1),'seconds',flush=True)
    assert len(out)==2 and out['ten_zero_contour_fixture.json']['count']==10 and out['union_7005_contour.json']['count']==5

if __name__=='__main__':main()
