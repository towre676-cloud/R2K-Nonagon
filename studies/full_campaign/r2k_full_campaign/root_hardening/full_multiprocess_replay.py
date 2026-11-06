"""Fresh analytic replay; parallelizes independent segments, never proof agents.
Every enclosure comes from r2kverify_hardened, with independent contour checks.
Atomic segment checkpoints permit a complete replay to survive interruption.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing,json,time,sys,hashlib
from flint import arb,acb,ctx
import r2kverify_hardened as h
HERE=Path(__file__).resolve().parent

def worker(job):
    ctx.prec=job['prec'];h.FUNCTION=job['function'];h.LOGPI=None
    s=job['segment'];A=acb(*s['A']);D=acb(*s['D']);cc=acb(h.fx(job['center']));a=acb(h.fx(job['halfwidth']))
    t=time.time()
    try:
        mine=h.segment_enclosure(A,D,h.fx(s['x0']),h.fx(s['x1']),h.fx(s['R2']),job['ks'],cc,a,job['K'],'dft',job['NDFT'])
        rec=[h.ball_from(c) for c in s['contrib']]
        h.need(all(r.overlaps(m) for r,m in zip(rec,mine)),'recorded ball disagrees with independent enclosure')
        return {'status':'PASS','id':job['id'],'name':job['name'],'seconds':time.time()-t,
                'enclosures':[[z.real.str(40,radius=True,more=True),z.imag.str(40,radius=True,more=True)] for z in mine]}
    except h.Reject as e:return {'status':'REJECTED','id':job['id'],'name':job['name'],'reason':str(e)}

def save(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2));tmp.replace(path)

def main():
    paths=[Path(p) for p in sys.argv[1:] if not p.startswith('--')] or [HERE/'certs/V_e12a.json',HERE/'certs/V_e12b.json']
    specs={};jobs=[]
    vh=hashlib.sha256((HERE/'r2kverify_hardened.py').read_bytes()).hexdigest()
    for p in paths:
        c=json.loads(p.read_text());name=p.stem;light=h.verify(c,False);assert light['status']=='CERTIFIED'
        assert c['method']['NDFT']>=c['method']['K']+2
        path,_,_,_=h.expected_path(c);_,groups=h.tiling(c,path)
        segments=[s for group in groups.values() for s in group]
        contract={'certificate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'verifier_sha256':vh,'algorithm':'independent-segment-replay/1'}
        check=HERE/(name+'_segment_checkpoint.json');done={}
        if '--resume' in sys.argv and check.exists():
            old=json.loads(check.read_text());assert old['contract']==contract;done={int(k):v for k,v in old['segments'].items()}
        specs[name]={'cert':c,'contract':contract,'done':done,'count':len(segments),'path':check,'start':time.time()}
        for i,s in enumerate(segments):
            if i not in done:jobs.append({'name':name,'id':i,'segment':s,'ks':c['moments_k'],'function':c['function'],
                   'center':c['chebyshev_center'],'halfwidth':c['chebyshev_halfwidth'],'prec':c['method']['prec_bits'],
                   'K':c['method']['K'],'NDFT':c['method']['NDFT']})
    with ProcessPoolExecutor(max_workers=6,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(worker,j) for j in jobs]
        for future in as_completed(futures):
            r=future.result();assert r['status']=='PASS',r
            spec=specs[r['name']];spec['done'][r['id']]=r
            save(spec['path'],{'contract':spec['contract'],'segments':spec['done']})
            if len(spec['done'])%8==0:print(r['name'],len(spec['done']),'/',spec['count'],'segments freshly recomputed',flush=True)
    for name,spec in specs.items():
        c=spec['cert'];ctx.prec=c['method']['prec_bits'];totals=[acb(0) for _ in c['moments_k']]
        assert len(spec['done'])==spec['count']
        for i in range(spec['count']):totals=[a+acb(*b) for a,b in zip(totals,spec['done'][i]['enclosures'])]
        result=h.decide(c,totals,c['moments_k']);result.update(status='CERTIFIED',mode='full',
               segments_checked=spec['count'],segments_recomputed=spec['count'],contract=spec['contract'],
               elapsed_seconds=time.time()-spec['start'],scope=c['rectangle'],function=c['function'])
        save(HERE/(name+'_parallel_replay.json'),result);print(name,'CERTIFIED',result['count'],result.get('inertia'),flush=True)

if __name__=='__main__':main()
