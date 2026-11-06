"""Check coverage/integrity; optionally recompute every direct-Arb tile.

Static checks inspect stored enclosures but do not prove their provenance.
--recompute reruns the independent direct engine; Arb remains trusted.
"""
import argparse
import copy
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path
from direct_engine import tile,require,sha,arb,acb,ctx,deserialize_exact

HERE=Path(__file__).resolve().parent


def load():
    p=json.loads((HERE/'protocol.json').read_text())
    ledger=json.loads((HERE/'campaign_results.json').read_text())
    blobs={r['certificate']:(HERE/r['certificate']).read_bytes() for r in ledger['leaves']}
    return p,ledger,blobs


def validate(p,ledger,blobs):
    ctx.prec=p['settings']['precision_bits']
    engine=hashlib.sha256((HERE/'direct_engine.py').read_bytes()).hexdigest()
    require(p['engine_sha256']==engine,'engine identity mismatch')
    require(ledger['engine_sha256']==engine,'ledger engine identity mismatch')
    require(ledger['protocol_sha256']==sha(p),'protocol identity mismatch')
    require(ledger['exact_range']==p['range'],'range claim differs from protocol')
    require(ledger['b']==p['b'],'height scope differs from protocol')
    low,high=map(Fraction,p['range']);last=low;total=0
    for row in ledger['leaves']:
        lo,hi=Fraction(row['lo']),Fraction(row['hi'])
        require(lo==last and lo<hi,'coverage gap, overlap, or invalid interval');last=hi
        raw=blobs[row['certificate']]
        require(hashlib.sha256(raw).hexdigest()==row['certificate_sha256'],'certificate identity mismatch')
        c=json.loads(raw)
        require(c['source_sha256']==engine,'certificate engine identity mismatch')
        require(c['function']==p['function'],'different analytic function')
        geometry={k:row[k] for k in ['lo','hi','b']}
        require(c['geometry']==geometry,'geometry mismatch')
        require(c['input_sha256']==sha({'function':p['function'],'geometry':geometry,'settings':p['settings']}),'input identity mismatch')
        require(c['precision_bits']==p['settings']['precision_bits'] and c['series_cap']>=2,'numerical context mismatch')
        require(c['status']==row['status'] and c['count']==row['count'],'ledger/certificate outcome mismatch')
        require(c['status'] in ['certified_empty','certified_real_simple'],'unresolved or negative region')
        require(c['boundary']['proved_zero_free'],'boundary not proved zero-free')
        for side in range(4):
            end=Fraction(0)
            for rec in c['boundary']['segments']:
                if rec['side']!=side:continue
                require(Fraction(rec['tau0'])==end,'boundary coverage gap')
                end=Fraction(rec['tau1'])
                require(end>Fraction(rec['tau0']) and end<=1,'boundary partition order')
                z=deserialize_exact(rec['function_value'])
                require(z.is_finite() and not z.contains(0),'stored boundary ball does not exclude zero')
            require(end==1,'boundary partition incomplete')
        moment=c['moments'][0];real=arb(moment['real']);imag=arb(moment['imag']);n=c['count']
        require(real.contains(n) and (real-n).abs_upper()<arb('0.25'),'integer count not contained and isolated')
        require(imag.contains(0) and imag.abs_upper()<arb('0.25'),'count not enclosed as real')
        if n:
            require(len(c['leading_minors'])==n and all(arb(x)>0 for x in c['leading_minors']),'positive minors not enclosed')
        total+=n
    require(last==high,'range incomplete')
    require(total==ledger['total_zero_count'],'total count mismatch')
    return {'ledger_consistent':True,'continuous_exact_range':p['range'],'b':p['b'],
            'leaf_tiles':len(ledger['leaves']),'total_zero_count':total,
            'scope':'static integrity, coverage and stored-enclosure consistency; not an independent proof of the saved interval evaluations'}


def attacks(p,ledger,blobs):
    out=[]
    for name in ['missing_tile','duplicate_tile','wrong_total_count','wrong_source_hash','forged_count_ball','missing_boundary_segment','zero_containing_boundary_ball']:
        q,l,b=copy.deepcopy(p),copy.deepcopy(ledger),copy.deepcopy(blobs)
        if name=='missing_tile':l['leaves'].pop(0)
        if name=='duplicate_tile':l['leaves'].insert(1,copy.deepcopy(l['leaves'][0]))
        if name=='wrong_total_count':l['total_zero_count']+=1
        if name=='wrong_source_hash':q['engine_sha256']='0'*64
        if name in ['forged_count_ball','missing_boundary_segment','zero_containing_boundary_ball']:
            row=l['leaves'][0];c=json.loads(b[row['certificate']])
            if name=='forged_count_ball':c['moments'][0]['real']='999'
            if name=='missing_boundary_segment':c['boundary']['segments'].pop(0)
            if name=='zero_containing_boundary_ball':
                z={'mid_mantissa':'0','mid_exponent':0,'rad_mantissa':'1','rad_exponent':0}
                c['boundary']['segments'][0]['function_value']={'real':'0 +/- 1','imag':'0 +/- 1','exact_real':z,'exact_imag':z}
            raw=json.dumps(c,indent=2).encode();b[row['certificate']]=raw
            row['certificate_sha256']=hashlib.sha256(raw).hexdigest()
        try:validate(q,l,b)
        except (ValueError,KeyError) as e:out.append({'attack':name,'rejected':True,'reason':str(e)})
        else:raise RuntimeError('tamper accepted '+name)
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--recompute',action='store_true');args=parser.parse_args()
    p,l,b=load();result=validate(p,l,b);result['tamper_checks']=attacks(p,l,b)
    if args.recompute:
        started=time.perf_counter();rows=[]
        for old in l['leaves']:
            geometry={k:old[k] for k in ['lo','hi','b']}
            c=tile(p['function'],geometry,p['settings'])
            same=c['status']==old['status'] and c.get('count')==old['count']
            row={'id':old['id'],'geometry':geometry,'status':c['status'],'count':c.get('count'),
                 'agrees':same,'reason':c.get('reason'),'seconds':c['seconds']}
            rows.append(row);print(json.dumps(row),flush=True)
            require(same,'direct replay disagrees with stored campaign')
        result['fresh_replay']={'tiles':rows,'all_agree':all(r['agrees'] for r in rows),'seconds':time.perf_counter()-started,
            'trust_boundary':'same independent direct engine and Arb rerun; no archive Taylor/DFT kernel involved'}
    (HERE/'replay_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
