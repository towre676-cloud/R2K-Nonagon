"""Frozen finite-range campaign using the independent direct contour engine."""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

from direct_engine import tile,sha,canonical,require,acb,arb

HERE=Path(__file__).resolve().parent


def write(path,obj):
    path.write_text(json.dumps(obj,indent=2))


def settings():
    return {'precision_bits':100,'tolerance_bits':65,'boundary_max_depth':18,
            'boundary_eval_limit':4096,'integral_eval_limit':12000,
            'integral_depth_limit':20,'max_count':6}


def protocol():
    return {'format':'nonagon-independent-direct-campaign-v1',
        'function':{'kind':'riemann_xi_t'},'range':['10','100'],'initial_width':'2','b':'3/5',
        'max_split_depth':6,'settings':settings(),
        'tile_placement':'fixed rational grid, no reference-zero input',
        'acceptance':'zero-free full boundary, integer contained+isolated, all positive Hermite leading minors or empty',
        'reference_policy':'no root-list use until all campaign certificates and coverage ledger are written',
        'engine_sha256':hashlib.sha256((HERE/'direct_engine.py').read_bytes()).hexdigest()}


def controls(p):
    specs=[('real_simple',{'kind':'polynomial','coefficients_ascending':['-1/4','0','1']},'certified_real_simple'),
        ('off_axis_pair',{'kind':'polynomial','coefficients_ascending':['1/16','0','1']},'certified_not_all_real_simple'),
        ('repeated_root',{'kind':'polynomial','coefficients_ascending':['0','0','1']},'unresolved'),
        ('boundary_root',{'kind':'polynomial','coefficients_ascending':['-1','1']},'unresolved'),
        ('close_real_pair',{'kind':'polynomial','coefficients_ascending':['0','-1/1000000','1']},'certified_real_simple')]
    out=[]
    for name,spec,expected in specs:
        c=tile(spec,{'lo':'-1','hi':'1','b':'1/2'},p['settings'])
        require(c['status']==expected, 'control status mismatch '+name)
        require(c.get('all_real_simple') is not True or name in ['real_simple','close_real_pair'],
                'negative control falsely certified')
        c['fixture_name']=name;c['expected_status']=expected
        write(HERE/'controls'/(name+'.json'),c)
        out.append({'name':name,'status':c['status'],'count':c.get('count'),
                    'reason':c.get('reason'),'seconds':c['seconds']})
    write(HERE/'controls_results.json',{'passed':len(out),'controls':out})
    return out


def campaign(p):
    attempts=[];leaves=[];started=time.perf_counter()
    def solve(lo,hi,depth,parent=None):
        ident='tile_%04d'%len(attempts)
        geometry={'lo':str(lo),'hi':str(hi),'b':p['b']}
        certificate=tile(p['function'],geometry,p['settings'])
        path=HERE/'tiles'/(ident+'.json');write(path,certificate)
        record={'id':ident,'parent':parent,'depth':depth,'lo':str(lo),'hi':str(hi),
                'b':p['b'],'status':certificate['status'],'count':certificate.get('count'),
                'reason':certificate.get('reason'),'seconds':certificate['seconds'],
                'certificate':str(path.relative_to(HERE)),
                'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'boundary_value_calls':certificate['boundary_value_calls'],
                'logarithmic_derivative_calls':certificate['logarithmic_derivative_calls']}
        attempts.append(record)
        write(HERE/'coverage_in_progress.json',{'protocol_sha256':sha(p),'attempts':attempts,'leaves':leaves})
        print(json.dumps(record),flush=True)
        accepted=certificate['status'] in ['certified_empty','certified_real_simple']
        if accepted or depth>=p['max_split_depth'] or certificate['status']=='certified_not_all_real_simple':
            leaves.append(record);return
        mid=(lo+hi)/2
        solve(lo,mid,depth+1,ident);solve(mid,hi,depth+1,ident)
    low,high=map(Fraction,p['range']);width=Fraction(p['initial_width']);left=low
    while left<high:
        right=min(left+width,high);solve(left,right,0);left=right
    sorted_leaves=sorted(leaves,key=lambda x:Fraction(x['lo']))
    end=low
    for r in sorted_leaves:
        require(Fraction(r['lo'])==end,'campaign coverage gap or overlap')
        end=Fraction(r['hi'])
    require(end==high,'campaign range incomplete')
    unresolved=[r for r in sorted_leaves if r['status'] not in ['certified_empty','certified_real_simple']]
    result={'protocol_sha256':sha(p),'engine_sha256':p['engine_sha256'],
            'exact_range':p['range'],'b':p['b'],
            's_plane_real_range':['-1/10','11/10'],
            'continuous_coverage':True,'all_regions_certified_real_simple_or_empty':not unresolved,
            'total_zero_count':sum(r['count'] for r in sorted_leaves) if not unresolved else None,
            'attempts':attempts,'leaves':sorted_leaves,'unresolved_leaves':unresolved,
            'seconds':time.perf_counter()-started,
            'interpretation':'local finite-height verification on the stated continuous range; not a new zero-height record or global RH proof'}
    write(HERE/'campaign_results.json',result)
    # This diagnostic is deliberately executed only after ledger completion.
    if not unresolved:
        n=result['total_zero_count'];refs=[v.imag for v in acb.zeta_zeros(1,n+2)]
        inside=[v for v in refs if v>arb(str(low)) and v<arb(str(high))]
        ref={'policy':'post-campaign diagnostic; no values used for placement/integration',
             'reference_method':'Arb zeta_zeros','count_inside_range':len(inside),
             'agrees_with_contour_count':len(inside)==n,
             'first_inside':inside[0].str(30,radius=True) if inside else None,
             'last_inside':inside[-1].str(30,radius=True) if inside else None}
        write(HERE/'reference_check.json',ref)
    return result


def lehmer(p,timeout):
    geometry={k:str(Fraction.from_float(v)) for k,v in [('lo',7004.95),('hi',7005.22),('b',.45)]}
    request={'function':p['function'],'geometry':geometry,'settings':p['settings']}
    write(HERE/'lehmer_request.json',request)
    started=time.perf_counter()
    try:
        process=subprocess.run([sys.executable,str(HERE/'worker.py'),str(HERE/'lehmer_request.json'),str(HERE/'lehmer_direct.json')],
                               capture_output=True,text=True,timeout=timeout,check=False)
        if process.returncode:
            result={'completed':False,'status':'worker_error','returncode':process.returncode,'stderr':process.stderr}
        else:
            c=json.loads((HERE/'lehmer_direct.json').read_text())
            result={'completed':True,'status':c['status'],'count':c.get('count'),
                    'reason':c.get('reason'),'gap_from_determinant':c.get('gap_from_determinant'),
                    'source_sha256':c['source_sha256']}
            # A reference is diagnostic after direct contour computation.
            if c.get('gap_from_determinant'):
                ref=acb.zeta_zero(6710).imag-acb.zeta_zero(6709).imag
                result['reference_gap']=ref.str(35,radius=True)
                result['gap_overlaps_reference']=bool(arb(c['gap_from_determinant']).overlaps(ref))
    except subprocess.TimeoutExpired as e:
        result={'completed':False,'status':'timeout','budget_seconds':timeout,
                'limitation':'direct derivative/interval integration at this height did not finish within the fixed wall-time budget'}
    result['seconds']=time.perf_counter()-started
    result['request_sha256']=sha(request)
    write(HERE/'lehmer_results.json',result);return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--lehmer-timeout',type=int,default=60);parser.add_argument('--skip-lehmer',action='store_true');args=parser.parse_args()
    for d in ['controls','tiles']:(HERE/d).mkdir(exist_ok=True)
    p=protocol();write(HERE/'protocol.json',p)
    # The protocol and engine hash are frozen before controls/campaign execute.
    write(HERE/'protocol_frozen.json',{'protocol_sha256':sha(p),'engine_sha256':p['engine_sha256']})
    c=controls(p);r=campaign(p)
    lp={'completed':False,'status':'not_requested'} if args.skip_lehmer else lehmer(p,args.lehmer_timeout)
    environment={'python':sys.version,'platform':platform.platform(),'python_flint':'0.9.0','precision_bits':100,'series_cap':8}
    result={'protocol_sha256':sha(p),'environment':environment,
            'controls_passed':len(c),'campaign_complete':r['all_regions_certified_real_simple_or_empty'],
            'campaign_range':r['exact_range'],'campaign_zero_count':r['total_zero_count'],
            'campaign_attempts':len(r['attempts']),'campaign_leaves':len(r['leaves']),
            'campaign_seconds':r['seconds'],'lehmer':lp,
            'trust_boundary':'independently implemented direct contour method; trusted Arb function balls and adaptive integrator, no archive Taylor/DFT code'}
    write(HERE/'results.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
