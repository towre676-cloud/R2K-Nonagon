#!/usr/bin/env python3
"""External-input controls for the finite-load planner.

The mathematical derivative checker remains independent of planner.py. This
separate interface check deliberately calls its public recommendation function.
"""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('planner_under_interface_test',HERE/'planner.py')
planner=importlib.util.module_from_spec(spec);spec.loader.exec_module(planner)


def fixture():
    h=F(5,18); rows=[{'nominal_loads':{},'normalized_response_interval':[str(h),str(h)]}]
    for site,b,q in [(1,F(4,9),F(1,81)),(2,F(11,18),F(1,324))]:
        for t in [F(1),F(2)]:
            value=h-t*q/(1+t*b)
            rows.append({'nominal_loads':{str(site):str(t)},'normalized_response_interval':[str(value),str(value)]})
    return {'schema':'nonagon.finite-load-pilot.v1','network_order':3,
            'input_and_readout_site_zero_based':0,'pilot_base_load':'1',
            'hardware':{'conductance_scale_S':'1/1000','minimum_shunt_resistance_ohm':'1',
                        'maximum_shunt_resistance_ohm':'100000','simultaneous_shunt_channels':2,
                        'candidate_base_loads':['1/4','1','2','4'],'voltage_error_V':'1/1000000000',
                        'drift_error_V':'0','current_A':'1/1000','current_relative_error':'0',
                        'shunt_relative_error':'0','relay_series_max_ohm':'0'},'readings':rows}


def set_path(doc,path,value):
    at=doc
    for key in path[:-1]:at=at[key]
    at[path[-1]]=value


def run():
    base=fixture();result=planner.recommend(base)
    if result['status']!='recommended':raise RuntimeError('Exact valid pilot did not recommend')
    if F(result['selected_base_load'])!=2:raise RuntimeError('Finite candidate choice disagrees with independent exact objective')
    singles={r['site_zero_based']:r for r in result['pilot_single_enclosures']}
    for site,b,q in [(1,F(4,9),F(1,81)),(2,F(11,18),F(1,324))]:
        if list(map(F,singles[site]['b_interval']))!=[b,b]:raise RuntimeError('Exact pilot b altered')
        if list(map(F,singles[site]['q_interval']))!=[q,q]:raise RuntimeError('Exact pilot q altered')
    controls=[
        ('boolean-order',['network_order'],True),('float-order',['network_order'],3.0),
        ('boolean-port',['input_and_readout_site_zero_based'],False),('float-port',['input_and_readout_site_zero_based'],0.0),
        ('wrong-port',['input_and_readout_site_zero_based'],1),('float-channel-count',['hardware','simultaneous_shunt_channels'],2.0),
        ('missing-fallback',['hardware','candidate_base_loads'],['2','4']),('empty-candidates',['hardware','candidate_base_loads'],[]),
        ('boolean-candidate',['hardware','candidate_base_loads'],['1',True]),
        ('negative-voltage-radius',['hardware','voltage_error_V'],'-1/1000000000'),
        ('negative-drift-radius',['hardware','drift_error_V'],'-1/10000000000'),
        ('negative-current',['hardware','current_A'],'-1'),('unbounded-current',['hardware','current_relative_error'],'1'),
        ('negative-shunt-radius',['hardware','shunt_relative_error'],'-1/1000000'),
        ('negative-relay',['hardware','relay_series_max_ohm'],'-1/1000'),
        ('invalid-resistance-range',['hardware','minimum_shunt_resistance_ohm'],'200000'),
        ('double-load-resistance-violation',['hardware','minimum_shunt_resistance_ohm'],'200'),
        ('reversed-response-interval',['readings',0,'normalized_response_interval'],['2','1']),
        ('boolean-response',['readings',0,'normalized_response_interval'],[True,True]),
        ('float-response',['readings',0,'normalized_response_interval'],[0.2,0.3]),
        ('extra-zero-load',['readings',1,'nominal_loads'],{'1':'1','99':'0'}),
        ('extra-pair-load',['readings',1,'nominal_loads'],{'1':'1','2':'1'}),
    ]
    cases=[]
    for name,path,value in controls:
        changed=deepcopy(base);set_path(changed,path,value)
        if name=='negative-voltage-radius':changed['hardware']['drift_error_V']='1/100000000'
        try:planner.recommend(changed)
        except (ValueError,TypeError,KeyError,OverflowError):cases.append({'case':name,'rejected':True})
        else:cases.append({'case':name,'rejected':False})
    for name,modify in [
        ('duplicate-reading',lambda d:d['readings'].append(deepcopy(d['readings'][0]))),
        ('missing-reading',lambda d:d['readings'].pop()),
        ('hidden-operator',lambda d:d.update({'K':[['4','-1','0'],['-1','3','-1'],['0','-1','2']]})),
    ]:
        changed=deepcopy(base);modify(changed)
        try:planner.recommend(changed)
        except (ValueError,TypeError,KeyError,OverflowError):cases.append({'case':name,'rejected':True})
        else:cases.append({'case':name,'rejected':False})
    failed=[r['case'] for r in cases if not r['rejected']]
    return {'schema':'nonagon.load-planner-input-controls.v1','status':'passed' if not failed else 'failed',
            'valid_exact_pilot':'passed','expected_candidate':'2','controls':cases,'rejected':sum(r['rejected'] for r in cases),
            'control_count':len(cases),'failed_controls':failed,
            'source_sha256':{'planner.py':hashlib.sha256((HERE/'planner.py').read_bytes()).hexdigest(),
                             Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=HERE/'planner_input_controls.json');args=ap.parse_args()
    result=run();args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if result['status']!='passed':raise SystemExit(1)
