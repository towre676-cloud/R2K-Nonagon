"""Independent synthetic-CSV and calibration-corner tests for bench_import.

This exercises the importer under test. Circuit voltages and physical envelopes
are recomputed by verify_engineering.gaussian, not by the tomography solver.
No physical instrument is connected.
"""
from __future__ import annotations
from copy import deepcopy
import csv
from fractions import Fraction as F
import json
from pathlib import Path

import bench_import as target
from verify_engineering import Checks, q, response, within, check_interval_result

HERE=Path(__file__).resolve().parent
K=[[F(4),F(-1),F(0)],[F(-1),F(3),F(-1)],[F(0),F(-1),F(2)]]


def make_rows(path,uncertain=False):
    rows=target.template(3,path)
    evidence=[]
    for row in rows:
        row.update(current_relative_error='0' if not uncertain else '1/100000',
                   voltage_absolute_error_V='1/1000000000000' if not uncertain else '1/10000000',
                   drift_absolute_error_V='0' if not uncertain else '1/100000000',
                   polarity_offset_change_bound_V='0',settled='true',
                   calibration_id='synthetic-reference-calibration',repeat_id='1',temperature_C='23')
        nominal={};actual={}
        for channel in ('A','B'):
            pad=row['shunt_'+channel+'_pad']
            if not pad:continue
            node=int(pad[5:])-1;res=q(row['shunt_'+channel+'_resistance_ohm'])
            row['shunt_'+channel+'_absolute_error_S']='0' if not uncertain else str((1/res)*F(1,10000))
            row['shunt_'+channel+'_relay_series_ohm']='0' if not uncertain else '1/20'
            row['shunt_'+channel+'_relay_error_ohm']='0' if not uncertain else '1/100'
            nominal[node]=1/(res*F(1,1000))
            resistance=res if not uncertain else res/F(1+F(1,20000))+F(1,20)
            actual[node]=1/(resistance*F(1,1000))
        current=F(1,1000) if not uncertain else F(1,1000)*(1+F(1,200000))
        truth=response(K,actual)
        signal=truth*current/F(1,1000)
        offset=F(1,10000)
        # Constant thermal EMF cancels through reversed current polarity.
        row['voltage_positive_V']=str(signal+offset)
        row['voltage_negative_V']=str(-signal+offset)
        evidence.append({'configuration':row['configuration'],'nominal':{str(j):str(v) for j,v in nominal.items()},
                         'actual':{str(j):str(v) for j,v in actual.items()},'current':str(current),
                         'normalized_truth':str(truth),'constant_emf_V':str(offset)})
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=target.FIELDS);w.writeheader();w.writerows(rows)
    return rows,evidence


def main():
    checks=Checks();report={'schema':'nonagon.bench-import.controls.v1','physical_hardware_measurements':0}
    accepted=[]
    for uncertain in (False,True):
        name='uncertain' if uncertain else 'exact_reference'
        path=HERE/f'bench_example_{name}.csv'
        rows,evidence=make_rows(path,uncertain)
        result=target.import_rows(3,rows,gap=F(1,4))
        for row,truth in zip(rows,evidence):
            reading=target.response_interval(row)
            checks.require(reading.contains(q(truth['normalized_truth'])),'reverse normalization contains actual circuit response')
            for channel in ('A','B'):
                pad=row['shunt_'+channel+'_pad']
                if pad:
                    node=int(pad[5:])-1
                    checks.require(target.load_interval(row,channel,F(1,1000)).contains(q(truth['actual'][str(node)])),
                                   'relay and conductance envelope contains actual shunt')
        check_interval_result(checks,K,result,'synthetic CSV '+name)
        checks.require(result['bench_import']['raw_voltage_samples']==12,'reversal samples explicitly counted')
        checks.require(result['bench_import']['scalar_configurations']==6,'canonical scalar settings count')
        accepted.append({'fixture':name,'status':result['status'],'csv':path.name,'evidence':evidence,
                         'result':result})
    rows,_=make_rows(HERE/'bench_example_exact_reference.csv',False)
    rejected=[]
    mutations=[]
    def edit(label,index,key,value):
        d=deepcopy(rows);d[index][key]=value;mutations.append((label,d))
    edit('missing calibration',0,'calibration_id','')
    edit('unsettled reading',0,'settled','false')
    edit('wrong actuator site',0,'source_pad','node_2')
    edit('wrong readout site',0,'sense_pad','node_2')
    edit('negative voltage error',0,'voltage_absolute_error_V','-1')
    edit('current envelope reaches zero',0,'current_relative_error','1')
    edit('negative drift bound',0,'drift_absolute_error_V','-1')
    edit('negative polarity-change bound',0,'polarity_offset_change_bound_V','-1')
    edit('unsupported polarity mode',0,'polarity_mode','uncalibrated')
    edit('missing negative-polarity reading',0,'voltage_negative_V','')
    edit('forbidden shunt at source',1,'shunt_A_pad','node_1')
    edit('unknown shunt site',1,'shunt_A_pad','node_99')
    edit('noncanonical load strength',1,'shunt_A_resistance_ohm','750')
    edit('negative relay bound',1,'shunt_A_relay_error_ohm','-1')
    edit('conductance calibration reaches zero',1,'shunt_A_absolute_error_S','1')
    edit('conductance units changed',1,'conductance_scale_S','1/2000')
    edit('missing configuration label',1,'configuration','')
    edit('same label different route',1,'configuration',rows[0]['configuration'])
    d=deepcopy(rows);d[-1]['shunt_B_pad']=d[-1]['shunt_A_pad'];mutations.append(('duplicate pair shunt pad',d))
    d=deepcopy(rows);d.pop();mutations.append(('missing canonical configuration',d))
    d=deepcopy(rows);d.append(deepcopy(rows[0]));d[-1]['voltage_positive_V']='1';mutations.append(('inconsistent repeated baseline',d))
    d=deepcopy(rows);d[0]['shunt_A_resistance_ohm']='1000';mutations.append(('unrouted resistor setting',d))
    for label,data in mutations:
        try:target.import_rows(3,data,gap=F(1,4))
        except (ValueError,KeyError) as e:rejected.append({'control':label,'exception':type(e).__name__})
        else:raise ValueError('Expected importer rejection: '+label)
        checks.require(True,'importer rejected '+label)
    report.update(status='passed',checks=checks.count,accepted=accepted,rejected=rejected)
    (HERE/'bench_import_controls.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','checks':checks.count,'accepted':[x['status'] for x in accepted],
                      'rejected':len(rejected),'physical_hardware_measurements':0}))


if __name__=='__main__':main()
