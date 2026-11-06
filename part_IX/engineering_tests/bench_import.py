"""Import calibrated DC bench observations; never sends hardware commands.

The import checks file completeness, routing, settling and declared numerical
envelopes. It does not establish metrological calibration or physical model
validity from the readings. Conductance inference is conditional on those inputs.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction as F
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
CORE=HERE.parent if (HERE.parent/'bounded_intervention.py').exists() else HERE.parent/'part8_work/part_VIII'
sys.path.insert(0,str(CORE))
from bounded_intervention import Interval, infer_bounds, canonical_plan, load_key

FIELDS=['configuration','source_pad','sense_pad','current_A','current_relative_error',
        'conductance_scale_S','polarity_mode','readout_V','voltage_positive_V','voltage_negative_V',
        'voltage_absolute_error_V','polarity_offset_change_bound_V','drift_absolute_error_V',
        'shunt_A_pad','shunt_A_resistance_ohm','shunt_A_absolute_error_S',
        'shunt_A_relay_series_ohm','shunt_A_relay_error_ohm',
        'shunt_B_pad','shunt_B_resistance_ohm','shunt_B_absolute_error_S',
        'shunt_B_relay_series_ohm','shunt_B_relay_error_ohm',
        'settled','calibration_id','repeat_id','temperature_C']


def number(row,key,positive=False):
    text=row.get(key,'').strip()
    if not text:raise ValueError(f'Missing calibration or reading: {key}')
    try:value=F(text)
    except (ValueError,ZeroDivisionError) as e:raise ValueError(f'Invalid rational: {key}') from e
    if positive and value<=0:raise ValueError(f'{key} must be positive')
    return value


def response_interval(row):
    current=number(row,'current_A',True);rel=number(row,'current_relative_error')
    scale=number(row,'conductance_scale_S',True)
    err=number(row,'voltage_absolute_error_V');drift=number(row,'drift_absolute_error_V')
    if not 0<=rel<1 or err<0 or drift<0:raise ValueError('Invalid current/voltage uncertainty')
    mode=row.get('polarity_mode','').strip()
    if mode=='single':
        value=number(row,'readout_V');radius=err+drift;den=current
    elif mode=='reversed':
        value=number(row,'voltage_positive_V')-number(row,'voltage_negative_V')
        offset=number(row,'polarity_offset_change_bound_V')
        if offset<0:raise ValueError('Negative polarity offset-change bound')
        radius=2*err+offset+2*drift;den=2*current
    else:raise ValueError('polarity_mode must be single or reversed')
    if value-radius<=0:raise ValueError('Normalized response is not separated from zero')
    return Interval((value-radius)*scale/(den*(1+rel)),
                    (value+radius)*scale/(den*(1-rel)))


def load_interval(row,channel,scale):
    resistance=number(row,f'shunt_{channel}_resistance_ohm',True)
    conductance_error=number(row,f'shunt_{channel}_absolute_error_S')
    relay=number(row,f'shunt_{channel}_relay_series_ohm')
    relay_error=number(row,f'shunt_{channel}_relay_error_ohm')
    if conductance_error<0 or relay<0 or relay_error<0 or resistance+relay-relay_error<=0:
        raise ValueError('Invalid shunt/relay uncertainty')
    # Calibration error is an absolute conductance error before adding series
    # relay resistance. This transformation encloses the physical shunt path.
    base=Interval(1/resistance-conductance_error,1/resistance+conductance_error)
    if base.lo<=0:raise ValueError('Shunt calibration interval reaches zero conductance')
    series=Interval(relay-relay_error,relay+relay_error)
    return (base.reciprocal()+series).reciprocal()/scale


def template(n,path,mode='reversed'):
    rows=[]
    for index,loads in enumerate(canonical_plan(n)):
        row=dict.fromkeys(FIELDS,'')
        row.update(configuration='baseline' if not loads else 'load_'+str(index),
            source_pad='node_1',sense_pad='node_1',current_A='1/1000',conductance_scale_S='1/1000',
            polarity_mode=mode)
        for channel,(j,t) in zip(('A','B'),sorted(loads.items())):
            row[f'shunt_{channel}_pad']=f'node_{j+1}'
            row[f'shunt_{channel}_resistance_ohm']=str(1000/t)
        rows.append(row)
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader();writer.writerows(rows)
    return rows


def import_rows(n,rows,gap=None):
    plans=canonical_plan(n);expected={load_key(p):p for p in plans}
    observed={};relative_load_error=F(0);calibration_ids=set();input_readings=[]
    scale=None;configuration_labels={}
    for count,row in enumerate(rows,1):
        if row.get('source_pad')!='node_1' or row.get('sense_pad')!='node_1':
            raise ValueError('This protocol fixes source and readout at node_1')
        if row.get('settled','').strip().lower()!='true':raise ValueError('Unsettled acquisition')
        calibration=row.get('calibration_id','').strip()
        if not calibration:raise ValueError('A calibration record identifier is required')
        calibration_ids.add(calibration)
        this_scale=number(row,'conductance_scale_S',True)
        if scale is None:scale=this_scale
        if this_scale!=scale:raise ValueError('Conductance unit changed within acquisition')
        nominal={};measured_loads=[]
        for channel in ('A','B'):
            pad=row.get(f'shunt_{channel}_pad','').strip()
            if not pad:
                if row.get(f'shunt_{channel}_resistance_ohm','').strip():
                    raise ValueError('Shunt resistance supplied without a routed pad')
                continue
            if not pad.startswith('node_') or not pad[5:].isdigit():raise ValueError('Unknown routed pad')
            node=int(pad[5:])-1
            if not 1<=node<n or node in nominal:raise ValueError('Duplicate or forbidden shunt pad')
            resistance=number(row,f'shunt_{channel}_resistance_ohm',True)
            target=1/(resistance*scale)
            if target not in (F(1),F(2)):raise ValueError('Noncanonical nominal shunt setting')
            nominal[node]=target
            physical=load_interval(row,channel,scale)
            measured_loads.append({'site_zero_based':node,'nominal':str(target),'actual_interval':physical.json()})
            relative_load_error=max(relative_load_error,abs(physical.lo/target-1),abs(physical.hi/target-1))
        key=load_key(nominal)
        if key not in expected:raise ValueError('Unexpected pair setting; only base-load pairs are canonical')
        label=row.get('configuration','').strip()
        if not label:raise ValueError('Missing configuration identifier')
        if label in configuration_labels and configuration_labels[label]!=key:
            raise ValueError('A configuration identifier refers to different physical routes')
        configuration_labels[label]=key
        reading=response_interval(row)
        if key in observed:
            old=observed[key];lo=max(old.lo,reading.lo);hi=min(old.hi,reading.hi)
            if lo>hi:raise ValueError('Repeated configuration intervals disagree; check drift/calibration')
            observed[key]=Interval(lo,hi)
        else:observed[key]=reading
        input_readings.append({'row':count,'configuration':label,'normalized_interval':reading.json(),
                               'loads':measured_loads,'polarity_mode':row['polarity_mode']})
    if set(observed)!=set(expected):raise ValueError('Incomplete canonical acquisition')
    if relative_load_error>=1:raise ValueError('Shunt uncertainty exceeds the supported positive-load envelope')
    result=infer_bounds(n,observed,gap=gap,relative_load_error=relative_load_error)
    result['bench_import']={'schema':'nonagon.bench-import.v1','physical_measurement_rows':len(rows),
        'scalar_configurations':len(expected),'raw_voltage_samples':sum(2 if r['polarity_mode']=='reversed' else 1 for r in rows),
        'calibration_records':sorted(calibration_ids),'conductance_scale_S':str(scale),
        'relative_shunt_envelope':str(relative_load_error),'readings':input_readings,
        'physical_execution':'imported data only; origin and calibration are external declarations',
        'interpretation':'conditional bounds for a fixed connected grounded reciprocal resistor network'}
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--nodes',type=int,required=True)
    ap.add_argument('--template',type=Path);ap.add_argument('--input',type=Path)
    ap.add_argument('--output',type=Path);ap.add_argument('--gap')
    args=ap.parse_args()
    if args.template:
        args.template.parent.mkdir(parents=True,exist_ok=True);rows=template(args.nodes,args.template)
        print(json.dumps({'template':str(args.template),'configurations':len(rows),'voltage_samples_with_reversal':2*len(rows)}));return
    if not args.input or not args.output:ap.error('--input and --output are required for import')
    with args.input.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    result=import_rows(args.nodes,rows,F(args.gap) if args.gap else None)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'output':str(args.output)}))


if __name__=='__main__':main()
