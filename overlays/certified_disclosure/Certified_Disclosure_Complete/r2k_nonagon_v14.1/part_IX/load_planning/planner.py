#!/usr/bin/env python3
"""Select a finite, globally shared t from calibrated pilot scalar intervals.

This is an acquisition recommendation. The objective is a local first-order
sensitivity diagnostic; only the production nonlinear interval gate can certify
the final matrix. No hidden operator, inferred topology or oracle is accepted.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

CORE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CORE))
from bounded_intervention import Interval, Unresolved, load_key, canonical_plan


def q(value):
    if isinstance(value, (float, bool)):
        raise ValueError('Use rational strings or integers, never floats or booleans')
    return F(value)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def calibration_radius(hardware, t):
    """Symmetric relative envelope for either independent shunt up to 2t."""
    e = q(hardware['shunt_relative_error'])
    relay = q(hardware['relay_series_max_ohm'])
    scale = q(hardware['conductance_scale_S'])
    require(0 <= e < 1 and relay >= 0 and scale > 0,
            'Invalid shunt, relay or unit-scale contract')
    smallest_ratio = (1-e)/(1+relay*2*t*scale*(1-e))
    return max(e, 1-smallest_ratio)


def _positive_box(box, name):
    if box.lo <= 0:
        raise Unresolved(name+' is not separated from zero')
    return box


def recommend(document):
    require(isinstance(document,dict) and set(document)=={
        'schema','network_order','input_and_readout_site_zero_based','pilot_base_load','hardware','readings'},
        'Unknown or missing pilot document fields')
    require(document.get('schema') == 'nonagon.finite-load-pilot.v1',
            'Unsupported pilot schema')
    n = document.get('network_order')
    require(type(n) is int and 2 <= n <= 100, 'Known labeled order 2..100 required')
    require(type(document.get('input_and_readout_site_zero_based')) is int
            and document['input_and_readout_site_zero_based'] == 0,
            'The calibrated pilot must use the fixed source/readout port 0')
    t0 = q(document['pilot_base_load'])
    require(t0 > 0, 'Positive pilot base load required')
    hardware = document['hardware']
    scale = q(hardware['conductance_scale_S'])
    minimum_r = q(hardware['minimum_shunt_resistance_ohm'])
    maximum_r = q(hardware['maximum_shunt_resistance_ohm'])
    require(0 < minimum_r <= maximum_r and scale > 0, 'Invalid resistance range')
    require(type(hardware.get('simultaneous_shunt_channels')) is int
            and hardware['simultaneous_shunt_channels'] == 2,
            'This canonical plan requires two independently routed shunt channels')
    candidates = sorted(set(q(x) for x in hardware['candidate_base_loads']))
    require(bool(candidates) and t0 in candidates, 'Candidates must include the pilot fallback')
    for t in candidates:
        require(t > 0 and minimum_r <= 1/(2*t*scale)
                and 1/(t*scale) <= maximum_r,
                'A candidate violates the declared resistance range at t or 2t')
    expected = {(): None}
    for j in range(1,n):
        expected[load_key({j:t0})] = None
        expected[load_key({j:2*t0})] = None
    observations = {}
    for row in document['readings']:
        require(isinstance(row,dict) and set(row)=={'nominal_loads','normalized_response_interval'},
                'Unknown or missing pilot reading fields')
        loads=row['nominal_loads']
        require(isinstance(loads,dict) and len(loads)<=1,'Pilot uses one hidden-site shunt at a time')
        clean={}
        for site,value in loads.items():
            require(type(site) is int or (isinstance(site,str) and site.isdecimal()
                    and str(int(site))==site),'Site indices must be canonical integers')
            j=int(site)
            require(1<=j<n,'Pilot load site is outside the labeled hidden-site set')
            strength=q(value)
            require(strength in (t0,2*t0),'Pilot strengths must be t0 or 2t0')
            clean[j]=strength
        key = load_key(clean)
        require(key in expected and key not in observations,
                'Pilot contains extra, repeated or unrecognized load configuration')
        values = row['normalized_response_interval']
        require(isinstance(values,list) and len(values)==2, 'Each reading needs a scalar interval')
        observations[key] = Interval(q(values[0]), q(values[1]))
    require(observations.keys() == expected.keys(), 'Complete baseline and two-load pilot required')
    base = observations[()]
    require(base.lo > 0, 'Pilot baseline cannot certify a positive response')
    voltage_error=q(hardware['voltage_error_V']);drift_error=q(hardware['drift_error_V'])
    require(voltage_error>=0 and drift_error>=0,'Each readout/drift radius must be nonnegative')
    read_error = voltage_error + drift_error
    current = q(hardware['current_A'])
    current_relative = q(hardware['current_relative_error'])
    require(read_error >= 0 and current > 0 and 0 <= current_relative < 1,
            'Invalid current/readout contract')
    # Every positive grounding load decreases h. Thus the pilot baseline upper
    # bound controls multiplicative source uncertainty for every future read.
    eta = read_error*scale/current + current_relative*base.hi
    common = {'schema':'nonagon.finite-load-recommendation.v1',
              'network_order':n, 'pilot_base_load':str(t0),
              'canonical_acquisition_count':n*(n+1)//2,
              'pilot_scalar_count':2*n-1,
              'total_scalar_count_with_pilot':2*n-1+n*(n+1)//2,
              'maximum_simultaneous_shunts':2,
              'decision_scope':'local first-order q sensitivity; final nonlinear interval gate remains authoritative',
              'error_model':'independent bounded raw h0,h(t),h(2t) errors with common measured baseline',
              'raw_response_error_bound':str(eta),
              'hardware':hardware}
    def finish(result):
        selected=q(result['selected_base_load'])
        result['nominal_configuration_sequence']=[
            {'loads':{str(j):str(v) for j,v in loads.items()}}
            for loads in canonical_plan(n,selected,2*selected)]
        return result
    singles = []
    try:
        e0 = calibration_radius(hardware, t0)
        u1 = Interval(t0*(1-e0), t0*(1+e0))
        u2 = Interval(2*t0*(1-e0), 2*t0*(1+e0))
        separation = _positive_box(u2-u1, 'Pilot load calibration separation')
        for j in range(1,n):
            d1 = _positive_box(base-observations[load_key({j:t0})], 'First pilot drop')
            d2 = _positive_box(base-observations[load_key({j:2*t0})], 'Second pilot drop')
            y1,y2 = u1/d1,u2/d2
            slope = (y2-y1)/separation
            intercept = _positive_box(y1-u1*slope, 'Pilot single-load intercept')
            qb = _positive_box(1/intercept, 'Pilot q')
            bb = (slope/intercept).positive_class_intersection('Pilot b')
            singles.append({'site_zero_based':j,'q_interval':qb.json(),'b_interval':bb.json()})
    except Unresolved as exc:
        return finish(dict(common, status='fallback', selected_base_load=str(t0),
                    reason=str(exc), pilot_single_enclosures=singles, candidates=[]))
    scores = []
    for t in candidates:
        e = calibration_radius(hardware,t)
        rows=[]
        for single in singles:
            qlo = q(single['q_interval'][0]); bhi = q(single['b_interval'][1])
            raw = eta*(4/t+8*bhi+4*t*bhi*bhi)/qlo
            # At nominal t,2t, differentiating the extracted q with respect to
            # the two actual load amplitudes gives 2q/t and -q/(2t).
            # Two independently bounded relative load errors therefore add 3e.
            rows.append({'site_zero_based':single['site_zero_based'],
                         'raw_reading_local_relative_q_bound':str(raw),
                         'load_calibration_local_relative_q_bound':str(3*e),
                         'total_local_relative_q_diagnostic':str(raw+3*e)})
        score=max(q(r['total_local_relative_q_diagnostic']) for r in rows)
        scores.append({'base_load':str(t),'worst_site_local_diagnostic':str(score),
                       'relative_load_envelope':str(e),'sites':rows})
    # Fixed deterministic tie-break prefers the existing pilot setting, then
    # lower t. This matters in an ideal exact campaign where all scores are zero.
    best=min(scores,key=lambda r:(q(r['worst_site_local_diagnostic']),
                                q(r['base_load']) != t0,q(r['base_load'])))
    return finish(dict(common,status='recommended',selected_base_load=best['base_load'],
                pilot_single_enclosures=singles,candidates=scores,
                selection='minimize the largest pilot-enclosed local q diagnostic over the finite hardware list'))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('pilot',type=Path);ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    result=recommend(json.loads(args.pilot.read_text()))
    out=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(out)
    else:print(out,end='')


if __name__=='__main__':main()
