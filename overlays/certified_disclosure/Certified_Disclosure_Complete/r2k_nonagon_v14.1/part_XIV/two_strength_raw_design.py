"""Equal-reading-budget comparison of propagated and raw-cancelled inversion."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time

import numpy as np

from two_strength_cert import (Interval, load_grid, pair_delta, certify_raw_pair,
                               verify_raw_certificate)


def run():
    samples,med,G00,meta,sha=load_grid('case118')
    out=dict(schema='nonagon.raw_cancellation_design.v1',case='case118',pairs=len(samples),
             input_sha256=sha,reading_budget='four single-node and two pair response differences; shared baseline',
             scope='positive calibrated port-column signs from the connected passive M-matrix',rows=[],fixtures=[])
    started=time.perf_counter()
    for c in (F(1,10),F(1),F(10)):
        theta=c/med
        for sigma in (F(1,10**8),F(1,10**6),F(1,10**4)):
            eps=2*sigma*abs(G00);old=[];new=[];ratios=[];rejected=0
            for k,(i,j,u,v,A,B,x) in enumerate(samples):
                def s(t,col,diag):return Interval.ball(t*col*col/(1+t*diag),eps)
                d1=Interval.ball(pair_delta(theta,u,v,A,B,x)+eps*F((17*k%9)-4,4),eps)
                d2=Interval.ball(pair_delta(2*theta,u,v,A,B,x)+eps*F((23*k%9)-4,4),eps)
                try:
                    cert=certify_raw_pair(theta,s(theta,u,A),s(2*theta,u,A),s(theta,v,B),s(2*theta,v,B),d1,d2)
                    X=Interval(*map(F,cert['enclosure']));Y=Interval(*map(F,cert['propagated']))
                    assert X.contains(x) and Y.contains(x) and X.width<=Y.width
                    old.append(float(Y.width/(2*abs(x))));new.append(float(X.width/(2*abs(x))));ratios.append(float(Y.width/X.width))
                    if k==0 and sigma==F(1,10**6):
                        assert verify_raw_certificate(cert)
                        out['fixtures'].append(dict(strength_scale=str(c),true_Gij=str(x),certificate=cert))
                except (ValueError,ZeroDivisionError):rejected+=1
            row=dict(strength_scale=str(c),sigma=str(sigma),certified=len(new),rejected=rejected,
                     truth_covered=len(new),median_radius_before=float(np.median(old)) if old else None,
                     median_radius_after=float(np.median(new)) if new else None,
                     p90_radius_before=float(np.quantile(old,.9)) if old else None,
                     p90_radius_after=float(np.quantile(new,.9)) if new else None,
                     median_per_pair_tightening=float(np.median(ratios)) if ratios else None,
                     within_one_percent_before=sum(r<=.01 for r in old),within_one_percent_after=sum(r<=.01 for r in new))
            out['rows'].append(row);print(row,flush=True)
            Path(__file__).with_name('two_strength_raw_design_results.json').write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.perf_counter()-started
    out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(__file__).with_name('two_strength_raw_design_results.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    run()
