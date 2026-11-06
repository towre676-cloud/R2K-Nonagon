"""Four-strength acquisition and exact certificate fusion on IEEE 118.

This experiment spends four pair readings instead of two. It measures how much
the additional responses tighten certified intervals, with no truth-based
selection. All intersections include the baseline two-strength certificate.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time

import numpy as np

from two_strength_cert import (Interval, certify_pair, certify_many, load_grid,
                               pair_delta, reconstruct_single,
                               verify_certificate, verify_many_certificate)


def run():
    samples,med,G00,meta,sha=load_grid('case118')
    scales=[F(1,10),F(1),F(2),F(10)]
    ts=[c/med for c in scales]
    out=dict(schema='nonagon.four_strength_design.v1',case='case118',pairs=len(samples),
             scales=list(map(str,scales)),pair_readings_baseline=2,pair_readings_extended=4,
             input_sha256=sha,arithmetic='exact rational interval endpoints; outward dyadic contraction',
             truth_use='external inclusion and relative-radius scoring only',rows=[],fixtures=[])
    started=time.perf_counter()
    for sigma in (F(1,10**8),F(1,10**6),F(1,10**4)):
        eps=2*sigma*abs(G00)
        for mode in ('calibrated','raw-single'):
            radii2=[];radii4=[];ratios=[];accepted=0;rejected=0
            for k,(i,j,u,v,A,B,x) in enumerate(samples):
                boxes=[Interval.ball(pair_delta(t,u,v,A,B,x)+eps*F(((k*17+ell*23)%9)-4,4),eps)
                       for ell,t in enumerate(ts)]
                inp=[u,v,A,B]
                try:
                    if mode=='raw-single':
                        def singles(g,d):
                            a=ts[1]*g*g/(1+ts[1]*d)
                            b=ts[2]*g*g/(1+ts[2]*d)
                            return reconstruct_single(ts[1],Interval.ball(a,eps),Interval.ball(b,eps))
                        Ai,ui=singles(u,A);Bj,vj=singles(v,B)
                        inp=[ui,vj,Ai,Bj]
                    X2=certify_pair(ts[1],ts[2],*inp,boxes[1],boxes[2],raw_interval=True)
                    X4=certify_many(ts,boxes,*inp,raw_interval=True)
                    assert X2.contains(x) and X4.contains(x) and X4.width<=X2.width
                    accepted+=1;radii2.append(float(X2.width/(2*abs(x))));radii4.append(float(X4.width/(2*abs(x))));ratios.append(float(X2.width/X4.width))
                    if k==0 and sigma==F(1,10**6):
                        c2=certify_pair(ts[1],ts[2],*inp,boxes[1],boxes[2])
                        c4=certify_many(ts,boxes,*inp)
                        assert verify_certificate(c2) and verify_many_certificate(c4)
                        out['fixtures'].append(dict(mode=mode,buses=[i,j],true_Gij=str(x),two=c2,four=c4))
                except (ZeroDivisionError,ValueError):
                    rejected+=1
            row=dict(sigma=str(sigma),calibration=mode,accepted=accepted,rejected=rejected,
                     truth_covered=accepted,median_radius_two=float(np.median(radii2)) if radii2 else None,
                     median_radius_four=float(np.median(radii4)) if radii4 else None,
                     median_per_pair_tightening=float(np.median(ratios)) if ratios else None,
                     p90_radius_two=float(np.quantile(radii2,.9)) if radii2 else None,
                     p90_radius_four=float(np.quantile(radii4,.9)) if radii4 else None)
            out['rows'].append(row)
            print(row,flush=True)
            # Save every completed row so publication code can consume progress.
            Path(__file__).with_name('two_strength_design_results.json').write_text(json.dumps(out,indent=2)+'\n')
    out['seconds']=time.perf_counter()-started
    out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(__file__).with_name('two_strength_design_results.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    run()
